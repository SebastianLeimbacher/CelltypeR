# compare between groups

# plotting function 

proportionplots <- function(seu, seu.var, seu.lable, groups = "Sample", my_colours= "default"){
  sample.lables <- as.data.frame(table(seu.var, seu.lable))
  sample.lables$Freq <- as.double(sample.lables$Freq)
  check <- length(unique(sample.lables$seu.lable))
  print(paste("Number of colours needed",check, sep = ""))
  print(paste("Number of colours entered ", length(my_colours), sep = ""))
  if (length(my_colours) < check){
    # bar chart of with percent 
    print("Default ggplot colours used")
    print(ggplot(sample.lables, aes(x = seu.var,y=Freq ,fill = seu.lable)) + 
            geom_bar(position= "fill", stat = "identity") + 
            scale_y_continuous(labels = scales::percent_format()) + 
            theme_classic() + 
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))  
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }else {
    # bar chart of with percent 
    print("Custome colours used.")
    print(ggplot(sample.lables, aes(x = seu.var,y=Freq ,fill = seu.lable)) + 
            geom_bar(position= "fill", stat = "identity") + 
            scale_y_continuous(labels = scales::percent_format()) + 
            scale_fill_manual(values = my_colours)+
            theme_classic() + 
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))  
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }
  
}

# run plotting function on list of variables 

plotproportions <- function(seu, var.list, xgroup, varnames, my_colours = 'default'){
  for (i in seq_along(var.list)) {
    var <- var.list[[i]]
    proportionplots(seu = seu, seu.var = var, seu.lable = xgroup, groups = varnames[i], my_colours = my_colours)
  }
}

# var.list is a list defining seurat metadata slots
# xgroup is the cell types or x axis grouping variable 
# varnames is a character vector of labels for the x axis of the plots
# Example
# plotproportions(seu.q, var.list = var.list, xgroup = seu.q$cell.types, varnames = varnames)
# can input costume colours
# later --- add in define sample order 

### dotplots and heatmaps mean expression by group

plotmean <- function(plot_type = 'heatmap',seu, group, markers, var_names, slot = 'scale.data',
                     xlab = 'Cell Types', ylab = 'Antibodies', var1order, var2order){
  express.by.cluster <- as.data.frame(AverageExpression(seu, features = markers, 
                                                        group.by = group, slot = 'scale.data'))
  express.by.cluster <- as.data.frame(scale(express.by.cluster))
  if(length(var_names) < 0){
    col.names <- colnames(express.by.cluster) 
  }else{
    col.names <- var_names
  }
  names(express.by.cluster) <- col.names
  AB <- row.names(express.by.cluster)
  ex.data <- cbind(AB,express.by.cluster)
  # must make a data table because melt is 
  dt <- data.table(ex.data)
  dt.long<- melt(dt, id.vars = "AB")
  # select dot plot or heat map
  if(plot_type == 'heatmap'){
    # Heat map
    print(ggplot(dt.long, aes(x = variable, y = AB)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient(low="blue", high="pink", na.value = "grey") +
      labs(x=xlab, y= ylab) +
      theme_bw() + theme(axis.text.x=element_text(size=12, angle=90, hjust=0.95,vjust=0.2),
                         axis.text.y=element_text(size=12),
                         plot.title=element_text(size=12)))
  } else if(plot_type== 'dotplot'){
    a <- DotPlot(seu, features = AB, group.by = group)
    pct.exp <- as.data.frame(a$data) %>% select(features.plot, id, pct.exp)
    # add the mean expression and percent cells
    # rename columns to match
    colnames(dt.long) <- c("Markers","id","expression")
    colnames(pct.exp) <- c("Markers","id","proportion")
    df.exp.pct <- merge(dt.long, pct.exp, by = c("Markers", "id"))
    data <- df.exp.pct %>% mutate(Cell.type = factor(id, levels = var1order))
    data <- data %>% mutate(Marker = factor(Markers, levels = var2order))
    print(ggplot(data = data, aes(x = Marker, y = Cell.type, color = expression, size = proportion)) +
      geom_point() +
      scale_color_gradient(low = "grey", high = "blue") + ylab("Cell Phenotypes") 
      + xlab("Antibodies") + RotatedAxis())
  }
  else {
    print("not a valid plotting option")
  }
}
  
# compare groups functions

# prep_for_stats
# run_stats
# make plots



##############################################################################################
# prep_for_stats
# input seurat object
# arguments input all variable that might be compared and where to find these variables
# Selects expression data from Seurat object organized by designated variables

require(data.table)

Prep_for_stats <- function(seu, marker_list, variables, marker_name = 'Marker'){
  # create a dataframe with all the expresson data
  df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
  # rename columns 
  colnames(df) <- marker_list
  # add the different grouping
  for (i in seq_along(variables)) {
    meta_data <- as.data.frame(seu@meta.data[, variables[i]])
    colnames(meta_data) <- variables[i]
    df <- cbind(df, meta_data)
  }
  # need the antibodies to also be a variable
  # convert to data table to use the melt function in data.table instead of reshape2
  dt <- as.data.table(df)
  molten <- melt(dt, measure.vars = marker_list, variable.name = marker_name)
  return(molten)
}




##############################################################################################
# run_stats
# takes in dataframe from prep_for_stats
# arguments, variables, level1, level2
# runs 2way anova and posthoc tests
# outputs tables in a list
# anovas
# posthoc tukey's main 
# tukey's interaction effects filtered matching sets of interest 
# tables  filtered by significance


### for using sample means 

# group_cols is a vector with the columns to get the means from 

runstats <- function(input_df, group_cols = c("Sample", "CellType","Marker"),
                     stat_type = "ANOVA",dp1, dp2 = NULL, use.means = TRUE){
  if(use.means = TRUE){
    get_means <- function(df, group_cols, value_col) {
      df_means <- df %>%
        group_by(across(all_of(group_cols))) %>%
        mutate(mean_value = mean(value)) %>%
        distinct(across(all_of(group_cols)), mean_value, .keep_all = TRUE) %>%
        select(-value)
    } else {
      df_means <- input_df
    }
     
  }
}






##############################################################################################

# make plots
# input seurat object?
# makes box plots by two variables



