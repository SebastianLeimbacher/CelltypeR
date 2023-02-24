# compare between groups

# plotting function

proportionplots <- function(seu, seu.var, seu.lable, groups = "Sample", my_colours){
  sample.lables <- as.data.frame(table(seu.var, seu.lable))
  sample.lables$Freq <- as.double(sample.lables$Freq)
  if (length(my_colours)== 0){
    # bar chart of with percent 
    print(ggplot(sample.lables, aes(x = seu.var,y=Freq ,fill = seu.lable)) + 
            geom_bar(position= "fill", stat = "identity") + 
            scale_y_continuous(labels = scales::percent_format()) + 
            theme_classic() + 
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))  
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }else {
    # bar chart of with percent 
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
