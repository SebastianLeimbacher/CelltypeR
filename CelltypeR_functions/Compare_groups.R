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

Prep_for_stats <- function(seu, marker_list, variables){
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
  molten <- melt(dt, measure.vars = marker_list, variable.name = 'Antibody')
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

##############################################################################################

# make plots
# input seurat object?
# makes box plots by two variables


