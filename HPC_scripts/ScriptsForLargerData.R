# Time course analysis

# load necessary libraries 
library(Seurat)
library(dplyr) 
library(ggplot2)
library(CelltypeR)
library(flowCore)

# # cluster large object

output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/TimeCourseAIW/Analysis/"
# read in the object
seu <- readRDS(paste(output_path,"SeuTimeCourseAll.RDS",sep=""))

seu <- get_clusters(seu, method = "louvain",
                    df_input = df.input,
                    k = 80,
                    resolution = c(0.5,0.8,1,1.5),
                    pcdim = 1:10,
                    plots = TRUE,
                    save_plots = FALSE)

# save with the clusters
saveRDS(seu, paste(output_path,"SeuTimeCourseAll.RDS",sep=""))


