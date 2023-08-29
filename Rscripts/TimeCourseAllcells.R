# load libraries
require(Seurat)
require("ggplot2")
require("tidyverse")
require("CelltypeR")

# read in the seurat object 
seu <- readRDS("SeuTimeCourseAll.RDS")
seu <- get_clusters(seu, method = "louvain",
                    k = 80,
                    resolution = c(0.2,0.6,1),
                    plots = TRUE,
                    save_plots = output_path
)

saveRDS(seu, paste(output_path, "TimeCourseAllcells.RDS", sep = ""))

# Look at expression patterns

AB <- c("TH","CD24","CD56","CD29","CD15","CD184","CD133","SSEA4","CD44","CD49f","CD140a")
# the cluster labels will match the active ident

# save for each resolution UMAPs and a heatmap

Idents(seu) <- "RNA_snn_res.1"
pdf(paste(output_path,"AllCellsExpressionRes1.pdf",sep = ""))
# UMAPs
for (i in AB) {
  print(FeaturePlot(seu, features = i, min.cutoff = 'q1', max.cutoff = 'q97', label = TRUE))
}
# heatmap
clustnum <- length(levels(seu))-1
print(plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.1', markers = AB, 
         var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
         ylab = "Markers"))
dev.off()

# resolution 0.6
Idents(seu) <- "RNA_snn_res.0.6"
pdf(paste(output_path,"AllCellsExpressionRes06.pdf",sep = ""))
# UMAPs
for (i in AB) {
  print(FeaturePlot(seu, features = i, min.cutoff = 'q1', max.cutoff = 'q97', label = TRUE))
}
# heatmap
clustnum <- length(levels(seu))-1
print(plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.0.6', markers = AB, 
               var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
               ylab = "Markers"))
dev.off()

# resolution 0.2
Idents(seu) <- "RNA_snn_res.0.2"
pdf(paste(output_path,"AllCellsExpressionRes02.pdf",sep = ""))
# UMAPs
for (i in AB) {
  print(FeaturePlot(seu, features = i, min.cutoff = 'q1', max.cutoff = 'q97', label = TRUE))
}
# heatmap
clustnum <- length(levels(seu))-1
print(plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.0.2', markers = AB, 
               var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
               ylab = "Markers"))
dev.off()

### Calculate and visualize correlation labels

reference_data2 <- read.csv(reference_path2)
# the reference matrix need to be in the format cell types as rows and markers as columns
# there is a column X with the cell type names
df1 <- reference_data2
rownames(df1) <- df1$X # add row names (these are the markers)
df1 <- df1 %>% select("Astrocyte","Endothelial","Epithelial","Neuron",
                    "NeuronDA","NPC","NPCDA","opc","Oligos","RG",
                    "RadialGliaDiv","StemCell")

colnames(df1) <- c("Astrocyte","Endothelial","Epithelial","Neuron",
                  "NeuronDA","NPC","NPCDA","OPC","Oligo","Radial Glia",
                  "Radial Glia Dividing","StemCell")
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/TimeCourseAIW/Analysis/"
input_df <- read.csv(paste(output_path,"FlowDFretroAllcellsTimecourse.csv"))

# The function will only work when the markers in the test data is are the same as in the reference data
# Both dataframes need to have the markers/antibodies as columns. 
# Calculate the correlations 
cor <- find_correlation(test = input_df, 
                        reference = df1, 
                        min_corr = 0.35, 
                        min_diff = 0.05)

write.csv(cor, paste(output_path, "AllcellsCor_thesh35.csv"),sep = "")
# creates a dataframe with cor1 cor2 and predicted cell type label

cor2 <- find_correlation(test = input_df, 
                         reference = df1, 
                         min_corr = 0.55, 
                         min_diff = 0.05)
write.csv(cor2, paste(output_path, "AllcellsCor_thesh55.csv"),sep = "")

cor3 <- find_correlation(test = input_df, 
                         reference = df1, 
                         min_corr = 0.05, 
                         min_diff = 0.05)
write.csv(cor3, paste(output_path, "AllcellsCor_thesh05.csv"),sep = "")

#### visualize the correlation results
pdf(paste(output_path, "CorrelationPlots.pdf",sep = ""))
print(plot_corr(cor, threshold = 0.35, min_cells = 400))
print(plot_corr(cor2, threshold = 0.55, min_cells = 400))
print(plot_corr(cor3, threshold = 0.05, min_cells = 400))
dev.off()

#### add the correlation results into the seurat object
# add the correlation predictions to the meta data
seu <- AddMetaData(object=seu, metadata=cor$cell.label, col.name = 'cor.labels3505')

# add the correlation predictions to the meta data
seu <- AddMetaData(object=seu, metadata=cor2$cell.label, col.name = 'cor.labels5505')

# add the correlation predictions to the meta data
seu <- AddMetaData(object=seu, metadata=cor3$cell.label, col.name = 'cor.labels0505')

### visualize the correlation labels 

pdf(paste(output_path, "CorrelationPlotsUMAPs.pdf",sep = ""), width = 20, height = 6)
print(DimPlot(seu, group.by = 'cor.labels3505', label = TRUE))
print(DimPlot(seu, group.by = 'cor.labels5505', label = TRUE))
print(DimPlot(seu, group.by = 'cor.labels0505', label = TRUE))
dev.off()

### add the cor labels per cluster - will just label the highest cluster number
cor.ann <- get_annotation(seu, seu.cluster = seu$RNA_snn_res.1, 
                          seu.label = seu$cor.labels3505, top_n = 3, 
                          filter_out = c("Unknown","unknown","Mixed","unassigned","Unassigned"), 
                          Label = "CAM3505")
seu <- annotate(seu,annotations = cor.ann$CAM3505, 
                to_label = "RNA_snn_res.1", 
                annotation_name = "CAM35")


cor.ann <- get_annotation(seu, seu.cluster = seu$RNA_snn_res.1, 
                          seu.label = seu$cor.labels5505, top_n = 3, 
                          filter_out = c("Unknown","unknown","Mixed","unassigned","Unassigned"), 
                          Label = "CAM5505")
seu <- annotate(seu,annotations = cor.ann$CAM5505, 
                to_label = "RNA_snn_res.1", 
                annotation_name = "CAM55")

cor.ann <- get_annotation(seu, seu.cluster = seu$RNA_snn_res.1, 
                          seu.label = seu$cor.labels0505, top_n = 3, 
                          filter_out = c("Unknown","unknown","Mixed","unassigned","Unassigned"), 
                          Label = "CAM0505")
seu <- annotate(seu,annotations = cor.ann$CAM0505, 
                to_label = "RNA_snn_res.1", 
                annotation_name = "CAM05")


# save the Correlation annotations
saveRDS(seu, paste(output_path, "TimeCourseAllcells.RDS", sep = ""))

########## Predict with Random Forest Classifier already trained



######## Predict with Seurat label transfer



