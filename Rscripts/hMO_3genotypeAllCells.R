# this script has the processing for the full data set of 9 hMOs from 3 iPSC lines
# Shown in Figure 4

# load libraries
require("Seurat")
require("ggplot2")
require("tidyverse")
require("CelltypeR")

setwd("~/Documents/Data/FlowCytometry/PhenoID/Analysis/NatMethodJuneSubmission/")
getwd()
# Read in and align all the data


AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")

# Cluster
seu <- get_clusters(seu, method = "louvain",
                    k = 80,
                    resolution = c(0.25,0.5,1),
                    plots = TRUE,
                    save_plots = output_path
)


# annotate clusters

# Look at expression patterns

# the cluster labels will match the active ident

# save for each resolution UMAPs and a heatmap

# resolution 0.5
Idents(seu) <- "RNA_snn_res.0.5"
pdf("All9MBOExpressionRes05.pdf")
# UMAPs
for (i in AB) {
  print(FeaturePlot(seu, features = i, min.cutoff = 'q1', max.cutoff = 'q97', 
                    label = TRUE, raster = FALSE))
}
# heatmap
clustnum <- length(levels(seu))-1
#pdf("HeatmapRes05expression.pdf")
print(plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.0.5', markers = AB, 
               var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
               ylab = "Markers"))
dev.off()


# 
# resolution 0.25
Idents(seu) <- "RNA_snn_res.0.25"
pdf("ExpressionRes025.pdf")
# UMAPs
for (i in AB) {
  print(FeaturePlot(seu, features = i, min.cutoff = 'q1', max.cutoff = 'q97', label = TRUE, raster = FALSE))
}
# heatmap
clustnum <- length(levels(seu))-1
print(plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.0.25', markers = AB, 
               var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
               ylab = "Markers"))
dev.off()


# save plots for cluster resolution 1
pdf("All9MBO_UMAPFeatures.pdf", width = 24,height = 20)
FeaturePlot(seu, features = AB, min.cutoff = 'q1', 
            max.cutoff = 'q97', label = TRUE, raster = FALSE,
            cols = 4)
dev.off()

# Clusters UMAPs
pdf("All9MBO_UMAPclustersRes1.pdf", width = 6,height = 4)
DimPlot(seu.q, raster = FALSE, group.by = "RNA_snn_res.1", label = TRUE)
dev.off()

# mean expression heatmap by cluster
clustnum <- length(levels(seu))-1
print(clustnum)
pdf("All9MBO_HMclustersRes1.pdf", width = 6,height = 4)
plotmean(plot_type = 'heatmap',seu = seu, group = 'RNA_snn_res.1', markers = AB, 
         var_names = 0:clustnum, slot = 'scale.data', xlab = "Cluster",
         ylab = "Markers")
dev.off()


### Calculate and visualize correlation labels

input_df <- read.csv(test_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs/retrotransformed_flowset.csv")

reference_data <- read.csv(reference_path <- "/Users/rhalenathomas/GITHUB/CelltypeR/ExampleOuts/FinalReferenceMatrix.csv")
# the reference matrix need to be in the format cell types as rows and markers as columns
# there is a column X with the cell type names
df1 <- reference_data
rownames(df1) <- df1$X # add row names (these are the markers)
df1 <- df1 %>% select(-"X") # remove the column with the marker names
colnames(df1) <- c("Astrocytes","Endothelial","Epithelial","Neurons",
                   "NPC","OPC","Oligo","RadialGlia","StemCell")
df.ref <- as.data.frame(t(df1))
df.ref$X <- rownames(df.ref)

# the reference data frame need columns as markers and rows as cell types.
# The cell type names column X is needed in the reference data.  The index column X is needed in the reference data
# extra columns are fine

cor <- find_correlation(test = input_df, 
                        reference = df.ref, 
                        min_corr = 0.35, 
                        min_diff = 0.005)

write.csv(cor,"AllcellsCor_thesh35005.csv")
# creates a dataframe with cor1 cor2 and predicted cell type label

cor2 <- find_correlation(test = input_df, 
                         reference = df.ref, 
                         min_corr = 0.55, 
                         min_diff = 0.005)
write.csv(cor2,"AllcellsCor_thesh55.csv")

cor3 <- find_correlation(test = input_df, 
                         reference = df1, 
                         min_corr = 0.05, 
                         min_diff = 0.005)
write.csv(cor3,"AllcellsCor_thesh05.csv")

#### visualize the correlation results
pdf("AllCellsCorrelationPlots.pdf")
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



