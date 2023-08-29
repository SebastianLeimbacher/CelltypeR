# training Random Forest Classifiers

# Using too few features (mytry), too low of starting nodes (tree complexity), too few nodes (max nodes) and too few trees can all lead to over fitting.  The num_folds variable provides a cross validation within the parameter search to help avoid overfitting. 
# mytry is for the number of features to be tested as input
# here we have a max of 13 putting all the features doesn't work to create a RFM because all the trees will be the same.
# Fewer nodes and trees can lead to a more accurate model in the training/test but an overfit model when applied to the rest of the dataset
# the cross validation feature will also help avoid overfitting as conditions that work in the various iterations will be selected.

# Whatever categories/cell types we define are what the model can predict.  No intermidiate cell types are predicted.
# RFM test does not output the second most probably annotation for each cell.

# To train a Random Forest model from the 9 MBO sample and test that on the new time course data with the new panel. We need to only select the Markers/Features that overlap.
panel1 <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
panel2 <- c("TH","CD24","CD56","CD29","CD15","CD184","CD133","SSEA4","CD44","CD49f","CD140a")
overlap <- intersect(panel1,panel2)
print(overlap)

# now there are only 8 Markers (features)

# Labelled data object
require(Seurat)
require(CelltypeR)
seu.r <- readRDS("~/Documents/Data/FlowCytometry/PhenoID/Analysis/NatMethodJuneSubmission/Figure4/All9MOannaoteAug.RDS")
# check metadata to see what the cell type annotations to use
print(colnames(seu.r@meta.data))

# Annotations in Figure 4 is under Celltypes
print(table(seu.r$Celltypes))
# some cell types like Stem cell like and OPC have very few cells and will not be predicted well.
# Down sample with the seurat function to improve cell type balance


Idents(seu.r) <- "Celltypes"
seu.sub <- subset(seu.r, downsample = 3000)
rf <- RFM_train(seurate_object = seu.sub,
                markers = overlap, 
                annotations = seu.sub$Celltypes,
                num_folds = 3,
                downsample = 8000,
                seed = 222,
                mytry = c(4:6),
                maxnodes = c(16:18),
                trees = c(1800),
                cores_to_use = 4)

# kappa on the training data is 0.83 but only 0.69 on the test data.  
# "Best parameter mytry 4 and max node 17 and tree number 1800 ."

# try some other conditions and include more cells
rf2 <- RFM_train(seurate_object = seu.sub,
                markers = overlap, 
                annotations = seu.sub$Celltypes,
                num_folds = 3,
                downsample = 20000,
                seed = 222,
                mytry = c(4:5),
                maxnodes = c(17:18),
                trees = c(1800,2000),
                cores_to_use = 4)

#"Best parameter mytry 4 and max node 18 and tree number 2000 ."
# kappa on training is 0.84 and test is 0.69
setwd("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/NatMethodJuneSubmission/RFM/")

saveRDS(rf2, "RFM_all9hMOsOverlapAB.RDS")


# train a model with all full panel and final cell types for annotating a new dataset.

# try some other conditions and include more cells
rf <- RFM_train(seurate_object = seu.sub,
                 markers = panel1, 
                 annotations = seu.sub$Celltypes,
                 num_folds = 3,
                 downsample = 30000,
                 seed = 222,
                 mytry = c(4:7),
                 maxnodes = c(18:20),
                 trees = c(1800,2000),
                 cores_to_use = 4)

saveRDS(rf, "RFM_all9hMOAB.RDS")
# training kappa 0.959 and test kappa 0.85

# train a model with the new antibody panel

class(rf)

unique(seu.r$Celltypes)

####### train the model from the subset of the all time course
seu.t <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/TimeCourseAIW/Analysis/SeuratSubTimeline.RDS")
colnames(seu.t@meta.data)
unique(seu.t$Celltypes2)

