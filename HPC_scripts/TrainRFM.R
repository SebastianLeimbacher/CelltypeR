R_LIB_PATH='/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2'
.libPaths(R_LIB_PATH)

# load libraries
require(Seurat)
require("ggplot2")
require("tidyverse")
require("randomForest")
require(caret)
require(doParallel)
require(foreach)


# run function - if CelltypeR library is installed this is not needed
RFM_train <- function(seurate_object,
                          markers, annotations = seu$CellTypes,
                          split = c(0.5, 0.5),
                          downsample = 'none',
                          seed = 222,
                          mytry = c(6:12),
                          maxnodes = c(12:25),
                          trees = c(250, 500, 1000, 2000),
                          start_node = 12,
                          cores_to_use = 4) {

  seu <- seurate_object
  all_features <- seu@assays[["RNA"]]@counts@Dimnames[[1]]
  df <- as.data.frame(t(as.data.frame(GetAssayData(seu, slot = 'scale.data'))))
  colnames(df) <- all_features
  df <- df[, markers]
  ann <- annotations
  df.l <- cbind(df, labels = ann)
  if (downsample == 'none') {
    df1 <- as.data.frame(df.l)
  } else {
    df1 <- as.data.frame(sample_n(df.l, downsample))
  }
  set.seed(seed)
  ind <- sample(2, nrow(df1), replace = TRUE, prob = split)
  train <- na.omit(df1[ind == 1,])
  test <- na.omit(df1[ind == 2,])
  trControl <- trainControl(method = "cv", number = 10, search = "grid")

  # Searching for the best mytry and maxnode
  registerDoParallel(cores = cores_to_use)
  results_list <- list()
  for (mytry_val in mytry) {
    for (maxnode_val in maxnodes) {
      for (tree_val in trees) {
        tryCatch({
          print(paste("Testing mytry:", mytry_val, ", maxnode:", maxnode_val, ", tree number:", tree_val))
          tuneGrid <- expand.grid(.mtry = mytry_val)
          rf_result <- train(labels ~ .,
                             data = train,
                             method = "rf",
                             metric = "Kappa",
                             tuneGrid = tuneGrid,
                             trControl = trControl,
                             importance = TRUE,
                             nodesize = maxnode_val,
                             ntree = tree_val)
          results_list[[paste("mytry_", mytry_val, "_maxnode_", maxnode_val, "_tree_", tree_val)]] <- rf_result
        }, error = function(e) {
          print(paste("Error for mytry:", mytry_val, ", maxnode:", maxnode_val, ", tree number:", tree_val))
          print(e)
        })
      }
    }
  }

  stopImplicitCluster()
  print("Finished testing mytry, maxnode, and tree number.")

  # Find the best mytry, maxnode, and tree number
  best_combination <- which.max(sapply(results_list, function(res) max(res$results$Kappa)))
  best_param <- names(results_list)[best_combination]
  print(paste("Best parameter combination:", best_param))

  # Extract the best mytry, maxnode, and tree number values
  param_values <- unlist(strsplit(gsub(".*mytry_ (\\d+) _maxnode_ (\\d+) _tree_ (\\d+)", "\\1 \\2 \\3", best_param), " "))
  best_mytry <- as.numeric(param_values[1])
  best_maxnode <- as.numeric(param_values[2])
  best_tree <- as.numeric(param_values[3])

  print(paste("Best parameter mytry", best_mytry, "and max node", best_maxnode, "and tree number", best_tree, "."))


  final_rf <- randomForest(labels ~ .,
                           data = train,
                           mtry = best_mytry,
                           nodesize = best_maxnode,
                           ntree = best_tree)

  prediction.train <- predict(final_rf, train)
  prediction.test <- predict(final_rf, test)

  print("Predict training data")
  print(confusionMatrix(prediction.train, train$labels))
  print("Predict test data")
  print(confusionMatrix(prediction.test, test$labels))

  return(final_rf)
}

# load in the reference object
print("Reading in reference seurat object labelled in June 2023 from the 9000 subset of cells")
seu.r <- readRDS("/lustre03/project/6070393/rhalena/Test_celltypeR/Input_Data/references/Seu9000lablesJune23.RDS")

print(dim(seu.r))
# call function to train model
print("Calling RFM training function.  Down sample to 36000 cells.")
markers <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
# mytry is for the number of features to be tested as input
# here we have a max of 13
rf <- RFM_train(seurate_object = seu.r, 
                    annotations = seu.r$labels,
                    markers = markers,
                      split = c(0.8,0.2),
                      downsample = 36000,
                      seed = 222,
                      mytry = c(3:7),
                      maxnodes = c(12:18),
                      trees = c(300,800,1200),
                      start_node = 10)


saveRDS(rf,"/lustre03/project/6070393/rhalena/Test_celltypeR/Input_Data/references/RFMtrainedFromJune2023lablesAug8.RDS")
