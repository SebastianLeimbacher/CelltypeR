library(dplyr) 
library(kit) 
library(reshape2) 
library(Seurat)

# Annotate clusters functions:

# see_features                (Shuming/Rhalena *** new)
# CAM                         (Shuming)
# RFM_train                   (Rhalena)
# RFM_predict                 (Rhalena)
# seurat_transfer             (Rhalena)
# cluster_annotate            (Rhalena)



##############################################################################################

# see_features
# takes in a seurat object with clustering and UMAP run in explore_param or otherwise
# creates featuremaps for each antibody
# creates a UMAP with cluster numbers labelled
# creates a heatmap split by cluster - takes an argument for which cluster resolution is desired
# takes in a features list








##############################################################################################

# CAM
# input correlation matrix
# input seurat object
#Correlation assignment method, 
# predicts cell types based on each cells correlation to the matrix types
# creates plots and tables of the prediction outputs. 
# takes arguement for "unknown" threshold (default 0.4)  and "double-label" thresholod(0.05)


find_correlation <- function(test, 
                             reference, 
                             min_corr = 0.1, 
                             min_diff = 0.05) {
  
  #1. process and scale test and reference matrix
  
  #find intersect between reference and test, ignore case, 
  #change reference spelling, order to test's
  testc <- vector()
  refc <- vector()
  for (i in colnames(test)) {
    for (j in colnames(reference)) {
      if (tolower(i) == tolower(j)) {
        testc <- c(testc, i)
        refc <- c(refc, j)
      }
    }
  } 
  
  reference <- reference[, refc] #select markers + X in reference
  colnames(reference) <- testc #change reference's spelling to test's
  test <- test[, testc] #select markers + X in test
  markers <- colnames(select_if(test, is.numeric)) #a list of markers (without X)
  test[, markers] <- scale(test[, markers]) #z score test (without X) 
  reference[, markers] <- scale(reference[, markers])  # z score the reference matrix 
  
  
  #2. find best and second correlation and cell type
  result <- vector()

  #the loop that will find the best and second best correlation and cell types
  for (i in 1:nrow(test)) {
    corr_ls <- vector() #list of correlation between the reference cell types and each sample
    for (j in 1:nrow(reference)) {
      corr <- cor(as.numeric(test[i,markers]), as.numeric(reference[j,markers])) # pearson by default and we use default
      corr_ls <- c(corr_ls, corr)
    }
    
    top <- topn(corr_ls, 2L) #return the index of the best 2
    result <- c(result, 
                test[i, 'X'], #col 1: cell sample
                corr_ls[top[1]], #col 2: 1st correlation
                reference[top[1], 'X'], #col 3: 1st best cell type
                corr_ls[top[2]], #col 4: 2nd correlation 
                reference[top[2], 'X'], #col 5: 2nd best cell type
                ifelse(corr_ls[top[1]] < min_corr, "unknown",
                       ifelse(corr_ls[top[1]] - corr_ls[top[2]] < min_diff, 
                              paste(reference[top[1], 1], 
                                    reference[top[2], 1], sep = "-"), 
                              reference[top[1], 1]))) #col 6: assigned cell type
    # if best corr < min_corr, assign unknown cell type
    # if best corr - second best corr < min diff, assign combined cell type 
    # else, assign best cell type
  }
  #convert the result list to a df  
  cdf <- data.frame(matrix(result, ncol=6, byrow = TRUE))
  colnames(cdf) <- c("X", "cor.1", "best.cell.type",
                     "cor.2", "second.cell.type", "cell.label")
  cdf$cor.1 <- as.numeric(cdf$cor.1)
  cdf$cor.2 <- as.numeric(cdf$cor.2)
  return(cdf)
}


plot_corr <- function(df) {
  # filter to get frequency table and save as csv
  df.f <- df %>% select(cell.label)
  freq.table <- as.data.frame(table(df.f))
  df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> 100)
  
  # plot the frequencies
  plot1 <- 
    ggplot(df.filter, aes(x = reorder(
      cell.label, cell.label, function(x) - length(x)),
      fill = cell.label)) + geom_bar() + theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) + 
    xlab('Assigned cell type') + 
    ylab('number of cell') + 
    labs(fill='Cell Types')
  # print(plot1)
  
  # violin plot of best correlation/cell type
  plot2 <- 
    ggplot(df, aes(x = best.cell.type, y = cor.1)) +
    geom_violin()+ 
    ylim(-0.1, 1)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("correlation coefficient") + 
    xlab("Cell type with max correlation coefficient")
  # print(plot2)
  
  df.melt <- melt(df) #reformat to long df
  
  # plot the best and second best correlation together
  plot3 <- 
    ggplot(df.melt, aes(x = cell.label, y = value ))+ 
    geom_boxplot()+ ylim(-0.1, 1)+theme_classic()+
    theme(axis.text.x = element_text(angle = 90))+ 
    ylab("correlation coefficient") + 
    xlab("Cell type label")
  # print(plot3)
  
  # plot the best and second best correlation separated on the same graph
  plot4 <- 
    ggplot(df.melt, aes(x = best.cell.type, y = value, fill = variable))+ 
    geom_boxplot()+ 
    ylim(-0.25, 1)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_fill_manual(values = c("#4E84C4", "#52854C")) + 
    ylab("correlation coefficient") + xlab("Cell type")
  # the second best correlation is so low it was removed from STEM with axis limit -0.1 and even -1
  
  # print(plot4)
  
  # down sample
  set.seed(64)
  df.downsample <- sample_n(df, 1000)
  df.melt.down <- melt(df.downsample)
  
  # # reformat the table to work with the before after plot
  # # y is the measurement in df.melt = value
  # # x is before after in df.melt = variable
  # # class another variable - in the example this is different shapes - for us this is best cell type
  # # might use facet to split the cell type - needs to be a factor
  # # id is the individual id this is the X column
  plot5 <- 
    ggplot(df.melt.down, aes(x = variable, y = value,colour = variable, group = X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + 
    geom_point()+ 
    scale_color_manual(values = c("#4E84C4", "#52854C")) + 
    ylim(-0.25, 0.95) +
    facet_wrap(~(as.factor(best.cell.type))) +
    theme(legend.position = "none") +
    ylab("Correlation Coefficient") +
    xlab("")
  # print(plot5)
  
  double.cells <- df[grep("-", df$cell.label),]
  df.melt.double <- melt(double.cells)
  
  # # this will be an excellent visualization but I need to subset only the double labels, then I can plot more cells and see more clearly.
  plot6 <- 
    ggplot(df.melt.double, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + 
    geom_point()+ 
    scale_color_manual(values = c("#4E84C4", "#52854C")) + 
    ylim(-0.15,0.8) +
    facet_wrap(~(as.factor(cell.label))) +
    ylab("Correlation Coefficient") +
    xlab("")
  # print(plot6)
  
  return(list(freq.table, plot1, plot2, plot3, plot4, plot5, plot6))
}


##############################################################################################
# RFM_train
# Input annotated FC dataset 
# Random Forest Model internally optimizing parameters and saving the best model. 
# requires randomForst, caret, data.table

# annotation needs to be the seurat object and meta data slot
require(data.table)
require(randomForest)
require(caret)

RFM_train <- function(seurate_object, 
                             AB_list, annotations = seu$CellTypes,
                      split = c(0.5,0.5),
                      downsample = 'none',
                      seed = 222,
                      mytry = c(1:10),
                      maxnodes = c(12: 25),
                      trees = c(250, 500, 1000,2000),
                      start_node = 15){
  # set up the data
  seu <- seurate_object
  AB <- AB_list
  df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
  nodes <- maxnodes # for renaming later
  # add antibody/marker names
  colnames(df) <- AB
  # add the annotations
  ann <- annotations
  df.l <- cbind(df, lables = ann)
  # downsample option
  if(downsample == 'none'){
    df1 <- df.l
  } else {
    # down sample - randomly selecting X rows from the dataframe
    df1 <- sample_n(df.l, downsample)
  }
  # split in to training and test
  set.seed(seed)
  ind <- sample(2, nrow(df1), replace = TRUE, prob = split) # prop is the proportions
  # split the data into train and test
  train <- df1[ind==1,]
  test <- df1[ind==2,]
  # need to start with a basic model to optimize
  trControl <- trainControl(method = "cv", number = 10, search ="grid")
  
  # test different parameters
  # mtry: number of features to draw - default is the squareroot of the number of columns
  mtry.range <- mytry
  print("optimize number of features to draw")
  tuneGrid <- expand.grid(.mtry = mtry.range)
  rf_mtry <- train(lables~., 
                   data=train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE,
                   nodesize = start_node,
                   ntree = 300)
  print(paste("Best number of features to draw is ",rf_mtry$bestTune$mtry, sep = ""))
  print(paste("Max accuracy is ", max(rf_mtry$results$Accuracy), sep=""))
  best_mtry <- rf_mtry$bestTune$mtry 
  
  # best node size
  print("Searching for best node size")
  store_maxnode <- list()
  tuneGrid <- expand.grid(.mtry = best_mtry)
  for (maxnodes in maxnodes) {
    rf_maxnode <- train(lables~.,
                        data = train,
                        method = "rf",
                        metric = "Accuracy",
                        tuneGrid = tuneGrid,
                        trControl = trControl,
                        importance = TRUE,
                        nodesize = start_node,
                        maxnodes = maxnodes,
                        ntree = 300)
    current_iteration <- toString(maxnodes)
    store_maxnode[[current_iteration]] <- rf_maxnode
  }
  results_node <- resamples(store_maxnode)
  results_sum <- summary(results_node)
  # summarize results makes another list object with accuracy
  results_sum$statistics$Accuracy
  # find the max accuracy and the max nodes with the max accuracy
  sum.df <- data.frame(results_sum$values)
  max_accuracy <- max(sum.df.2)
  #print(paste("Max accuracy is ",max_accuracy))
  sum.df.2 <- select(sum.df, -contains("Kappa"))
  colnames(sum.df.2) <- nodes
  for(i in (1:10)){
    #print(colnames(sum.df.2)[which.max(sum.df.2[i,])])
    #print(max(sum.df.2[i,]))
    max.temp <- max(sum.df.2[i,])
    if(max.temp == max_accuracy){
      max_node <-colnames(sum.df.2)[which.max(sum.df.2[i,])]
      #print(paste("Max node is ", max_node, sep=""))
    } 
  }
  print(paste("The max accuracy is ",max_accuracy, ". For the maxnodes of ",
              max_node, sep = ""))
  
  # now search of the best number of trees
  print("searching for best number of trees")
  store_maxtrees <- list()
  for (ntree in trees) {
    rf_maxtrees <- train(lables~.,
                         data = train,
                         method = "rf",
                         metric = "Accuracy",
                         tuneGrid = tuneGrid,
                         trControl = trControl,
                         importance = TRUE,
                         nodesize = 25,
                         maxnodes = as.numeric(max_node),
                         ntree = ntree)
    key <- toString(ntree)
    store_maxtrees[[key]] <- rf_maxtrees
  }
  # get the best number of trees
  results_tree <- resamples(store_maxtrees)
  results_sum <- (summary(results_tree))
  # make a data frame with the summary results
  sum.df <- data.frame(results_sum$values)
  #print(paste("Max accuracy is ",max_accuracy))
  sum.df.2 <- select(sum.df, -contains("Kappa"))
  colnames(sum.df.2) <- trees
  max_accuracy <- max(sum.df.2)
  for(i in (1:10)){
    #print(colnames(sum.df.2)[which.max(sum.df.2[i,])])
    #print(max(sum.df.2[i,]))
    max.temp <- max(sum.df.2[i,])
    if(max.temp == max_accuracy){
      best_tree <-colnames(sum.df.2)[which.max(sum.df.2[i,])]
      print(paste("Best trees  ", best_tree, sep=""))
    } 
  }
  print(paste("The max accuracy is ",max_accuracy, ". For ",
              best_tree," trees", sep = ""))
  
  
  # now run model with the selected conditions
  rf <- randomForest(lables~.,
                     train,
                     mtry = best_mtry,
                     nodesize = start_node,
                     ntree = as.numeric(best_tree),
                     maxnodes = as.numeric(max_node))  
  
  # check the results of the true model
  prediction.train <-predict(rf, train)
  prediction.test <-predict(rf, test)
  
  print("predict training data")
  print(confusionMatrix(prediction.train, train$lables))
  print("predict test data")
  print(confusionMatrix(prediction.test, test$lables))
  
  # return the model - can save after  
  return(rf)
}
  





##############################################################################################
# RFM_predict
# take the saved trained RFM
# takes seurat data object
# predicts cell types
# creates table and umap
# outputs prediction vector


RFM_predict <- function(seu, rf){
  # prepare a data object to be test data
  df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
  colnames(df) <- AB
  # run the predictions
  rfm.pred <- as.data.frame(predict(rf,df))
  
}
  


##############################################################################################
# seurat_transfer
# takes in a reference seurat object
# follows seurat workflow, finds anchors, predicts labels
# default of no threshold, a threshold can be applied
# outputs prediction vector

seurat_predict <- function(seu.q, seu.r, ref_id = 'Labels', down.sample = 500, markers){
  Idents(seu.r) <- ref_id 
  seu.r <- subset(seu.r, downsample = down.sample)
  # find anchors
  anchors <- FindTransferAnchors(reference = seu.r, 
                                 query = seu.q, features = markers,
                                 reference.reduction = "pca", 
                                 dim= 1:10) 
  predictions <- TransferData(anchorset = anchors, refdata = seu.r$subgroups, dims = 1:10)
  seu.q <- AddMetaData(seu.q, metadata = predictions$predicted.id, col.name = 'seu.pred')
  
}





##############################################################################################

# plot output of predictions

# takes in a seurat object with the labels added 
# makes a dataframe with the count of predicted labels for each cluster
# input seurat object with the predicted labels in the meta data
# input the clusters meta data slot to be labels
# input the meta data slot with the labels (correlation, random forest, seurat predicted)


plot_lab_clust <- function(seu, seu.cluster, seu.labels){
  t.lables <- as.data.frame(table(seu.clusters, seu.labels))
  t.lables$Freq <- as.double(t.lables$Freq)
  colnames(t.lables) <- c("Cluster","Lable","Freq")
  print(ggplot(t.lables, aes(y = Freq, x = Cluster, fill = Lable)) + geom_bar(position = "stack", stat= "identity"))
  t.lab.known <- t.lables %>% filter(!Lable == "unknown")
  print(ggplot(t.lab.known, aes(y = Freq, x = Cluster, fill = Lable)) + geom_bar(position = "stack", stat= "identity"))
}





## examples usage

# input arguments 
#seu.clusters <- seu.t$RNA_snn_res.0.8   # cluster resolution to label
#seu.lable <- seu.t$cor.labels   # predicted labels per cell meta data slot

#test.cltable <- plot_lab_clust(seu, seu.clusters, seu.lable)

#plot_lab_clust(seu, seu.clusters, seu.lable)

##############################################################################################
# get annotation predicted for each cluster
# a function to get a dataframe for annotation

# takes in a seurat object with clusters and predicted labels
# takes the cluster meta data slot, the predicted labels slot
# take number of different top predictions per cluster to view

# later add in option if > X% of cells are unknown keep unknown
# will return a dataframe with label for each cluster number

get_annotation <- function(seu, seu.cluster, seu.label, top_n = 3, 
                           ignore_unknown = FALSE, Label = "Label"){
  t.lables <- as.data.frame(table(seu.cluster, seu.label))
  t.lables$Freq <- as.double(t.lables$Freq)
  colnames(t.lables) <- c("Cluster", Label,"Freq")
  top.labs <- t.lables  %>% group_by(Cluster)  %>% top_n(top_n, Freq)
  sort.tops <- top.labs %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
  print(sort.tops)
  # now filter out the unknown if desired to get a label for each cluster
  if(ignore_unknown == TRUE){
    t.lab.known <- t.lables %>% filter(!Label == "unknown")
    top.lab <- t.lab.known  %>% group_by(Cluster)  %>% top_n(1, Freq)
    sort.tops <- top.lab %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
    #print(sort.tops)
  }
  else{
    top.lab <- t.lables  %>% group_by(Cluster)  %>% top_n(1, Freq)
    sort.tops <- top.lab %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
    #print(sort.tops)
  }
  vec <- sort.tops[,Label]
  vec <- paste(vec, collapse = ",")
  print("Annotations in order of clusters starting at 0")
  return(sort.tops %>% select(-"Freq"))
}

## example
# input arguments 
#seu.clusters <- seu.t$RNA_snn_res.0.8   # cluster resolution to label
#seu.lable <- seu.t$cor.labels   # predicted labels per cell meta data slot
# top number of cell type labels. This will only be for the printed table
# the top most frequently label will be used for cluster annotation
# ignore unknow is an option to filter out unknown cells and then the next more frequently
# predicted cell type will be used. 

#ann_cor <- get_annotation(seu.t, seu.clusters, seu.lable, top_n = 3, ignore_unknown = FALSE)



##############################################################################################

# cluster_annotate
# input dfs of cluster annotation, at least 1 is required
# dfs were created with get_annotation function 
# manual annotation can be df as well.  Make in excel read in csv or create 
# type manually, be careful to match syntax to the other predictions
# input seurat data object
# finds most common label per cluster and applies that label


cluster_annotate <- function(seu, ann.list, 
                             annotation_name, 
                             to_label){
  # not all annotations methods will always present 
  # easiest to make the input a list
  df.merge <- Reduce(function(x,y) merge(x,y, by="Cluster", all = TRUE),
                     ann.list)
  # convert all to lowercase 
  df2 <- lapply(df.merge, function(x) {tolower(as.character(x))})
  # back into a dataframe
  df3 <- as.data.frame(do.call(cbind, df2))
  
  df3$consensus <- apply(df3, 1, function(row) {
    if (is.data.frame(row)) {
      word_counts <- table(unlist(row[,2:ncol(row)]))
    } else if (is.matrix(row)) {
      word_counts <- table(unlist(row))
    } else {
      word_counts <- table(unlist(c(row)))
    }
    names(word_counts)[which.max(word_counts)]
  })
  
  # now reorder by cluster number
  df3$Cluster <- as.integer(as.character(df3$Cluster))
  df.sort <- df3  %>% arrange(Cluster) 
  # get a vector 
  print(df.sort %>% select("Cluster","consensus"))
  
  # now add the annotations into the seurat
  # use the simple annnotation function
  seu <- annotate(seu, annotations = df.sort$consensus, to_label,
           annotation_name)
}




annotate <- function(seu, annotations, to_label, annotation_name = "CellType"){
  Idents(seu) <- to_label
  names(annotations) <- levels(seu)
  seu <- RenameIdents(seu, annotations)
  seu <- AddMetaData(object=seu, metadata=Idents(seu), col.name = annotation_name)
  
}


