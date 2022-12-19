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




# # prepare data
# seu.r<- readRDS(pathway_to_data)
# AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
# 
# 
# df <- transpose(as.data.frame(GetAssayData(seu.r,slot = 'scale.data')))
# dim(df)
# 
# 
# 
# mtry.range <- c(1:10)
# 
# set.seed(48)
# 
# tuneGrid <- expand.grid(.mtry = mtry.range)
# rf_mtry <- train(lables~., 
#                  data=train,
#                  method = "rf",
#                  metric = "Accuracy",
#                  tuneGrid = tuneGrid,
#                  trControl = trControl,
#                  importance = TRUE,
#                  nodesize = 15,
#                  ntree = 300)
# rf_mtry$bestTune$mtry
# max(rf_mtry$results$Accuracy)
# 
# best_mtry <- rf_mtry$bestTune$mtry 
# best_mtry
# 
# store_maxnode <- list()
# tuneGrid <- expand.grid(.mtry = best_mtry)
# for (maxnodes in c(12: 25)) {
#   set.seed(1234)
#   rf_maxnode <- train(lables~.,
#                       data = train,
#                       method = "rf",
#                       metric = "Accuracy",
#                       tuneGrid = tuneGrid,
#                       trControl = trControl,
#                       importance = TRUE,
#                       nodesize = 14,
#                       maxnodes = maxnodes,
#                       ntree = 300)
#   current_iteration <- toString(maxnodes)
#   store_maxnode[[current_iteration]] <- rf_maxnode
# }
# results_mtry <- resamples(store_maxnode)
# summary(results_mtry)






##############################################################################################
# RFM_predict
# take the saved trained RFM
# takes seurat data object
# predicts cell types
# creates table and umap
# outputs prediction vector





##############################################################################################
# seurat_transfer
# takes in a reference seurat object
# follows seurat workflow, finds anchors, predicts labels
# default of no threshold, a threshold can be applied
# outputs prediction vector









##############################################################################################

# cluster_annotate
# input vectors of cluster annotation, at least 1 is required, prediction vector from visualization
# type manually, be careful to match syntax to the other predictions
# input seurat data object
# finds most common label per cluster and applies that label