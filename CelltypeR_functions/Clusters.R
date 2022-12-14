# Create cell type clusters functions

# explore_param        (Shuming)
#Intrinsic stats       (Shuming)
# # clust_stability    (Shuming)


library(clusterSim) #new package for dbi
library(FlowSOM)
library(flowCore)
library(cluster) #for silhouette score
library(fpc) #for calinhara
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree) #for clustree plot
library(Rphenograph) #for phenograph



# install.packages('flexclust') 
library(flexclust)#for adjusted rand index
library(ggplot2) #for the plot function

input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/Old/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
df <- read.csv(input_path)
# subsample <- sample(nrow(df), 3000) #subsample
# df <- df[subsample, ]
df2 <- df %>% dplyr::select(c("AQP4", "CD24", "CD44", "CD184", "CD15",
                              "HepaCAM", "CD29", "CD56", "O4", "CD140a",
                              "CD133", "GLAST", "CD71"))

tm <- t(df2)
rownames(tm) <- colnames(df2)
colnames(tm) <- rownames(df2)
s <- CreateSeuratObject(tm)
s <- AddMetaData(object=s, metadata=df$Batch, col.name = 'Batch')
AB <- colnames(df2) # save antibody names for feature plotting later
s <- ScaleData(s) # add to scale data slot
print(DoHeatmap(s, group.by = "Batch", features = AB)) # check the data
s <- RunPCA(s, features = AB, npcs = 12, approx = FALSE)

#*shuming notes end


####################################################################################################



# #*testing:
# test1 <- explore_param(input = s, 
#               cluster_method = c("phenograph"), 
#               for.phenograph.and.louvain.k = c(20, 25, 30),
#               save.plot = TRUE, 
#               output_path = "/Users/shumingli/Desktop/nov4/")
# test2 <- explore_param(input = s, 
#                       cluster_method = c("flowsom"), 
#                       for.flowsom.k = c(3, 4, 5),
#                       save.plot = TRUE, 
#                       output_path = "/Users/shumingli/Desktop/nov4/")
# test3 <- explore_param(input = s, 
#                        cluster_method = c("louvain"), 
#                        for.phenograph.and.louvain.k = c(20, 25, 30),
#                        for.louvain.resolution = c(0.1,0.2,0.3),
#                        save.plot = TRUE, 
#                        output_path = "/Users/shumingli/Desktop/nov4/")

test$phenograph[[1]]

#*end testing

# explore_param
# reads in csv files with flow cytometry experiments or seurat data object depending on arguments
# runs FlowSom, Phenograph, Louvain (Seurat) 
# input arguments - cluster method, for flowsom k, for phenograph and louvain k paramater, for louvain resolution
# select outputs, generates - UMAPS, heatmaps, clustree  : save to folder or put in global environ
# creates list object to run stats

explore_param <- function(input, 
                          cluster_method, #takes a list
                          for.flowsom.k = NULL, 
                          for.phenograph.and.louvain.k = NULL, 
                          for.louvain.resolution = NULL,
                          save.plot = FALSE, 
                          save.stats = TRUE,
                          output_path = NULL) {
  # if output_path is not given, set output_path to working directory 
  if(is.null(output_path)) {
    output_path <- getwd()
    
  }
  if (save.plot || save.stats) { 
    dir.create(file.path(output_path), showWarnings = FALSE)
    setwd(file.path(output_path))
    fs <- NULL
    pg <- NULL
    lv <- NULL
  }

  
  if ("flowsom" %in% cluster_method) {
    fs <- flowsom(input = input,
            for.flowsom.k = for.flowsom.k,
            save.stats = save.stats,
            save.plot = save.plot,
            output_path = output_path)
  }
  if ("phenograph" %in% cluster_method) {
    pg <- phenograph(seu = input,
          for.phenograph.and.louvain.k = for.phenograph.and.louvain.k,
          output_path = output_path,
          save.stats = save.stats,
          save.plot = save.plot) 
  }
  if ("louvain" %in% cluster_method) {
    lv <- louvain(seu = input, #seu object
            kn = for.phenograph.and.louvain.k,
            resolutions = for.louvain.resolution,
            save.plot = save.plot, #option to save the graphs  
            save.stats = save.stats, #option to save stats list  
            output_path #only required when save is TRUE
    )
  }
  if (save.plot || save.stats) {
    return(list(flowsom = fs, phenograph = pg, louvain = lv))
  }
}


#helper functions


# test4 <- flowsom(input = s,
#                 for.flowsom.k = c(3, 4, 5),
#                 save.plot = TRUE,
#                 output_path = "/Users/shumingli/Desktop/nov4/")
# input <- s
# for.flowsom.k <- c(3, 4, 5)


#helper functions #2:
flowsom <- function(input, #csv
                    for.flowsom.k,
                    save.stats = TRUE,
                    save.plot = FALSE,
                    output_path) {
  clust_method <- "flowsom"
  # create the flowframe. If reading in a csv convert to flowset
  frame <- new("flowFrame", exprs = as.matrix(df2)) #convert input to flowframe
  fs <- ReadInput(frame) #convert flowframe to flowsom object
  fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because I cleaned the df about before making flowset
  fs <- BuildMST(fs) # build minimum spanning tree
  
  
  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  m <- t(as.matrix(GetAssayData(object = seu, slot = "counts"))) 
  row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
  dis <- daisy(m[row_n, ])  #dissimilarities distance
  # statsl <- vector() #stats list returned
  
  
  #store seurat objects and stats:
  cnl <- vector() #store names of clusters, names of seul 
  seul <- vector() #list of seu objects returned
  stats_ls <- vector() #create a list to store all stats
  
  #subsample for silhouette score
  m <- as.matrix(df2) # create a matrix for later
  row_n <- sample(1:nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m))) 
  dis <- dist(m[row_n, ]) 
  
  kn = round(sqrt(dim(df2)[1]))
  seu <- FindNeighbors(seu, dims = 1:12, k = kn)
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn)
  
  # save feature plots UMAP
  if(save.plot) {
    pdf(paste(clust_method, "UMAPfeatures_kn", kn, ".pdf", sep = ""), 
        width = 20, height = 10)
    print(FeaturePlot(seu, features = AB,
                      slot = 'scale.data',
                      min.cutoff = 'q1',
                      max.cutoff ='99',
                      label.size = 1)+
            theme(plot.title = element_text(size = 0.1)))
    dev.off()
  }
  
  for (i in for.flowsom.k){
    
    flowSOMcluster <- metaClustering_consensus(fs$map$codes, k = i, seed=42)
    print("flowsom test1")
    clust_name = paste('FlowSom.k.', i, sep="")
    
    # add the cluster ID into seurat object to visualize
    seu <- AddMetaData(object = seu,
                       metadata = flowSOMcluster[fs$map$mapping[, 1]],
                       col.name = paste('FlowSom.k.', i, sep="")) #clustering name
    number.clusters <- length(unique(flowSOMcluster[fs$map$mapping[,1]]))
    print("flowsom test2")
    # save feature plots of this UMAP
    if (save.plot) {
      png(paste(clust_method, "UMAPclusters_k", i, ".png", sep = ""))
      print(DimPlot(seu, reduction = "umap", repel = TRUE, label = TRUE, 
                    group.by = clust_name)) # will automatically group by active ident
      dev.off()
      print("flowsom test3")
      # heatmap
      heatmap_name = paste("Heatmapclusters_k",i,".png",sep="")
      png(paste(output_path, clust_method, heatmap_name, sep = ""), 
          width = 600, height = 500)
      print(DoHeatmap(seu, features = AB,group.by = clust_name))
      dev.off()
      print("flowsom test4")
    }
    
    # add stats
    stats_ls <- c(stats_ls, 
                  kn, #kn 
                  i, #krange
                  number.clusters, #number of clusters
                  mean(silhouette(flowSOMcluster[fs$map$mapping[, 1]][row_n], dis)[, 3]), #silhouette score
                  calinhara(m, flowSOMcluster[fs$map$mapping[, 1]], cn = i), #Calinski-Harabasz index
                  index.DB(df2, as.numeric(flowSOMcluster[fs$map$mapping[, 1]]))$DB) # Davies–Bouldin index
    print("flowsom test5")
    cnl <- c(cnl, paste("clusters_krange_", i, "_", sep = "")) #update column name lists
    seul <- c(seul, seu) #rds list of seu objects to be returned
  }
  print("flowsom test6")
  #make statsl into matrix
  stats_ls <- matrix(stats_ls, ncol = 6, byrow = TRUE)
  colnames(stats_ls) <- c("kn", "krange", "nc", "si", "ch", "db")
  print("flowsom test7")
  if (save.stats) {
    saveRDS(stats_ls, paste(output_path, clust_method, 'statslist.Rds', sep = ""))
  }
  
  
  print("flowsom test8")
  # save feature plots of this UMAP
  if (save.plot) {
    # make clustree plot
    pdf(paste(output_path, clust_method, 'Clustree.pdf',sep = ""), width = 8, height = 8)
    print(clustree(seu, prefix ='FlowSom.k.'))
    dev.off()
    print("flowsom test9")
    # save the UMAP with cell types
    pdf(paste(output_path, clust_method,'UMAPcelltype.pdf',sep=""),width =8, height = 6)
    print(DimPlot(seu,group.by = 'Batch'))
    dev.off()
    print("flowsom test10")
  }
  
  if (save.stats) { # save the stats
    stats_ls <- matrix(stats_ls, ncol = 6, byrow = TRUE)
    colnames(stats_ls) <- c("kn", "resolution", "nc", "si", "ch", "db")
    saveRDS(stats_ls, paste(clust_method, 'statslist.Rds', sep = ""))
    print("flowsom test11")
  } 
  
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(clust_method, 'seul.Rds', sep = ""))
  print("flowsom test12")
  if (save.stats) {
    # stats_plot(data.frame(stats_ls), output_path, clust_method)
    return(list(stats_ls, seul))} 
  else {
    return(seul)
    }
}

# #*test:
# test <- phenograph(seu = s, 
#         for.phenograph.and.louvain.k = c(25, 30, 35),
#         save.plot = TRUE,
#         output_path = "/Users/shumingli/Desktop/nov4/")
# seu <- s
# for.phenograph.and.louvain.k = c(25, 30, 35)
# save.plot = TRUE
# save.stats = TRUE
# output_path = "/Users/shumingli/Desktop/nov4/"

#helper functions #3:
phenograph <- function(seu,
                       for.phenograph.and.louvain.k,
                       output_path,
                       save.stats = TRUE,
                       save.plot = FALSE) {
  
  clust_method <- "phenograph"
  # if save plot or stats but output_path is not given, give warning


  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  m <- t(as.matrix(GetAssayData(object = seu, slot = "counts"))) 
  row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
  dis <- daisy(m[row_n, ])  #dissimilarities distance

  
  #store seurat objects:
  cnl <- vector() #store names of clusters, names of seul 
  seul <- vector() #list of seu objects returned
  stats_ls <- vector() #stats list returned
  
  
  # we also only need to plot the features once
  # file name
  UMAP_name = paste("UMAPfeatures_kn", for.phenograph.and.louvain.k, ".pdf", sep = "") #*weird name 271, check later
  # print(UMAP_name) #testing
  
  # save feature plots UMAP
  pdf(paste(output_path, clust_method, UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(seu, features = AB,slot = 'scale.data', 
                    min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ 
          theme(plot.title = element_text(size = 0.1)))
  dev.off()
  print("test7")
  # we also want to see the batch on the UMAP
  pdf(paste(output_path, clust_method, UMAP_name,sep = ""), width =8, height = 6)
  print(DimPlot(seu, group.by = 'Batch'))
  dev.off()
  print("test8")
  
  # #subsampling for silhouette score, n=1000, can make n bigger if needed
  # row_n <- sample(1:nrow(m), 1000)
  # dis <- dist(m[row_n,])
  
  print("test9")
  
  ############################# loop to explore parameters ########################################
  # kn = c(25,50,75,100,125,150,175,200,225,250,300,350,400,450,500)
  # kn = c(25,50,75,100,125,150,175,200,225,250,275,300)
  # larger kn fewer clusters in general but not always
  #kn = c(50,500)
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  
  for (i in for.phenograph.and.louvain.k){
    # kn_umap = round(sqrt(dim(df2)[1]))
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i) #*was outside of the loop with a calculated kn of 271, i moved it in, check later
    print("test5")
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i) #* same here, check later
    print("test6")
    # save feature plots of this UMAP
    # just for testing print
    
    ### run phenograph clustering
    Rphenograph_out_flow <- Rphenograph(m, k = i)
    phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
    print("test10")
    clust_name = paste('Pheno.kn.',i,sep="")
    # add the cluster ID into seurat object to visualize
    seu <- AddMetaData(object=seu, phenocluster, col.name = clust_name)
    print("test11")
    number.clusters <- length(unique(phenocluster))
    print("test12")
    ### make umap
    
    if (save.plot) {
      UMAP_name = paste("UMAPclusters_kn",i,".pdf",sep="")
      print(UMAP_name) #testing
      pdf(paste(output_path, clust_method, UMAP_name,sep=""),width =20, height = 10)
      # save UMAP grouped
      print(DimPlot(seu, reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
      dev.off()
      # heatmap
      heatmap_name = paste("Heatmapclusters_kn",i,".pdf",sep="")
      # testing
      pdf(paste(output_path,clust_method,heatmap_name,sep=""),width =25, height = 10)
      print(DoHeatmap(seu, features = AB, group.by = clust_name)) #doesn't work when the sample is smaller
      dev.off()
    }
    print("test13")
    #### add stats
    # "kn", "nc","si", "ch", "db"
    # get the cluster indexes
    # phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
    print("test14")
    stats_ls <- c(stats_ls, 
                  i, #kn
                  number.clusters, #number of clusters
                  mean(silhouette(as.numeric(phenocluster[row_n]),dis)[, 3]), #silhouette score
                  calinhara(m,phenocluster,cn=i), #Calinski-Harabasz index
                  index.DB(df2, as.numeric(phenocluster))$DB # Davies–Bouldin index
    )
    print("test15")
    cnl <- c(cnl, paste("clusters_kn_", i, sep = "")) #update column name lists
    seul <- c(seul, seu) #rds list of seu objects to be returned
  }
  
  if (save.plot) {
    # make clustree plot
    pdf(paste(output_path, clust_method,'Clustree.pdf',sep=""),width =15, height = 10)
    print(clustree(seu, prefix ='Pheno.kn.'))
    dev.off()
    print("test16")
  }
  
  if (save.stats) { # save the stats
    stats_ls <- matrix(stats_ls, ncol = 5, byrow = TRUE)
    colnames(stats_ls) <- c("kn", "nc","si", "ch", "db")
    saveRDS(stats_ls, paste(clust_method, 'statslist.Rds', sep = ""))
    print("test17")
  }
  
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(clust_method, 'seul.Rds', sep = ""))
  
  if (save.stats) {return(list(stats_ls, seul))} else {return(seul)}
}


#helper functions #4:
louvain <- function(seu, #seu object
                    kn,
                    resolutions,
                    save.plot = FALSE, #option to save the graphs  
                    save.stats = TRUE, #option to save stats list  
                    output_path #only required when save is TRUE
) {
  clust_method <- "louvain"
  # if save plot or stats but output_path is not given, give warning
  if ((save.plot && is.null(output_path)) || (save.stats && is.null(output_path))) { 
    stop("output_path is required if save.plot or save.stats is true")
  }  
  
  if (save.plot || save.stats) { #one of them true
    dir.create(file.path(output_path), showWarnings = FALSE)
    setwd(file.path(output_path))
  } 
  if (save.stats) {
    #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
    m <- t(as.matrix(GetAssayData(object = seu, slot = "counts"))) 
    row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
    dis <- daisy(m[row_n, ])  #dissimilarities distance
    statsl <- vector() #stats list returned
  } 
  #store seurat objects:
  cnl <- vector() #store names of clusters, names of seul 
  seul <- vector() #list of seu objects returned
  
  
  # In the loop
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  for (i in kn){
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    
    # save feature plots of this UMAP
    if (save.plot) {
      pdf(paste(clust_method, "UMAPfeatures_kn", i, ".pdf", sep = ""),
          width =20, height = 10)
      print(FeaturePlot(seu,
                        features = AB,
                        slot = 'scale.data',
                        min.cutoff = 'q1',
                        max.cutoff ='99',
                        label.size = 1)+
              theme(plot.title = element_text(size = 0.1)))
      dev.off()
      
      # look at batches
      pdf(paste(output_path, clust_method, "UMAPbatches_kn", i, ".pdf", sep = ""),
          width = 20, height = 10)
      print(DimPlot(seu, group.by = 'Batch', label.size = 1))
      dev.off()
    }
    
    for (j in resolutions) {
      seu <- FindClusters(seu, resolution = j)
      louvainCluster <- seu@meta.data$seurat_clusters
      
      #stats
      if (save.stats) {
        statsl <- c(statsl,
                    i, # kn
                    j, # resolution
                    length(unique(louvainCluster))) # number of clusters (nc
        
        if (length(unique(louvainCluster)) == 1)  {#skip the ones with only 1 cluster
          statsl <- c(statsl, rep(NA, 3))
        } else { 
          statsl <- c(statsl,
                      mean(silhouette(as.numeric(louvainCluster[row_n]), dis)[, 3]), # silhouette score, summary(dis)[[4]] = subsample number
                      calinhara(m, louvainCluster, cn = i),# Calinski-Harabasz index
                      index.DB(df2, as.numeric(louvainCluster))$DB) # Davies–Bouldin index
        }
      } 
      if (save.plot) {# make UMAP grouped plots
        pdf(paste(clust_method, "UMAPclusters_kn", i, "_res_", j, ".pdf", sep = ""),
            width = 15, height = 10)
        print(DimPlot(seu, reduction = "umap", repel = TRUE, label = TRUE)) # will automatically group by active ident
        dev.off()
        
        # heatmap
        pdf(paste(clust_method, "Heatmapclusters_kn", i, "_res_", j, ".pdf", sep = ""),
            width = 15, height = 10)
        print(DoHeatmap(seu, features = AB, size = 10)+
                theme(text = element_text(size = 30)))
        dev.off()
      }
      cnl <- c(cnl, paste("clusters_kn_", i, "_res_", j, sep = "")) #update column name lists
      seul <- c(seul, seu) #rds list of seu objects to be returned
    }
    #outside of resolution loop, in kn loop now:
    if (save.plot & (length(resolutions) > 1)) {
      # run clustree
      pdf(paste(clust_method, "kn", i, 'Clustree.pdf', sep = ""),
          width = 15, height = 10)
      print(clustree(seu, prefix ='RNA_snn_res.'))
      dev.off()
    }
    # # save seurat object
    # if (!save.plot) {
    #   saveRDS(seu, paste(input_name,
    #                      clust_method,
    #                      "SeuratObject", i, ".Rds", sep = ""))
    # }
  }
  #outside of kn loop
  if (save.stats) { # save the stats
    statsl <- matrix(stats_ls, ncol = 6, byrow = TRUE)
    colnames(statsl) <- c("kn", "resolution", "nc", "si", "ch", "db")
    saveRDS(statsl, paste(clust_method, 'statslist.Rds', sep = ""))
  } 
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(clust_method, 'seul.Rds', sep = ""))
  
  if (save.stats) {return(list(statsl, seul))} else {return(seul)}
}



####################################################################################################


#Intrinsic stats
# takes in a dataframe or list of statistics
# plots intrinsic statistics from pararmater explorations


stats_plot <- function(stats_ls,
                       output_path,
                       clust_method) {
  # drop any rows containing NAs, they contain NAs because some kn x res
  #give number of clusters of 1 (there is no split), and you can't run
  #internal stats on them
  stats_ls <- stats_ls[complete.cases(stats_ls), ]
  
  if (clust_method == "louvain") {
    ##silhouette score: ranges from -1  to 1
    ##-1: bad clusters  0: neutral, indifferent  1: good clusters
    
    #x axis = number of cluster
    siplot1 <- ggplot(stats_ls, aes(x = nc, y = si, label = resolution)) +
      geom_line(aes(group = kn, color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)),
                check_overlap = TRUE,
                position = position_jitter(width = 0.2),
                size = 3) +
      labs(color = "kn", title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    #x axis = kn
    siplot2 <- ggplot(stats_ls, aes(kn, si)) +
      geom_point(aes(colour = factor(resolution), group = (resolution))) +
      geom_line(aes(colour = factor(resolution), group = (resolution)), size = 0.2) +
      labs(title = "Silhouette Scores",
           x = "kn",
           y = "Average Silhouette Scores",
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))
    
    ##Calinski-Harabasz index:
    ## the highest value is the optimal number of clusters
    
    #x axis = number of cluster
    chplot1 <- ggplot(stats_ls, aes(x = nc, y = ch, label = resolution)) +
      geom_line(aes(group = kn,color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)),
                check_overlap = TRUE,
                position = position_jitter(width = 0.2),
                size = 3) +
      labs(color = "kn",
           title = "Calinski-Harabasz Index",
           x = "Number of Clusters",
           y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    #x axis = kn
    chplot2 <- ggplot(stats_ls, aes(kn, ch)) +
      geom_point(aes(colour = factor(resolution), group = factor(resolution))) +
      geom_line(aes(colour = factor(resolution), group = factor(resolution)), size = 0.2) +
      labs(title = "Calinski-Harabasz Index",
           x = "kn",
           y = "Calinski-Harabasz Index",
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ## #Davies–Bouldin index: minimum score is zero
    ## #the lowest value is the optimal number of clusters
    
    #x axis = number of cluster
    dbplot1 <- ggplot(stats_ls,
                      aes(x = nc, y = db, label = resolution)) +
      geom_line(aes(group = kn,color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)),
                check_overlap = TRUE,
                position = position_jitter(width = 0.2), size = 3) +
      labs(color = "kn", title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    #x axis = kn
    dbplot2 <- ggplot(stats_ls, aes(kn, db)) +
      geom_point(aes(colour = factor(resolution), group = factor(resolution))) +
      geom_line(aes(colour = factor(resolution), group = factor(resolution)),
                size = 0.2) +
      labs(title = "Davies-Bouldin index",
           x = "kn",
           y = "Davies-Bouldin index",
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))
    
  } else if (clust_method == "flowsom") {
    
    siplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = krange, y = si)) +
      geom_line(aes(x = krange, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "krange", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    siplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = si)) +
      geom_line(aes(x = nc, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = krange, y = ch)) +
      geom_line(aes(x = krange, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "krange", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = ch)) +
      geom_line(aes(x = nc, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "Number of Clusters", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    dbplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = krange, y = db)) +
      geom_line(aes(x = krange, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "krange", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    dbplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = db)) +
      geom_line(aes(x = nc, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
    
  } else if (clust_method == "phenograph") {
    siplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = kn, y = si)) +
      geom_line(aes(x = kn, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "kn", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    siplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = si)) +
      geom_line(aes(x = nc, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = kn, y = ch)) +
      geom_line(aes(x = kn, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "kn", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = ch)) +
      geom_line(aes(x = nc, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "Number of Clusters", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    dbplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = kn, y = db)) +
      geom_line(aes(x = kn, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "kn", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    dbplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = nc, y = db)) +
      geom_line(aes(x = nc, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
  } else (
    warning("clust_method is incorrect. ")
  )
  
  pdf(paste(output_path, clust_method, "_stats_plot.pdf", sep = ""))
  print(siplot1)
  print(siplot2)
  print(chplot1)
  print(chplot2)
  print(dbplot1)
  print(dbplot2)
  dev.off()
  return(list(siplot1, siplot2, chplot1, chplot2, dbplot1, dbplot2))
}


####################################################################################################

# clust_stability
# select cluster method and one pararmeter to vary (resolutions 0.1,0.3,0.5,0.8) 
# runs n iterations randomizing starting point - default 100
# calculates the RandIndex (RI) between each iteration of clustering for each resolution
# calculates the mean and standard deviation of the number of clusters and the RI for each resoloution
# outputs a table and plot of the results

#*input is s for this one
# #*test:
# df <- s
# resolutions <- c(0.1, 0.2)
# kn <- c(20, 25)
# n <- 3
# rdf <- Rand_index(df, resolutions, kn, n)


Rand_index <- function(input,
                       resolutions,
                       kn,
                       n = 100, #number of iterations
                       output_path = NULL  #if null, will not save seul, ril, ncl, rdf
) {
  
  #list of random integers
  rn_ls <- round(runif(n, min = 0, max = 100000), 0)
  
  #final df with mean sd rand index and nc
  rdf <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(rdf) <- c('kn', 'resolution', 
                     'meanri', 'medianri', 'sdri', 
                     'meannc', 'mediannc', 'sdnc')
  
  #helper functions:
  hp1 <- function(x) { #store n repeats of the same ixj clustering
    seu <- FindClusters(seu, random.seed = rn_ls[x], resolution = j)
    return(Idents(object = seu))
  }
  
  #find rand index between every 2 elements in a list, no repetition 
  #(ex: 12, 21), no same number (ex: 11, 22):
  hp2 <- function(x) { 
    #1. this only caluclates ari, no subsampling needed:
    ri<-randIndex(table(as.numeric(seul[, x[1]]), 
                        as.numeric(seul[, x[2]])))
    
    ##2. this calculates ari and ri, but it's slowers:
    # ri <- comPart(unlist(seu[[paste("repeat_", x[1], sep = "")]]), 
    # unlist(seu[[paste("repeat_", x[2], sep = "")]]), type=c("ARI","RI"))
    
    ##3. this is the fossol function, toooo slow, and has to subsample:
    # ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", x[1], sep = "")]), 
    # as.numeric(seu@meta.data[row_n, paste("repeat_", x[2], sep = "")]))
    return(ri)
  }
  
  #find number of clusters
  hp3 <- function(x) { return((length(unique(x)))) }
  
  #main loops:
  for(i in kn) {
    #shuming note: why dims=12 here? 
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i) 
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    for (j in resolutions) {
      seul <- sapply(1:n, hp1) #list of n repeats of clustering (Ident(seu) object)
      ncl <- apply(seul, 2, hp3) #list of number of clustering 
      ril <- apply(t(combn(1:n, 2)), 1, hp2) #list of rand index
      
      if (!is.null(output_path)) {
        saveRDS(seul,paste(output_path, "seu_ls_kn", i, "_j", j, ".Rds",sep=""))
        saveRDS(ncl,paste(output_path, "nc_ls_kn", i, "_j", j, ".Rds",sep=""))
        saveRDS(ncl,paste(output_path, "ri_ls_kn", i, "_j", j, ".Rds",sep=""))
      }
      
      rdf <-rbind(rdf, list(
        kn = i, resolution = j,
        meanri = mean(ril), medianri = median(ril), sdri = sd(ril), 
        meannc = mean(ncl), mediannc = median(ncl), sdnc = sd(ncl)))
    }
  }
  p <- plot_randindex(rdf, c('pink', 'violet'), c(0.7, 1))
  return(list(list = rdf, figure = p))
}


#helper function:
plot_randindex <- function (
    rdf,
    cp = c("orange", "violet"),
    view = c(0, 1) #zoom in x axis, this does not changes scales, just the viewable sections
) {
  
  s <- (view[2]-view[1])/max(rdf$meannc+rdf$sdnc)
  
  p <- rdf %>% 
    ggplot(aes(x = resolution)) +
    geom_line(aes(y = meanri), color = cp[1]) +
    geom_point(aes(y = meanri), color = cp[1], size=1)+ 
    geom_errorbar(aes(ymin=meanri-sdri, ymax=meanri+sdri), color = cp[1], width=.01,
                  position=position_dodge(.9))+
    geom_line(aes(y = meannc*s+view[1]), color = cp[2]) +
    geom_point(data = rdf, mapping = aes(x = resolution, y = meannc*s+view[1]), color = cp[2])+ 
    geom_errorbar(rdf, mapping = aes(x=resolution, ymin=((meannc-sdnc)*s+view[1]), ymax=((meannc+sdnc)*s+view[1])), width=.01,
                  position=position_dodge(.9), color = cp[2])+
    scale_y_continuous(limits= view, name="Mean Rand Index",
                       sec.axis = sec_axis(~ . /s - view[1]/s, 
                                           name = "Mean Number of Clusters"))+
    theme(axis.text.y  = element_text(color = cp[1]),
          axis.title.y = element_text(color=cp[1]),
          axis.text.y.right =  element_text(color = cp[2]),
          axis.title.y.right = element_text(color=cp[2]),
          plot.title = element_text(hjust = 0.5, size=10))+
    ggtitle("Plot of Mean Rand Index and \n Mean Number of Clusters")
  
  return(p)
}


