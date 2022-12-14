# Create cell type clusters functions

# explore_param        (Shuming)
#Intrinsic stats       (Shuming)
# # clust_stability    (Shuming)

library(clusterSim) #for dbi
library(FlowSOM)
library(flowCore)
library(cluster) #for silhouette score
library(fpc) #for calinhara
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree) #for clustree plot
library(Rphenograph) #for phenograph
library(flexclust)#for adjusted rand index
library(ggplot2) #for the plot function

####################################################################################################



# explore_param
# reads in csv files with flow cytometry experiments or seurat data object depending on arguments
# runs FlowSom, Phenograph, Louvain (Seurat) 
# input arguments - cluster method, for flowsom k, for phenograph and louvain k paramater, for louvain resolution
# select outputs, generates - UMAPS, heatmaps, clustree  : save to folder or put in global environ
# creates list object to run stats


# input is a seurat object
explore_param <- function(input, #for phenograph and louvain only
                          cluster_method, #takes a list
                          df_input, #if input "flowsom", this is needed
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
                  df_input = df_input,
            for.flowsom.k = for.flowsom.k,
            save.stats = save.stats,
            save.plot = save.plot,
            output_path = output_path)
  }
  if ("phenograph" %in% cluster_method) {
    pg <- phenograph(input = input,
          for.phenograph.and.louvain.k = for.phenograph.and.louvain.k,
          output_path = output_path,
          save.stats = save.stats,
          save.plot = save.plot) 
  }
  if ("louvain" %in% cluster_method) {
    lv <- louvain(input = input, #seu object
            for.phenograph.and.louvain.k = for.phenograph.and.louvain.k,
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


#helper functions #1:
flowsom <- function(input, #seurat
                    df_input, #the processed df2 file before being converted to seurat
                    for.flowsom.k,
                    save.stats = TRUE,
                    save.plot = FALSE,
                    output_path) {
  
  clust_method <- "flowsom"
  # create the flowframe. If reading in a csv convert to flowset
  frame <- new("flowFrame", exprs = as.matrix(df_input)) #convert input to flowframe
  
  fs <- ReadInput(frame) #convert flowframe to flowsom object
  fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because I cleaned the df about before making flowset
  fs <- BuildMST(fs) # build minimum spanning tree
  
  
  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
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
  input <- FindNeighbors(input, dims = 1:12, k = kn)
  input <- RunUMAP(input, dims = 1:12, n.neighbors = kn)
  
  # save feature plots UMAP
  if(save.plot) {
    pdf(paste(output_path, clust_method, "UMAPfeatures_kn", kn, ".pdf", sep = ""), 
        width = 20, height = 10)
    print(FeaturePlot(input, features = rownames(input),
                      slot = 'scale.data',
                      min.cutoff = 'q1',
                      max.cutoff ='99',
                      label.size = 1)+
            theme(plot.title = element_text(size = 0.1)))
    dev.off()
  }
  
  for (i in for.flowsom.k){
    
    flowSOMcluster <- metaClustering_consensus(fs$map$codes, k = i, seed=42)
    
    clust_name = paste('FlowSom.k.', i, sep="")
    
    # add the cluster ID into seurat object to visualize
    input <- AddMetaData(object = input,
                       metadata = flowSOMcluster[fs$map$mapping[, 1]],
                       col.name = paste('FlowSom.k.', i, sep="")) #clustering name
    number.clusters <- length(unique(flowSOMcluster[fs$map$mapping[,1]]))
    
    # save feature plots of this UMAP
    if (save.plot) {
      png(paste(output_path, clust_method, "UMAPclusters_k", i, ".png", sep = ""))
      print(DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE, 
                    group.by = clust_name)) # will automatically group by active ident
      dev.off()
      
      # heatmap
      heatmap_name = paste("Heatmapclusters_k",i,".png",sep="")
      png(paste(output_path, clust_method, heatmap_name, sep = ""), 
          width = 600, height = 500)
      print(DoHeatmap(input, features = rownames(input), group.by = clust_name))
      dev.off()
     
    }
   
    # add stats
    stats_ls <- c(stats_ls, 
                  kn, #kn 
                  i, #krange
                  number.clusters, #number of clusters
                  mean(silhouette(flowSOMcluster[fs$map$mapping[, 1]][row_n], dis)[, 3]), #silhouette score
                  calinhara(m, flowSOMcluster[fs$map$mapping[, 1]], cn = i), #Calinski-Harabasz index
                  index.DB(df2, as.numeric(flowSOMcluster[fs$map$mapping[, 1]]))$DB) # Davies–Bouldin index
    
    cnl <- c(cnl, paste("clusters_krange_", i, "_", sep = "")) #update column name lists
    seul <- c(seul, input) #rds list of seu objects to be returned
  }
 
  #make statsl into matrix
  stats_ls <- matrix(stats_ls, ncol = 6, byrow = TRUE)
  colnames(stats_ls) <- c("kn", "krange", "nc", "si", "ch", "db")
 
  if (save.stats) {
    saveRDS(stats_ls, paste(output_path, clust_method, 'statslist.Rds', sep = ""))
  }
  
  # save feature plots of this UMAP
  if (save.plot) {
    if(length(for.flowsom.k) > 2) {
      # make clustree plot
      pdf(paste(output_path, clust_method, 'Clustree.pdf',sep = ""), width = 8, height = 8)
      print(clustree(input, prefix ='FlowSom.k.'))
      dev.off()
    }
    # save the UMAP with cell types
    pdf(paste(output_path, clust_method,'UMAPcelltype.pdf',sep=""),width =8, height = 6)
    print(DimPlot(input, group.by = 'Batch'))
    dev.off()
  }
  
  if (save.stats) { # save the stats
    stats_ls <- matrix(stats_ls, ncol = 6, byrow = TRUE)
    colnames(stats_ls) <- c("kn", "resolution", "nc", "si", "ch", "db")
    saveRDS(stats_ls, paste(output_path, clust_method, 'statslist.Rds', sep = ""))
  } 
  
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(output_path, clust_method, 'seul.Rds', sep = ""))
  if (save.stats) {
    # stats_plot(data.frame(stats_ls), output_path, clust_method)
    return(list(stats_ls, seul))} 
  else {
    return(seul)
    }
}


#helper functions #3:
phenograph <- function(input,
                       for.phenograph.and.louvain.k,
                       output_path,
                       save.stats = TRUE,
                       save.plot = FALSE) {
  
  clust_method <- "phenograph"
  # if save plot or stats but output_path is not given, give warning


  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
  row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
  dis <- daisy(m[row_n, ])  #dissimilarities distance

  
  #store seurat objects:
  cnl <- vector() #store names of clusters, names of seul 
  seul <- vector() #list of seu objects returned
  stats_ls <- vector() #stats list returned
  
  
  # we also only need to plot the features once
  # file name
  UMAP_name = paste("UMAPfeatures_kn", for.phenograph.and.louvain.k, ".pdf", sep = "") #*weird name 271, check later
  
  # save feature plots UMAP
  pdf(paste(output_path, clust_method, UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(input, features = rownames(input), slot = 'scale.data', 
                    min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ 
          theme(plot.title = element_text(size = 0.1)))
  dev.off()
  
  # we also want to see the batch on the UMAP
  pdf(paste(output_path, clust_method, UMAP_name,sep = ""), width =8, height = 6)
  print(DimPlot(input, group.by = 'Batch'))
  dev.off()
 
  
  for (i in for.phenograph.and.louvain.k){
    input <- FindNeighbors(input, dims = 1:12, k.param = i)
    input <- RunUMAP(input, dims = 1:12, n.neighbors = i) 

    ### run phenograph clustering
    Rphenograph_out_flow <- Rphenograph(m, k = i)
    phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
    clust_name = paste('Pheno.kn.',i,sep="")
    # add the cluster ID into seurat object to visualize
    input <- AddMetaData(object = input, phenocluster, col.name = clust_name)
    number.clusters <- length(unique(phenocluster))
    ### make umap
    
    if (save.plot) {
      UMAP_name = paste("UMAPclusters_kn",i,".pdf",sep="")
      pdf(paste(output_path, clust_method, UMAP_name,sep=""),width =20, height = 10)
      # save UMAP grouped
      print(DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
      dev.off()
      # heatmap
      pdf(paste(output_path,clust_method, paste("Heatmapclusters_kn",i,".pdf", sep=""),sep=""),width =25, height = 10)
      print(DoHeatmap(input, features = rownames(input), group.by = clust_name)) #doesn't work when the sample is smaller
      dev.off()
    }

    stats_ls <- c(stats_ls, 
                  i, #kn
                  number.clusters, #number of clusters
                  mean(silhouette(as.numeric(phenocluster[row_n]),dis)[, 3]), #silhouette score
                  calinhara(m,phenocluster,cn=i), #Calinski-Harabasz index
                  index.DB(df2, as.numeric(phenocluster))$DB # Davies–Bouldin index
    )
    
    cnl <- c(cnl, paste("clusters_kn_", i, sep = "")) #update column name lists
    seul <- c(seul, input) #rds list of seu objects to be returned
  }
  
  if (save.plot && (length(for.flowsom.k) > 2)) {
    # make clustree plot
    if(length(for.flowsom.k) > 2) 
    pdf(paste(output_path, clust_method,'Clustree.pdf',sep=""),width =15, height = 10)
    print(clustree(input, prefix ='Pheno.kn.'))
    dev.off()
  }
  
  if (save.stats) { # save the stats
    stats_ls <- matrix(stats_ls, ncol = 5, byrow = TRUE)
    colnames(stats_ls) <- c("kn", "nc","si", "ch", "db")
    saveRDS(stats_ls, paste(output_path, clust_method, 'statslist.Rds', sep = ""))
  }
  
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(output_path, clust_method, 'seul.Rds', sep = ""))
  if (save.stats) {return(list(stats_ls, seul))} else {return(seul)}
}


#helper functions #4:
louvain <- function(input, #seu object
                    for.phenograph.and.louvain.k,
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
    m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
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
  for (i in for.phenograph.and.louvain.k){
    input <- FindNeighbors(input, dims = 1:12, k.param = i)
    input <- RunUMAP(input, dims = 1:12, n.neighbors = i)
    
    # save feature plots of this UMAP
    if (save.plot) {
      pdf(paste(output_path, clust_method, "UMAPfeatures_kn", i, ".pdf", sep = ""),
          width =20, height = 10)
      print(FeaturePlot(input,
                        features = rownames(input),
                        slot = 'scale.data',
                        min.cutoff = 'q1',
                        max.cutoff ='99',
                        label.size = 1)+
              theme(plot.title = element_text(size = 0.1)))
      dev.off()
      
      # look at batches
      pdf(paste(output_path, clust_method, "UMAPbatches_kn", i, ".pdf", sep = ""),
          width = 20, height = 10)
      print(DimPlot(input, group.by = 'Batch', label.size = 1))
      dev.off()
    }
    
    for (j in resolutions) {
      input <- FindClusters(input, resolution = j)
      louvainCluster <- input@meta.data$seurat_clusters
      
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
        pdf(paste(output_path, clust_method, "UMAPclusters_kn", i, "_res_", j, ".pdf", sep = ""),
            width = 15, height = 10)
        print(DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE)) # will automatically group by active ident
        dev.off()
        
        # heatmap
        pdf(paste(output_path, clust_method, "Heatmapclusters_kn", i, "_res_", j, ".pdf", sep = ""),
            width = 15, height = 10)
        print(DoHeatmap(input, features = rownames(input), size = 10)+
                theme(text = element_text(size = 30)))
        dev.off()
      }
      cnl <- c(cnl, paste("clusters_kn_", i, "_res_", j, sep = "")) #update column name lists
      seul <- c(seul, input) #rds list of seu objects to be returned
    }
    #outside of resolution loop, in kn loop now:
    if (save.plot & (length(resolutions) > 1)) {
      # run clustree
      pdf(paste(output_path, clust_method, "kn", i, 'Clustree.pdf', sep = ""),
          width = 15, height = 10)
      print(clustree(input, prefix ='RNA_snn_res.'))
      dev.off()
    }

  }
  #outside of kn loop
  if (save.stats) { # save the stats
    stats_ls <- matrix(statsl, ncol = 6, byrow = TRUE)
    colnames(stats_ls) <- c("kn", "resolution", "nc", "si", "ch", "db")
    saveRDS(stats_ls, paste(output_path, clust_method, 'statslist.Rds', sep = ""))
  } 
  #return seu lists:
  names(seul) <- cnl
  saveRDS(seul, paste(output_path, clust_method, 'seul.Rds', sep = ""))
  
  if (save.stats) {return(list(stats_ls, seul))} else {return(seul)}
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
    input <- FindClusters(input, random.seed = rn_ls[x], resolution = j)
    return(Idents(object = input))
  }
  
  #find rand index between every 2 elements in a list, no repetition 
  #(ex: 1-2, 2-1) or same number (ex: 1-1, 2-2):
  hp2 <- function(x) { 
    ri<-randIndex(table(as.numeric(seul[, x[1]]), 
                        as.numeric(seul[, x[2]])))
    return(ri)
  }
  
  #find number of clusters
  hp3 <- function(x) { return((length(unique(x)))) }
  
  #main loops:
  for(i in kn) {
    #shuming note: why dims=12 here? 
    input <- FindNeighbors(input, dims = 1:12, k.param = i) 
    input <- RunUMAP(input, dims = 1:12, n.neighbors = i)
   
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

