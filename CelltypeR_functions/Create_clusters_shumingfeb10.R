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

#JAN2923EDIT:im going to make it so that you can only run 1 clustering at a time, 
#if the graph could not be run, i will report an error message and skip it


# input is a seurat object
explore_param <- function(input, #for phenograph and louvain only
                          cluster_method, #take 1 cluster method
                          df_input, #needed  if input "flowsom"
                          flow_k = NULL, #k for flowsom
                          pheno_lou_kn = NULL, #kn for phenograph or louvain
                          lou_resolution = NULL, #resolution for louvain
                          run.plot = FALSE, #print if run
                          run.stats = TRUE, #print and return if run
                          save_to = NULL) { #need it if save plots or stats
  
  
  #call helper functions depending on the cluster_method requested
  if (cluster_method == "louvain") {
    cl <- louvain(input = input, #seu object
                  df_input = df_input,
                  pheno_lou_kn = pheno_lou_kn,
                  resolutions = lou_resolution,
                  run.plot = run.plot,
                  run.stats = run.stats,
                  save_to = save_to) 
  }
  if (cluster_method == "phenograph") {
    cl <- phenograph(input = input,
                     df_input = df_input,
                     pheno_lou_kn = pheno_lou_kn,
                     run.stats = run.stats,
                     run.plot = run.plot,
                     save_to = save_to) 
  }
  if (cluster_method == "flowsom") {
    cl <- flowsom(input = input,
                  df_input = df_input,
                  flow_k = flow_k,
                  run.stats = run.stats,
                  run.plot = run.plot,
                  save_to = save_to)
  }

  if (run.plot || run.stats) {
    return(list(cluster = cl))
  }
}



#helper functions #1:
flowsom <- function(input, #seurat
                    df_input, #the processed df2 file before being converted to seurat
                    flow_k,
                    run.stats = TRUE,
                    run.plot = FALSE,
                    save_to) {

  clust_method <- "flowsom" #for naming the graphs
  
  # create the flowframe. If reading in a csv convert to flowset
  #only select numerical columns, exclude columns like X-column
  frame <- new("flowFrame", exprs = as.matrix(df_input %>% select(where(is.numeric)))) #convert input to flowframe
  fs <- ReadInput(frame) #convert flowframe to flowsom object
  fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because X column is removed
  fs <- BuildMST(fs) # build minimum spanning tree
  
  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  if (run.stats) {
    m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
    row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
    dis <- daisy(m[row_n, ])  #dissimilarities distance
  }
  
  # #store seurat objects and stats:
  stats_ls <- vector() #create a list to store all stats
  
  # kn = round(sqrt(dim(df2)[1]))
  input <- FindNeighbors(input, dims = 1:12)
  input <- RunUMAP(input, dims = 1:12)
  
  # save feature plots UMAP
  if(run.plot) {
    p1 <- FeaturePlot(input, features = rownames(input),
                      slot = 'scale.data',
                      min.cutoff = 'q1',
                      max.cutoff ='99',
                      label.size = 1)+
      theme(plot.title = element_text(size = 0.1))
    print(p1)
    
    # see sample on the UMAP
    p2 <- DimPlot(input, group.by = 'Sample')
    print(p2)  
    
    if (!is.null(save_to)) {
      png(filename=paste(save_to, clust_method, "UMAPfeatures.png", sep = ""), width = 5000, height = 3000, res = 300)
      print(p1)
      dev.off()
      
      png(filename=paste(save_to, clust_method, "UMAP_sample.png", sep = ""), width = 5000, height = 3000, res = 300)
      print(p2)
      dev.off()
    }
  }

  for (i in flow_k){
    flowSOMcluster <- metaClustering_consensus(fs$map$codes, k = i, seed=42)

    clust_name = paste('FlowSom.k.', i, sep="")
    
    # add the cluster ID into seurat object to visualize
    input <- AddMetaData(object = input,
                       metadata = flowSOMcluster[fs$map$mapping[, 1]],
                       col.name = paste('FlowSom.k.', i, sep="")) #clustering name
    number.clusters <- length(unique(flowSOMcluster[fs$map$mapping[,1]]))
    
    # save feature plots of this UMAP
    if (run.plot) {
      p3 <- DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE, 
                   group.by = clust_name) # will automatically group by active ident
      p4 <- DoHeatmap(input, features = rownames(input), group.by = clust_name) # heatmap
  
      print(p3)
      print(p4)
      
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, "UMAPclusters_k", i, ".png", sep = ""), width = 2000, height = 1200, res = 300)
        print(p3)
        dev.off()
        
        png(filename = paste(save_to, clust_method, "_Heatmap_clusters_k",i,".png", sep = ""), width = 5000, height = 3000, res = 300)
        print(p4)
        dev.off()
      }
      
      suppressWarnings(rm(p3, p4))
    }
   
    if (run.stats) {
      # add stats
      stats_ls <- c(
        stats_ls, 
        # kn, i, #krange
        number.clusters, #number of clusters
        mean(silhouette(flowSOMcluster[fs$map$mapping[, 1]][row_n], dis)[, 3]), #silhouette score
        calinhara(m, flowSOMcluster[fs$map$mapping[, 1]], cn = i), #Calinski-Harabasz index
        index.DB(x=df_input %>% select(where(is.numeric)), 
                 cl=as.numeric(flowSOMcluster[fs$map$mapping[, 1]]))$DB) # Davies–Bouldin index
    }
  }
 
  if (run.stats) {
    #make statsl into matrix then dataframe
    stats_df <- as.data.frame(matrix(stats_ls, ncol = 4, byrow = TRUE))
    
    colnames(stats_df) <- c(
      # "kn", "krange", 
      "number.cluster", "silhouette.score", "calinski.harabasz", "davies.bouldin")
    
    print(stats_df)
    write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
  }
 
  # save feature plots of this UMAP
  if (run.plot) {
    if(length(flow_k) > 2) {
      # make clustree plot
      p5 <- clustree(input, prefix ='FlowSom.k.')
      print(p5)
      
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, 'Clustree.png',sep = ""), width = 1200, height = 2000, res = 300)
        print(p5)
        dev.off()
      }
    }
    # save the UMAP with cell types
    p6 <- DimPlot(input, group.by = 'Sample')
    print(p6)
    
    if (!is.null(save_to)) {
      png(filename=paste(save_to, clust_method,'UMAPcelltype.png',sep=""), width = 2000, height = 1200, res = 300)
      print(p6)
      dev.off()
    }
  }
  
  if (run.stats && !is.null(save_to)) {
    return(stats_df)
  }
}


#helper functions #2:
phenograph <- function(input,
                       df_input,
                       pheno_lou_kn,
                       run.stats = TRUE,
                       run.plot = FALSE,
                       save_to = NULL) {
  
  clust_method <- "phenograph"

  #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
  if (run.stats) {
    m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
    row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
    dis <- daisy(m[row_n, ])  #dissimilarities distance
  }

  #create a list to store stats:
  stats_ls <- vector() #stats list returned

  # print /& save feature plots UMAP
  if(run.plot) {
    p1 <- FeaturePlot(input, features = rownames(input), slot = 'scale.data', 
                      min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ 
      theme(plot.title = element_text(size = 0.1))
    print(p1)
    if (!is.null(save_to)) {
      png(filename=paste(save_to, clust_method, "UMAPfeatures", pheno_lou_kn, ".png", sep = ""), width = 5000, height = 3000, res = 300)
      print(p1)
      dev.off()
    }
    
    #sample on the UMAP
    p2 <- DimPlot(input, group.by = 'Sample')
    print(p2)
    if (!is.null(save_to)) {
      png(filename=paste(save_to, clust_method, "UMAP_dim_by_sample.png",sep = ""), width = 2000, height = 1200, res = 300)
      print(p2)
      dev.off()
    }
  }
  

  for (i in pheno_lou_kn){
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
    
    if (run.plot) {
      # heatmap
      p3 <- DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)
      print(p3)
      
      # save UMAP grouped
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, "_UMAPclusters_kn",i,".png", sep=""), width =2000, height = 1200, res=300)
        print(p3)
        dev.off()
      }
      suppressWarnings(rm(p3))
    }
    
    stats_ls <- c(stats_ls, 
                  i, #kn
                  number.clusters, #number of clusters
                  mean(silhouette(as.numeric(phenocluster[row_n]),dis)[, 3]), #silhouette score
                  calinhara(m,phenocluster,cn=i), #Calinski-Harabasz index
                  index.DB(x=df_input %>% select(where(is.numeric)),cl=as.numeric(phenocluster))$DB # Davies–Bouldin index
    )
  }
  
  if (run.plot && (length(pheno_lou_kn) > 2)) {
    # make clustree plot
    if(length(pheno_lou_kn) > 2) {
      p4 <- clustree(input, prefix ='Pheno.kn.')
      print(4)
      
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method,'Clustree.png',sep=""), width = 1000, height = 2000, res=300)
        print(clustree(input, prefix ='Pheno.kn.'))
        dev.off()
      }
    }
  }
  
  if (run.stats) { # save the stats
    stats_df <- as.data.frame(matrix(stats_ls, ncol = 5, byrow = TRUE))
    colnames(stats_df) <- c("kn",  "number.cluster", "silhouette.score", "calinski.harabasz", "davies.bouldin")
    write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
  }
  
  if (run.stats && !is.null(save_to)) {
    return(stats_df)
  }

}


#helper functions #3:
louvain <- function(input, #seu object
                    df_input,
                    pheno_lou_kn,
                    resolutions,
                    run.plot = FALSE, #option to save the graphs  
                    run.stats = TRUE, #option to save stats list  
                    save_to #only required when save is TRUE
) {
  clust_method <- "louvain"
  if (run.stats) {
    #subsampling for silhouette score, n=9000, if less than 9000 cells, use the full size 
    m <- t(as.matrix(GetAssayData(object = input, slot = "counts"))) 
    row_n <- sample(1 : nrow(m), ifelse(nrow(m) > 9000, 9000, nrow(m)))
    dis <- daisy(m[row_n, ])  #dissimilarities distance
    statsl <- vector() #stats list returned
  } 
  
  # In the loop
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  for (i in pheno_lou_kn){
    input <- FindNeighbors(input, dims = 1:12, k.param = i)
    input <- RunUMAP(input, dims = 1:12, n.neighbors = i)
    
    # save feature plots of this UMAP
    if (run.plot) {
      p1 <- FeaturePlot(input,
                        features = rownames(input),
                        slot = 'scale.data',
                        min.cutoff = 'q1',
                        max.cutoff ='99',
                        label.size = 1)+
        theme(plot.title = element_text(size = 0.1))
      
      p2 <- DimPlot(input, group.by = 'Sample', label.size = 1)
      
      print(p1)
      print(p2)
      
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, "UMAPfeatures_kn", i, ".pdf", sep = ""), width = 2000, height = 1200, res = 300)
        print(p1)
        dev.off()
        # look at batches
        png(filename=paste(save_to, clust_method, "UMAP_Sample_kn", i, ".pdf", sep = ""), width = 2000, height = 1200, res = 300)
        print(p2)
        dev.off()
      }
      suppressWarnings(rm(p1, p2))
    }
    
    for (j in resolutions) {
      input <- FindClusters(input, resolution = j)
      louvainCluster <- input@meta.data$seurat_clusters
      
      #stats
      if (run.stats) {
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
                      index.DB(x=df_input %>% select(where(is.numeric)),cl=as.numeric(louvainCluster))$DB) # Davies–Bouldin index
        }
      } 
      if (run.plot) {# make UMAP grouped plots
        p3 <- DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE)
        p4 <- DoHeatmap(input, features = rownames(input), size = 10)+
          theme(text = element_text(size = 30)) # heatmap
        print(p3)
        print(p4)
        
        if (!is.null(save_to)) {
          png(filename=paste(save_to, clust_method, "UMAPclusters_kn", 
                                 i, "_res_", j, ".pdf", sep = ""), 
              width = 2000, height = 1200, res = 300)
          print(p3)
          dev.off()
          
          png(filename=paste(save_to, clust_method, "Heatmapclusters_kn", 
                             i, "_res_", j, ".pdf", sep = ""), 
              width = 2000, height = 1200, res = 300)
          print(p4)
          dev.off()
        }
        suppressWarnings(rm(p3, p4))
      }
    }
    #outside of resolution loop, in kn loop now:
    if (run.plot & (length(resolutions) > 1)) {
      # run clustree
      p5 <- clustree(input, prefix ='RNA_snn_res.')
      print(p5)
      
      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, "kn", i, "_res_", j, 'Clustree.pdf', sep = ""), 
            width = 1200, height = 2000, res = 300)
        print(p5)
        dev.off()
      }
      suppressWarnings(rm(p5))
    }
  }
  #outside of kn loop
  if (run.stats) { # save the stats
    stats_df <- as.data.frame(matrix(statsl, ncol = 6, byrow = TRUE))
    colnames(stats_df) <- c("kn", "resolution", "number.cluster", "silhouette.score", "calinski.harabasz", "davies.bouldin")
    write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
  } 
  
  if (run.stats && !is.null(save_to)) {
    return(stats_df)
  }
}



####################################################################################################

#Intrinsic stats
# takes in a dataframe or list of statistics
# plots intrinsic statistics from pararmater explorations

# c("kn", "resolution", "number.cluster", "silhouette.score", "calinski.harabasz", "davies.bouldin")

stats_plot <- function(stats_ls,
                       save_to,
                       clust_method) {
  # drop any rows containing NAs, they contain NAs because some kn x res
  #give number of clusters of 1 (there is no split), and you can't run
  #internal stats on them
  stats_ls <- stats_ls[complete.cases(stats_ls), ]
  
  if (clust_method == "louvain") {
    ##silhouette score: ranges from -1  to 1
    ##-1: bad clusters  0: neutral, indifferent  1: good clusters
    
    #x axis = number of cluster
    siplot1 <- ggplot(stats_ls, aes(x = number.cluster, y = silhouette.score, label = resolution)) +
      geom_line(aes(group = kn, color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)),
                check_overlap = TRUE,
                position = position_jitter(width = 0.2),
                size = 3) +
      labs(color = "kn", title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    #x axis = kn
    siplot2 <- ggplot(stats_ls, aes(kn, silhouette.score)) +
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
    chplot1 <- ggplot(stats_ls, aes(x = number.cluster, y = calinski.harabasz, label = resolution)) +
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
    chplot2 <- ggplot(stats_ls, aes(kn, calinski.harabasz)) +
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
                      aes(x = number.cluster, y = db, label = resolution)) +
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
      geom_point(aes(x = krange, y = silhouette.score)) +
      geom_line(aes(x = krange, y = silhouette.score), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "krange", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    siplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = number.cluster, y = silhouette.score)) +
      geom_line(aes(x = number.cluster, y = silhouette.score), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = krange, y = calinski.harabasz)) +
      geom_line(aes(x = krange, y = calinski.harabasz), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "krange", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = number.cluster, y = calinski.harabasz)) +
      geom_line(aes(x = number.cluster, y = calinski.harabasz), size = 0.1) +
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
      geom_point(aes(x = number.cluster, y = db)) +
      geom_line(aes(x = number.cluster, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
    
  } else if (clust_method == "phenograph") {
    siplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = kn, y = silhouette.score)) +
      geom_line(aes(x = kn, y = silhouette.score), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "kn", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    siplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = number.cluster, y = silhouette.score)) +
      geom_line(aes(x = number.cluster, y = silhouette.score), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot1 <- ggplot(stats_ls) +
      geom_point(aes(x = kn, y = calinski.harabasz)) +
      geom_line(aes(x = kn, y = calinski.harabasz), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "kn", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
    
    chplot2 <- ggplot(stats_ls) +
      geom_point(aes(x = number.cluster, y = calinski.harabasz)) +
      geom_line(aes(x = number.cluster, y = calinski.harabasz), size = 0.1) +
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
      geom_point(aes(x = number.cluster, y = db)) +
      geom_line(aes(x = number.cluster, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))
  } else (
    warning("clust_method is incorrect. ")
  )
  
  pdf(paste(save_to, clust_method, "_stats_plot.pdf", sep = ""))
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
                       save_to = NULL  #if null, will not save seul, ril, ncl, rdf
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
      if (!is.null(save_to)) {
        saveRDS(seul,paste(save_to, "seu_ls_kn", i, "_j", j, ".Rds",sep=""))
        saveRDS(ncl,paste(save_to, "nc_ls_kn", i, "_j", j, ".Rds",sep=""))
        saveRDS(ncl,paste(save_to, "ri_ls_kn", i, "_j", j, ".Rds",sep=""))
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

