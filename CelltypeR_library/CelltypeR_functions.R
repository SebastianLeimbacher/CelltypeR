# CellTypeR functions



####### Preprocessing functions to create seurat object from Flow cytometry multichannel data ####

# fsc_to_df
# reads in fsc files creates a flowset file
# filtering if wanted on each input file
# users can rename samples after


fsc_to_fs <- function(input_path, downsample = 'none'){ 
  flowset = read.flowSet(path=input_path,transformation = FALSE ,
                         emptyValue = FALSE,truncate_max_range = FALSE, 
                         package="flowCore")
  fsApply(flowset,
          function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(flowset))]})
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  if (downsample == 'none') { 
    return(flowset)
  } else if (downsample == 'min') {
    desired_size <- min(fsApply(flowset,function(x){nrow(x@exprs)}))
    sf <- sampleFilter(filterId = "SizeFilter", size =desired_size) #Creates a "filter" object to subset
    flowset <- fsApply(flowset,function(x){Subset(x,sf,truncate_max_range = FALSE)}) #apply the filter on every flowframe
    return(flowset)
  } else {
    desired_size <- downsample
    set.seed(42) #Set a seed for reproducibility 
    sf <- sampleFilter(filterId = "SizeFilter", size =desired_size) #Creates a "filter" object to subset
    flowset <- fsApply(flowset,function(x){Subset(x,sf,truncate_max_range = FALSE)}) #apply the filter on every flowframe
    return(flowset)
  }
}






##############################################################################################
# harmonize
# transform, align, retro transform flowset object and save a csv file

# will be called by transform function
inversebiexponentialTransform <- function(flowset,a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0){
  copy_flowset=flowset[seq(along=flowset)] #makes a copy of the input flowset
  for (i in 1:length(copy_flowset)){ #loop though index to get each flowframe
    copy_flowset[[i]]@exprs=a*exp(b*(copy_flowset[[i]]@exprs-w))-c*exp(-d*(copy_flowset[[i]]@exprs-w))+f
  }
  return(copy_flowset)
}


#### 
### process to different levels - 
### if alignment isn't desired there is no need to transform
## this function can be used to see different transformations
## processing = 'retro' will biex transform, align, retro transform
## processing = 'align' will biex transform, align
## processing = 'biexp' will biex transform



harmonize <-  function(flowset, processing = 'retro', 
                       two_peaks = c(9:length(colnames(transformed_flowset))),
                       one_peak = c(1:8), threshold = 0.01) {
  # biexp transform conditions
  biexp  <- biexponentialTransform("biexp transform",a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0)
  # run biexp transform function 
  transformed_flowset <- transform(flowset, transformList(colnames(flowset), biexp))
  # if we want to see the biexp tranform
  if (processing == 'biexp') {  
    return(transformed_flowset)
  } else if (processing == 'align' ) {
    normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[two_peaks],max.lms = 2,peak.density.thr = threshold) #Detects and align 2 peaks on the marker 3,5,6,9...14. 
    expbe_norm2=normtr$flowset
    normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[one_peak],max.lms = 1,peak.density.thr = 0.05)#Detects and align 1 peak 
    aligned_transformed_flowset=normtr$flowset
    return(aligned_transformed_flowset)
  } else{
    normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[two_peaks],max.lms = 2,peak.density.thr = threshold) #Detects and align 2 peaks on the marker 3,5,6,9...14. 
    expbe_norm2=normtr$flowset
    normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[one_peak],max.lms = 1,peak.density.thr = 0.05)#Detects and align 1 peak 
    aligned_transformed_flowset=normtr$flowset
    retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset) 
    return(retrotransformed_flowset)
  }
}



#### plotting function for visulization ####
plotdensity_flowset <- function(flowset){ 
  ggplot(melt(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)})), 
         aes(x=value,y=L1,fill=L1)) + geom_density_ridges(alpha=.4,verbose=FALSE) +
    facet_wrap(~variable)+theme_light()} 






###### saving the csv 


rename_markers<-function(flowset){#Defines a function to use marker names 
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  return(copy_flowset)
}



flowset_to_csv=function(flowset, output_path, save.csv = FALSE){  
  list_of_flowframes=fsApply(rename_markers(flowset),function(x){as.data.frame(x@exprs)})#Makes a list of dataframes
  list_names=names(list_of_flowframes)#extract names for not calling names() function at each loop
  for (index in seq_along(list_of_flowframes)){ #Iterates along the index for adding sample names
    list_of_flowframes[[index]]$Sample = list_names[index]
    colnames(list_of_flowframes[[index]]) = colnames(list_of_flowframes[[1]]) 
  }
  ds=list.rbind(list_of_flowframes)#binds every fcs file in a single dataframe
  ds$cell=as.factor(unlist(lapply(as.list(c(1:length(flowset))),function(x){c(1:nrow(flowset[[x]]@exprs))})))#add cell IDs - cell count per sample 
  if (save.csv){
    write.csv(ds,file=paste0(output_path,deparse(substitute(flowset)),".csv"))#save the R data for further usage
  } else {
    return(ds)
  }
  
}

# add in an argument to save or return


##############################################################################################
# df_to_seurat
# creates a seurat object from the expression matrix and adds in meta data



make_seu <- function(df, AB_vector){
  df2 <- df %>% dplyr::select(all_of(AB_vector))
  m <- as.matrix(df2)
  tm <- t(df2)
  rownames(tm) <- colnames(df2)
  colnames(tm) <- rownames(df2)
  seu <- CreateSeuratObject(tm)
  seu <- AddMetaData(object=seu, metadata=df$Sample, col.name = 'Sample')
  #seu <- NormalizeData(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = AB_vector)
}




###### Functions to explore cluster parameters and create clusters in seurat object ############


# explore_param
# reads in csv files with flow cytometry experiments or seurat data object depending on arguments
# runs FlowSom, Phenograph, Louvain (Seurat) 
# input arguments - cluster method, for flowsom k, for phenograph and louvain k paramater, for louvain resolution
# select outputs, generates - UMAPS, heatmaps, clustree  : save to folder or put in global environ
# creates list object to run stats


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



# function #1:
flowsom <- function(input, #seurat
                    df_input, #the processed df2 file before being converted to seurat
                    flow_k,
                    run.stats = TRUE,
                    run.plot = TRUE,
                    save_to = NULL) {
  
  clust_method <- "flowsom" #for naming the graphs
  
  # create the flowframe. If reading in a csv convert to flowset
  #only select numerical columns, exclude columns like X-column
  frame <- new("flowFrame", exprs = as.matrix(df_input[, sapply(df_input, is.numeric)])) #convert input to flowframe, remove non-numerical col
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
        index.DB(x = df_input[, sapply(df_input, is.numeric)], 
                 cl=as.numeric(flowSOMcluster[fs$map$mapping[, 1]]))$DB)# Davies–Bouldin index 
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
  
  # plots clustree result
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
                       run.plot = TRUE,
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
  
  # here is the clustering
  
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
      p4 <- DoHeatmap(input, features = rownames(input), group.by = clust_name) # heatmap
      
      # UMAP
      p3 <- DimPlot(input, reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)
      print(p3)
      print(p4)
      
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
                  index.DB(x = df_input[, sapply(df_input, is.numeric)], cl=as.numeric(phenocluster))$DB) # Davies–Bouldin index
  }
  
  if (run.plot && (length(pheno_lou_kn) > 2)) {
    # make clustree plot
    if(length(pheno_lou_kn) > 2) {
      p4 <- clustree(input, prefix ='Pheno.kn.')
      print(p4)
      
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


#function #3 - seurat louvain clustering:
louvain <- function(input, #seu object
                    df_input,
                    pheno_lou_kn,
                    resolutions,
                    run.plot = TRUE, #option to save the graphs  
                    run.stats = TRUE, #option to save stats list  
                    save_to = NULL #only required when save is TRUE
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
    print(paste("finding neighbours for kn ",i,sep=""))
    input <- FindNeighbors(input, dims = 1:12, k.param = i)
    input <- RunUMAP(input, dims = 1:12, n.neighbors = i)
    for (j in resolutions) {
      input <- FindClusters(input, resolution = j)
      print(paste("completed louvain kn ",i, "resolution ",j))
      louvainCluster <- input$seurat_clusters
      print(length(louvainCluster))
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
                      index.DB(x = df_input[, sapply(df_input, is.numeric)], cl = as.numeric(louvainCluster))$DB) # Davies–Bouldin index
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






# clust_stability
# select cluster method and one pararmeter to vary (resolutions 0.1,0.3,0.5,0.8) 
# runs n iterations randomizing starting point - default 100
# calculates the RandIndex (RI) between each iteration of clustering for each resolution
# calculates the mean and standard deviation of the number of clusters and the RI for each resoloution
# outputs a table and plot of the results


clust_stability <- function(input,
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
  print(p)
  return(list(list = rdf, figure = p))
}


# Calculate RAND INDEX
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


# choose cluster conditions to make the seurat object
# method can be "louvain", "phenograph", "flowsom"

get_clusters <- function(seu, method = "louvain",
                         df_input = NULL, #needed  if input "flowsom"
                         k = 60, #k for flowsom or kn for Phenograph and Seurat Louvain
                         resolution = 0.8,
                         plots = TRUE,
                         save_plots = FALSE) {
  if(method == "louvain"){
    print("method is Louvain")
    seu <- FindNeighbors(seu, dims = 1:12, k.param = k)
    # must take one less than the number of antibodies 
    seu <- FindClusters(seu, resolution = resolution)
    group_name <- "seurat_clusters"
  }
  else if(method == "phenograph"){
    # phenograph clustering 
    # add to R object
    print("method is phenograph")
  }
  else if(method == "flowsom"){
    # cluster flowsom 
    print("method is flowsom")
    frame <- new("flowFrame", exprs = as.matrix(df_input[, sapply(df_input, is.numeric)])) #convert input to flowframe, remove non-numerical col
    fs <- ReadInput(frame) #convert flowframe to flowsom object
    fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because X column is removed
    fs <- BuildMST(fs)
    # get the meta data
    flowSOMcluster <- metaClustering_consensus(fs$map$codes, k = k, seed=42)
    
    # add the cluster ID into seurat object to visualize
    seu <- AddMetaData(object = seu,
                       metadata = flowSOMcluster[fs$map$mapping[, 1]],
                       col.name = paste('FlowSom.k.', k, sep="")) 
    group_name <- paste('FlowSom.k.', k, sep="")
  }
  else{ 
    print("select a valid clustering method: 'louvain','phenograph','flowsom' ")
  }
  # make the UMAP for all object
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = k, min.dist = 0.4,
                 spread = 1.5)
  if(plots){
    
    p <- DimPlot(seu, group.by = group_name)
    print(p)
    p1 <- FeaturePlot(seu, features = rownames(seu),
                      slot = 'scale.data',
                      min.cutoff = 'q1',
                      max.cutoff ='99',
                      label.size = 1) +
      theme(plot.title = element_text(size = 0.1))
    
  }
  if(save_plots){
    png(paste(save_plots, "UMAPfeatureplots.png", sep = ""))
    p1
    dev.off()
    png(paste(save_plots, "UMAPclusters.png", sep = ""))
    p
    dev.off()
  }
  return(seu)
}


########## Functions for annotating clusters ##################################


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
  print("Searching for best max node size")
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
  sum.df.2 <- data.frame(results_sum$values)
  max_accuracy <- max(sum.df.2)
  #print(paste("Max accuracy is ",max_accuracy))
  sum.df.2 <- select(sum.df.2, -contains("Kappa"))
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
                         nodesize = start_node,
                         maxnodes = as.numeric(max_node),
                         ntree = ntree)
    key <- toString(ntree)
    store_maxtrees[[key]] <- rf_maxtrees
  }
  # get the best number of trees
  results_tree <- resamples(store_maxtrees)
  results_sum <- (summary(results_tree))
  # make a data frame with the summary results
  sum.df.2 <- data.frame(results_sum$values)
  #print(paste("Max accuracy is ",max_accuracy))
  sum.df.2 <- select(sum.df.2, -contains("Kappa"))
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


plot_lab_clust <- function(seu, seu.cluster, seu.labels, filter_out = c("unknown","Unknown","Mixed")){
  t.lables <- as.data.frame(table(seu.cluster, seu.labels))
  t.lables$Freq <- as.double(t.lables$Freq)
  colnames(t.lables) <- c("Cluster","Label","Freq")
  print(ggplot(t.lables, aes(y = Freq, x = Cluster, fill = Label)) + geom_bar(position = "stack", stat= "identity"))
  # now filter out the unknown and mixed or whatever labels are over powering
  t.lab.known <- t.lables %>% dplyr::filter(!Label %in% filter_out)
  print(ggplot(t.lab.known, aes(y = Freq, x = Cluster, fill = Label)) + 
          geom_bar(position = "stack", stat= "identity"))
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
                           filter_out = 'none', Label = "Label"){
  t.lables <- as.data.frame(table(seu.cluster, seu.label))
  t.lables$Freq <- as.double(t.lables$Freq)
  colnames(t.lables) <- c("Cluster", "Label","Freq")
  top.labs <- t.lables  %>% group_by(Cluster) %>% top_n(top_n, Freq)
  sort.tops <- top.labs %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
  print(sort.tops)
  # now filter out the unknown if desired to get a label for each cluster
  if(length(filter_out) == 0){
    #(all(filter_out %in% c("","FALSE","none"))){
    print("not filtering")
    #top.lab <- t.lables  %>% group_by(Cluster)  %>% top_n(1, Freq)
    top.lab <- t.lables %>%
      group_by(Cluster) %>%
      mutate(rank = dense_rank(desc(Freq))) %>%
      filter(rank == 1) %>%
      select(-rank)
    sort.tops <- top.lab %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
  }
  else{
    print("filtering")
    t.lab.known <- t.lables %>% dplyr::filter(!(Label %in% filter_out))
    top.lab <- t.lab.known  %>% group_by(Cluster) %>% top_n(1, Freq)
    #top.lab <- t.lab.known %>%
    #  group_by(Cluster) %>%
    #  mutate(rank = dense_rank(desc(Freq))) %>%
    # filter(rank == 1) %>%
    #  select(-rank)
    sort.tops <- top.lab %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster) 
  }
  ann.df <- sort.tops %>% select(-"Freq")
  colnames(ann.df) <- c("Cluster", Label)
  return(ann.df)
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
  dfsort <- df3  %>% arrange(Cluster) 
  # get a vector 
  dfcon <- dfsort %>% dplyr::select("Cluster","consensus")
  print(dfcon)
  
  # now add the annotations into the seurat
  # use the simple annnotation function
  seu <- annotate(seu, annotations = dfcon$consensus, to_label,
                  annotation_name)
}




annotate <- function(seu, annotations, to_label, annotation_name = "CellType"){
  Idents(seu) <- to_label
  names(annotations) <- levels(seu)
  seu <- RenameIdents(seu, annotations)
  seu <- AddMetaData(object=seu, metadata=Idents(seu), col.name = annotation_name)
  
}


######### functions to compare groups in an annotated seurat object ################


# plotting function 

proportionplots <- function(seu, seu.var, seu.lable, groups = "Sample", my_colours= "default"){
  sample.lables <- as.data.frame(table(seu.var, seu.lable))
  sample.lables$Freq <- as.double(sample.lables$Freq)
  check <- length(unique(sample.lables$seu.lable))
  print(paste("Number of colours needed",check, sep = ""))
  print(paste("Number of colours entered ", length(my_colours), sep = ""))
  if (length(my_colours) < check){
    # bar chart of with percent 
    print("Default ggplot colours used")
    print(ggplot(sample.lables, aes(x = seu.var,y=Freq ,fill = seu.lable)) + 
            geom_bar(position= "fill", stat = "identity") + 
            scale_y_continuous(labels = scales::percent_format()) + 
            theme_classic() + 
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))  
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }else {
    # bar chart of with percent 
    print("Custome colours used.")
    print(ggplot(sample.lables, aes(x = seu.var,y=Freq ,fill = seu.lable)) + 
            geom_bar(position= "fill", stat = "identity") + 
            scale_y_continuous(labels = scales::percent_format()) + 
            scale_fill_manual(values = my_colours)+
            theme_classic() + 
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))  
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }
  
}

# run plotting function on list of variables 

plotproportions <- function(seu, var.list, xgroup, varnames, my_colours = 'default'){
  for (i in seq_along(var.list)) {
    var <- var.list[[i]]
    proportionplots(seu = seu, seu.var = var, seu.lable = xgroup, groups = varnames[i], my_colours = my_colours)
  }
}

# var.list is a list defining seurat metadata slots
# xgroup is the cell types or x axis grouping variable 
# varnames is a character vector of labels for the x axis of the plots
# Example
# plotproportions(seu.q, var.list = var.list, xgroup = seu.q$cell.types, varnames = varnames)
# can input costume colours
# later --- add in define sample order 

### dotplots and heatmaps mean expression by group

plotmean <- function(plot_type = 'heatmap',seu, group, markers, var_names, slot = 'scale.data',
                     xlab = 'Cell Types', ylab = 'Antibodies', var1order, var2order){
  express.by.cluster <- as.data.frame(AverageExpression(seu, features = markers, 
                                                        group.by = group, slot = 'scale.data'))
  express.by.cluster <- as.data.frame(scale(express.by.cluster))
  if(length(var_names) < 0){
    col.names <- colnames(express.by.cluster) 
  }else{
    col.names <- var_names
  }
  names(express.by.cluster) <- col.names
  AB <- row.names(express.by.cluster)
  ex.data <- cbind(AB,express.by.cluster)
  # must make a data table because melt is 
  dt <- data.table(ex.data)
  dt.long<- melt(dt, id.vars = "AB")
  # select dot plot or heat map
  if(plot_type == 'heatmap'){
    # Heat map
    print(ggplot(dt.long, aes(x = variable, y = AB)) + 
            geom_raster(aes(fill=value)) + 
            scale_fill_gradient(low="blue", high="pink", na.value = "grey") +
            labs(x=xlab, y= ylab) +
            theme_bw() + theme(axis.text.x=element_text(size=12, angle=90, hjust=0.95,vjust=0.2),
                               axis.text.y=element_text(size=12),
                               plot.title=element_text(size=12)))
  } else if(plot_type== 'dotplot'){
    a <- DotPlot(seu, features = AB, group.by = group)
    pct.exp <- as.data.frame(a$data) %>% select(features.plot, id, pct.exp)
    # add the mean expression and percent cells
    # rename columns to match
    colnames(dt.long) <- c("Markers","id","expression")
    colnames(pct.exp) <- c("Markers","id","proportion")
    df.exp.pct <- merge(dt.long, pct.exp, by = c("Markers", "id"))
    data <- df.exp.pct %>% mutate(Cell.type = factor(id, levels = var1order))
    data <- data %>% mutate(Marker = factor(Markers, levels = var2order))
    print(ggplot(data = data, aes(x = Marker, y = Cell.type, color = expression, size = proportion)) +
            geom_point() +
            scale_color_gradient(low = "grey", high = "blue") + ylab("Cell Phenotypes") 
          + xlab("Antibodies") + RotatedAxis())
  }
  else {
    print("not a valid plotting option")
  }
}

# compare groups functions

# prep_for_stats
# run_stats
# make plots



##############################################################################################
# prep_for_stats
# input seurat object
# arguments input all variable that might be compared and where to find these variables
# Selects expression data from Seurat object organized by designated variables

require(data.table)

Prep_for_stats <- function(seu, marker_list, variables, marker_name = 'Marker'){
  # create a dataframe with all the expresson data
  df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
  # rename columns 
  colnames(df) <- marker_list
  # add the different grouping
  for (i in seq_along(variables)) {
    meta_data <- as.data.frame(seu@meta.data[, variables[i]])
    colnames(meta_data) <- variables[i]
    df <- cbind(df, meta_data)
  }
  # need the antibodies to also be a variable
  # convert to data table to use the melt function in data.table instead of reshape2
  dt <- as.data.table(df)
  molten <- melt(dt, measure.vars = marker_list, variable.name = marker_name)
  return(molten)
}




##############################################################################################
# run_stats
# takes in dataframe from prep_for_stats
# arguments, variables, level1, level2
# runs 2way anova and posthoc tests
# outputs tables in a list
# anovas
# posthoc tukey's main 
# tukey's interaction effects filtered matching sets of interest 
# tables  filtered by significance


### for using sample means 


# group_cols is a vector with the columns to get the means from 
# where id1 is the independent variable to compare (dependent variable is the expression)
# id2 is for 2 way anova and should be Marker or Celltype, 
# but other options are possible





## try getting results DF possible list of DF
## new
library(dplyr)
run_stats <- function(input_df, group_cols = c("Sample", "CellType", "Marker"),
                      value_col = "value",
                      stat_type = "ANOVA",
                      id1,
                      id2 = NULL,
                      use_means = TRUE,
                      loop_by = "CellType") {
  aov.l <- list() # for ANOVA output
  tuk.l <- list() # for TUKEY output
  if (use_means) {
    df_means <- get_means(input_df, group_cols, value_col)
    
  } else {
    df_means <- input_df
    names(df_means)[names(df_means) == 'value'] <- 'expression'
  }
  
  if (stat_type == "ANOVA") {
    if (loop_by == "CellType") {
      var_list <- unique(df_means$CellType)
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(CellType == i)
        one_way <- aov(expression ~ df[[id1]], data = df)
        output <- summary(one_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(one_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(Celltype = names(aov.l)[i],
                   IndependentVariable = id1,
                   Fvalue = aov.l[[i]][[1]][["F value"]][1],
                   pValue = aov.l[[i]][[1]][["Pr(>F)"]][1],
                   Df = aov.l[[i]][[1]][["Df"]][1],
                   SumSqEffect = aov.l[[i]][[1]][["Sum Sq"]][1],
                   MeanSqEffect = aov.l[[i]][[1]][["Mean Sq"]][1],
                   SumSqError = aov.l[[i]][[1]][["Sum Sq"]][2],
                   MeanSqError = aov.l[[i]][[1]][["Mean Sq"]][2]
        )
      }
      ))  
      print(paste("ANOVA results for each cell type comparing ", id1))
      tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
        rownames <- row.names(tuk.l[[i]][[1]])
        data.frame(Celltype = names(tuk.l)[i],
                   IndependentVariable = id1,
                   Contrast = rownames,
                   tuk.l[[i]][[1]], row.names = NULL)
      }))
      print(paste("TukeyHSD results for each cell type comparing ", id1)) 
    } else if (loop_by == "Marker") {
      var_list <- unique(df_means$Marker)
      for (i in var_list) {
        df <- df_means %>% filter(Marker == i)
        one_way <- aov(expression ~ df[[id1]], data = df)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(one_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(Marker = names(aov.l)[i],
                   IndependentVariable = id1,
                   Fvalue = aov.l[[i]][[1]][["F value"]][1],
                   pValue = aov.l[[i]][[1]][["Pr(>F)"]][1],
                   Df = aov.l[[i]][[1]][["Df"]][1],
                   SumSqEffect = aov.l[[i]][[1]][["Sum Sq"]][1],
                   MeanSqEffect = aov.l[[i]][[1]][["Mean Sq"]][1],
                   SumSqError = aov.l[[i]][[1]][["Sum Sq"]][2],
                   MeanSqError = aov.l[[i]][[1]][["Mean Sq"]][2]
        )
      }
      )) 
      print(paste("ANOVA results for each marker comparing ", id1))
      tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
        rownames <- row.names(tuk.l[[i]][[1]])
        data.frame(Marker = names(tuk.l)[i],
                   IndependentVariable = id1,
                   Contrast = rownames,
                   tuk.l[[i]][[1]], row.names = NULL)}))
      print(paste("Tukey results for each marker comparing ", id1))
    } else {
      one_way <- aov(expression ~ df_means[[id1]], data = df_means)
      aov.l <- summary(one_way)
      # now the posthoc test
      tukey <- TukeyHSD(one_way)
      tuk.l <- tukey
      aov_df <- data.frame(
        IndependentVariable = id1,
        Fvalue = aov.l[[1]][["F value"]][1],
        pValue = aov.l[[1]][["Pr(>F)"]][1],
        Df = aov.l[[1]][["Df"]][1],
        SumSqEffect = aov.l[[1]][["Sum Sq"]][1],
        MeanSqEffect = aov.l[[1]][["Mean Sq"]][1],
        SumSqError = aov.l[[1]][["Sum Sq"]][2],
        MeanSqError = aov.l[[1]][["Mean Sq"]][2]
      )
      print("ANOVA results for all celltypes and markers combined") 
      rownames <- row.names(tuk.l[[1]])
      tuk_df  <- data.frame(
        IndependentVariable = id1,
        Contrast = rownames,
        tuk.l[[1]], row.names = NULL)
      print("TukeyHSD results for all celltypes and markers combined")
    }
    output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_df)
  }else if (stat_type == "ANOVA2") {
    if (loop_by == "CellType") {
      var_list <- unique(df_means$CellType)
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(CellType == i)
        formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
        two_way <- aov(formula, data = df)
        output <- summary(two_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(two_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      # make a dataframe of the 2 way anova results
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(
          rownames <- row.names(aov.l[[i]][[1]]),
          data.frame(Celltype = names(aov.l)[i],
                     Contrast = rownames,
                     aov.l[[i]][[1]], row.names = NULL
                     
          )
        )
      }))
      # get ride of the residual which will confuse things and the extra row added
      aov_df <- aov_df[2:8]
      # fix the column names
      colnames(aov_df) <- c("Celltype","Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      print(paste("ANOVA 2way results for each cell type comparing ", id1,
                  "and ", id2))
      # make dataframes from each of the tukey results and put them into a list
      tuk_summary.l <- list()
      
      for(J in 1:length(tuk.l[[1]])){
        tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
          rownames <- row.names(tuk.l[[i]][[J]])
          data.frame(Celltype = names(tuk.l)[i],
                     Contrast = names(tuk.l[[i]])[J],
                     Subgroups = rownames,
                     tuk.l[[i]][[J]], row.names = NULL)
        }))
        # put the dataframe in a list 
        tuk_summary.l[[as.character(names(tuk.l[[i]])[J])]] <- tuk_df
      }
      # filter the interaction dataframe to have only df1 contrasts
      dft <- tuk_summary.l[[3]]
      filtered_df <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>% 
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe for only id2 pairs 
      tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
      # filter the interaction dataframe
      filtered_df2 <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>% 
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe to have matching id1 or id2 contrasts
      tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      
      print(paste("TukeyHSD results for each cell type comparing ", id1,
                  "and ", id2)) 
    } else if (loop_by == "Marker") {
      var_list <- unique(df_means$Marker)
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(Marker == i)
        formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
        two_way <- aov(formula, data = df)
        output <- summary(two_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(two_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      # make a dataframe of the 2 way anova results
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(
          rownames <- row.names(aov.l[[i]][[1]]),
          data.frame(Marker = names(aov.l)[i],
                     Contrast = rownames,
                     aov.l[[i]][[1]], row.names = NULL
                     
          )
        )
      }))
      # get ride of the residual which will confuse things and the extra row added
      aov_df <- aov_df[2:8]
      # fix the column names
      colnames(aov_df) <- c("Marker","Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      print(paste("ANOVA 2way results for each Marker comparing ", id1,
                  "and ", id2))
      # make dataframes from each of the tukey results and put them into a list
      tuk_summary.l <- list()
      
      for(J in 1:length(tuk.l[[1]])){
        tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
          rownames <- row.names(tuk.l[[i]][[J]])
          data.frame(Marker = names(tuk.l)[i],
                     Contrast = names(tuk.l[[i]])[J],
                     Subgroups = rownames,
                     tuk.l[[i]][[J]], row.names = NULL)
        }))
        # put the dataframe in a list 
        tuk_summary.l[[as.character(names(tuk.l[[i]])[J])]] <- tuk_df
      }
      # filter the interaction dataframe to have only df1 contrasts
      dft <- tuk_summary.l[[3]]
      filtered_df <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>% 
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe for only id2 pairs 
      tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
      # filter the interaction dataframe
      filtered_df2 <- dft %>%
        filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>% 
                 apply(2, function(x) all(x == x[1])))
      # add the filtered dataframe to have matching id1 or id2 contrasts
      tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      
      print(paste("TukeyHSD results for each Marker comparing ", id1,
                  "and ", id2)) 
    } else {
      df <- df_means
      formula <- as.formula(paste0("expression ~ ", id1, "*", id2))
      two_way <- aov(formula, data = df)
      aov.l <- summary(two_way)
      # make the df of the ANOVA two way results
      rownames <- row.names(aov.l[[1]])
      aov_df <- data.frame(
        Contrast = rownames,
        aov.l[[1]], row.names = NULL
        
      )
      # get ride of the residual which will confuse things and the extra row added
      # fix the column names
      colnames(aov_df) <- c("Contrast","Df","SumSq","MeanSq",
                            "Fvalue","Pvalue")
      # remove extra spaces
      aov_df <- aov_df %>% mutate_all(trimws)
      # remove residuals
      aov_df <- aov_df[aov_df$Contrast != "Residuals", ]
      
      # now the posthoc test
      tukey <- TukeyHSD(two_way)
      tuk.l <- tukey
      print("ANOVA results for all celltypes and markers combined") 
      # get the Tukey outputs
      for(J in 1:length(tuk.l)){
        rownames <- row.names(tuk.l[[J]])
        tuk_df <- data.frame(
          Contrast = names(tuk.l[J]),
          Subgroups = rownames,
          tuk.l[[J]], row.names = NULL)
        # put the dataframe in a list 
        tuk_summary.l[[as.character(names(tuk.l)[J])]] <- tuk_df
        dft <- tuk_summary.l[[3]]
        filtered_df <- dft %>%
          filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(2, 4)) %>% 
                   apply(2, function(x) all(x == x[1])))
        # add the filtered dataframe for only id2 pairs 
        tuk_summary.l[[paste("Interactions_",id1)]] <- filtered_df
        # filter the interaction dataframe
        filtered_df2 <- dft %>%
          filter(sapply(strsplit(Subgroups, "[:-]"), "[", c(1, 3)) %>% 
                   apply(2, function(x) all(x == x[1])))
        # add the filtered dataframe to have matching id1 or id2 contrasts
        tuk_summary.l[[paste("Interactions_",id2)]] <- filtered_df2
      }
      
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)
    }
    return(output_list)
  }
}
  
  # 
  
  # original get means function
  
  get_means <- function(df, group_cols, value_col) {
    df_means <- df %>%
      group_by(across(all_of(group_cols))) %>%
      mutate(expression = mean(value)) %>%
      distinct(across(all_of(group_cols)), expression, .keep_all = TRUE) %>%
      select(-value)
    return(df_means)
  }
  

