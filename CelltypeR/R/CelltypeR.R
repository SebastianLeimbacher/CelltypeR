# CellTypeR functions

# note: warnings occur for packages maptools rgdal and rgeos
# rgeos is phenograph, ggplot2 dependency
# maptools is a ggplot2 dependency


####### Preprocessing functions to create seurat object from Flow cytometry multichannel data ####

# reads in fsc files creates a flowset file
# filtering if wanted on each input file
# users can rename samples after


#' Read in fsc (Flow cytometry) files and creates a flowset file with Channel names
#'
#' This function to take in a folder/file path of fsc flow cytometer file selecting the FSC-A
#' value for each channel.  Each sample file in the folder will be a slot in the flowset object.
#' The samples will be names by the file names and can be renamed. The channel names will automatically
#' be added from within the flow data. If only a subset of cells is desired downsampling is an option.

#' @export
#' @examples
#' fsc_to_fs(input_folder_fsc, downsample = "none")
#' @import flowCore
#' @importFrom flowCore read.flowSet fsApply

fsc_to_fs <- function(input, downsample = 'none'){
  flowset = read.flowSet(path=input_path,transformation = FALSE ,
                         emptyValue = FALSE,truncate_max_range = FALSE,
                         package="flowCore")
  fsApply(flowset,
          function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(flowset))]})
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){gsub(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs) <- unlist(lapply(marker.names, function(x){sapply(strsplit(x,"_"),head,1)}))
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
# transform, align, retro transform flowset object

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


#' Transform and align samples in flowset object
#'
#' This function to takes in the flowset object created from a fsc files and transformed the data
#' If 'retro' (default) is selected the samples will be biexponentially transformed, aligned and then reverse transformed
#' If 'biexp' is selected the samples will only be biexponentially transformed
#' If 'align' is selected the samples will be transformed and aligned but not reverse transformed
#' The two_peaks and one_peak argument refer to the biexp transformed data and are the indexes of which measures have
#' one or two peaks. At least one measurement of two peaks must be included. You will receive an error if there is only one peak in a measurement defined as having two peaks.

#' @export
#' @examples
#' harmonize(flowset, processing = 'retro', two_peaks = c(10:20), one_peak = c(1:9), theshold = 0.01)
#'
#' @import flowCore
#' @importFrom flowCore read.flowSet
#' @importFrom flowStats gaussNorm


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
#' Create a density plot from a flowset object
#'
#' Takes in a flowset object. Plot density plots for each channel in multichannel fsc
#' Flow Cytometry data. Creates stacked plots where each row is a sample in the
#' set object.

#' @export
#' @import ggridges
#' @importFrom ggridges geom_density_ridges

plotdensity_flowset <- function(flowset, nsample_include = "all") {
  # Create the data frame
  df <- melt(lapply(as.list(flowset@frames), function(x) { x = as.data.frame(x@exprs) }))

  # Count the number of unique samples
  nsamples <- length(unique(df$L1))
  if(nsample_include != "all") {
  # Check if there are more than the number to include samples
  if(nsamples > nsample_include) {
    # Select the first 9 unique samples
    samples_to_plot <- df %>%
      distinct(L1) %>%
      head(nsample_include)
    # Filter the data frame to include only the selected samples
    df <- df %>%
      filter(L1 %in% samples_to_plot$L1)
  }
  }
  # make the plot
  ggplot(df, aes(x = value, y = L1)) +
    geom_density_ridges(aes(fill = L1), alpha = 0.4, verbose = FALSE) +
    facet_wrap(~variable) +
    theme_light() +
    guides(fill = FALSE)+
    theme(axis.text.y = element_blank())  # Remove y-axis labels
}



###### function called by fsc_to_df

rename_markers <- function(flowset) {
  copy_flowset <- flowset[seq(along = flowset)]
  for (i in 1:length(copy_flowset)) {
    marker.names <- copy_flowset[[i]]@parameters@data$desc
    marker.names <- gsub("-", "_", marker.names, fixed = TRUE)
    colnames(copy_flowset[[i]]@exprs) <- unlist(lapply(marker.names, function(x) sapply(strsplit(x, "_"), head, 1)))
  }
  return(copy_flowset)
}


#' Create a dataframe from a flowset object and option to save csv
#'
#' Takes in a flowset object and creates a dataframe. By default no output path is added and a data frame is returned.
#' The name of each sample will be added as the column "Sample"
#' into the final dataframe. To save the dataframe directly as a csv set save.csv = filepath/tosave/
#' The filename will be the name of the flowset object input into the function.
#'
#' @export
#' @examples
#' flowset_to_csv(flowset)
#' flowset_to_csv(flowset, output_path = "path/to/location/", save.csv = TRUE)
#' @importFrom flowCore fsApply
flowset_to_csv=function(flowset, output_path, save.csv = FALSE){
  list_of_flowframes=fsApply(rename_markers(flowset),function(x){as.data.frame(x@exprs)})#Makes a list of dataframes
  list_names=names(list_of_flowframes)#extract names for not calling names() function at each loop
  for (index in seq_along(list_of_flowframes)){ #Iterates along the index for adding sample names
    list_of_flowframes[[index]]$Sample = list_names[index]
    colnames(list_of_flowframes[[index]]) = colnames(list_of_flowframes[[1]])
  }
  ds=combine_flowframes(list_of_flowframes)#binds every fcs file in a single dataframe
  ds$cell=as.factor(unlist(lapply(as.list(c(1:length(flowset))),function(x){c(1:nrow(flowset[[x]]@exprs))})))#add cell IDs - cell count per sample
  if (save.csv){
    write.csv(ds,file=paste0(output_path,deparse(substitute(flowset)),".csv"))#save the R data for further usage
  }
    return(ds)
  }


# alternate to list.rbind
combine_flowframes <- function(list_of_flowframes) {
  # Use do.call and rbind to combine the flow frames
  combined_flowframe <- do.call(rbind, list_of_flowframes)

  return(combined_flowframe)
}

##############################################################################################
# df_to_seurat
# creates a seurat object from the expression matrix and adds in meta data

#' Creates a seurat object from the Flow Cytometry expression data frame.
#'
#' Takes in a dataframe created from Flow Cytometry data using the fsc_to_fs and fsc_to_csv
#' functions. This function creates a Seurat object from the Flow Cytometry expression data
#' in a data frame format where cells are rows and Markers are columns.  A column indicating
#' the starting fsc file "Sample" is required and added as meta data into the Seurat object.

#' @export
#' @examples
#' make_seu(df = flow_dataframe, AB_vector = markers_names)
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @importFrom dplyr select

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

#' Explore different clustering methods and parameters in single cell data
#'
#' Takes in a Seurat object. Choose a clustering method or vector of clustering methods:
#' cluster_method = c("louvain","flowsom","phenograph"). The flowsom method requires the
#' original data frame of the Flow Cytometry expression data. For flows k is the number of
#' final clusters. To explore different cluster numbers flow_k = c(1,5,10,15). For louvain and
#' phenograph (phenograph uses Louvain clustering with a Jaccard index) the kn value is the
#' starting value for the snn nearest neighbour network detection: pheno_lou_kn = c(20, 60, 200).
#' For louvain network detection run with the Seurat function a vector of resolutions can be
#' input: lou_resolution = (0,0.2,0.5,1.2). If run.plot = TRUE the results will be shown in
#' using the clustree function.  Putting a starting resolution of 0 will make the visualization
#' better. Setting run.plot = TRUE will also produce UMAPs and heatmaps for each cluster
#' condition. If the option run.stats = TRUE internal statistics will be calculated. Data tables will be
#' produced in a list. If you wish to save the plots generated and the list of statistic set
#' save_to = filepath/tosaveto/ . To run intrinsic statistics. Silhouette score ranges from -1 to 1.
#' A value near -1 indicates a poor quality of the clusters, a value near 0 indicates the neutrality and
#' insignificance of the clusters, and a value near 1 indicates a good quality of the clusters
#' (Rousseeuw, 1987). Calinski-Harabasz index does not have a range, but a higher index indicates a
#' better cluster. The number of clusters with the highest Calinski-Harabasz index is the optimal number
#' of clusters (Calinski & Harabasz, 1974). Davies-Bouldin index has no upper bound but the minimum index
#' is 0. A lower index indicates better clusters (Davies & Bouldin, 1979).
#'
#' @examples
#' explore_param(seuratObj, cluster_method = "flowsom",df_input = flowset_dataframe, flow_k = 10)
#' explore_param(seuratObj, cluster_method = "louvain",df_input = flowset_dataframe, flow_k = NULL,
#' pheno_lou_kn = 100, lou_resolution = c(0,0.5,1), pcdim = 1:10)
#' @export
#' @import Seurat dplyr FlowSOM cluster fpc clusterSim clustree
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @importFrom dplyr select
#' @importFrom FlowSOM BuildSOM BuildMST ReadInput
#' @importFrom cluster daisy
#' @importFrom fpc calinhara
#' @importFrom clusterSim index.DB
#' @importFrom clustree clustree
#' @importFrom igraph membership

# input is a seurat object
explore_param <- function(input, #for phenograph and louvain only
                          cluster_method, #take 1 cluster method
                          df_input, #required for stats
                          flow_k = NULL, #k for flowsom
                          pheno_lou_kn = NULL, #kn for phenograph or louvain
                          lou_resolution = NULL, #resolution for louvain
                          pcdim = 1:10,
                          run.plot = TRUE, #print if run
                          run.stats = TRUE, #print and return if run
                          save_to = NULL) { #need it if save plots or stats


  #call helper functions depending on the cluster_method requested
  if (cluster_method == "louvain") {
    cl <- louvain(input = input, #seu object
                  df_input = df_input,
                  pheno_lou_kn = pheno_lou_kn,
                  resolutions = lou_resolution,
                  ldim = pcdim,
                  run.plot = run.plot,
                  run.stats = run.stats,
                  save_to = save_to)
  }
  if (cluster_method == "phenograph") {
    cl <- phenograph(input = input,
                     df_input = df_input,
                     pheno_lou_kn = pheno_lou_kn,
                     pdim = pcdim,
                     run.stats = run.stats,
                     run.plot = run.plot,
                     save_to = save_to)
  }
  if (cluster_method == "flowsom") {
    cl <- flowsom(input = input,
                  df_input = df_input,
                  flow_k = flow_k,
                  fdim = pcdim,
                  run.stats = run.stats,
                  run.plot = run.plot,
                  save_to = save_to)
  }

  if (run.plot || run.stats) {
    return(list(cluster = cl))
  }
}

# function #1 flowsom run but the main function above:
flowsom <- function(input, #seurat
                    df_input, #the processed df2 file before being converted to seurat
                    flow_k,
                    fdim = pcdim,
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
  #input <- FindNeighbors(input, dims = dim)
  input <- RunUMAP(input, dims = fdim)

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
                    group.by = clust_name, raster = FALSE) # will automatically group by active ident
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
    if(!is.null(save_to)){
      write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
    }

  }

  # plots clustree result
  if (run.plot) {
    if(length(flow_k) > 2) {
      # make clustree plot
      p5 <- clustree::clustree(input, prefix ='FlowSom.k.')
      print(p5)
      print("Clustree")

      if (!is.null(save_to)) {
        png(filename=paste(save_to, clust_method, 'Clustree.png',sep = ""), width = 1200, height = 2000, res = 300)
        print(p5)
        dev.off()
      }
    }
  }

  if (run.stats) {
    return(stats_df)
  }
}


#helper functions #2:
phenograph <- function(input,
                       df_input,
                       pheno_lou_kn,
                       pdim = pcdim,
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
    input <- FindNeighbors(input, dims = pdim, k.param = i, reduction = "pca")
    input <- RunUMAP(input, dims = pdim, n.neighbors = i)

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
      p3 <- DimPlot(input, reduction = "umap", repel = TRUE,
                    label = TRUE, group.by = clust_name, raster = FALSE)
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
    if(!is.null(save_to)){
      write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
    }
  }

  if (run.stats) {
    return(stats_df)
  }

}


#function #3 - seurat louvain clustering:
louvain <- function(input, #seu object
                    df_input,
                    pheno_lou_kn,
                    resolutions,
                    ldim = pcdim,
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

    input <- FindNeighbors(input, dims = ldim, k.param = i, reduction = "pca")
    input <- RunUMAP(input, dims = ldim, n.neighbors = i, reduction = "pca")
    print(paste("finding neighbours for kn ",i,sep=""))
    for (j in resolutions) {
      DefaultAssay(input) <- "RNA"
      input <- FindClusters(input, resolution = j, reduction = "pca")
      print(paste("completed louvain kn ",i, "resolution ",j))
      louvainCluster <- input$seurat_clusters
      print(length(louvainCluster))
      #stats
      if (run.stats) {
        statsl <- c(statsl,
                    i, # kn
                    j, # resolution
                    length(unique(louvainCluster))) # number of clusters (nc)

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
        p4 <- DoHeatmap(input, features = rownames(input), size = 10) +
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
    if (!is.null(save_to)) {
      write.csv(stats_df, paste(save_to, clust_method, '_stats.csv', sep = ""), row.names = F)
    }
  }

  if (run.stats) {
    return(stats_df)
  }
}


# clust_stability


#' Test the stability or reproducibility of clusters across kn or resolutions
#'
#' For Louvain clustering in a seurat object test calculate the RAND Index. For n iterations
#' each cluser identity for each cell is compared with the other iterations of clustering.
#' Select one pararmeter to vary resolutions = c(0.1,0.3,0.5,0.8) or kn = c(20,40,60,80).
#' The function also calculates the mean and standard deviation of the number of clusters
#' for each resolution or kn. The function returns an outputs and plot of the results.
#' @export
#' @importFrom Seurat FindNeighbors FindClusters RunUMAP
#' @importFrom flexclust randIndex

clust_stability <- function(input,
                            resolutions,
                            kn,
                            pcdim = 1:10,
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
    input <- FindNeighbors(input, dims = pcdim, k.param = i)
    input <- RunUMAP(input, dims = pcdim, n.neighbors = i)

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


#' Plot cluster stability results from output data
#'
#' Requires the rdf object output from the function clust_stability. The colours and
#' viewable area of the plot can be altered. cp takes a vector of two colours. view takes
#' vector of for the x axis of two coordinates.

#' @export
#' @import ggplot2


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

#' Get clusters and add to seurat object
#'
#' This function to in a seurat object and performs the desired method of clustering.
#' ("louvain","phenograph","flowsom"). With the already optomized conditions. The cluster
#' indexes are added to the seurat object. A heat map of the expression of the a Marker
#' list by cluster features plot of a UMAP for each marker are plotted plots = TRUE.
#' Plots will be saved if an input filepath is given.

#' @export
#' @import Seurat ggplot2 FlowSOM Rphenograph
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom ggplot2 theme
#' @importFrom FlowSOM ReadInput
get_clusters <- function(seu, method = "louvain",
                         df_input = NULL, #needed  if input "flowsom"
                         k = 60, #k for flowsom or kn for Phenograph and Seurat Louvain
                         resolution = 0.8,
                         pcdim = 1:10,
                         plots = TRUE,
                         save_plots = FALSE) {
  # make the UMAP for all object
  seu <- RunUMAP(seu, dims = pcdim, n.neighbors = k, min.dist = 0.4,
                 spread = 1.5)
  print("UMAP run")
  if(method == "louvain"){
    seu <- FindNeighbors(seu, dims = pcdim, k.param = k, reduction = "pca")
    # must take one less than the number of antibodies
    seu <- FindClusters(seu, resolution = resolution)
    group_name <- "seurat_clusters"
    print("method is Louvain")
  }
  else if(method == "phenograph"){
    # phenograph clustering
    m <- t(as.matrix(GetAssayData(object = input, slot = "counts")))
    Rphenograph_out_flow <- Rphenograph(m, k = k)
    phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
    clust_name = paste('Pheno.kn.',k,sep="")
    # add the cluster ID into seurat object
    seu <- AddMetaData(object = seu, phenocluster, col.name = clust_name)
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
  if(plots){

    p <- DimPlot(seu, group.by = group_name)
    print(p)
    p1 <- FeaturePlot(seu, features = rownames(seu),
                      slot = 'scale.data',
                      min.cutoff = 'q1',
                      max.cutoff ='99',
                      label.size = 1)
    print(p1)
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

# input correlation matrix
# input seurat object
#Correlation assignment method,
# predicts cell types based on each cells correlation to the matrix types
# creates plots and tables of the prediction outputs.
# takes argument for "unassigned" threshold (default 0.4)  and "double-label" thresholod(0.05)


#' CAM predict cell types in Flow Cytometry data by correlation to an existing reference matrix
#'
#' The function requires an reference matix already created with the same markers
#' (antibodies) as used in the Flow Cytometry expression data frame. Correlation of
#' expression across markers for each cell type in the reference matrix and each cell in the
#' sample data is calculated. A cutoff for the R value (pearson correlation coefficient) is set,
#' where an R value must be above this threshold for a cell type to be assigned.  Cells with only
#' correlations below the threshold are set as "unknown". The threshold is defined by "min_corr".
#' The best and second best R value is recorded for each cell. If the difference between R max 1
#' and R max 2 is less than or equal to "min_diff" that cell is assigned as a joint cell type. For
#' example "neuron-npc". Test is a data frame of the expression for each marker and cell. Reference
#' is the reference matrix, cell type by marker.

#' @export
#' @importFrom dplyr filter
#' @importFrom kit topn
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
    # topn is from the package kit
    top <- topn(corr_ls, 2) #return the index of the best 2
    result <- c(result,
                test[i, 'X'], #col 1: cell sample
                corr_ls[top[1]], #col 2: 1st correlation
                reference[top[1], 'X'], #col 3: 1st best cell type
                corr_ls[top[2]], #col 4: 2nd correlation
                reference[top[2], 'X'], #col 5: 2nd best cell type
                ifelse(corr_ls[top[1]] < min_corr, "unassigned",
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


#' Plot correlation assignment model (CAM) results
#'
#' The function requires the data frame output from the find_correlation function
#' Filter out cell type annotation with fewer than 'min_cells' for the plot
#' The threshold with draw a line on the plot.  This can be set to the R theshold used in the find_correlation function.

#' @export
#' @examples
#' plot_corr(cor_df, threshold = 0.3, min_cells = 300)
#' @importFrom dplyr filter group_by select
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt

plot_corr <- function(df, threshold = 0, min_cells = 100) {
  # filter to get frequency table and save as csv
  df.f <- df %>% select(cell.label)
  freq.table <- as.data.frame(table(df.f))
  df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> min_cells)
  # plot the frequencies
  plot1 <- ggplot(df.filter, aes(x = reorder(
    cell.label, cell.label, function(x) -length(x)),
    fill = cell.label)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust=0.99,vjust=0.5),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black")) +

    scale_y_continuous(expand = c(0, 0)) + # take the space away below the bars
    xlab('Assigned cell type') +
    ylab('number of cell') +
    labs(fill = 'Cell Types')
  # print(plot1)

  # plot without the unassigned
  df.filter2 <- df %>%
    dplyr::filter(!cell.label %in% c("Unassigned", "unassigned")) %>%
    group_by(cell.label) %>%
    dplyr::filter(n() > min_cells)


  plot1b <-
    ggplot(df.filter2, aes(x = reorder(
      cell.label, cell.label, function(x) - length(x)),
      fill = cell.label)) + geom_bar() + theme_classic() +
    theme(axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust=0.99,vjust=0.5),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black")) +

    scale_y_continuous(expand = c(0, 0)) + # take the space away below the bars
    xlab('Assigned cell type') +
    ylab('number of cell') +
    labs(fill='Cell Types')

  # violin plot of best correlation/cell type with threshold shown or not

  plot2 <- ggplot(df.filter, aes(x = best.cell.type, y = cor.1, fill = best.cell.type)) +
    geom_violin(trim = FALSE) +
    ylim(0, 1) +
    theme_classic() +
    theme(text = element_text(size = 14), axis.text.x = element_text(angle = 90, size = 12)) +
    ylab("correlation coefficient") +
    xlab("Cell type with max correlation coefficient") +
    geom_hline(yintercept = threshold) +
    guides(fill = guide_legend(title = "Cell Type"))

  df.melt <- melt(df) #reformat to long df
  # melt from reshape2 will depreciate
  #df.melt <- pivot_longer(cor1, cols = "cor.1" )

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
    scale_fill_manual(values = c("#4E84C4", "#52854C","purple","orange")) +
    ylab("correlation coefficient") + xlab("Cell type")
  # the second best correlation is so low it was removed from STEM with axis limit -0.1 and even -1
  # down sample
  set.seed(64)
  n_samples <- min(1000, nrow(df))
  df.downsample <- sample_n(df, n_samples)
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
    scale_color_manual(values = c("#4E84C4", "#52854C","purple","orange")) +
    ylim(-0.25, 0.95) +
    facet_wrap(~(as.factor(best.cell.type))) +
    theme(legend.position = "none") +
    ylab("Correlation Coefficient") +
    xlab("")
  # print(plot5)

  double.cells <- df[grep("-", df$cell.label),]
  df.melt.double <- melt(double.cells)

  # # this will be an excellent visualization but I need to subset only the double labels,
  #then I can plot more cells and see more clearly.
  if (length(unique(df.melt.double$cell.label)) > 1) {
    plot6 <- ggplot(df.melt.double, aes(x = variable, y = value, colour = variable, group = X)) +
      geom_line(show.legend = FALSE, size = 0.1, color = "black") +
      geom_point() +
      scale_color_manual(values = c("#4E84C4", "#52854C","purple","orange")) +
      #ylim(-0.15, 0.8) +
      facet_wrap(~ as.factor(cell.label)) +
      ylab("Correlation Coefficient") +
      xlab("")
  } else {
    # If there are no valid values for cell.label, create an empty plot
    plot6 <- ggplot() +
      theme_void()
    print("length of double cells is less than 1")
  }

  return(list(freq.table, plot1, plot1b, plot2, plot3, plot4, plot5,plot6))
}


##############################################################################################

#' Train a Random Forest Classifier to predict cell types in a Seurat object
#'
#'Requires a Seurat object for the labels to be predicted. AB_list is the list of markers(antibodies)
#'or selected features(genes) to input into the training model. An annotated Seurat object of matching
#'data is required. The location of the annotations metadata needs to be designated:
#'annotated_seurat_object$MetadataAnnotations_slot. If the number of cells is high the function may not run.  If this is the case
#'use downsample to reduce the number of cells.  The starting data will be randomly selected from to
#'reduce the total size. This data will then be split in 80/20 training to test data.
#'A random seed can be set to repeat the training with different random starts and random splitting of
#'data. The tuning parameters for the random forst model can be optomized. The best values are
#'selected for the final model from mytry, maxnodes, trees
#'mtry: number of features to draw at once.
#'maxnodes: the terminal number of nodes trees in the forest can have.
#'num_folds is for cross validation to avoid overfitting. The validation data is split from within the training data.
#'Markers is the vector of the markers to select features from. If a RFM is need with only some of
#'the markers in a reference data set selecting only those markers can be done with this argument.
# after the best conditions are selected using the parameter search and cross validation these are used on the full training data.
# cores_to_use is for parallel processing. Make sure not to exceed the number of cores you have available.

#' @export
#' @import data.table randomForest caret
#' @importFrom randomForest randomForest
#' @importFrom caret trainControl train
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @import foreach
RFM_train <- function(seurate_object,
                       markers,
                       annotations = seu$CellTypes,
                       num_folds = 3,
                       downsample = 'none',
                       seed = 222,
                       mytry = c(6:9),
                       maxnodes = c(12:15),
                       trees = c(1000, 1500),
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

  # Split data into training (80%) and test (20%) sets
  sample_indices <- sample(nrow(df1), size = round(0.8 * nrow(df1)))
  train <- df1[sample_indices, ]
  test <- df1[-sample_indices, ]

  trControl <- trainControl(method = "cv", number = num_folds, search = "grid")

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


##############################################################################################


#' Predict cell types in a Seurat object with a trained Random Forest Model
#'
#'Predict cell types in a Seurat object using a trained Random Forest Model.
#'A trained model is required. Returns the Seurat object with the predictions
#'added as the metadata slot 'rfm.labels'.

#' @export
#' @importFrom stats predict
RFM_predict <- function(seu, rf){
  # prepare a data object to be test data
  df <- t(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
  # run the predictions
  rfm.pred <- as.data.frame(predict(rf,df))
  colnames(rfm.pred) <- "Prediction"
  return(rfm.pred)
}


##############################################################################################
# seurat_transfer
# takes in a reference seurat object
# follows seurat workflow, finds anchors, predicts labels
# default of no threshold, a threshold can be applied
# outputs prediction vector

#' Predict cell types in a Seurat object with Seurat label transfer
#'
#'Predict cell types in a Seurat object using a reference Seurat object.
#'An annotated Seurat object is required. The list of marker or features is required.
#'Predictions of cell type for each cell are made using anchors. To increase speed use down.sample
#'to down sample the reference by the cell types. Markers is a list of the flow cytometery antibodies,
#'a vector of selected features(genes) is possible. The predictions will be added to the
#'Seurat object as 'seu.pred'.
#' @export
#' @importFrom Seurat AddMetaData FindTransferAnchors TransferData
seurat_predict <- function(seu.q, seu.r, ref_id = 'labels', seu_meta = seu.r$labels,
                           down.sample = 500, markers, kw = 50, kfilter = 100, dims =10){
  Idents(seu.r) <- ref_id
  seu.r <- subset(seu.r, downsample = down.sample)
  # find anchors
  anchors <- FindTransferAnchors(reference = seu.r,
                                 query = seu.q, features = markers,
                                 reference.reduction = "pca",
                                 dim= 1:dims,
                                 npcs = dims,
                                 k.filter = kfilter,
                                 max.features = length(markers), verbose = FALSE)
  n_anchors <- dim(as.data.frame(anchors@anchors))[1]
  # Ensure kw is always an even number
  if(n_anchors < kw) {
    kw <- n_anchors - 1
    if(kw %% 2 == 0){
      kw <- kw -1
    }
  }
  predictions <- TransferData(anchorset = anchors, refdata = seu_meta,
                              dims = 1:dims, k.weight = kw, query = seu.q)
  seu.q <- AddMetaData(seu.q, metadata = predictions$predicted.id,
                       col.name = 'seu.pred')
}



##############################################################################################

# plot output of predictions

# takes in a seurat object with the labels added
# makes a dataframe with the count of predicted labels for each cluster
# input seurat object with the predicted labels in the meta data
# input the clusters meta data slot to be labels
# input the meta data slot with the labels (correlation, random forest, seurat predicted)


#' Creates a UMAP and barchart to visualize cell type predictions
#'
#'Takes a Seurat object with predicted cell type meta data slot(s) and cluster. Creates
#'a UMAP with the predictions labelled and a bar chart with the proportion of predicted
#'cell types in each cluster. There is an option to ignore the selected cell types using
#'filter_out. The followin inputs are required:
#'seu = Seurat_object
#'seu.cluster = seu$RNA_snn_res.0.5 (the metadata slot with the clusters you want to annotate)
#'seu.labels = seu$seu.pred (the metadata slot with the prediciton or labels to visualize)
#'filter_out = c("unknown") enter the annotation to filter out or set
#'filter_out = "none" to leave all annotations.
#' @export
#' @import dplyr ggplot2
#' @importFrom dplyr filter
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


#' Predict the cell type annotation for each cluster
#'
#'This function creates a data frame with the predicted cell type for each cluster. The
#'Seurat object with clusters and predicted labels must be input. The cluster meta data slot
#'to label and the predicted labels slot must both be input. Sometime many cells are predicted
#'as unknown or another unwanted type. The option filter_out = c("unknown","mixed") can be used to
#'ignore these cell types if they are the most frequent cell type in a given cluster and annotate
#'using the second most frequent cell type.
#'seu.cluster = seu$RNA_snn_res.0.5
#'seu.label = seu$seu.pred or seu$rfm.pred
#'top_n = 2 (this is the number of cell types to visualize)
#'Only the most frequent cell type for each cluster is used for the cluster annotation.

#' @export
#' @import dplyr
#' @importFrom dplyr filter group_by mutate
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
    sort.tops <- top.lab %>% as.data.frame() %>% arrange(desc(Freq))  %>% arrange(Cluster)
  }
  ann.df <- sort.tops %>% select(-"Freq")
  colnames(ann.df) <- c("Cluster", Label)
  return(ann.df)
}





##############################################################################################





#' Add the consensus of annotations to a Seurat object from a list
#'
#'This function takes in a list of data frames with different annotation options.
#'Use at least 2 options.  If you simple want to add annotations use "annotate",
#'which is called by this function.The data frames can be created in R, read in
#'from a csv of txt or created with the "get_annotations" function. All names are
#'changed to lower case letters.  Cell type names must be exactly the same to match.
#'For example "astrocytes" and "Astrocytes" will be considered the same, but
#'"astrocytes" and "astrocyte" will not. The data slot to annotate must be indicated.
#'To label the data slot with the consensus annotation set annotation_name. For example
#'annotation_name = "CellTypeCon".

#' @export
#' @import dplyr
#' @importFrom dplyr select

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
  df3$consensus <- apply(df3[, -1], 1, function(row) {
    if (is.data.frame(row)) {
      word_counts <- table(unlist(row))
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
  print(dfsort)
  dfcon <- dfsort %>% dplyr::select("Cluster","consensus")
  print(dfcon)

  # now add the annotations into the seurat
  # use the simple annnotation function
  seu <- annotate(seu, annotations = dfcon$consensus, to_label,
                  annotation_name)
}


#' Add the annotations to a Seurat object
#'
#'This function takes in a column from a data frame or vector of cell type
#'annotations. For example annotations = df$CellTypes or annotations = c("Neurons",
#'"astrocytes","radial glia").The data slot to annotate must be indicated.
#'To label the data slot with the consensus annotation set annotation_name.
#'Default annotation_name = "CellType".
#'Example: seu <- annotate(seu, annotations = df$CellTypes, to_label = seu$clusters, annotation_name = "CellTypes")

#' @export
#' @import Seurat
#' @importFrom Seurat AddMetaData

annotate <- function(seu, annotations, to_label, annotation_name = "CellType"){
  Idents(seu) <- to_label
  names(annotations) <- levels(seu)
  seu <- RenameIdents(seu, annotations)
  seu <- AddMetaData(object=seu, metadata=Idents(seu), col.name = annotation_name)

}

####### function to make a table with the most frequent cell type per cluster
#' This function takes in a list of annotation data frames
#' The outputs of get_annotations or a manually created data frame
#' The colnames must be "Cluster" and a different column name for the prediction method.
#' @export
#' @import dplyr
#' @importFrom dplyr select
annotate_df <- function(ann.list){
  # not all annotations methods will always present
  # easiest to make the input a list
  df.merge <- Reduce(function(x,y) merge(x,y, by="Cluster", all = TRUE),
                     ann.list)
  # convert all to lowercase
  df2 <- lapply(df.merge, function(x) {tolower(as.character(x))})
  # back into a dataframe
  df3 <- as.data.frame(do.call(cbind, df2))
  df3$consensus <- apply(df3[, -1], 1, function(row) {
    if (is.data.frame(row)) {
      word_counts <- table(unlist(row))
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
  return(dfsort)
}




######### functions to compare groups in an annotated seurat object ################

#' Plot barcharts of proportion of cell types grouped by a list of variables
#'
#'This function takes a Seurat object and creates a boxplot each variable in the list.

#' @export

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
            scale_y_continuous(labels = scales::percent_format(),
                               expand = c(0,0)) +
            scale_fill_manual(values = my_colours)+
            theme_classic() +
            theme(text = element_text(size=15),
                  axis.text.x = element_text(angle=90, hjust=1))
          + xlab(groups) + ylab('Percent of Cell type') + RotatedAxis())
  }

}


# function called above

plotproportions <- function(seu, var.list, xgroup, varnames, my_colours = 'default'){
  for (i in seq_along(var.list)) {
    var <- var.list[[i]]
    proportionplots(seu = seu, seu.var = var, seu.lable = xgroup, groups = varnames[i], my_colours = my_colours)
  }
}


#' Plot mean values by group in dotplot or heatmap
#'
#'This function takes a Seurat object and creates a heatmap or dotplot.
#'This function takes a Seurat object and a list of variables to plot. One
#'var_names is a vector with the levels of the seurat groups to be plotted
#'group is the cell types or x axis grouping variable
#'cluster_order is a character vector that can be used to reorder cell type annotations
#'marker_order is the order to plot the markers


#' @export
#' @import data.table
plotmean <- function(plot_type = 'heatmap',seu, group, markers, var_names, slot = 'scale.data',
                     xlab = 'Cell Types', ylab = 'Markers',
                     cluster_order=NULL, marker_order=NULL, low_colour = "grey",
                     mid_colour = "white", high_colour = "red"){
  express.by.cluster <- as.data.frame(AverageExpression(seu, features = markers,
                                                        group.by = group, slot = 'scale.data'))
  express.by.cluster <- as.data.frame(scale(express.by.cluster))
  if(length(var_names) == length(express.by.cluster)){
    col.names <- var_names
  }else{
    col.names <- colnames(express.by.cluster)
  }
  # Check if cluster_order is NULL and set it to var_names if it is
  if (is.null(cluster_order)) {
    cluster_order <- var_names
  }
  if (is.null(marker_order)) {
    marker_order <- markers
  }
  colnames(express.by.cluster) <- col.names
  AB <- row.names(express.by.cluster)
  ex.data <- cbind(AB,express.by.cluster)
  dt <- data.table(ex.data)
  dt.long <- melt(dt, id.vars = "AB")

  if (plot_type == 'heatmap') {

    print(ggplot(dt.long, aes(x = factor(variable, levels = cluster_order),
                              y = factor(AB, levels = marker_order))) +
            geom_raster(aes(fill = value)) +
            scale_fill_gradient2(low = low_colour, mid = mid_colour,
                                 high = high_colour, na.value = "grey") +
            labs(x = xlab, y = ylab) +
            theme_bw() + theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.95, vjust = 0.2),
                               axis.text.y = element_text(size = 12),
                               plot.title = element_text(size = 12)))
  } else if (plot_type == 'dotplot') {
    a <- DotPlot(seu, features = AB, group.by = group)
    pct.exp <- as.data.frame(a$data) %>% select(features.plot, id, pct.exp)
    colnames(dt.long) <- c("Markers", "id", "expression")
    colnames(pct.exp) <- c("Markers", "id", "proportion")
    df.exp.pct <- merge(dt.long, pct.exp, by = c("Markers", "id"))
    data <- df.exp.pct %>% mutate(Cell.type = factor(id, levels = cluster_order))
    data <- data %>% mutate(Marker = factor(Markers, levels = marker_order))
    print(ggplot(data = data, aes(x = Marker, y = Cell.type, color = expression, size = proportion)) +
            geom_point() +
            scale_color_gradient2(low = low_colour, mid = mid_colour, high = high_colour) +
            ylab(ylab) + xlab(xlab) + RotatedAxis() +
            theme_classic() +
            theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.95, vjust = 0.2),
                  axis.text.y = element_text(size = 12),
                  plot.title = element_text(size = 12))
            )
  } else {
    print("not a valid plotting option")
  }
}


##############################################################################################


#' Prepare a data frame for statistics
#'
#'This function takes a Seurat object and creates a data frame with the expression
#'values for each cell for each marker.  Columns are created for all the indicated
#'variables.  All the variables must exist as meta data in the Seurat object.
#'
#' @export
#' @import data.table
#' @importFrom data.table melt

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

#' Run one-way or two-way ANOVAs and Tukey's HSD test
#'
#'This function takes a data frame prepared by "Prep_for_stats" and outputs
#'a list of ANOVA and TUKEY results. For two-way ANOVAs, stat_type = "ANOVA2". The
#'interaction effect is calculated.
#' group_cols is a vector with the columns to get the means from if use_means = NULL
#' the means will not be calculated and n = number of cells instead of samples per group.
#'id1 is the independent variable to compare (dependent variable is the expression)
#' id2 is for 2 way anova and should be Marker or Celltype but other options are possible.
#' value_col is the column name with the expression values.  This only needs to be
#' changed if the input data frame wasn't created with the "Prep_for_stats" function.

#' @export

run_stats <- function(input_df, group_cols = c("Sample", "CellType", "Marker","Genotype"),
                      value_col = "value",
                      stat_type = "ANOVA",
                      id1 = "Genotype",
                      id2 = NULL,
                      use_means = TRUE,
                      loop_by = "none") {
  # checking if the input is correct to the function
  if (!all(group_cols %in% colnames(input_df))) {
    stop("Some group columns are not present in the input data frame.")
  }

  if (!value_col %in% colnames(input_df)) {
    stop("The value column is not present in the input data frame.")
  }
  # make empty list for the results
  aov.l <- list() # for ANOVA output
  tuk.l <- list() # for TUKEY output
  # we can either use each cell as an n or the replicates as n
  # if get means is true we will use a replicate as an n and need to calculate means
  # get_means is a function below
  if (use_means){
    df_means <- get_means(input_df, group_cols, value_col)
    print("Mean values per group have been calculated")
  } else {
    df_means <- input_df
    names(df_means)[names(df_means) == 'value'] <- 'expression'
  }
  # we can select one-way or two way anova
  # here we are selecting one-way anova
  # we have only one independent variable
  if (stat_type == "ANOVA") {
    print("Running one way ANOVA")
    # if we want all the variables that are not the independed variable to be merged loop by is none
    if(loop_by== "none"){
      one_way <- aov(expression ~ df_means[[id1]], data = df_means)
      aov.l <- summary(one_way)
      # now run the posthoc test
      tukey <- TukeyHSD(one_way)
      # add the results to the tukey list
      tuk.l <- tukey
      # create a data frame for the anova
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
      # create a data frame for the tukey tests
      rownames <- row.names(tuk.l[[1]])
      tuk_df  <- data.frame(
        IndependentVariable = id1,
        Contrast = rownames,
        tuk.l[[1]], row.names = NULL)
      print("TukeyHSD results for all celltypes and markers combined")
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_df)
    }

    # if each cell type or another variable is to be evaluated independently then we loop over that variable defined by loop by
    else{
      var_list <- unique(df_means[[loop_by]])
      # to store outputs and format in a readable way
      for (i in var_list) {
        # from rlang library function !!syn sets the variable loop_by to be a character
        df <- df_means %>% filter(!!sym(loop_by) == i)
        one_way <- aov(expression ~ df[[id1]], data = df)
        output <- summary(one_way)
        aov.l[[as.character(i)]] <- output # Append output to list
        # now the posthoc test
        tukey <- TukeyHSD(one_way)
        tuk.l[[as.character(i)]] <- tukey
      }
      # create a data frame for each anova output and put it in the list
      aov_df <- do.call(rbind, lapply(seq_along(aov.l), function(i) {
        data.frame(LoopVariable = names(aov.l)[i],
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
      print(paste("ANOVA results for each ", loop_by," comparing ", id1))
      # create a data frame for each tukey post hoc test  and put it in the list
      tuk_df <- do.call(rbind, lapply(seq_along(tuk.l), function(i){
        rownames <- row.names(tuk.l[[i]][[1]])
        data.frame(LoopVariable = names(tuk.l)[i],
                   IndependentVariable = id1,
                   Contrast = rownames,
                   tuk.l[[i]][[1]], row.names = NULL)
      }))
      print(paste("TukeyHSD results for each ", loop_by," comparing ", id1))
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_df)
    }
    # If two way anova is select two independent variables are compared
  }else if (stat_type == "ANOVA2") {
    # run the 2 way anova
    print("Running 2-way ANOVA")
    if (loop_by == "none"){
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
        # put the outputs into a list
        output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)
      }
      # now if you are running a loop it will loop by the selected variable
    }else {
      # get the vector of entries to loop by
      var_list <- unique(df_means[[loop_by]])
      # to store outputs and format in a readable way
      for (i in var_list) {
        df <- df_means %>% filter(!!sym(loop_by) == i)
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
          data.frame(LoopVariable = names(aov.l)[i],
                     Contrast = rownames,
                     aov.l[[i]][[1]], row.names = NULL

          )
        )
      }))
      # get ride of the residual which will confuse things and the extra row added
      aov_df <- aov_df[2:8]
      # fix the column names
      colnames(aov_df) <- c("LoopVariable","Contrast","Df","SumSq","MeanSq",
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
      output_list <- list(ANOVA = aov_df,TukeyHSD = tuk_summary.l)

    }
  } # end of option of for 2 way anova
  return(output_list)
}


# original get means function

get_means <- function(df, group_cols, value_col) {
  df_means <- df %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(expression = mean(value)) %>%
    distinct(across(all_of(group_cols)), expression, .keep_all = TRUE) %>%
    select(-value)
  return(df_means)
}


