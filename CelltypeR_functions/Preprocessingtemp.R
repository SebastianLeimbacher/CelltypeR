# Data preprocessing functions
# fsc_to_fs                        (Rhalena)
# harmonize                        (Rhalena)
# flowset_to_csv
# df_to_seurat                     (Rhalena)



##### other required libraries - check each if actually needed #####
require("flowCore") #Used for reading the data
require("ggplot2")
require("flowStats") #Alignment functions
require("scales") #scale colour intensity for visualization
require("dplyr")
require("reshape2") # for plotting density function
require("ggridges") # for plotting density function
require("stringr")
require("rlist")
require("Seurat")




##############################################################################################
# fsc_to_df
# reads in fsc files creates a flowset file
# filtering if wanted on each input file
# users can rename samples after


fsc_to_df <- function(input_path, downsample = 'none'){ 
  flowset = read.flowSet(path=input_path,transformation = FALSE ,
               emptyValue = FALSE,truncate_max_range = FALSE, 
               package="flowCore")
  fsApply(flowset,
          function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(flowset))]})
  return(flowset)
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
  df2 <- df %>% dplyr::select(AB_vector)
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








