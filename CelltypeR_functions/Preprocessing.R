# Data preprocessing functions
# fsc_to_df                        (Rhalena)
# harmonize                        (Rhalena)
# df_to_seurat                     (Rhalena)


##############################################################################################
# fsc_to_df
# reads in fsc files in a folder to create and R list and save a dataframe
# sample names are input by user

fsc_to_df <- function(input_path){ 
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
  return(copy_flowset)
}





##############################################################################################
# harmonize
# transform, align, retro transform flowset object and save a csv file


inversebiexponentialTransform <- function(flowset,a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0){
  copy_flowset=flowset[seq(along=flowset)] #makes a copy of the input flowset
  for (i in 1:length(copy_flowset)){ #loop though index to get each flowframe
    copy_flowset[[i]]@exprs=a*exp(b*(copy_flowset[[i]]@exprs-w))-c*exp(-d*(copy_flowset[[i]]@exprs-w))+f
  }
  return(copy_flowset)
}

biexp  <- biexponentialTransform("biexp transform",a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0) #creates the transformation object (the values make it equivalent to arcsinh transform)
transformed_flowset <- transform(flowset, transformList(colnames(flowset), biexp)) #Apply the transformation

if(write_fcs_files==TRUE){#Check if user set the conditional for writing several fcs files
  write.flowSet(transformed_flowset,outdir=paste0(output_path,"transformed_flowset"))#writes the flowset
  transformed_flowset=read.flowSet(path=paste0(output_path,"transformed_flowset"),phenoData = "annotation.txt")
}


normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[c(3,5:6,9:length(colnames(transformed_flowset)))],max.lms = 2,peak.density.thr = 0.01) #Detects and align 2 peaks on the marker 3,5,6,9...14. 
expbe_norm2=normtr$flowset
normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[c(4,7:8)],max.lms = 1,peak.density.thr = 0.05)#Detects and align 1 peak 
aligned_transformed_flowset=normtr$flowset
retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset) #apply the function for cancelling biexp transform 


transformed_flowset <- transform(flowset, transformList(colnames(flowset), biexp))

# a function that takes in a subset filter if desired
# aligns if desired
# saves the csv file

# needs to have if statements
# needs to get an input of the column names that should be detected

harmonize <- function(flowset) {
  normtr <- gaussNorm(flowset = transformed_flowset,channel.names = colnames(transformed_flowset)[c(3:length(colnames(transformed_flowset)))], 
                      max.lms = 2, peak.density.thr = 0.01, peak.distance.thr = 0.5) #Detects and align 2 peaks (max.lms) for all markers except FSC and SSC (columns 3 to number of markers) the threshold are data-dependant and may have to differ from one analysis to another.
  aligned_transformed_flowset <- normtr$flowset #Extract the flowset from the result of the alignment
  retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset)
}


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



flowset_to_csv=function(flowset, output_path){  
  list_of_flowframes=fsApply(rename_markers(flowset),function(x){as.data.frame(x@exprs)})#Makes a list of dataframes
  list_names=names(list_of_flowframes)#extract names for not calling names() function at each loop
  for (index in seq_along(list_of_flowframes)){ #Iterates along the index for adding sample names
    list_of_flowframes[[index]]$Sample = list_names[index]
    colnames(list_of_flowframes[[index]]) = colnames(list_of_flowframes[[1]]) 
  }
  ds=list.rbind(list_of_flowframes)#binds every fcs file in a single dataframe
  ds$cell=as.factor(unlist(lapply(as.list(c(1:length(flowset))),function(x){c(1:nrow(flowset[[x]]@exprs))})))#add cell IDs - cell count per sample 
  write.csv(ds,file=paste0(output_path,deparse(substitute(flowset)),".csv"))#save the R data for further usage
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
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu)
}








