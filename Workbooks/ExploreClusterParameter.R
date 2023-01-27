# explore cluster parameters
# read in the seurat object and dataframe
# read in the antibody vector


output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/testingLibrary/exp_clusters/"
df <- read.csv(paste(output_path,"df9000.csv", sep = ""))
seu <- readRDS(paste(output_path,"seu9000.RDS", sep = ""))
AB <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29",
        "CD56", "O4","CD140a","CD133","GLAST","CD71")

df2 <- df %>% dplyr::select(AB) # need to add this line into the main explore param function 

param.test <- explore_param(input = seu, 
                       cluster_method = c("flowsom", "phenograph", "louvain"),
                       df_input = df2,
                       for.flowsom.k = c(3,5,10,15,20), 
                       for.phenograph.and.louvain.k = c(20,60,100),
                       for.louvain.resolution = c(0.1,0.5,0.8,1.2),
                       save.stats = TRUE, 
                       save.plot = TRUE, 
                       output_path = output_path)