# CelltypeR

An R library with functions for a work flow to use flow cytometry to quantify cell types within a complex tissue.

R library "CelltypeR"

Workbooks:
- Workbooks using the CelltypeR workflow and functions for each data set analysis
- Workbooks with all code used to generate figures.

FlowCytometry_Data: 
- raw flow cytometry data organized by data set

ExampleOuts:
- reference matrix for correlation predicitons
- trained random forest classifier for predictions
- statistics outputs

Rscripts: 
- Analysis using CelltypeR for larger amounts of data run outside workbooks

CelltypeR: 
- this folder contains all library documents
- All functions are in the file CelltypeR.R

If you use this library please site: 
Thomas, Rhalena A., et al. "CelltypeR: A flow cytometry pipeline to annotate, characterize and isolate single cells from brain organoids." bioRxiv (2022): 2022-11.
doi: https://doi.org/10.1101/2022.11.11.516066


# To install the CelltypeR library

Install devtools package

```
install.packages("devtools")
# load the library
library(devtools)
```

Install CelltypeR library using devtools

```
devtools::install_github("RhalenaThomas/CelltypeR/CelltypeR")
```
Load CelltypeR library

```
library("CelltypeR")
```

These are dependencies that you may not have installed: Rphenograph, FlowSOM, flowCore 
Note: Rphenograph or FlowSOM are only used if you cluster with these algorithms instead of the Seurat implementation of Louvain network detection.  
```
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("JinmiaoChenLab/Rphenograph")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FlowSOM")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")

```
Note: the current version of CelltypeR is not compatible with Seurat Version 5. Please be sure to install version 4 or earlier. The package will be updated to accommodate different data structures for clustering and plotting functions soon. 


# Help and Contributions

If you are encountering an error or need help please open an issue. https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue
Contributions are welcome please contact Rhalena Thomas. 

# Citation

If you use this R package or any of the code or data in this repository in your research please cite:
https://www.biorxiv.org/content/10.1101/2022.11.11.516066v3

