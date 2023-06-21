# CelltypeR

An R library with functions for a work flow to use flow cytometry to quantify cell types within a complex tissue.

R library "CelltypeR"

Workbooks:
- analysis steps and usage of CelltypeR functions.
- Figures organized by figure.

Data folder:
- raw flow cytometry data (sample.fsc)
- reference matrix for correlation predicitons
- trained random forest classifier for predictions

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

There are dependencies that you may not have installed: Rphenograph, FlowSOM, flowCore 

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


