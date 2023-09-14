# updating the R library

# open the R project


# set working directory to the R library (it will already be set by the project)

setwd("/Users/rhalenathomas/GITHUB/CelltypeR/CelltypeR")

# update documentation
library(roxygen2)
roxygen2::roxygenise()  # this seems to load the package

#
roxygenise()



# update the NAMESPACE and DESCRIPTION

devtools::build()
devtools::install()


# to check package versions

packageVersion("packagename")  # use the quotes
