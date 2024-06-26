% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CelltypeR.R
\name{RFM_train}
\alias{RFM_train}
\title{Train a Random Forest Classifier to predict cell types in a Seurat object}
\usage{
RFM_train(
  seurate_object,
  markers,
  annotations = seu$CellTypes,
  num_folds = 3,
  downsample = "none",
  seed = 222,
  mytry = c(6:9),
  maxnodes = c(12:15),
  trees = c(1000, 1500),
  cores_to_use = 4
)
}
\description{
Requires a Seurat object for the labels to be predicted. AB_list is the list of markers(antibodies)
or selected features(genes) to input into the training model. An annotated Seurat object of matching
data is required. The location of the annotations metadata needs to be designated:
annotated_seurat_object$MetadataAnnotations_slot. If the number of cells is high the function may not run.  If this is the case
use downsample to reduce the number of cells.  The starting data will be randomly selected from to
reduce the total size. This data will then be split in 80/20 training to test data.
A random seed can be set to repeat the training with different random starts and random splitting of
data. The tuning parameters for the random forst model can be optomized. The best values are
selected for the final model from mytry, maxnodes, trees
mtry: number of features to draw at once.
maxnodes: the terminal number of nodes trees in the forest can have.
num_folds is for cross validation to avoid overfitting. The validation data is split from within the training data.
Markers is the vector of the markers to select features from. If a RFM is need with only some of
the markers in a reference data set selecting only those markers can be done with this argument.
}
