rm(list = ls())
source("dataset.R")
source("MixModel.R")
source("ModelOptions.R")
source("ModelLearner.R")



mixData <- MyData$new()
mixData$setData("data/generated_data_1.txt")
G <- 3; # number of clusters
K <- 3; # number of regimes (polynomial regression components)
p <- 1; # degree of the polynomials
q <- 1; # order of the logistic regression (by default 1 for contiguous segmentation)
mixModel <- MixModel(mixData,G,K,p,q)

n_tries <- 20
max_iter <- 1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE
modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$free)

####
# EM Algorithm
####
solution <- EM(mixModel, modelOptions)
mixParamSolution <- solution[[1]]
mixStatsSolution <- solution[[2]]

# show the results
mixStatsSolution$showDataClusterSegmentation(mixModel, mixParamSolution)

####
# CEM Algorithm
####
#solution <- CEM(mixModel, modelOptions)
#mixParamSolution <- solution[[1]]
#mixStatsSolution <- solution[[2]]
#mixStatsSolution$showDataClusterSegmentation(mixModel, mixParamSolution)

