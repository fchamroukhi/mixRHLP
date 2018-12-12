rm(list = ls())
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/ModelLearner.R")

mixData <- MyData$new()
mixData$setData("data/generated_data_1.txt")
G <- 3; # nombre de clusters
K <- 3; #nombre de regimes
p <- 1; #dimension de beta (ordre de reg polynomiale)
q <- 1; #dimension de w (ordre de reg logistique)
mixModel <- MixModel(mixData,G,K,p,q)

n_tries=1
max_iter=1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE
modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$free)

solution <- EM(mixModel, modelOptions)
mixParamSolution <- solution[[1]]
mixStatsSolution <- solution[[2]]
