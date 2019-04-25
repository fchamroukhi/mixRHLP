"
Curve clustering with the MixFRHLP model and the EM (or a CEM-like) algorithm

%%%% based on matlab code by Faicel Chamroukhi (2011)%%%%%%%

 When using this code please cite the following papers : The two first
 ones concern the model and its use in clusterng and the two last ones
 concern the model and its use in discrimination, and a review.

 @article{Chamroukhi-RHLP-2009,
 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
 	Journal = {Neural Networks},
 	Number = {5-6},
 	Pages = {593--602},
 	Publisher = {Elsevier Science Ltd.},
 	Title = {Time series modeling by a regression approach based on a latent process},
 	Volume = {22},
 	Year = {2009}
 }

 @article{Chamroukhi-MixRHLP-2011,
 	Author = {Sam{\'e}, A. and Chamroukhi, F. and Govaert, G{\'e}rard and Aknin, P.},
 	Issue = 4,
 	Journal = {Advances in Data Analysis and Classification},
 	Pages = {301--321},
 	Publisher = {Springer Berlin / Heidelberg},
 	Title = {Model-based clustering and segmentation of time series with changes in regime},
 	Volume = 5,
 	Year = {2011}
 }

 @article{Chamroukhi-FDA-2018,
 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
 	Note = {DOI: 10.1002/widm.1298.},
 	Volume = {},
 	Title = {Model-Based Clustering and Classification of Functional Data},
 	Year = {2018},
 	Month = {Dec}}

 @article{Chamroukhi-RHLP-FLDA,
 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
 	Journal = {Neurocomputing},
 	Number = {7-9},
 	Pages = {1210--1221},
 	Title = {A hidden process regression model for functional data description. Application to curve discrimination},
 	Volume = {73},
 	Year = {2010}
 }

 @article{Chamroukhi-FMDA-2013,
 	Author = {Chamroukhi, F. and Glotin, H. and Sam{\'e}, A.},
 	Journal = {Neurocomputing},
 	Pages = {153-163},
 	Title = {Model-based functional mixture discriminant analysis with hidden process regression for curve classification},
 	Volume = {112},
 	Year = {2013}
 }
"


rm(list = ls())
source("R/FData.R")
source("R/ModelMixRHLP.R")
source("R/enums.R")
source("R/ModelLearner.R")

load("data/simulatedTimeSeries.RData")
fData <- FData$new()
fData$setData(X, Y)


# setting the model
G <- 3; # number of clusters
K <- 3; # number of regimes (polynomial regression components)
p <- 1; # degree of the polynomials
q <- 1; # order of the logistic regression (by default 1 for contiguous segmentation)
variance_type <- variance_types$hetereskedastic
modelMixRHLP <- ModelMixRHLP(fData, G, K, p, q, variance_type)

n_tries <- 1
max_iter <- 1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE

####
# EM Algorithm
####
# 1. running the em algorithm giving mixModel (data, and model itself) and the model options
solution <- EM(modelMixRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans)
# show the results
solution$plot()



####
# CEM Algorithm
####
#solution <- CEM(modelMixRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans)
# show the results
#solution$plot()
