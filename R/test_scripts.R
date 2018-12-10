rm(list = ls())
source("R/Phi.R")
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/MixParam.R")
source("R/MixStats.R")
source("R/ModelLearner.R")
source("R/enums.R")
source("R/utils.R")


testData <- function(){
  m <- MyData$new()
  m$setData("data/generated_data_1.txt")
  #m$generateExampleData()
  return(m)
}
m <- testData()

phiTest <- function(){
  mixData <- MyData$new()
  mixData$setData("data/generated_data_1.txt")
  G <- 3; # nombre de clusters
  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  mix <- MixModel(mixData,G,K,p,q)

  phi <- Phi$new()
  phi$setPhiN(m$t,mix$p,mix$q, mix$n)
  return(phi)
}
phi <- phiTest()

testModelCreation <- function(){
  mixData <- MyData$new()
  mixData$setData("data/generated_data_1.txt")
  G <- 3; # nombre de clusters
  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  mix <- MixModel(mixData,G,K,p,q)

  n_tries=1
  max_iter=1000
  threshold <- 1e-5
  verbose <- TRUE
  verbose_IRLS <- TRUE
  init_kmeans <- TRUE
  mixOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$common)

  phi <- Phi$new()
  phi$setPhiN(mixData$t,mix$p,mix$q, mix$n)

  param <- MixParam(mix, mixOptions)
  param$initParam(mix, phi, mixOptions, try_algo = 1)

  return(param)
}
param<-testModelCreation()

testModelOptions <- function(){
  n_tries=1
  max_iter=1000
  threshold <- 1e-5
  verbose <- TRUE
  verbose_IRLS <- TRUE
  init_kmeans <- TRUE
  options <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$common)

  return(options)
}
o <- testModelOptions()

testEnum <- function(){
  variance_type <- enum(free, common)
  v_type <- variance_type$common
}
testEnum()


test_MAP <- function(){
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
  verbose_IRLS <- TRUE
  init_kmeans <- TRUE
  mixOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$common)

  stmix <- MixStats(mixModel, mixOptions)
  h_ig <- read.csv("data/tests/h_ig.csv", header = FALSE)
  stmix$h_ig <- as.matrix(h_ig)
  m <- stmix$MAP()
  klas <- m[1]
  Z <- m[2]
}
test_MAP()

test_modellogit <- function(){
  Wg <- array(0,dim=c(2, 2, 3))
  Wg[,,1] <- rand(2, 2)
  Wg[,,2] <- rand(2, 2)
  Wg[,,3] <- rand(2, 2)


  mixData <- MyData$new()
  mixData$setData("data/generated_data_1.txt")
  G <- 3; # nombre de clusters
  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  mix <- MixModel(mixData,G,K,p,q)

  m <- MyData$new()
  m$setData("data/generated_data_1.txt")

  phi <- Phi$new()
  phi$setPhiN(m$t,mix$p,mix$q, mix$n)

  problik <- modele_logit(Wg[,,1], phi$phiW)
  print(dim(problik[[1]]))
}

test_modellogit()




testKmeans <- function(){
  mixData <- MyData$new()
  mixData$setData("data/generated_data_1.txt")
  kmeans(mixData$X, iter.max = 400, nstart = 20, trace=FALSE)
}


