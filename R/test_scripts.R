rm(list = ls())
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/MixParam.R")
source("R/MixStats.R")
source("R/ModelLearner.R")
source("R/enums.R")


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
  mixData <- testData()
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
  options <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_types$common)


  param <- MixParam(mix, options)
  param$initParam(mix)

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
  stmix <- MixStats$new()
  h_ig <- read.csv("data/tests/h_ig.csv", header = FALSE)
  stmix$h_ig <- as.matrix(h_ig)
  m <- stmix$MAP()
  klas <- m[1]
  Z <- m[2]
}
test_MAP()

