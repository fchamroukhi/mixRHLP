source("dataset.R")
source("utils.R")

MixModel <- setRefClass(
  "MixModel",
  contains = "MyData",
  # Define the fields
  fields = list(G="numeric", #  number of clusters
                K="numeric", # number of regimes (polynomial regression components)
                p="numeric", # degree of the polynomials
                q="numeric" # order of the logistic regression (by default 1 for contiguous segmentation)
                ),
  methods = list(

  )
)

MixModel<-function(mixData, G,K,p,q){
  new("MixModel",X=mixData$X, XR=mixData$XR, t=mixData$t, n=mixData$n, m=mixData$m, G=G,K=K, p=p, q=q)
}


