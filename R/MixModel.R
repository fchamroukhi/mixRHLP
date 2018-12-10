source("R/dataset.R")
source("R/utils.R")

MixModel <- setRefClass(
  "MixModel",
  contains = "MyData",
  # Define the fields
  fields = list(G="numeric", # nombre de clusters
                K="numeric", # nombre de regimes
                p="numeric", # dimension de beta (ordre de reg polynomiale)
                q="numeric" # dimension de w (ordre de reg logistique)
                ),
  methods = list(

  )
)

MixModel<-function(mixData, G,K,p,q){
  new("MixModel",X=mixData$X, XR=mixData$XR, t=mixData$t, n=mixData$n, m=mixData$m, G=G,K=K, p=p, q=q)
}


