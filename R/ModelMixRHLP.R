ModelMixRHLP <- setRefClass(
  "ModelMixRHLP",
  contains = "FData",
  # Define the fields
  fields = list(
    G = "numeric",
    # number of clusters
    K = "numeric",
    # number of regimes
    p = "numeric",
    # dimension of beta (order of polynomial regression)
    q = "numeric",
    # dimension of w (order of logistic regression)
    variance_type = "numeric",
    nu = "numeric" # degree of freedom
  )
)

ModelMixRHLP <- function(fData, G, K, p, q, variance_type) {
  if (variance_type == variance_types$homoskedastic) {
    nu <<- (G - 1) + G * ((q + 1) * (K - 1) + K * (p + 1) + 1)
  }
  else{
    nu <<- (G - 1) + G * ((q + 1) * (K - 1) + K * (p +  1) + K)
  }

  new(
    "ModelMixRHLP",
    Y = fData$Y,
    X = fData$X,
    m = fData$m,
    n = fData$n,
    vecY = fData$vecY,
    G = G,
    K = K,
    p = p,
    q = q,
    variance_type = variance_type,
    nu = nu
  )
}
