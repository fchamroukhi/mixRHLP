ModelOptions <- setRefClass(
  "ModelOptions",
  fields = list(
    n_tries = "numeric",
    max_iter = "numeric",
    threshold = "numeric",
    verbose = "logical",
    verbose_IRLS = "logical",
    init_kmeans = "logical",
    variance_type = "numeric"
  )
)

ModelOptions<-function(n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans, variance_type){
  new("ModelOptions",n_tries=n_tries, max_iter=max_iter, threshold=threshold, verbose=verbose, verbose_IRLS=verbose_IRLS, init_kmeans=init_kmeans, variance_type=variance_type)
}



