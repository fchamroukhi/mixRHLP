MixStats <- setRefClass(
  "MixStats",
  fields = list(
    h_ig="matrix", #post probabilities
    c_ig="matrix",
    # pi_jgk = "matrix",
    Ex_g = "matrix",
    log_lik="numeric",
    com_loglik="numeric",
    BIC="numeric",
    ICL="numeric",
    AIC="numeric",
    log_alphag_fg_xij="matrix",
    polynomials="array",
    weighted_polynomials="array"
  ),
  methods=list(
    MAP = function(){
      "
      calcule une partition d'un echantillon par la regle du Maximum A Posteriori à partir des probabilites a posteriori
      Entrees : post_probas , Matrice de dimensions [n x K] des probabibiltes a posteriori (matrice de la partition floue)
           n : taille de l'echantillon
           K : nombres de classes
           klas(i) = arg   max (post_probas(i,k)) , for all i=1,...,n
                         1<=k<=K
                   = arg   max  p(zi=k|xi;theta)
                         1<=k<=K
                   = arg   max  p(zi=k;theta)p(xi|zi=k;theta)/sum{l=1}^{K}p(zi=l;theta) p(xi|zi=l;theta)
                         1<=k<=K
      Sorties : classes : vecteur collones contenant les classe (1:K)
           Z : Matrice de dimension [nxK] de la partition dure : ses elements sont zik, avec zik=1 si xi
           appartient à la classe k (au sens du MAP) et zero sinon.
      "
      N <- nrow(h_ig)
      K <- ncol(h_ig)
      ikmax <- max.col(h_ig)
      ikmax <- matrix(ikmax, ncol = 1)
      Z = ikmax%*%ones(1,K) == ones(N,1)%*%(1:K)
      klas = ones(N,1)
      for (k in 1:K){
        klas[Z[,k]==1]=k
      }
      return(list(klas, Z))
    }
  )
)


MixStats<-function(mixModel, options){
  h_ig <- matrix(NA,mixModel$n, mixModel$G)
  c_ig <- matrix(NA,mixModel$n, mixModel$G)
  Ex_g <- matrix(NA,mixModel$n, mixModel$G)
  log_lik <- -Inf
  com_loglik <- -Inf
  BIC <- -Inf
  ICL <- -Inf
  AIC <- -Inf
  log_alphag_fg_xij <- matrix(NA, mixModel$n, mixModel$G)
  polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))
  weighted_polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))

  new("MixStats", h_ig=h_ig, c_ig=c_ig, Ex_g=Ex_g, log_lik=log_lik, com_loglik=com_loglik, BIC=BIC, ICL=ICL, AIC=AIC, log_alphag_fg_xij=log_alphag_fg_xij, polynomials=polynomials, weighted_polynomials=weighted_polynomials)
}
