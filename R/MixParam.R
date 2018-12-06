source("R/enums.R")



MixParam <- setRefClass(
  "MixParam",
  fields = list(
    Wg = "array",
    betag = "array",
    sigmag = "matrix",
    pi_jgk = "array",
    alpha_g = "matrix"
  ),
  methods = list(
    init_hlp <- function(mixModel,phiW,try_algo){
      nm <- mixModel$m*mixModel$n
      if  (try_algo == 1){
        for (g in (1:mixModel$G)){
          pi_jgk[,,g] <<- modele_logit(Wg[,,g],phiW);
        }
      }
      else{
        for (g in (1:mixModel$G)){
          Wg[,,g] <<- rand(mixModel$q+1, mixModel$K-1);#initialisation aléatoire du vercteur param�tre du IRLS
          pi_jgk[,,g] <<- modele_logit(Wg[,,g],phiW);
        }
      }
    },

    initParam = function(mixModel){
      alpha_g <<- 1/(mixModel$G*ones(mixModel$G,1))

    }
  )
)

MixParam<-function(mixModel, options){
  #mixModel <- mixModel
  Wg <- array(0,dim=c(mixModel$q+1, mixModel$K-1, mixModel$G))
  betag <- array(NA, dim=c(mixModel$q+1, mixModel$K, mixModel$G))
  if (options$variance_type == variance_types$common){
    sigmag <- matrix(NA, mixModel$G)
  }
  else{
    sigmag <- matrix(NA, mixModel$G, mixModel$K)
  }
  pi_jgk <- array(0, dim=c(mixModel$m*mixModel$n, mixModel$K, mixModel$G))
  alpha_g <- matrix(NA, mixModel$G)
  new("MixParam", Wg=Wg, betag=betag, sigmag=sigmag, pi_jgk=pi_jgk, alpha_g=alpha_g)#, mixModel = mixModel)
}
