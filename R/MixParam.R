source("R/enums.R")
source("R/utils.R")

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
    init_hlp = function(mixModel, phiW, try_algo){
      nm <- mixModel$m * mixModel$n
      if  (try_algo == 1){
        for (g in (1:mixModel$G)){
          problik <- modele_logit(Wg[,,g],phiW)
          pi_jgk[,,g] <<- problik[[1]]
        }
      }
      else{
        for (g in (1:mixModel$G)){
          Wg[,,g] <<- rand(mixModel$q+1, mixModel$K-1);#initialisation aléatoire du vercteur paramètre du IRLS
          problik <- modele_logit(Wg[,,g], phiW)
          pi_jgk[,,g] <<- problik[[1]]
        }
      }
    },

    initParam = function(mixModel, phi, mixOptions, try_algo){
      alpha_g <<- 1/(mixModel$G * ones(mixModel$G, 1))
      init_hlp(mixModel, phi$phiW, try_algo) # setting Wg and pi_jgk
      if (mixOptions$init_kmeans){
        # run k means
        kmeans_res <- kmeans(mixModel$X, iter.max = 400, centers=mixModel$G, nstart=20, trace=FALSE)
        klas <- kmeans_res$cluster
        for (g in 1:mixModel$G){
          Xg <- mixModel$X[klas==g,]
          initRegressionParam(Xg, g, mixModel$K, mixModel$p, phi$phiBeta, mixOptions$variance_type, try_algo)
        }
      }
      else{
        ind <- sample(mixModel$n)
        for (g in 1:mixModel$G){
          if (g<G){
            Xg <- mixModel$X[ind[(g-1)*round(mixModel$n/mixModel$G) +1 : g*round(mixModel$n/mixModel$G)],]
          }
          else{
            Xg = mixModel$X[ind[(g-1)*round(n/G) +1 : length(mixModel$X)],]
          }
          initRegressionParam(Xg, g, mixModel$K, mixModel$p, phi$phiBeta, mixOptions$variance_type, try_algo)
        }
      }
    },

    initRegressionParam = function(Xg, g, K, p, phiBeta, variance_type, try_algo){
       n <- nrow(Xg)
       m <- ncol(Xg)

       if (try_algo==1){
          # decoupage de l'echantillon (signal) en K segments
          zi <- round(m/K)-1

          beta_k <- matrix(NA, p+1, K)
          sigma <- c()

          for (k in 1:K){
            i <- (k-1)*zi+1
            j <- k*zi
            Xij <- Xg[,i:j]
            Xij <- matrix(t(Xij), ncol = 1)
            phi_ij <- phiBeta[i:j,]
            Phi_ij <- repmat(phi_ij,n,1)

            bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%Xij
            beta_k[,k] <- bk

            if (variance_type == variance_types$common){
              sigma <- var(Xij)
            }
            else{
              mk <- j-i+1 #length(Xij);
              z <- Xij - Phi_ij %*% bk;
              sk <- t(z) %*% z/(n*mk);
              sigma[k] <- sk;
            }
          }
       }
       else{ # initialisation aléatoire
         Lmin <- round(m/(K+1)) #nbr pts min dans un segments
         tk_init <- zeros(1,K+1)
         K_1 <- K
         for (k in 2:K) {
           K_1 <- K_1-1;
           temp <- tk_init[k-1] + Lmin : m - K_1*Lmin;
           ind <- randperm(length(temp));
           tk_init[k] <- temp(ind(1))
         }
         tk_init[K+1] <- m

         beta_k <- matrix(NA, p+1, K)
         sigma <- c()
         for (k in 1:K){
           i <- tk_init[k] + 1
           j <- tk_init[k+1]
           Xij <- Xg[,i:j]
           Xij <- matrix(t(Xij), ncol = 1)
           phi_ij <- phiBeta[i:j,]
           Phi_ij <- repmat(phi_ij,n,1)


           bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%Xij
           beta_k[,k] <- bk

           if (variance_type == variance_types$common){
             sigma <- var(Xij)
           }
           else{
             mk <- j-i+1 #length(Xij);
             z <- Xij - Phi_ij %*% bk;
             sk <- t(z) %*% z/(n*mk);
             sigma[k] <- sk;
           }
         }
       }

       betag[,,g] <<- beta_k
       if (variance_type == variance_types$common){
         sigmag[g] <<- sigma
       }
       else{
         sigmag[,g] <<- sigma
       }
    },


    MStep = function(mixModel, mixStats, phi, mixOptions){
      alpha_g <<- t(colSums(mixStats$h_ig))/mixModel$n
      for (g in 1:mixModel$G){
        temp <- repmat(mixStats$h_ig[,g], 1, mixModel$m) # [m x n]
        cluster_weights <- matrix(t(temp), mixModel$m*mixModel$n, 1) # cluster_weights(:)% [mn x 1]
        tauijk <- mixStats$tau_ijgk[,,g] #[(nxm) x K]

        if (mixOptions$variance_type == variance_types$common){
          s <- 0
        }
        else{
          sigma_gk <- zeros(mixModel$K, 1)
        }

        beta_gk <- matrix(NA, mixModel$p+1, mixModel$K)
        for (k in 1:mixModel$K){
          segments_weights <- tauijk[,k] # poids du kieme segment   pour le cluster g
          # poids pour avoir K segments floues du gieme cluster flou
          phigk <- (sqrt(cluster_weights * segments_weights) %*% ones(1,mixModel$p+1))*phi$phiBeta #[(n*m)*(p+1)]
          Xgk <- sqrt(cluster_weights * segments_weights) * mixModel$XR

          # maximization w.r.t beta_gk: Weighted least squares
          beta_gk[,k] <- solve(t(phigk)%*%phigk + .Machine$double.eps*diag(p+1))%*%t(phigk)%*%Xgk # Maximization w.r.t betagk

          #    the same as
          #                 W_gk = diag(cluster_weights.*segment_weights);
          #                 beta_gk(:,k) = inv(phiBeta'*W_gk*phiBeta)*phiBeta'*W_gk*X;
          #   Maximization w.r.t au sigma_gk :
          if (mixOptions$variance_type == variance_types$common){
            sk <- sum((Xgk - phigk %*% beta_gk[,k])^2)
            sk <- s+sk
            sigma_gk <- s/sum((cluster_weights%*%ones(1,mixModel$K))*tauijk)
          }
          else{
            sigma_gk[k] <- sum((Xgk-phigk%*%beta_gk[,k])^2)/(sum(cluster_weights*segments_weights));
          }
        }
        betag[,,g] <<- beta_gk
        sigmag[,g] <<- sigma_gk;

        # Maximization w.r.t W
        #  IRLS : Regression logistique multinomiale pondérée par cluster
        # setting of Wg[,,g] and pi_jgk
        Wg_init <- Wg[,,g]
        res_irls <- IRLS_MixFRHLP(cluster_weights, tauijk, phi$phiW, Wg_init, mixOptions$verbose_IRLS)
        Wg[,,g] <<- res_irls[[1]]
        piik <- res_irls[[2]]
        pi_jgk[,,g] <<- repmat(piik[1:mixModel$m,], mixModel$n, 1)

      }

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
    sigmag <- matrix(NA, mixModel$K, mixModel$G)
  }
  pi_jgk <- array(0, dim=c(mixModel$m*mixModel$n, mixModel$K, mixModel$G))
  alpha_g <- matrix(NA, mixModel$G)
  new("MixParam", Wg=Wg, betag=betag, sigmag=sigmag, pi_jgk=pi_jgk, alpha_g=alpha_g)#, mixModel = mixModel)
}
