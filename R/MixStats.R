MixStats <- setRefClass(
  "MixStats",
  fields = list(
    h_ig="matrix", #post probabilities
    c_ig="matrix",
    klas="matrix",
    # pi_jgk = "matrix",
    Ex_g = "matrix",
    log_lik="numeric",
    com_loglik="numeric",
    stored_loglik = "list",
    stored_com_loglik = "list",
    tau_ijgk = "array", # segments post probabilities
    log_tau_ijgk = "array",
    BIC="numeric",
    ICL="numeric",
    AIC="numeric",
    cpu_time = "numeric",
    log_fg_xij="matrix",
    log_alphag_fg_xij="matrix",
    polynomials="array",
    weighted_polynomials="array"
  ),
  methods=list(
    showDataClusterSegmentation = function(mixModel, mixParamSolution){

      t <- seq(0,mixModel$m-1)
      colors = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')
      colors_cluster_means = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')

      par(mfrow=c(1,1))
      for (g in 1:mixModel$G){
        cluster_g = mixModel$X[klas==g,]

        if (g==1){
          plot(t, t(cluster_g)[,1], type = 'l', col=colors[g], xlab='Time', ylab = 'y')
        }
        else{
          lines(t, t(cluster_g)[,1], type = 'l', col=colors[g])
        }

        for (i in 1:nrow(cluster_g)){
          #print(colors[i])
          lines(t, t(cluster_g)[,i], col=colors[g])
        }

      }

      for(g in 1:mixModel$G){
        lines(t,Ex_g[,g], col=colors_cluster_means[g], lwd=5)
      }


      ##########################################################################"
      ##########################################################################
      ##########################################################################

      for (g in 1:mixModel$G){
        par(mfrow=c(2,1))
        cluster_g = mixModel$X[klas==g,]
        plot(t, t(cluster_g)[,1], type = 'l', col=colors[g], ylab = 'y')
        for (k in 1:K){
          lines(t,polynomials[,k,g], lty =2 , col="black",lwd=1)
        }
        lines(t,Ex_g[,g], col=colors_cluster_means[g], lwd=5)

        for (k in 1:mixModel$K){
          if (k==1){
            plot(t,mixParamSolution$pi_jgk[1:mixModel$m,k,g], xlab='Time', ylab = 'Logistic proportions', type = "l", lwd=2, col=colors_cluster_means[k])
          }
          else{
            lines(t,mixParamSolution$pi_jgk[1:mixModel$m,k,g] , lwd=2, col=colors_cluster_means[k])
          }
        }
      }


    },

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
      c_ig <<- ikmax%*%ones(1,K) == ones(N,1)%*%(1:K)
      klas <<- ones(N,1)
      for (k in 1:K){
        klas[c_ig[,k]==1] <<- k
      }
    },

    #######
    # compute the final solution stats
    #######
    computeStats = function(mixModel, mixParam, phi, cpu_time_all){
      for (g in 1:mixModel$G){
        polynomials[,,g] <<- phi$phiBeta[1:mixModel$m, ] %*% mixParam$betag[,,g]
        if (K!=1 && G!=1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk[,,g] * polynomials[,,g]
          Ex_g[,g] <<- rowSums(weighted_polynomials[,,g])
        }
        else if (K==1 && G!=1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk[,g] * polynomials[,,g]
          Ex_g[,g] <<- weighted_polynomials[,,g]
        }
        else if (K!=1 && G==1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk * polynomials[,,g]
          Ex_g[,g] <<- matrix(rowSums(weighted_polynomials[,,g]))
        }
        else{ #(K==1 && G==1)
          weighted_polynomials[,,g] <<- mixParam$pi_jgk * polynomials[,,g]
          Ex_g[,g] <<- weighted_polynomials[,,g]
        }
      }

      Ex_g <<- matrix(Ex_g, nrow = mixModel$m)
      cpu_time <<- mean(cpu_time_all)
      Psi <- c(as.vector(mixParam$alpha_g), as.vector(mixParam$Wg), as.vector(mixParam$betag), as.vector(mixParam$sigmag))
      nu <- length(Psi)
      BIC <<- log_lik - (nu*log(mixModel$n)/2)
      AIC <<- log_lik - nu

      cig_log_alphag_fg_xij <- (c_ig)*(log_alphag_fg_xij);
      com_loglik <<- sum(rowSums(cig_log_alphag_fg_xij));

      ICL <<- com_loglik - nu*log(mixModel$n)/2;

    },

    CStep = function(reg_irls){
      #CStep
      h_ig <<- lognormalize(log_alphag_fg_xij)

      MAP() # setting klas and c_ig

      # Compute the optimized criterion
      cig_log_alphag_fg_xij <- (c_ig)*log_alphag_fg_xij
      com_loglik <<- sum(cig_log_alphag_fg_xij) +  reg_irls
    },

    #######
    # EStep
    #######

    EStep = function(mixModel, mixParam, phi, variance_type){
      for (g in 1:mixModel$G){
        alpha_g <- mixParam$alpha_g[g]
        beta_g <- mixParam$betag[,,g]
        Wg <- mixParam$Wg[,,g]
        pi_jgk <- mixParam$pi_jgk[,,g]
        if (!is.matrix(beta_g)){
          beta_g <- matrix(beta_g)
          pi_jgk <- matrix(pi_jgk)
        }
        log_pijgk_fgk_xij <- zeros(mixModel$n*mixModel$m, mixModel$K)
        for (k in 1:mixModel$K){
          beta_gk <- beta_g[,k]
          if (variance_type == variance_types$common){
            sgk <- mixParam$sigmag[g]
          }
          else{
            sgk <- mixParam$sigmag[k,g]
          }
          z <- ((mixModel$XR - phi$phiBeta %*% beta_gk)^2)/sgk
          log_pijgk_fgk_xij[,k] <- log(pi_jgk[,k]) - 0.5 * (log(2*pi) + log(sgk)) - 0.5 * z # pdf cond à c_i = g et z_i = k de xij
        }

        log_pijgk_fgk_xij <- pmin(log_pijgk_fgk_xij, log(.Machine$double.xmax))
        log_pijgk_fgk_xij <- pmax(log_pijgk_fgk_xij, log(.Machine$double.xmin))

        pijgk_fgk_xij <- exp(log_pijgk_fgk_xij)
        sumk_pijgk_fgk_xij <- rowSums(pijgk_fgk_xij) # sum over k
        log_sumk_pijgk_fgk_xij <- log(sumk_pijgk_fgk_xij) # [n*m, 1]

        log_tau_ijgk[,,g] <<- log_pijgk_fgk_xij - log_sumk_pijgk_fgk_xij %*% ones(1,mixModel$K)
        tau_ijgk[,,g] <<- exp(lognormalize(log_tau_ijgk[,,g]))

        log_fg_xij[,g] <<- rowSums(t(matrix(log_sumk_pijgk_fgk_xij,mixModel$m,mixModel$n))) # [n x 1]:  sum over j=1,...,m: fg_xij = prod_j sum_k pi_{jgk} N(x_{ij},mu_{gk},s_{gk))
        log_alphag_fg_xij[,g] <<- log(alpha_g) + log_fg_xij[,g] # [nxg]
      }
      log_alphag_fg_xij <<- pmin(log_alphag_fg_xij, log(.Machine$double.xmax))
      log_alphag_fg_xij <<- pmax(log_alphag_fg_xij, log(.Machine$double.xmin))

      h_ig <<- exp(lognormalize(log_alphag_fg_xij))

      log_lik <<- sum(log(rowSums(exp(log_alphag_fg_xij))))
    }

  )
)


MixStats<-function(mixModel, options){
  h_ig <- matrix(NA,mixModel$n, mixModel$G)
  c_ig <- matrix(NA,mixModel$n, mixModel$G)
  klas <- matrix(NA, mixModel$n,1)
  Ex_g <- matrix(NA,mixModel$m, mixModel$G)
  log_lik <- -Inf
  com_loglik <- -Inf
  stored_loglik <- list()
  stored_com_loglik <- list()
  BIC <- -Inf
  ICL <- -Inf
  AIC <- -Inf
  cpu_time <- Inf
  log_fg_xij <- matrix(0, mixModel$n, mixModel$G)
  log_alphag_fg_xij <- matrix(0, mixModel$n, mixModel$G)
  polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))
  weighted_polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))
  tau_ijgk <- array(0, dim = c(mixModel$n*mixModel$m, mixModel$K, mixModel$G))
  log_tau_ijgk <- array(0, dim = c(mixModel$n*mixModel$m, mixModel$K, mixModel$G))

  new("MixStats", h_ig=h_ig, c_ig=c_ig, klas=klas, Ex_g=Ex_g, log_lik=log_lik, com_loglik=com_loglik, stored_loglik=stored_loglik, stored_com_loglik=stored_com_loglik, BIC=BIC, ICL=ICL, AIC=AIC, cpu_time=cpu_time,
      log_fg_xij=log_fg_xij, log_alphag_fg_xij=log_alphag_fg_xij, polynomials=polynomials, weighted_polynomials=weighted_polynomials, tau_ijgk=tau_ijgk, log_tau_ijgk=log_tau_ijgk)
}
