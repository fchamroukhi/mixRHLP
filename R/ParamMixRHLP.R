#' @export
ParamMixRHLP <- setRefClass(
  "ParamMixRHLP",
  fields = list(
    fData = "FData",
    phi = "list",

    G = "numeric",
    # number of clusters
    K = "numeric",
    # number of regimes
    p = "numeric",
    # dimension of beta (order of polynomial regression)
    q = "numeric",
    # dimension of w (order of logistic regression)
    variance_type = "numeric",
    nu = "numeric", # degree of freedom

    Wg = "array",
    # Wg = (Wg1,...,w_gK-1) parameters of the logistic process:
    # matrix of dimension [(q+1)x(K-1)] with q the order of logistic regression.
    betag = "array",
    # betag = (beta_g1,...,beta_gK) polynomial regression coefficient vectors: matrix of
    # dimension [(p+1)xK] p being the polynomial  degree.
    sigmag = "matrix",
    # sigma_g = (sigma_g1,...,sigma_gK) : the variances for the K regmies. vector of dimension [Kx1]
    pi_jgk = "array",
    # pi_jgk :logistic proportions for cluster g
    alpha_g = "matrix" #cluster weights
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), G = 1, K = 1, p = 2, q = 1, variance_type = 1) {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p, q = q, n = fData$n)

      G <<- G
      K <<- K
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == variance_types$homoskedastic) {
        nu <<- (G - 1) + G * ((q + 1) * (K - 1) + K * (p + 1) + 1)
      }
      else{
        nu <<- (G - 1) + G * ((q + 1) * (K - 1) + K * (p +  1) + K)
      }

      Wg <<- array(0, dim = c(q + 1, K - 1, G))
      betag <<- array(NA, dim = c(p + 1, K, G))
      if (variance_type == variance_types$homoskedastic) {
        sigmag <<- matrix(NA, G)
      }
      else{
        sigmag <<- matrix(NA, K, G)
      }
      pi_jgk <<- array(0, dim = c(fData$m * fData$n, K, G))
      alpha_g <<- matrix(NA, G)
    },

    init_hlp = function(try_algo) {
      "
        initialize the Hidden Logistic Process
      "
      nm <- fData$m * fData$n
      if (try_algo == 1) {
        for (g in (1:G)) {
          problik <- multinomialLogit(Wg[, , g], phi$Xw, ones(nrow(phi$Xw), ncol(Wg[, , g]) + 1), ones(nrow(phi$Xw), 1))
          pi_jgk[, , g] <<- problik$piik
        }
      }
      else{
        for (g in (1:G)) {
          Wg[, , g] <<- rand(q + 1, K - 1)
          # random initialization of parameter vector for IRLS
          problik <- multinomialLogit(Wg[, , g], phi$Xw, ones(nrow(phi$Xw), ncol(Wg[, , g]) + 1), ones(nrow(phi$Xw), 1))
          pi_jgk[, , g] <<- problik$piik
        }
      }
    },

    initParam = function(init_kmeans, try_algo) {
      # 1. Initialization of cluster weights
      alpha_g <<- 1 / (G * ones(G, 1))
      # 2. Initialization of the model parameters for each cluster: W (pi_jgk), betak and sigmak
      init_hlp(try_algo) # setting Wg and pi_jgk
      # betagk and sigmagk
      if (init_kmeans) {
        # run k means original R
        # kmeans_res <- kmeans(mixModel$X, iter.max = 400, centers=mixModel$G, nstart=20, trace=FALSE)
        # klas <- kmeans_res$cluster

        # run myKmeans
        kmeans_res <- kmeans(fData$Y, G, nbr_runs = 20, nbr_iter_max = 400, verbose = FALSE)
        klas <- kmeans_res$klas
        for (g in 1:G) {
          Xg <- fData$Y[klas == g, ]
          initRegressionParam(Xg, g, try_algo)
        }
      }
      else{
        ind <- sample(fData$n)
        for (g in 1:G) {
          if (g < G) {
            Xg <- fData$Y[ind[((g - 1) * round(fData$n / G) + 1):(g * round(fData$n / G))], ]
          }
          else{
            Xg <- fData$Y[ind[((g - 1) * round(fData$n / G) + 1):length(ind)], ]
          }
          initRegressionParam(Xg, g, try_algo)
        }
      }
    },

    initRegressionParam = function(Xg, g, try_algo) {
      "
        Initialize the Regresssion model with Hidden Logistic Process
      "
      n <- nrow(Xg)
      m <- ncol(Xg)
      if (try_algo == 1) {
        # cutting the sample (signal) in K segments
        zi <- round(m / K) - 1

        beta_k <- matrix(NA, p + 1, K)
        sigma <- c()

        for (k in 1:K) {
          i <- (k - 1) * zi + 1
          j <- k * zi
          Xij <- Xg[, i:j]
          Xij <- matrix(t(Xij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, ]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Xij
          beta_k[, k] <- bk

          if (variance_type == variance_types$homoskedastic) {
            sigma <- var(Xij)
          }
          else{
            mk <- j - i + 1 #length(Xij);
            z <- Xij - Phi_ij %*% bk

            sk <- t(z) %*% z / (n * mk)

            sigma[k] <- sk

          }
        }
      }
      else{
        # random initialization
        Lmin <- round(m / K) #nbr pts min into one segment
        tk_init <- zeros(1, K + 1)
        K_1 <- K
        for (k in 2:K) {
          K_1 <- K_1 - 1

          #temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))

          temp <-
            tk_init[k - 1] + Lmin:(m - (K_1 * Lmin) - tk_init[k - 1])

          ind <- sample(length(temp))

          tk_init[k] <- temp[ind[1]]
        }
        tk_init[K + 1] <- m
        beta_k <- matrix(NA, p + 1, K)
        sigma <- c()
        for (k in 1:K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]
          Xij <- Xg[, i:j]
          Xij <- matrix(t(Xij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, ]
          Phi_ij <- repmat(phi_ij, n, 1)


          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Xij
          beta_k[, k] <- bk

          if (variance_type == variance_types$homoskedastic) {
            sigma <- var(Xij)
          }
          else{
            mk <- j - i + 1 #length(Xij);
            z <- Xij - Phi_ij %*% bk

            sk <- t(z) %*% z / (n * mk)

            sigma[k] <- sk

          }
        }
      }

      betag[, , g] <<- beta_k
      if (variance_type == variance_types$homoskedastic) {
        sigmag[g] <<- sigma
      }
      else{
        sigmag[, g] <<- sigma
      }
    },

    CMStep = function(mixStats) {
      good_segmentation = TRUE
      #MStep for CEM algorithm
      alpha_g <<- t(colSums(mixStats$c_ig)) / fData$n
      # Maximization w.r.t betagk et sigmagk
      cluster_labels <- t(repmat(mixStats$klas, 1, fData$m)) # [m x n]
      cluster_labels <- as.vector(cluster_labels)

      for (g in 1:G) {
        Xg = fData$vecY[cluster_labels == g, ] # cluster g (found from a hard clustering)
        tauijk <- mixStats$tau_ijgk[cluster_labels == g, , g] #[(ng xm) x K]
        if (!is.matrix(tauijk)) {
          tauijk <- matrix(tauijk)
        }

        if (variance_type == variance_types$homoskedastic) {
          s <- 0
        }
        else{
          sigma_gk <- zeros(K, 1)
        }

        beta_gk <- matrix(NA, p + 1, K)

        for (k in 1:K) {
          segments_weights <- tauijk[, k] # poids du kieme segment   pour le cluster g
          # poids pour avoir K segments floues du gieme cluster flou
          phigk <- (sqrt(segments_weights) %*% ones(1, p + 1)) * phi$XBeta[cluster_labels ==
                                                                                   g, ] #[(ng*m)*(p+1)]
          Xgk <- sqrt(segments_weights) * Xg

          # maximization w.r.t beta_gk: Weighted least squares
          beta_gk[, k] <- solve(t(phigk) %*% phigk + .Machine$double.eps * diag(p + 1)) %*% t(phigk) %*% Xgk # Maximization w.r.t betagk

          #    the same as
          #                 W_gk = diag(cluster_weights.*segment_weights);
          #                 beta_gk(:,k) = inv(phiBeta'*W_gk*phiBeta)*phiBeta'*W_gk*X;
          #   Maximization w.r.t au sigma_gk :
          if (variance_type == variance_types$homoskedastic) {
            sk <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2)
            s <- s + sk
            sigma_gk <- s / sum(tauijk)
          }
          else{
            sigma_gk[k] <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2) / (sum(segments_weights))
            if ((sum(segments_weights) == 0)) {
              good_segmentation = FALSE
              return(list(0, good_segmentation))
            }
          }
        }

        betag[, , g] <<- beta_gk
        if (variance_type == variance_types$homoskedastic) {
          sigmag[g] <<- sigma_gk
        }
        else{
          sigmag[, g] <<- sigma_gk

        }

        # Maximization w.r.t W
        #  IRLS : Regression logistique multinomiale ponderee par cluster
        # setting of Wg[,,g] and pi_jgk
        Wg_init <- Wg[, , g]
        if (!is.matrix(Wg_init)) {
          Wg_init <- matrix(Wg_init)
        }

        res_irls <- IRLS(phi$Xw, tauijk, cluster_weights, Wg_init, verbose_IRLS)

        Wg[, , g] <<- res_irls$W
        piik <- res_irls$piik
        pi_jgk[, , g] <<- repmat(piik[1:fData$m, ], fData$n, 1)
        reg_irls <- res_irls$reg_irls
      }
      return(list(reg_irls, good_segmentation))
    },

    MStep = function(mixStats, verbose_IRLS) {
      alpha_g <<- t(colSums(mixStats$h_ig)) / fData$n
      for (g in 1:G) {
        temp <- repmat(mixStats$h_ig[, g], 1, fData$m) # [m x n]
        cluster_weights <- matrix(t(temp), fData$m * fData$n, 1) # cluster_weights(:) [mn x 1]
        tauijk <- mixStats$tau_ijgk[, , g] #[(nxm) x K]
        if (!is.matrix(tauijk)) {
          tauijk <- matrix(tauijk)
        }
        if (variance_type == variance_types$homoskedastic) {
          s <- 0
        }
        else{
          sigma_gk <- zeros(K, 1)
        }

        beta_gk <- matrix(NA, p + 1, K)


        for (k in 1:K) {
          # foreach (k = 1:mixModel$K, .combine = rbind) %dopar%{
          segments_weights <- tauijk[, k] # poids du kieme segment   pour le cluster g
          # poids pour avoir K segments floues du gieme cluster flou
          phigk <- (sqrt(cluster_weights * segments_weights) %*% ones(1, p + 1)) * phi$XBeta #[(n*m)*(p+1)]
          Xgk <- sqrt(cluster_weights * segments_weights) * fData$vecY

          # maximization w.r.t beta_gk: Weighted least squares
          beta_gk[, k] <- solve(t(phigk) %*% phigk + .Machine$double.eps * diag(p + 1)) %*% t(phigk) %*% Xgk # Maximization w.r.t betagk

          #    the same as
          #                 W_gk = diag(cluster_weights.*segment_weights);
          #                 beta_gk(:,k) = inv(phiBeta'*W_gk*phiBeta)*phiBeta'*W_gk*X;
          #   Maximization w.r.t au sigma_gk :
          if (variance_type == variance_types$homoskedastic) {
            sk <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2)
            s <- s + sk
            sigma_gk <- s / sum(colSums((cluster_weights %*% ones(1, K)) * tauijk))
          }
          else{
            sigma_gk[k] <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2) / (colSums(cluster_weights * segments_weights))
          }
        }


        betag[, , g] <<- beta_gk
        if (variance_type == variance_types$homoskedastic) {
          sigmag[g] <<- sigma_gk
        }
        else{
          sigmag[, g] <<- sigma_gk

        }


        # Maximization w.r.t W
        #  IRLS : Regression logistique multinomiale pond??r??e par cluster
        # setting of Wg[,,g] and pi_jgk
        Wg_init <- Wg[, , g]

        if (!is.matrix(Wg_init)) {
          Wg_init <- matrix(Wg_init)
        }

        res_irls <- IRLS(phi$Xw, tauijk, cluster_weights, Wg_init, verbose_IRLS)

        Wg[, , g] <<- res_irls$W
        piik <- res_irls$piik
        pi_jgk[, , g] <<- repmat(piik[1:fData$m, ], fData$n, 1)

      }
    }
  )
)
