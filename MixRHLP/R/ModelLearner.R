source("R/enums.R")
source("R/utils.R")
source("R/ParamMixRHLP.R")
source("R/StatMixRHLP.R")
source("R/FittedMixRHLP.R")

################################################
# The EM algorithm for the MixFRHLP model
################################################

EM <- function(modelMixRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans){
  phi <- designmatrix(modelMixRHLP$X, modelMixRHLP$p, modelMixRHLP$q, modelMixRHLP$n)

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf
  cpu_time_all <- c()

  while(try_EM < n_tries){
    try_EM <- try_EM+1
    message("EM try nr ",try_EM)
    time <- Sys.time()

    # Initialization
    mixParam <- ParamMixRHLP(modelMixRHLP)
    mixParam$initParam(modelMixRHLP, phi, init_kmeans, try_EM)

    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    mixStats <- StatMixRHLP(modelMixRHLP)

    while(!converge && (iter <= max_iter)){
      mixStats$EStep(modelMixRHLP, mixParam, phi)

      mixParam$MStep(modelMixRHLP, mixStats, phi, verbose_IRLS)
      # FIN EM

      iter <- iter + 1
      if (verbose){
        message("EM     : Iteration : ",iter, "  log-likelihood : "  , mixStats$log_lik)
      }
      if (prev_loglik - mixStats$log_lik > 1e-5){
        message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ",  mixStats$log_lik)
        top <- top + 1
        if (top > 20) break
      }

      # TEST OF CONVERGENCE
      converge <- abs((mixStats$log_lik - prev_loglik)/prev_loglik) <= threshold
      if (is.na(converge)) {converge <- FALSE}

      prev_loglik <- mixStats$log_lik
      mixStats$stored_loglik[iter] <- mixStats$log_lik
    }# End of the EM LOOP

    cpu_time_all[try_EM] <- Sys.time()-time


    if (mixStats$log_lik > best_loglik){
      mixStatsSolution <- mixStats$copy()
      mixParamSolution <- mixParam$copy()

      if (modelMixRHLP$K==1 && modelMixRHLP$G==1){
        mixParamSolution$pi_jgk <- matrix(mixParam$pi_jgk, nrow = modelMixRHLP$m, ncol = 1)
      }
      else{
        mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:modelMixRHLP$m,,]
      }

      best_loglik <- mixStats$log_lik
    }
    if (n_tries > 1){
      message("max value: ", mixStats$log_lik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  mixStatsSolution$MAP()

  if (n_tries > 1){
    message("max value: ", mixStatsSolution$log_lik)
  }


  mixStatsSolution$computeStats(modelMixRHLP, mixParamSolution, phi, cpu_time_all)

  return(FittedMixRHLP(modelMixRHLP, mixParamSolution, mixStatsSolution))
}


################################################
####                  CEM algorithm
################################################

CEM <- function(modelMixRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans){
  phi <- designmatrix(modelMixRHLP$X, modelMixRHLP$p, modelMixRHLP$q, modelMixRHLP$n)

  top <- 0
  try_CEM <- 0
  best_com_loglik <- -Inf
  cpu_time_all <- c()

  while(try_CEM < n_tries){
    try_CEM <- try_CEM+1
    message("CEM try nr ",try_CEM)
    time <- Sys.time()

    # Initialization
    mixParam <- ParamMixRHLP(modelMixRHLP)
    mixParam$initParam(modelMixRHLP, phi, init_kmeans, try_CEM)

    iter <- 0
    converge <- FALSE
    prev_com_loglik <- -Inf
    reg_irls <- 0
    mixStats <- StatMixRHLP(modelMixRHLP)

    while(!converge && (iter <= max_iter)){
      mixStats$EStep(modelMixRHLP, mixParam, phi)
      mixStats$CStep(reg_irls)
      res <- mixParam$CMStep(modelMixRHLP, mixStats, phi)
      reg_irls = res[[1]]
      good_segmentation = res[[2]]
      if (good_segmentation==FALSE){
        try_CEM <- try_CEM-1
        break # try one more time CEM
      }

      #

      iter <- iter + 1
      if (verbose){
        message("CEM     : Iteration : ",iter, "  complete log-likelihood : "  , mixStats$com_loglik)
      }
      if (prev_com_loglik - mixStats$com_loglik > 1e-5){
        message("!!!!! CEM complete log-likelihood is decreasing from ", prev_loglik, "to ",  mixStats$com_loglik)
        top <- top + 1
        if (top > 20) break
      }

      # TEST OF CONVERGENCE
      converge <- abs((mixStats$com_loglik - prev_com_loglik)/prev_com_loglik) <= threshold
      if (is.na(converge)) {converge <- FALSE}

      prev_com_loglik <- mixStats$com_loglik
      mixStats$stored_loglik[iter] <- mixStats$com_loglik
    }# End of the CEM LOOP

    cpu_time_all[try_CEM] <- Sys.time()-time


    if (mixStats$com_loglik > best_com_loglik){
      mixStatsSolution <- mixStats$copy()
      mixParamSolution <- mixParam$copy()

      if (modelMixRHLP$K==1 && modelMixRHLP$G==1){
        mixParamSolution$pi_jgk <- matrix(mixParam$pi_jgk, nrow = modelMixRHLP$m, ncol = 1)
      }
      else{
        mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:modelMixRHLP$m,,]
      }

      best_com_loglik <- mixStats$com_loglik
    }
    if (n_tries > 1){
      message("max value: ", mixStats$com_loglik)
    }
  }

  if (n_tries > 1){
    message("max value: ", mixStatsSolution$com_loglik)
  }


  mixStatsSolution$computeStats(modelMixRHLP, mixParamSolution, phi, cpu_time_all)

  return(FittedMixRHLP(modelMixRHLP, mixParamSolution, mixStatsSolution))
}
