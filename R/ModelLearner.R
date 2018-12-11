source("R/enums.R")
source("R/utils.R")
source("R/MixParam.R")
source("R/MixStats.R")
source("R/Phi.R")

EM <- function(mixModel, modelOptions){
  phi <- Phi$new()
  phi$setPhiN(mixModel$t,mixModel$p,mixModel$q, mixModel$n)

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf
  cpu_time_all <- list()

  while(try_EM < modelOptions$n_tries){
    try_EM <- try_EM+1
    message("EM try nr ",try_EM)
    time <- Sys.time()

    # Initializations
    mixParam <- MixParam(mixModel, modelOptions)
    mixParam$initParam(mixModel, phi, modelOptions, try_EM)

    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    mixStats <- MixStats(mixModel, modelOptions)

    while(!converge && (iter <= modelOptions$max_iter)){
      mixStats$EStep(mixModel, mixParam, phi, modelOptions$variance_type)

      mixParam$MStep(mixModel, mixStats, phi, modelOptions)

      # FIN EM

      iter <- iter + 1
      if (modelOptions$verbose){
        message("EM     : Iteration : ",iter, "  log-likelihood : "  , mixStats$log_lik)
      }
      if (prev_loglik - mixStats$log_lik > 1e-5){
        message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ",  mixStats$log_lik)
        top <- top + 1
        if (top > 20) break
      }

      # TEST OF CONVERGENCE
      converge <- abs((mixStats$log_lik - prev_loglik)/prev_loglik) <= modelOptions$threshold
      if (is.na(converge)) {converge <- FALSE}

      prev_loglik <- mixStats$log_lik
      mixStats$stored_loglik[iter] <- mixStats$log_lik
    }# FIN EM LOOP

    cpu_time_all[try_EM] <- Sys.time()-time

    # at this point we have computed param and mixStats that contains all the information

    if (mixStats$log_lik > best_loglik){
      mixStatsSolution <- mixStats$copy()
      mixParamSolution <- mixParam$copy()
      mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:mixModel$m,,]
      best_loglik <- mixStats$log_lik
    }
    if (modelOptions$n_tries > 1){
      message("max value: ", mixStats$log_lik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  mixStatsSolution$MAP()

  if (modelOptions$n_tries > 1){
    message("max value: ", mixStatsSolution$log_lik)
  }


  # FINISH computation of mixStatsSolution
  mixStatsSolution$computeStats(mixModel, mixParam, phi, cpu_time_all)

  return(list(mixParamSolution, mixStatsSolution))
}



