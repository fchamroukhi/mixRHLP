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
  cpu_time_all <- c()

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
      if (is.na(converge)) {converge <- FALSE} # basicly for the first iteration when prev_loglik is Inf

      prev_loglik <- mixStats$log_lik
      mixStats$stored_loglik[iter] <- mixStats$log_lik
    }# FIN EM LOOP

    cpu_time_all[try_EM] <- Sys.time()-time

    # at this point we have computed param and mixStats that contains all the information

    if (mixStats$log_lik > best_loglik){
      mixStatsSolution <- mixStats$copy()
      mixParamSolution <- mixParam$copy()

      if (mixModel$K==1 && mixModel$G==1){
        mixParamSolution$pi_jgk <- matrix(mixParam$pi_jgk, nrow = mixModel$m, ncol = 1)
      }
      else{
        mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:mixModel$m,,]
      }

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
  mixStatsSolution$computeStats(mixModel, mixParamSolution, phi, cpu_time_all)

  return(list(mixParamSolution, mixStatsSolution))
}


################################################
####                  CEM
################################################

CEM <- function(mixModel, modelOptions){
  phi <- Phi$new()
  phi$setPhiN(mixModel$t,mixModel$p,mixModel$q, mixModel$n)

  top <- 0
  try_EM <- 0
  best_com_loglik <- -Inf
  cpu_time_all <- c()

  while(try_EM < modelOptions$n_tries){
    try_EM <- try_EM+1
    message("CEM try nr ",try_EM)
    time <- Sys.time()

    # Initializations
    mixParam <- MixParam(mixModel, modelOptions)
    mixParam$initParam(mixModel, phi, modelOptions, try_EM)

    iter <- 0
    converge <- FALSE
    prev_com_loglik <- -Inf
    reg_irls <- 0
    mixStats <- MixStats(mixModel, modelOptions)

    while(!converge && (iter <= modelOptions$max_iter)){
      mixStats$EStep(mixModel, mixParam, phi, modelOptions$variance_type)
      mixStats$CStep(reg_irls)
      reg_irls <- mixParam$CMStep(mixModel, mixStats, phi, modelOptions)

      # FIN EM

      iter <- iter + 1
      if (modelOptions$verbose){
        message("CEM     : Iteration : ",iter, "  complete log-likelihood : "  , mixStats$com_loglik)
      }
      if (prev_com_loglik - mixStats$com_loglik > 1e-5){
        message("!!!!! CEM complete log-likelihood is decreasing from ", prev_loglik, "to ",  mixStats$com_loglik)
        top <- top + 1
        if (top > 20) break
      }

      # TEST OF CONVERGENCE
      converge <- abs((mixStats$com_loglik - prev_com_loglik)/prev_com_loglik) <= modelOptions$threshold
      if (is.na(converge)) {converge <- FALSE} # basicly for the first iteration when prev_com_loglik is Inf

      prev_com_loglik <- mixStats$com_loglik
      mixStats$stored_loglik[iter] <- mixStats$com_loglik
    }# FIN EM LOOP

    cpu_time_all[try_EM] <- Sys.time()-time

    # at this point we have computed param and mixStats that contains all the information

    if (mixStats$com_loglik > best_com_loglik){
      mixStatsSolution <- mixStats$copy()
      mixParamSolution <- mixParam$copy()

      if (mixModel$K==1 && mixModel$G==1){
        mixParamSolution$pi_jgk <- matrix(mixParam$pi_jgk, nrow = mixModel$m, ncol = 1)
      }
      else{
        mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:mixModel$m,,]
      }

      best_com_loglik <- mixStats$com_loglik
    }
    if (modelOptions$n_tries > 1){
      message("max value: ", mixStats$com_loglik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  mixStatsSolution$MAP()

  if (modelOptions$n_tries > 1){
    message("max value: ", mixStatsSolution$com_loglik)
  }


  # FINISH computation of mixStatsSolution
  mixStatsSolution$computeStats(mixModel, mixParamSolution, phi, cpu_time_all)

  return(list(mixParamSolution, mixStatsSolution))
}
