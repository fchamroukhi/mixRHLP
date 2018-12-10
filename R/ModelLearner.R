source("R/utils.R")


ModelLearner <- setRefClass(
  "ModelLearner",
  fields = list(
    modelOptions="ModelOptions",
    mixModel = "MixModel"
  ),

  methods = list(
    EM = function(){
      phi <- Phi$new()
      phi$setPhiN(mixModel$t,mixModel$p,mixModel$q, mixModel$n)

      top <- 0
      try_EM <- 0
      best_loglik <- 0
      cpu_time <- list()

      while(try_EM < modelOptions$n_tries){
        try_EM <- try_EM+1
        message("EM try nr "+try_EM)
        time <- Sys.time()

        # Initializations
        mixParam <- MixParam(mixModel, mixOptions)
        mixParam$initParam(mixModel, phi, mixOptions, try_EM)

        iter <- 0
        converge <- FALSE
        prev_loglik <- -Inf

        mixStats <- MixStats(mixModel, mixOptions)

        while(!converge && (iter <= modelOptions$max_iter)){

          mixStats$EStep(mixModel, mixParam, modelOptions$variance_type)

          mixStats$MStep()
          # FIN EM

          iter <- iter + 1
          if (mixOptions$verbose){
            message("EM     : Iteration : ",iter, "  log-likelihood : "  , mixStats$log_lik)
          }
          if (prev_loglik - mixStats$log_lik > 1e-5){
            message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ",  mixStats$log_lik)
            top <- top + 1
            if (top > 20) break
          }

          # TEST OF CONVERGENCE
          converge <- abs((mixStats$log_lik - prev_loglik)/prev_loglik) <= mixOptions$threshold
          prev_loglik <- mixStats$log_lik
          mixStats$stored_loglik[iter] <- mixStats$log_lik
        }# FIN EM LOOP

        cpu_time[try_EM] <- Sys.time()-time

        # at this point we have computed param and mixStats that contains all the information

        if (mixStats$log_lik > best_loglik){
          mixStatsSolution <- mixStats$copy()
          mixParamSolution <- param$copy()
          mixParamSolution$pi_jgk <- param$pi_jgk[1:mixModel$m,,]
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


      # todo: FINISH computation of mixStatsSolution

      return(list(mixParamSolution, mixStatsSolution))
    }

  )
)

ModelLearner<-function(myData, mixModel, mixStats, modelOptions){
  new("ModelLearner",myData=myData, mixModel=mixModel, mixStats=mixStats, modelOptions=modelOptions)
}



