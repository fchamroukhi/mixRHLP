#' emMixRHLP is used to fit a MixRHLP model.
#'
#' emMixRHLP is used to fit a MixRHLP model. The estimation method is performed by
#' the Expectation-Maximization algorithm.
#'
#' @details emMixRHLP function is based on the EM algorithm. This functions starts
#' with an initialization of the parameters done by the method `initParam` of
#' the class [ParamMixRHLP][ParamMixRHLP], then it alternates between a E-Step
#' (method of the class [StatMixRHLP][StatMixRHLP]) and a M-Step (method of the class
#' [ParamMixRHLP][ParamMixRHLP]) until convergence (until the absolute difference of
#' log-likelihood between two steps of the EM algorithm is less than the
#' `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates.
#' @param Y Matrix of size \eqn{(n, m)} representing \emph{n} functions of `X`
#' observed at points \eqn{1,\dots,m}.
#' @param G The number of clusters.
#' @param K The number of regimes (mixture components).
#' @param p The order of the polynomial regression.
#' @param q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @param variance_type Optional character indicating if the model is
#' "homoskedastic" or "heteroskedastic". By default the model is
#' "heteroskedastic".
#' @param n_tries Number of times EM algorithm will be launched.
#' The solution providing the highest log-likelihood will be returned.
#'
#' If `n_tries` > 1, then for the first pass, parameters are initialized
#' by uniformly segmenting the data into K segments, and for the next passes,
#' parameters are initialized by randomly segmenting the data into K contiguous
#'  segments.
#' @param max_iter The maximum number of iterations for the EM algorithm.
#' @param threshold A numeric value specifying the threshold for the relative
#'  difference of log-likelihood between two steps  of the EM as stopping
#'  criteria.
#' @param verbose A logical value indicating whether values of the
#' log-likelihood should be printed during EM iterations.
#' @param verbose_IRLS A logical value indicating whether values of the
#' criterion optimized by IRLS should be printed at each step of the EM
#' algorithm.
#' @param init_kmeans A logical value indicating wheather the initialization of
#' the model parameters is performed by kmeans.
#' @return EM returns an object of class [ModelMixRHLP][ModelMixRHLP].
#' @seealso [ModelMixRHLP], [ParamMixRHLP], [StatMixRHLP]
#' @export
emMixRHLP <-
  function(X, Y, G, K, p = 3, q = 1, variance_type = c("heteroskedastic", "homoskedastic"), n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans) {

    fData <- FData(X, Y)

    top <- 0
    try_EM <- 0
    best_loglik <- -Inf
    cpu_time_all <- c()

    while (try_EM < n_tries) {
      try_EM <- try_EM + 1
      message("EM try nr ", try_EM)
      time <- Sys.time()

      # Initialization
      variance_type <- match.arg(variance_type)
      mixParam <- ParamMixRHLP$new(fData = fData, G = G, K = K, p = p, q = q, variance_type = variance_type)
      mixParam$initParam(init_kmeans, try_EM)

      iter <- 0
      converge <- FALSE
      prev_loglik <- -Inf

      mixStats <- StatMixRHLP(mixParam)

      while (!converge && (iter <= max_iter)) {
        mixStats$EStep(mixParam)

        mixParam$MStep(mixStats, verbose_IRLS)
        # FIN EM

        iter <- iter + 1
        if (verbose) {
          message("EM     : Iteration : ", iter, "  log-likelihood : "  , mixStats$log_lik)
        }
        if (prev_loglik - mixStats$log_lik > 1e-5) {
          message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ", mixStats$log_lik)
          top <- top + 1
          if (top > 20)
            break
        }

        # TEST OF CONVERGENCE
        converge <- abs((mixStats$log_lik - prev_loglik) / prev_loglik) <= threshold
        if (is.na(converge)) {
          converge <- FALSE
        }

        prev_loglik <- mixStats$log_lik
        mixStats$stored_loglik[iter] <- mixStats$log_lik
      }# End of the EM LOOP

      cpu_time_all[try_EM] <- Sys.time() - time


      if (mixStats$log_lik > best_loglik) {
        mixStatsSolution <- mixStats$copy()
        mixParamSolution <- mixParam$copy()

        if (mixParam$K == 1 && mixParam$G == 1) {
          mixParamSolution$pi_jgk <-
            matrix(mixParam$pi_jgk, nrow = mixParam$fData$m, ncol = 1)
        }
        else{
          mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:mixParam$fData$m, , ]
        }

        best_loglik <- mixStats$log_lik
      }
      if (n_tries > 1) {
        message("max value: ", mixStats$log_lik)
      }
    }

    # Computation of c_ig the hard partition of the curves and klas
    mixStatsSolution$MAP()

    if (n_tries > 1) {
      message("max value: ", mixStatsSolution$log_lik)
    }


    mixStatsSolution$computeStats(mixParamSolution, cpu_time_all)

    return(ModelMixRHLP(paramMixRHLP = mixParamSolution, statMixRHLP = mixStatsSolution))
  }


#' cemMixRHLP is used to fit a MixRHLP model.
#'
#' cemMixRHLP is used to fit a MixRHLP model. The estimation method is performed by
#' the Expectation-Maximization algorithm.
#'
#' @details cemMixRHLP function is based on the CEM algorithm. This functions starts
#' with an initialization of the parameters done by the method `initParam` of
#' the class [ParamMixRHLP][ParamMixRHLP], then it alternates between a E-Step
#' (method of the class [StatMixRHLP][StatMixRHLP]), a C-Step and a M-Step (methods of the class
#' [ParamMixRHLP][ParamMixRHLP]) until convergence (until the absolute difference of
#' log-likelihood between two steps of the EM algorithm is less than the
#' `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates.
#' @param Y Matrix of size \eqn{(n, m)} representing \emph{n} functions of `X`
#' observed at points \eqn{1,\dots,m}.
#' @param G The number of clusters.
#' @param K The number of regimes (mixture components).
#' @param p The order of the polynomial regression.
#' @param q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @param variance_type Numeric indicating if the model is homoskedastic
#' (`variance_type` = 1) or heteroskedastic (`variance_type` = 2).
#' @param n_tries Number of times EM algorithm will be launched.
#' The solution providing the highest log-likelihood will be returned.
#'
#' If `n_tries` > 1, then for the first pass, parameters are initialized
#' by uniformly segmenting the data into K segments, and for the next passes,
#' parameters are initialized by randomly segmenting the data into K contiguous
#'  segments.
#' @param max_iter The maximum number of iterations for the EM algorithm.
#' @param threshold A numeric value specifying the threshold for the relative
#'  difference of log-likelihood between two steps  of the EM as stopping
#'  criteria.
#' @param verbose A logical value indicating whether values of the
#' log-likelihood should be printed during EM iterations.
#' @param verbose_IRLS A logical value indicating whether values of the
#' criterion optimized by IRLS should be printed at each step of the EM
#' algorithm.
#' @param init_kmeans A logical value indicating wheather the initialization of
#' the model parameters is performed by kmeans.
#' @return EM returns an object of class [ModelMixRHLP][ModelMixRHLP].
#' @seealso [ModelMixRHLP], [ParamMixRHLP], [StatMixRHLP]
#' @export
CEM <- function(X, Y, G, K, p, q = 1, variance_type = 2, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans) {

    fData <- FData(X, Y)

    top <- 0
    try_CEM <- 0
    best_com_loglik <- -Inf
    cpu_time_all <- c()

    while (try_CEM < n_tries) {
      try_CEM <- try_CEM + 1
      message("CEM try nr ", try_CEM)
      time <- Sys.time()

      # Initialization
      mixParam <- ParamMixRHLP$new(fData = fData, G = G, K = K, p = p, q = q, variance_type = variance_type)
      mixParam$initParam(init_kmeans, try_CEM)

      iter <- 0
      converge <- FALSE
      prev_com_loglik <- -Inf
      reg_irls <- 0
      mixStats <- StatMixRHLP(mixParam)

      while (!converge && (iter <= max_iter)) {
        mixStats$EStep(mixParam)
        mixStats$CStep(reg_irls)
        res <- mixParam$CMStep(mixStats)
        reg_irls = res[[1]]
        good_segmentation = res[[2]]
        if (good_segmentation == FALSE) {
          try_CEM <- try_CEM - 1
          break # try one more time CEM
        }

        #

        iter <- iter + 1
        if (verbose) {
          message("CEM     : Iteration : ", iter, "  complete log-likelihood : "  , mixStats$com_loglik)
        }
        if (prev_com_loglik - mixStats$com_loglik > 1e-5) {
          message("!!!!! CEM complete log-likelihood is decreasing from ", prev_com_loglik, "to ", mixStats$com_loglik)
          top <- top + 1
          if (top > 20)
            break
        }

        # TEST OF CONVERGENCE
        converge <- abs((mixStats$com_loglik - prev_com_loglik) / prev_com_loglik) <= threshold
        if (is.na(converge)) {
          converge <- FALSE
        }

        prev_com_loglik <- mixStats$com_loglik
        mixStats$stored_loglik[iter] <- mixStats$com_loglik
      }# End of the CEM LOOP

      cpu_time_all[try_CEM] <- Sys.time() - time


      if (mixStats$com_loglik > best_com_loglik) {
        mixStatsSolution <- mixStats$copy()
        mixParamSolution <- mixParam$copy()

        if (mixParam$K == 1 && mixParam$G == 1) {
          mixParamSolution$pi_jgk <- matrix(mixParam$pi_jgk, nrow = mixParam$fData$m, ncol = 1)
        }
        else{
          mixParamSolution$pi_jgk <- mixParam$pi_jgk[1:mixParam$fData$m, , ]
        }

        best_com_loglik <- mixStats$com_loglik
      }
      if (n_tries > 1) {
        message("max value: ", mixStats$com_loglik)
      }
    }

    if (n_tries > 1) {
      message("max value: ", mixStatsSolution$com_loglik)
    }


    mixStatsSolution$computeStats(mixParamSolution, cpu_time_all)

    return(ModelMixRHLP(paramMixRHLP = mixParamSolution, statMixRHLP = mixStatsSolution))
  }
