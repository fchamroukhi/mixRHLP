#' A Reference Class which represents a fitted MixRHLP model.
#'
#' ModelMixRHLP represents a [MixRHLP][ModelMixRHLP] model for which parameters have
#' been estimated.
#'
#' @usage NULL
#' @field paramMixRHLP A [ParamMixRHLP][ParamMixRHLP] object. It contains the estimated values of the parameters.
#' @field statMixRHLP A [StatMixRHLP][StatMixRHLP] object. It contains all the statistics associated to the MixRHLP model.
#' @seealso [ParamMixRHLP], [StatMixRHLP]
#' @export
ModelMixRHLP <- setRefClass(
  "ModelMixRHLP",
  fields = list(
    paramMixRHLP = "ParamMixRHLP",
    statMixRHLP = "StatMixRHLP"
  ),
  methods = list(
    plot = function(what = c("meancurve", "regressors", "loglikelihood")) {

      what <- match.arg(what, several.ok = TRUE)

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(modelMixRHLP$Y) - 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))), max(modelMixRHLP$Y) + 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))))

      colorsvector = rainbow(paramMixRHLP$G)

      if (any(what == "meancurve")) {
        # Cluster and means
        par(mfrow = c(ceiling(sqrt(paramMixRHLP$G + 1)), round(sqrt(paramMixRHLP$G + 1))), mai = c(0.6, 0.6, 0.5, 0.25), mgp = c(2, 1, 0))

        matplot(paramMixRHLP$fData$X, t(paramMixRHLP$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y")
        title(main = "Dataset")

        for (g in 1:paramMixRHLP$G) {
          cluster_g = as.matrix(paramMixRHLP$fData$Y[statMixRHLP$klas == g, ])

          if (length(cluster_g) != 0) {
            dim(cluster_g) <- c(sum(mixrhlp$statMixRHLP$klas == g), paramMixRHLP$fData$m)
            matplot(paramMixRHLP$fData$X, t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "x", ylab = "y")
            lines(paramMixRHLP$fData$X, as.matrix(statMixRHLP$Ex)[, g], col = "black", lty = "solid", lwd = 1.5)
            title(main = sprintf("Cluster %1.1i", g))
          }
        }
      }

      if (any(what == "regressors")) {
        par(mfrow = c(2, 1), mai = c(0.6, 0.8, 0.5, 0.5))
        for (g in 1:paramMixRHLP$G) {
          cluster_g = as.matrix(paramMixRHLP$fData$Y[statMixRHLP$klas == g, ])

          if (length(cluster_g) != 0) {
            dim(cluster_g) <- c(sum(mixrhlp$statMixRHLP$klas == g), paramMixRHLP$fData$m)

            matplot(paramMixRHLP$fData$X, t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "x", ylab = "y")

            # Polynomial regressors
            if (paramMixRHLP$K > 0) {
              for (k in 1:paramMixRHLP$K) {
                lines(paramMixRHLP$fData$X, statMixRHLP$polynomials[, k, g], col = "black", lty = "dotted", lwd = 1.5)
              }
            } else {
              lines(paramMixRHLP$fData$X, statMixRHLP$polynomials[, , g], col = "black", lty = "dotted", lwd = 1.5)
            }

            lines(paramMixRHLP$fData$X, as.matrix(statMixRHLP$Ex)[, g], col = "black", lty = "solid", lwd = 1.5)
            title(main = sprintf("Cluster %1.1i", g))

            matplot(paramMixRHLP$fData$X, statMixRHLP$pi_jgk[1:paramMixRHLP$fData$m, , g], type = "l", lty = "solid", xlab = "x", ylab = "Logistic proportions", ylim = c(0, 1))
          }
        }
      }

      if (any(what == "loglikelihood")) {
        par(mfrow = c(1, 1))
        plot.default(unlist(statMixRHLP$stored_loglik), type = "l", col = "blue", xlab = "Iteration", ylab = "Log Likelihood")
        title(main = "Log-Likelihood")
      }

    },

    summary = function() {
      digits = getOption("digits")

      title <- paste("Fitted mixRHLP model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("MixRHLP model with G = ", paramMixRHLP$G, ifelse(paramMixRHLP$G > 1, " clusters", " cluster"), " and K = ", paramMixRHLP$K, ifelse(paramMixRHLP$K > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = statMixRHLP$loglik, "nu" = paramMixRHLP$nu,
                        "AIC" = statMixRHLP$AIC, "BIC" = statMixRHLP$BIC, "ICL" = statMixRHLP$ICL,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table:")
      print(table(statMixRHLP$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(paramMixRHLP$alpha)
      colnames(pro) <- 1:paramMixRHLP$G
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (g in 1:paramMixRHLP$G) {
        cat(txt)
        cat("\nCluster ", g, " (G = ", g, "):\n", sep = "")

        cat("\nRegression coefficients:\n\n")
        if (paramMixRHLP$p > 0) {
          row.names = c("1", sapply(1:paramMixRHLP$p, function(x) paste0("X^", x)))
        } else {
          row.names = "1"
        }

        betas <- data.frame(paramMixRHLP$beta[, , g], row.names = row.names)
        colnames(betas) <- sapply(1:paramMixRHLP$K, function(x) paste0("Beta(K = ", x, ")"))
        print(betas, digits = digits)

        cat(paste0(ifelse(paramMixRHLP$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(paramMixRHLP$sigma2[, g]))
        if (paramMixRHLP$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:paramMixRHLP$K, function(x) paste0("Sigma2(K = ", x, ")"))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }

  )
)
