#' @export
ModelMixRHLP <- setRefClass(
  "ModelMixRHLP",
  fields = list(
    paramMixRHLP = "ParamMixRHLP",
    statMixRHLP = "StatMixRHLP"
  ),
  methods = list(
    plot = function() {

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(modelMixRHLP$Y) - 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))), max(modelMixRHLP$Y) + 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))))

      # Cluster and means
      colorsvector = rainbow(paramMixRHLP$G)
      par(mfrow = c(round(sqrt(paramMixRHLP$G + 1)), round(sqrt(paramMixRHLP$G + 1))), mai = c(0.6, 0.6, 0.5, 0.25), mgp = c(2, 1, 0))

      matplot(t(paramMixRHLP$fData$Y), type = "l", lty = "solid", col = "black", xlab = "Time", ylab = "y(t)")
      title(main = "Original time series")

      for (g in 1:paramMixRHLP$G) {
        cluster_g = paramMixRHLP$fData$Y[statMixRHLP$klas == g,]
        matplot(t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "Time", ylab = "y(t)")
        if (is.null(dim(statMixRHLP$Ex_g))) {
          lines(statMixRHLP$Ex_g, col = "black", lty = "solid", lwd = 1.5)
        } else {
          lines(statMixRHLP$Ex_g[, g], col = "black", lty = "solid", lwd = 1.5)
        }
        title(main = sprintf("Cluster %1.1i", g))
      }

      par(mfrow = c(2, 1), mai = c(0.6, 0.8, 0.5, 0.5))
      for (g in 1:paramMixRHLP$G) {

        cluster_g = paramMixRHLP$fData$Y[statMixRHLP$klas == g,]
        matplot(t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "Time", ylab = "y(t)")
        # Polynomial regressors
        for (i in 1:ncol(statMixRHLP$polynomials[, , g])) {
          lines(statMixRHLP$polynomials[, i, g], col = "black", lty = "dotted", lwd = 1.5)
        }
        if (is.null(dim(statMixRHLP$Ex_g))) {
          lines(statMixRHLP$Ex_g, col = "black", lty = "solid", lwd = 1.5)
        } else {
          lines(statMixRHLP$Ex_g[, g], col = "black", lty = "solid", lwd = 1.5)
        }
        title(main = sprintf("Cluster %1.1i", g))

        matplot(paramMixRHLP$pi_jgk[1:paramMixRHLP$fData$m, , g], type = "l", lty = "solid", xlab = "Time", ylab = "Logistic proportions")

      }

      par(mfrow = c(1, 1))
      plot.default(unlist(statMixRHLP$stored_loglik), type = "l", col = "blue", xlab = "Iteration", ylab = "Log Likelihood")
      title(main = "Log-Likelihood")

    }
  )
)

