FittedMixRHLP <- setRefClass(
  "FittedMixRHLP",
  fields = list(
    modelMixRHLP = "ModelMixRHLP",
    paramMixRHLP = "ParamMixRHLP",
    statMixRHLP = "StatMixRHLP"
  ),
  methods = list(
    plot = function() {

      #### cluster and means
      couleur = rainbow(modelMixRHLP$G)
      par(mfrow=c(round(sqrt(modelMixRHLP$G + 1)), round(sqrt(modelMixRHLP$G + 1))))

      matplot(t(modelMixRHLP$Y), col = "black", type = "l", lty = "solid", xlab = "Time", ylab = "y")
      title(main = "Initial datset")

      for (g in 1:modelMixRHLP$G) {
        cluster_g = modelMixRHLP$Y[statMixRHLP$klas == g,]
        matplot(t(cluster_g), col = couleur[g], type = "l", lty = "dotted", xlab = "Time", ylab = "y")
        if (is.null(dim(statMixRHLP$Ex_g))) {
          lines(statMixRHLP$Ex_g, lty = "solid", lwd = 2, col = "black")
        } else {
          lines(statMixRHLP$Ex_g[, g], lty = "solid", lwd = 2, col = "black")
        }
        title(main = sprintf("Cluster %1.1i", g))
      }

      par(mfrow = c(2, 1))
      for (g in 1:modelMixRHLP$G) {

        # First graph
        cluster_g = modelMixRHLP$Y[statMixRHLP$klas == g,]
        matplot(t(cluster_g), col = couleur[g], type = "l", lty = "dotted", xlab = "Time", ylab = "y")
        # polynomial regressors
        for (i in 1:ncol(statMixRHLP$polynomials[, , g])) {
          lines(statMixRHLP$polynomials[, i, g], lty = "dotted", lwd = 2, col = "black")
        }
        if (is.null(dim(statMixRHLP$Ex_g))) {
          lines(statMixRHLP$Ex_g, lty = "solid", lwd = 2, col = "black")
        } else {
          lines(statMixRHLP$Ex_g[, g], lty = "solid", lwd = 2, col = "black")
        }
        title(main = sprintf("Cluster %1.1i", g))

        matplot(paramMixRHLP$pi_jgk[1:modelMixRHLP$m, , g], type = "l", lty = "solid",
                xlab = "Time", ylab = "Logistic proportions") # ['\pi_{jk}( w^g) , g = ',int2str(g)])

      }

      par(mfrow = c(1, 1))

      plot.default(unlist(statMixRHLP$stored_loglik), type = "l", col = "blue", xlab = "Iteration", ylab = "Log Likelihood")
      title(main = "Log-Likelihood")

    }
  )
)

FittedMixRHLP <- function(modelMixRHLP, paramMixRHLP, statMixRHLP) {
  new("FittedMixRHLP", modelMixRHLP = modelMixRHLP, paramMixRHLP = paramMixRHLP, statMixRHLP = statMixRHLP)
}
