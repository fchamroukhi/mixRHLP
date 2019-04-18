FittedMixRHLP <- setRefClass(
  "FittedMixRHLP",
  fields = list(
    modelMixRHLP = "ModelMixRHLP",
    paramMixRHLP = "ParamMixRHLP",
    statMixRHLP = "StatMixRHLP"
  ),
  methods = list(
    plot = function() {
      t <- seq(0, modelMixRHLP$m - 1)
      colors = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')
      colors_cluster_means = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')

      par(mfrow=c(1,1))
      for (g in 1 : modelMixRHLP$G){
        cluster_g <- modelMixRHLP$Y[statMixRHLP$klas == g,]
        if (g==1){
          plot.default(t, t(cluster_g)[,1], type = 'l', col=colors[g], xlab='Time', ylab = 'y')
        }
        else{
          lines(t, t(cluster_g)[,1], type = 'l', col=colors[g])
        }
        for (i in 1 : nrow(cluster_g)){
          lines(t, t(cluster_g)[,i], col=colors[g])
        }

      }

      for(g in 1 : modelMixRHLP$G){
        lines(t, statMixRHLP$Ex_g[,g], col = colors_cluster_means[g], lwd = 5)
      }

      ##########################################################################
      for (g in 1 : modelMixRHLP$G){
        par(mfrow = c(2,1))
        cluster_g = modelMixRHLP$Y[statMixRHLP$klas == g,]
        plot.default(t, t(cluster_g)[, 1], type = 'l', col=colors[g], ylab = 'y')
        for (k in 1 : modelMixRHLP$K){
          lines(t, statMixRHLP$polynomials[, k, g], lty = 2 , col="black",lwd = 1)
        }
        lines(t, statMixRHLP$Ex_g[, g], col = colors_cluster_means[g], lwd = 5)

        for (k in 1 : modelMixRHLP$K){
          if (k == 1){
            plot.default(t, paramMixRHLP$pi_jgk[1 : modelMixRHLP$m, k, g], xlab = 'Time', ylab = 'Logistic proportions', type = "l", lwd=2, col=colors_cluster_means[k])
          }
          else{
            lines(t, paramMixRHLP$pi_jgk[1 : modelMixRHLP$m, k, g] , lwd = 2, col = colors_cluster_means[k])
          }
        }
      }

    }
  )
)

FittedMixRHLP <- function(modelMixRHLP, paramMixRHLP, statMixRHLP) {
  new("FittedMixRHLP", modelMixRHLP = modelMixRHLP, paramMixRHLP = paramMixRHLP, statMixRHLP = statMixRHLP)
}
