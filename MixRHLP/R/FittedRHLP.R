FittedRHLP <- setRefClass(
  "FittedRHLP",
  fields = list(
    modelRHLP = "ModelRHLP",
    paramRHLP = "ParamRHLP",
    statRHLP = "StatRHLP"
  ),
  methods = list(
    plot = function() {
      par(mfrow = c(2, 1))
      plot.default(modelRHLP$Y, type = "l", ylab = "y", xlab = "")
      title(main = "Time series, RHLP regimes and process probabilities")
      colors = rainbow(modelRHLP$K)
      for (k in 1:modelRHLP$K) {
        index <- (statRHLP$klas == k)
        polynomials <- statRHLP$polynomials[(statRHLP$klas == k), k]
        lines(statRHLP$polynomials[, k], lty = "dotted", lwd = 2, col = colors[k])
        lines(seq(1:modelRHLP$m)[index], polynomials, lwd = 2, col = colors[k])
      }

      plot.default(statRHLP$piik[, 1], type = "l", lwd = 2, col = colors[1], xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)))
      if (modelRHLP$K > 1) {
        for (k in 2:modelRHLP$K) {
          lines(statRHLP$piik[, k], type = "l", lwd = 2, col = colors[k])
        }
      }

      par(mfrow = c(2, 1))
      plot.default(modelRHLP$Y, type = "l", ylab = "y", xlab = "")
      title(main = "Time series, estimated RHLP model, and segmentation")

      tk <- which(diff(statRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], lty = "dotted", lwd = 2, col = "red")
      }
      plot.default(statRHLP$klas, type = "l", lwd = 2, col = "red", xlab = "", ylab = "Estimated class labels")
    }
  )
)

FittedRHLP <- function(modelRHLP, paramRHLP, statRHLP) {
  new("FittedRHLP", modelRHLP = modelRHLP, paramRHLP = paramRHLP, statRHLP = statRHLP)
}
