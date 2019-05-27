
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

<!-- badges: start -->

<!-- badges: end -->

R code for the **clustering** and **segmentation** of time series
(including with regime changes) by mixture of Hidden Logistic Processes
(MixRHLP) and the EM algorithm; i.e functional data clustering and
segmentation.

## Installation

You can install the development version of mixRHLP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixRHLP")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixRHLP", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("mixRHLP")
```

## Usage

``` r
library(mixRHLP)

data("simulatedtimeseries")
fData <- FData$new()
fData$setData(simulatedtimeseries$X, t(simulatedtimeseries[, 2:ncol(simulatedtimeseries)]))

G <- 3 # number of clusters
K <- 3 # number of regimes (polynomial regression components)
p <- 1 # degree of the polynomials
q <- 1 # order of the logistic regression (by default 1 for contiguous segmentation)
variance_type <- variance_types$heteroskedastic
modelMixRHLP <- ModelMixRHLP(fData, G, K, p, q, variance_type)

n_tries <- 1
max_iter <- 1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE

solution <- EM(modelMixRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS, init_kmeans)

solution$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-4.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-5.png" style="display: block; margin: auto;" />
