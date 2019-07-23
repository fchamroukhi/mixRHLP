
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

# Overview

R code for the **clustering** and **segmentation** of time series
(including with regime changes) by mixture of Hidden Logistic Processes
(MixRHLP) and the EM algorithm; i.e functional data clustering and
segmentation.

# Installation

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

# Usage

``` r
library(mixRHLP)
```

``` r
# Application to a toy data set
data("toydataset")

G <- 3 # Number of clusters
K <- 3 # Number of regimes (polynomial regression components)
p <- 1 # Degree of the polynomials
q <- 1 # Order of the logistic regression (by default 1 for contiguous segmentation)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

n_tries <- 1
max_iter <- 1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE

mixrhlp <- emMixRHLP(toydataset$x, t(as.matrix(toydataset[,2:ncol(toydataset)])), 
                     G, K, p, q, variance_type, init_kmeans, n_tries, max_iter, 
                     threshold, verbose, verbose_IRLS)
#> EM: Iteration : 1 || log-likelihood : -18129.8169520025
#> EM: Iteration : 2 || log-likelihood : -16642.732267463
#> EM: Iteration : 3 || log-likelihood : -16496.947898833
#> EM: Iteration : 4 || log-likelihood : -16391.6755568235
#> EM: Iteration : 5 || log-likelihood : -16308.151649539
#> EM: Iteration : 6 || log-likelihood : -16242.6749975019
#> EM: Iteration : 7 || log-likelihood : -16187.9951484578
#> EM: Iteration : 8 || log-likelihood : -16138.360050325
#> EM: Iteration : 9 || log-likelihood : -16092.9430959116
#> EM: Iteration : 10 || log-likelihood : -16053.588838999
#> EM: Iteration : 11 || log-likelihood : -16020.7365667916
#> EM: Iteration : 12 || log-likelihood : -15993.7513179937
#> EM: Iteration : 13 || log-likelihood : -15972.7088032469
#> EM: Iteration : 14 || log-likelihood : -15957.3889127412
#> EM: Iteration : 15 || log-likelihood : -15946.5663566082
#> EM: Iteration : 16 || log-likelihood : -15938.693534838
#> EM: Iteration : 17 || log-likelihood : -15932.584112949
#> EM: Iteration : 18 || log-likelihood : -15927.5299507605
#> EM: Iteration : 19 || log-likelihood : -15923.1499635319
#> EM: Iteration : 20 || log-likelihood : -15919.2392546398
#> EM: Iteration : 21 || log-likelihood : -15915.6795793534
#> EM: Iteration : 22 || log-likelihood : -15912.3944381959
#> EM: Iteration : 23 || log-likelihood : -15909.327585346
#> EM: Iteration : 24 || log-likelihood : -15906.4326405988
#> EM: Iteration : 25 || log-likelihood : -15903.6678636145
#> EM: Iteration : 26 || log-likelihood : -15900.9933370165
#> EM: Iteration : 27 || log-likelihood : -15898.3692402859
#> EM: Iteration : 28 || log-likelihood : -15895.7545341827
#> EM: Iteration : 29 || log-likelihood : -15893.1056775993
#> EM: Iteration : 30 || log-likelihood : -15890.3751610539
#> EM: Iteration : 31 || log-likelihood : -15887.5097378815
#> EM: Iteration : 32 || log-likelihood : -15884.4482946475
#> EM: Iteration : 33 || log-likelihood : -15881.1193453446
#> EM: Iteration : 34 || log-likelihood : -15877.4381561224
#> EM: Iteration : 35 || log-likelihood : -15873.3037170772
#> EM: Iteration : 36 || log-likelihood : -15868.595660791
#> EM: Iteration : 37 || log-likelihood : -15863.171868441
#> EM: Iteration : 38 || log-likelihood : -15856.8678694783
#> EM: Iteration : 39 || log-likelihood : -15849.5002500459
#> EM: Iteration : 40 || log-likelihood : -15840.8778843568
#> EM: Iteration : 41 || log-likelihood : -15830.8267303162
#> EM: Iteration : 42 || log-likelihood : -15819.2343887404
#> EM: Iteration : 43 || log-likelihood : -15806.11425583
#> EM: Iteration : 44 || log-likelihood : -15791.6651550126
#> EM: Iteration : 45 || log-likelihood : -15776.2575311116
#> EM: Iteration : 46 || log-likelihood : -15760.2525673176
#> EM: Iteration : 47 || log-likelihood : -15743.6600428386
#> EM: Iteration : 48 || log-likelihood : -15725.8494727209
#> EM: Iteration : 49 || log-likelihood : -15705.5392028324
#> EM: Iteration : 50 || log-likelihood : -15681.0330055801
#> EM: Iteration : 51 || log-likelihood : -15650.7058006772
#> EM: Iteration : 52 || log-likelihood : -15614.1891628978
#> EM: Iteration : 53 || log-likelihood : -15574.3209962234
#> EM: Iteration : 54 || log-likelihood : -15536.9561042095
#> EM: Iteration : 55 || log-likelihood : -15505.9888676546
#> EM: Iteration : 56 || log-likelihood : -15480.3479747868
#> EM: Iteration : 57 || log-likelihood : -15456.7432033066
#> EM: Iteration : 58 || log-likelihood : -15432.855894347
#> EM: Iteration : 59 || log-likelihood : -15408.4123139152
#> EM: Iteration : 60 || log-likelihood : -15384.7708355233
#> EM: Iteration : 61 || log-likelihood : -15363.3704926307
#> EM: Iteration : 62 || log-likelihood : -15344.3247788467
#> EM: Iteration : 63 || log-likelihood : -15326.444200793
#> EM: Iteration : 64 || log-likelihood : -15308.1502066517
#> EM: Iteration : 65 || log-likelihood : -15288.3650661699
#> EM: Iteration : 66 || log-likelihood : -15267.1380314858
#> EM: Iteration : 67 || log-likelihood : -15245.8151021308
#> EM: Iteration : 68 || log-likelihood : -15226.3007649639
#> EM: Iteration : 69 || log-likelihood : -15209.9671868432
#> EM: Iteration : 70 || log-likelihood : -15197.3697193674
#> EM: Iteration : 71 || log-likelihood : -15187.8845852548
#> EM: Iteration : 72 || log-likelihood : -15180.4065779427
#> EM: Iteration : 73 || log-likelihood : -15174.1897193241
#> EM: Iteration : 74 || log-likelihood : -15168.8680084075
#> EM: Iteration : 75 || log-likelihood : -15164.1615627415
#> EM: Iteration : 76 || log-likelihood : -15159.6679572457
#> EM: Iteration : 77 || log-likelihood : -15155.1488045656
#> EM: Iteration : 78 || log-likelihood : -15150.9231858137
#> EM: Iteration : 79 || log-likelihood : -15147.2212168192
#> EM: Iteration : 80 || log-likelihood : -15144.078942659
#> EM: Iteration : 81 || log-likelihood : -15141.3516305636
#> EM: Iteration : 82 || log-likelihood : -15138.8602529876
#> EM: Iteration : 83 || log-likelihood : -15136.5059345662
#> EM: Iteration : 84 || log-likelihood : -15134.2384537766
#> EM: Iteration : 85 || log-likelihood : -15132.0298589309
#> EM: Iteration : 86 || log-likelihood : -15129.8608706576
#> EM: Iteration : 87 || log-likelihood : -15127.7157936565
#> EM: Iteration : 88 || log-likelihood : -15125.5797196054
#> EM: Iteration : 89 || log-likelihood : -15123.4372146492
#> EM: Iteration : 90 || log-likelihood : -15121.2712280838
#> EM: Iteration : 91 || log-likelihood : -15119.0622569401
#> EM: Iteration : 92 || log-likelihood : -15116.7874031382
#> EM: Iteration : 93 || log-likelihood : -15114.4192658119
#> EM: Iteration : 94 || log-likelihood : -15111.9245293407
#> EM: Iteration : 95 || log-likelihood : -15109.262047444
#> EM: Iteration : 96 || log-likelihood : -15106.3802520661
#> EM: Iteration : 97 || log-likelihood : -15103.2137059945
#> EM: Iteration : 98 || log-likelihood : -15099.6787565231
#> EM: Iteration : 99 || log-likelihood : -15095.6664401258
#> EM: Iteration : 100 || log-likelihood : -15091.0341403017
#> EM: Iteration : 101 || log-likelihood : -15085.5952981967
#> EM: Iteration : 102 || log-likelihood : -15079.1100803411
#> EM: Iteration : 103 || log-likelihood : -15071.2863215881
#> EM: Iteration : 104 || log-likelihood : -15061.8155026615
#> EM: Iteration : 105 || log-likelihood : -15050.4931948422
#> EM: Iteration : 106 || log-likelihood : -15037.4728804542
#> EM: Iteration : 107 || log-likelihood : -15023.5663638262
#> EM: Iteration : 108 || log-likelihood : -15010.227713049
#> EM: Iteration : 109 || log-likelihood : -14998.9216243488
#> EM: Iteration : 110 || log-likelihood : -14990.3428946115
#> EM: Iteration : 111 || log-likelihood : -14984.2931646741
#> EM: Iteration : 112 || log-likelihood : -14980.0317050997
#> EM: Iteration : 113 || log-likelihood : -14976.7574542595
#> EM: Iteration : 114 || log-likelihood : -14973.9768267566
#> EM: Iteration : 115 || log-likelihood : -14971.5304235767
#> EM: Iteration : 116 || log-likelihood : -14969.3710026547
#> EM: Iteration : 117 || log-likelihood : -14967.3301314624
#> EM: Iteration : 118 || log-likelihood : -14965.1319732928
#> EM: Iteration : 119 || log-likelihood : -14962.818626259
#> EM: Iteration : 120 || log-likelihood : -14961.1657986148
#> EM: Iteration : 121 || log-likelihood : -14960.1001793804
#> EM: Iteration : 122 || log-likelihood : -14959.2029493404
#> EM: Iteration : 123 || log-likelihood : -14958.3643653619
#> EM: Iteration : 124 || log-likelihood : -14957.5579272948
#> EM: Iteration : 125 || log-likelihood : -14956.7769206505
#> EM: Iteration : 126 || log-likelihood : -14956.0220832192
#> EM: Iteration : 127 || log-likelihood : -14955.2990068376
#> EM: Iteration : 128 || log-likelihood : -14954.6080936987
#> EM: Iteration : 129 || log-likelihood : -14953.9546052572
#> EM: Iteration : 130 || log-likelihood : -14953.3424683065
#> EM: Iteration : 131 || log-likelihood : -14952.7742704947
#> EM: Iteration : 132 || log-likelihood : -14952.2512735504
#> EM: Iteration : 133 || log-likelihood : -14951.7732467988
#> EM: Iteration : 134 || log-likelihood : -14951.3384384815
#> EM: Iteration : 135 || log-likelihood : -14950.9439547413
#> EM: Iteration : 136 || log-likelihood : -14950.5860673359
#> EM: Iteration : 137 || log-likelihood : -14950.2605961901
#> EM: Iteration : 138 || log-likelihood : -14949.9632302133
#> EM: Iteration : 139 || log-likelihood : -14949.6897803656
#> EM: Iteration : 140 || log-likelihood : -14949.4363440458
#> EM: Iteration : 141 || log-likelihood : -14949.1993934329
#> EM: Iteration : 142 || log-likelihood : -14948.9758045711
#> EM: Iteration : 143 || log-likelihood : -14948.7628462595
#> EM: Iteration : 144 || log-likelihood : -14948.5581447387
#> EM: Iteration : 145 || log-likelihood : -14948.3596363733
#> EM: Iteration : 146 || log-likelihood : -14948.1655161518
#> EM: Iteration : 147 || log-likelihood : -14947.9741866833
#> EM: Iteration : 148 || log-likelihood : -14947.7842100466
#> EM: Iteration : 149 || log-likelihood : -14947.5942633197
#> EM: Iteration : 150 || log-likelihood : -14947.4030977377
#> EM: Iteration : 151 || log-likelihood : -14947.2095010109
#> EM: Iteration : 152 || log-likelihood : -14947.0122620331
#> EM: Iteration : 153 || log-likelihood : -14946.8101371804
#> EM: Iteration : 154 || log-likelihood : -14946.6018173877
#> EM: Iteration : 155 || log-likelihood : -14946.3858952193
#> EM: Iteration : 156 || log-likelihood : -14946.1608312027
#> EM: Iteration : 157 || log-likelihood : -14945.9249187549
#> EM: Iteration : 158 || log-likelihood : -14945.676247118
#> EM: Iteration : 159 || log-likelihood : -14945.4126618353
#> EM: Iteration : 160 || log-likelihood : -14945.1317224602
#> EM: Iteration : 161 || log-likelihood : -14944.8306573941
#> EM: Iteration : 162 || log-likelihood : -14944.5063160023
#> EM: Iteration : 163 || log-likelihood : -14944.1551184229
#> EM: Iteration : 164 || log-likelihood : -14943.7730037188
#> EM: Iteration : 165 || log-likelihood : -14943.355377134
#> EM: Iteration : 166 || log-likelihood : -14942.8970570836
#> EM: Iteration : 167 || log-likelihood : -14942.3922219831
#> EM: Iteration : 168 || log-likelihood : -14941.8343559995
#> EM: Iteration : 169 || log-likelihood : -14941.2161912546
#> EM: Iteration : 170 || log-likelihood : -14940.5296397031
#> EM: Iteration : 171 || log-likelihood : -14939.7657190993
#> EM: Iteration : 172 || log-likelihood : -14938.9144460343
#> EM: Iteration : 173 || log-likelihood : -14937.9647057519
#> EM: Iteration : 174 || log-likelihood : -14936.9040831122
#> EM: Iteration : 175 || log-likelihood : -14935.7186499891
#> EM: Iteration : 176 || log-likelihood : -14934.3927038884
#> EM: Iteration : 177 || log-likelihood : -14932.9084527435
#> EM: Iteration : 178 || log-likelihood : -14931.245639997
#> EM: Iteration : 179 || log-likelihood : -14929.3811026273
#> EM: Iteration : 180 || log-likelihood : -14927.2882537299
#> EM: Iteration : 181 || log-likelihood : -14924.9364821865
#> EM: Iteration : 182 || log-likelihood : -14922.2904675358
#> EM: Iteration : 183 || log-likelihood : -14919.3094231961
#> EM: Iteration : 184 || log-likelihood : -14915.9463144684
#> EM: Iteration : 185 || log-likelihood : -14912.1471647651
#> EM: Iteration : 186 || log-likelihood : -14907.8506901999
#> EM: Iteration : 187 || log-likelihood : -14902.9887290339
#> EM: Iteration : 188 || log-likelihood : -14897.4883102736
#> EM: Iteration : 189 || log-likelihood : -14891.27676833
#> EM: Iteration : 190 || log-likelihood : -14884.2919447409
#> EM: Iteration : 191 || log-likelihood : -14876.4995909623
#> EM: Iteration : 192 || log-likelihood : -14867.9179321727
#> EM: Iteration : 193 || log-likelihood : -14858.6442978196
#> EM: Iteration : 194 || log-likelihood : -14848.8804338117
#> EM: Iteration : 195 || log-likelihood : -14838.9872847758
#> EM: Iteration : 196 || log-likelihood : -14829.6292321768
#> EM: Iteration : 197 || log-likelihood : -14821.8717823403
#> EM: Iteration : 198 || log-likelihood : -14816.6461672058
#> EM: Iteration : 199 || log-likelihood : -14813.7497363742
#> EM: Iteration : 200 || log-likelihood : -14812.2267827519
#> EM: Iteration : 201 || log-likelihood : -14811.4198287137
#> EM: Iteration : 202 || log-likelihood : -14811.0049217051
#> EM: Iteration : 203 || log-likelihood : -14810.7960368513
#> EM: Iteration : 204 || log-likelihood : -14810.6883875777

mixrhlp$summary()
#> ------------------------
#> Fitted mixRHLP model
#> ------------------------
#> 
#> MixRHLP model with G = 3 clusters and K = 3 regimes:
#> 
#>  log-likelihood nu       AIC       BIC       ICL
#>       -14810.69 41 -14851.69 -14880.41 -14880.41
#> 
#> Clustering table (Number of curves in each clusters):
#> 
#>  1  2  3 
#> 10 10 10 
#> 
#> Mixing probabilities (cluster weights):
#>          1         2         3
#>  0.3333333 0.3333333 0.3333333
#> 
#> 
#> --------------------
#> Cluster 1 (G = 1):
#> 
#> Regression coefficients:
#> 
#>     Beta(K = 1) Beta(K = 2) Beta(K = 3)
#> 1     6.8902863   5.1134337  3.90153421
#> X^1   0.9265632  -0.3959402  0.08748466
#> 
#> Variances:
#> 
#>  Sigma2(K = 1) Sigma2(K = 2) Sigma2(K = 3)
#>       0.981915     0.9787717     0.9702211
#> 
#> --------------------
#> Cluster 2 (G = 2):
#> 
#> Regression coefficients:
#> 
#>     Beta(K = 1) Beta(K = 2) Beta(K = 3)
#> 1    4.96556671   6.7326717   4.8807183
#> X^1  0.08880479   0.4984443   0.1350271
#> 
#> Variances:
#> 
#>  Sigma2(K = 1) Sigma2(K = 2) Sigma2(K = 3)
#>      0.9559969       1.03849     0.9506928
#> 
#> --------------------
#> Cluster 3 (G = 3):
#> 
#> Regression coefficients:
#> 
#>     Beta(K = 1) Beta(K = 2) Beta(K = 3)
#> 1     6.3513369    4.214736   6.6536553
#> X^1  -0.2449377    0.839666   0.1024863
#> 
#> Variances:
#> 
#>  Sigma2(K = 1) Sigma2(K = 2) Sigma2(K = 3)
#>      0.9498285     0.9270384      1.001413

mixrhlp$plot()
```

<img src="man/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-6-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-6-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-6-4.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-6-5.png" style="display: block; margin: auto;" />
