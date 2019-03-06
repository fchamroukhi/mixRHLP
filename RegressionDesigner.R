RegressionDesigner <- setRefClass(
  "RegressionDesigner",
  fields = list(
    XBeta = "matrix",
    Xw = "matrix"
  ),

  methods = list(
    designmatrix_FRHLP = function(x,p,q=NULL){
      "
       constructs the desing matrix of a polynomial regression with degree p
       (for the regimes), and with degree q (for the logistic regression)
      "
      if (ncol(x) == 1){
        x<-t(x)
      }
      order_max <- p
      if (!is.null(q)){
        order_max <- max(p,q)
      }

      X <- matrix(NA, length(x), order_max+1)
      for (i in 0:(order_max)){
        X[,i+1] <- x^i # X = [1 x x^2 x^3 x^p;......;...]
      }

      XBeta <<- X[,1:(p+1)]; # design matrix for Beta (the polynomial regressions)
      if (!is.null(q)){
        Xw <<- X[,1:(q+1)]; # design matrix for w (the logistic regression)
      }
    },

    setPhi1 = function(x,p,q){
      ## for 1 curve
      designmatrix_FRHLP(x, p, q)
    },

    setPhiN = function(x,p,q,n){
      ## for 1 curve
      designmatrix_FRHLP(x, p, q)

      ##for n curves
      XBeta <<- repmat(XBeta, n, 1)
      Xw <<- repmat(Xw, n, 1)
    }
  )
)
