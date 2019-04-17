modele_logit <- function(Wg, phiW, Y = NULL, Gamma = NULL) {
  "
  calculates the pobabilities according to multinomial logistic model
  "
  if (!is.null(Y)) {
    n1 <- nrow(Y)
    K <- ncol(Y)
    if (!is.null(Gamma)) {
      Gamma <- Gamma %*% ones(1, K)
    }
    n2 <- nrow(phiW)
    q <- ncol(phiW)
    if (n1 == n2) {
      n <- n1
    }
    else{
      stop("Wg and Y must have the same number of lines")
    }
  }
  else{
    n <- nrow(phiW)
    q <- ncol(phiW)
  }

  if (!is.null(Y)) {
    if (ncol(Wg) == (K - 1)) {
      # W doesnt contain the null vector associated with the last class
      wK <- zeros(q, 1)
      Wg <- cbind(Wg, wK) #add the null vector wK for the last component probability
    }
    else{
      stop("Wg and Y must have the same number of lines")
    }
  }
  else{
    wK <- zeros(q, 1)
    Wg <- cbind(Wg, wK)
    q <- nrow(Wg)
    K <- ncol(Wg)
  }

  MW <- phiW %*% Wg
  maxm <- apply(MW, 1, max)
  MW <- MW - maxm %*% ones(1, K) # to avoid overfolow

  expMW <- exp(MW)
  if (ncol(expMW) == 1) {
    probas <- expMW / (expMW[, 1:K] %*% ones(1, K))
  }
  else{
    probas <- expMW / (rowSums(expMW[, 1:K]) %*% ones(1, K))
  }


  if (!is.null(Y)) {
    if (is.null(Gamma)) {
      loglik <- sum((Y * MW) - (Y * log(rowSums(expMW) %*% ones(1, K))))
    }
    else {
      loglik <-
        sum(((Gamma * Y) * MW) - ((Gamma * Y) * log(rowSums(expMW) %*% ones(1, K))))
    }

    if (is.nan(loglik)) {
      # to avoid numerical overflow since exp(XW=-746)=0 and exp(XW=710)=inf)
      MW <- phiW %*% Wg
      minm <- -745.1
      MW <- pmax(MW, minm)
      maxm <- 709.78
      MW <- pmin(MW, maxm)
      expMW <- exp(MW)

      eps <- .Machine$double.eps

      if (is.null(Gamma)) {
        loglik <- sum((Y * MW) - (Y * log(rowSums(expMW) %*% ones(1, K) + eps)))
      }
      else {
        loglik <- sum(((Gamma * Y) * MW) - ((Gamma * Y) * log(rowSums(expMW) %*% ones(1, K) + eps)))
      }
    }
    if (is.nan(loglik)) {
      stop("Probleme loglik NaN (!!!)")
    }
  }
  else{
    loglik <- c()
  }

  return(list(probas = probas, loglik = loglik))
}
