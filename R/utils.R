ones <- function(n,d,g=1){
  if (g==1){
    return(matrix(1,n,d))
  }
  else{
    return(array(1,dim=c(n,d,g)))
  }
}

zeros <- function(n,d,g=1){
  if (g==1){
    return(matrix(0,n,d))
  }
  else{
    return(array(0,dim=c(n,d,g)))
  }
}

rand <- function(n,d,g=1){
  if (g==1){
    return(matrix(runif(n*d), n,d))
  }
  else{
    return(array(runif(n*d),dim=c(n,d,g)))
  }
}

repmat <- function(M, n, d){
  return(kronecker(matrix(1, n, d), M))
}

drnorm <- function(n, d, mean, sd){
  A <- matrix(nrow = n, ncol = d)
  for (i in 1:d){
    A[,i] <- rnorm(n,mean,sd)
  }
  return(A)
}


lognormalize <- function(M){
  n <- nrow(M)
  d <- ncol(M)
  #to do: continue develop
}


modele_logit <- function(Wg, phiW, Y=NULL, Gamma=NULL){
  #W - Wg
  #M - phiW
  #Y - ?
  #Gamma - ?
  if (!is.null(Y)) {
    n1 <- nrow(Y)
    K <- ncol(Y)
    if (!is.null(Gamma)) {
      Gamma <- Gamma * ones(1,K)
    }
    n2 <- nrow(phiW)
    q <- ncol(phiW)
    if (n1==n2){
      n <- n1
    }
    else{
      stop("Wg and Y must have the same number of lines")
    }
  }
  else{
    n <- nrow(Wg)
    q <- nrow(Wg)
  }

  if (!is.null(Y)){
    if (ncol(Wg) == (K-1)){ # pas de vecteur nul dans W donc l'ajouter
      wK <- zeros(q,1)
      Wg <- cbind(Wg, wK)
    }
    else{
      stop("Wg and Y must have the same number of lines")
    }
  }
  else{
    wK <- zeros(q,1)
    Wg <- cbind(Wg, wK)
    q <- nrow(Wg)
    K <- ncol(Wg)
  }

  MW <- phiW %*% Wg
  maxm <- max.col(MW)
  MW <- MW - maxm %*% ones(1,K)

  expMW <- exp(MW)

  probas <- expMW / (apply(expMW[,1:K],1,sum) %*% ones(1,K))

  if (!is.null(Y)){
    if (is.null(Gamma)) {
      loglik <- sum(apply((Y*MW) - (Y*log(apply(expMW,1,sum)%*%ones(1,K))),1,sum))
    }
    else {
      loglik <- sum(apply((Gamma*(Y*MW)) - ((Gamma*Y)*log(apply(expMW,1,sum)%*%ones(1,K))),1,sum))
    }

    # todo: verify R computation of loglik gives nan
    if (is.nan(loglik)){
      MW <- phiW %*% Wg
      minm <- -745.1
      MW <- max(MW, minm)
      maxm <- 709.78
      MW <- min(MW, maxm)
      expMW <- exp(MW)

      eps <- 2^(-52)

      if (is.null(Gamma)) {
        loglik <- sum(apply((Y*MW) - (Y*log(apply(expMW,1,sum)%*%ones(1,K)+eps)),1,sum))
      }
      else {
        loglik <- sum(apply((Gamma*(Y*MW)) - ((Gamma*Y)*log(apply(expMW,1,sum)%*%ones(1,K)+eps)),1,sum))
      }
    }
    if (is.nan(loglik)){
      stop("Probleme loglik NaN (!!!)")
    }
  }
  else{
    loglik <- list()
  }

  return(list(probas, loglik))
}


