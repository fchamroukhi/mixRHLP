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
  a <- max.col(M)
  return(M-repmat(a + log(rowSums(exp(M - repmat(a,1,d)))), 1, d))
}

IRLS_MixFRHLP <- function(cluster_weights, tauijk, phiW, Wg_init=NULL, verbose_IRLS=FALSE){
  K <- ncol(tauijk)
  n <- nrow(phiW)
  q <- ncol(phiW)
  if (is.null(Wg_init)){
    Wg_init <- zeros(q,K-1)
  }
  lambda <- 1e-9
  I <- diag(q*(K-1))

  #Initialisation du IRLS (iter = 0)
  W_old <- Wg_init

  problik <- modele_logit(W_old, phiW, tauijk, cluster_weights)
  piik_old <- problik[[1]]
  loglik_old <- problik[[2]]

  loglik_old <- loglik_old - lambda * norm(as.vector(W_old),"2")^2
  iter <- 0
  converge <- FALSE
  max_iter <- 300
  LL <- list()
  if (verbose_IRLS){
    message("IRLS : Iteration ", iter, "Log-likehood : ", loglik_old)
  }

  while(!converge && (iter<max_iter)){
    #Hw_old matrice carree de dimensions hx x hx
    hx <- q*(K-1)
    Hw_old <- zeros(hx,hx)
    gw_old <- zeros(q,K-1)

    #Gradient
    for (k in 1:K-1){
      gwk <- cluster_weights*(tauijk[,k] - piik_old[,k])
      for (qq in 1:q){
        vq <- phiW[,qq]
        gw_old[qq,k] <- t(gwk) %*% vq
      }
    }
    gw_old <- matrix(gw_old, q*(K-1), 1)


    #Hessienne
    for (k in 1:K-1){
      for (ell in 1:K-1){
        delta_kl <- (k==ell)
        gwk <- cluster_weights * (piik_old[,k] * (ones(n,1)%*%delta_kl - piik_old[,ell]))
        Hkl <- zeros(q,q)
        for (qqa in 1:q){
          vqa <- phiW[,qqa]
          for (qqb in 1:q){
            vqb <- phiW[, qqb]
            hwk <- t(vqb) %*% (gwk * vqa)
            Hkl[qqa,qqb] <- hwk
          }
        }
        Hw_old[(k-1)*q+1 : k*q, (ell-1)*q+1 : ell*q] <- -Hkl
      }
    }

    # si a priori gaussien sur W (lambda ~ 1e-9)
    Hw_old <- Hw_old + lambda%*%I
    gw_old = gw_old - lambda*as.vector(W_old)

    # Newton Raphson : W(c+1) = W(c) - H(W(c))^(-1)g(W(c))
    w <- as.vector(W_old) - solve(Hw_old)%*%gw_old # [(q+1)x(K-1),1]
    W <- matrix(w,q,K-1) #[(q+1)*(K-1)]

    # mise a jour des probas et de la loglik
    problik <- modele_logit(W, phiW, tauijk, cluster_weights)
    piik <- problik[[1]]
    loglik <- problik[[2]]
    loglik <- loglik - lambda*(norm(as.vector(W_old),"2"))^2

    # Verifier si Qw1(w^(c+1),w^(c))> Qw1(w^(c),w^(c))
    #(adaptation) de Newton Raphson : W(c+1) = W(c) - pas*H(W)^(-1)*g(W)
    pas <- 1
    alpha <- 2

    while(loglik < loglik_old){
      pas <- pas/alpha # pas d'adaptation de l'algo Newton raphson
      #recalcul du parametre W et de la loglik
      #Hw_old = Hw_old + lambda*I;
      w <- as.vector(W_old) - pas * solve(Hw_old)%*%gw_old
      W = matrix(w,q,K-1)
      problik <- modele_logit(W, phiW, tauijk, cluster_weights)
      piik <- problik[[1]]
      loglik <- problik[[2]]
      loglik <- loglik - lambda*(norm(as.vector(W_old),"2"))^2
    }

    converge1 <- abs((loglik - loglik_old)/loglik_old) <= 1e-7
    converge2 <- abs(loglik - loglik_old) <= 1e-6

    converge <- converge1 | converge2

    piik_old <- piik
    W_old <- W
    iter <- iter + 1
    LL[iter] <- loglik_old
    loglik_old <- loglik
    if(verbose_IRLS){
      message("IRLS: Iteration ", iter, "Log-likelihood: ", loglik_old)
    }
  }

  if (converge){
    if (verbose_IRLS){
      message("IRLS : convergence  OK ; nbre d''iterations : ", iter)
    }
  }
  else{
    message("IRLS : pas de convergence (augmenter le nombre d''iterations > ", max_iter, ")")
  }

  reg_irls <- 0
  if (lambda!=0){
    reg_irls <- lambda * (norm(as.vector(W),"2"))^2
  }

  return(list(W, piik, reg_irls, LL, loglik))
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


