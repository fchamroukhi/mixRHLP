ones <- function(n,d){
  return(matrix(1,n,d))
}

repmat <- function(M, n, d){
  return(kronecker(matrix(1,n,d),M))
}



drnorm <- function(n,d,mean,sd){
  A=matrix(nrow = n, ncol = d)
  for (i in 1:d){
    A[,i]<-rnorm(n,mean,sd)
  }
  return(A)
}




