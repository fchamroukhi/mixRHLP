myKmeans = function(X, K, nbr_runs=20, nbr_iter_max=100, verbose = FALSE) {
  #
  #
  #   K-means algorithm
  #
  #   X: dataset with the dimension nxp
  #   K: the number of groups in the dataset X
  #   
  #   using euclidian distance
  #
  #   based on matlab code: Faicel CHAMROUKHI Septembre 2008 (update)
  ###########################################################################
  
  n = nrow(X)
  p = ncol(X)
  
  res = list()
  
  if (K == 1) {
    ## one class example
    klas_one = rep(1, n)
    global_mean = apply(X, 2, mean)
    dmin = rowSums((X - klas_one %*% t(global_mean)) ^ 2)
    res$muk = global_mean
    res$klas = klas_one
    res$err = sum(dmin)
    res$Zik = klas_one
    
    return(res)
    
  } else {
    nbr_run = 0
    best_solution = list()
    best_solution$err = Inf
    
    while (nbr_run < nbr_runs) {
      nbr_run = nbr_run + 1
      if (verbose) {
        cat(sprintf("Kmeans : run nr. %1i", nbr_run), "\n")
      }
      critere = vector()
      iter = 0
      partition = matrix(0, nrow = n, ncol = K) # initialisation de la partition
      Zik = matrix(0, nrow = n, ncol = K)
      converged = FALSE
      err = -Inf
      ## 1.Initialisation aleatoire des centres
      indices_aleatoires = sample(n)
      centres = X[indices_aleatoires[1:K], ]
      while ((iter <= nbr_iter_max) && !(converged)) {
        iter = iter + 1
        clas_vide = vector()
        old_centres = centres
        
        # Computing the euclidian distance
        dist_eucld = matrix(0, nrow = n, ncol = K) # distance between each datapoint and mean
        for (k in 1:K) {
          mk = centres[k,]
          dist_eucld[, k] = rowSums((X - rep(1, n) %*% t(mk)) ^ 2)
        }
        
        ## The classification step
        dmin = apply(dist_eucld, 1, min)
        klas = apply(dist_eucld, 1, which.min)
        Zik = ((klas %*% t(rep(1, K))) == rep(1, n) %*% t(1:K)) * 1
        
        ## Relocation step
        sigmak = array(data = 0, dim = c(p, p, K))
        pik = rep(0, K)
        for (m in 1:K) {
          ind_ck = which(klas == m)
          
          #if empty classes
          if (length(ind_ck) == 0) {
            clas_vide = rbind(clas_vide, m)
          } else {
            classek = X[ind_ck, ]
            # update the centers
            # handle the case where a class is only one point
            if (length(ind_ck) == 1) {
              centres[m,] = classek
            } else {
              centres[m,] = apply(classek, 2, mean) 
            }
          }
        }
        
        centres[clas_vide, ] = old_centres[clas_vide, ]
        
        ## 1. first stop criteria : 
        err2 = sum(rowSums(Zik * dist_eucld)) 
        crit1 = (abs(err2 - err)) / err < 1e-6
        
        #Handle NA value for crit1 at the first pass
        if (is.na(crit1)) {
          crit1 = FALSE
        }
        
        ## 2. second stop criteria :
        crit2 = max(abs(centres - old_centres)) < 1e-6
        
        if (any(crit1, crit2)) {
          converged = TRUE # convergence
        }
        
        err = err2
        if (verbose) {
          cat(sprintf("Kmeans : Iteration  %1.1i, Objective %2.2f", iter, err), "\n")
        }
        
      } # one run
      
      ##  solution for each run
      res$muk = centres
      res$klas = klas
      res$err = err
      res$Zik = Zik
      
      if (err < best_solution$err) {
        best_solution = res
      }
    } # End of Kmeans runs
    
    res = best_solution
    
    return(res)
  }
}