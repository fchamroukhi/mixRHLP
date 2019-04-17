IRLS <- function(tauijk, phiW, Wg_init = NULL, cluster_weights = NULL, verbose_IRLS = FALSE, piik_len = NULL) {
    "an efficient Iteratively Reweighted Least-Squares (IRLS) algorithm for estimating
    the parameters of a multinomial logistic regression model given the
    predictors X and a partition (hard or smooth) Tau into K>=2 segments,
    and a cluster weights Gamma (hard or smooth)
    %% References
    % Please cite the following papers for this code:
    %
    %
    % @INPROCEEDINGS{Chamroukhi-IJCNN-2009,
    %   AUTHOR =       {Chamroukhi, F. and Sam\'e,  A. and Govaert, G. and Aknin, P.},
    %   TITLE =        {A regression model with a hidden logistic process for feature extraction from time series},
    %   BOOKTITLE =    {International Joint Conference on Neural Networks (IJCNN)},
    %   YEAR =         {2009},
    %   month = {June},
    %   pages = {489--496},
    %   Address = {Atlanta, GA},
    %  url = {https://chamroukhi.users.lmno.cnrs.fr/papers/chamroukhi_ijcnn2009.pdf}
    % }
    %
    % @article{chamroukhi_et_al_NN2009,
    %  	Address = {Oxford, UK, UK},
    % 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
    % 	Date-Added = {2014-10-22 20:08:41 +0000},
    % 	Date-Modified = {2014-10-22 20:08:41 +0000},
    % 	Journal = {Neural Networks},
    % 	Number = {5-6},
    % 	Pages = {593--602},
    % 	Publisher = {Elsevier Science Ltd.},
    % 	Title = {Time series modeling by a regression approach based on a latent process},
    % 	Volume = {22},
    % 	Year = {2009},
    % 	url  = {https://chamroukhi.users.lmno.cnrs.fr/papers/Chamroukhi_Neural_Networks_2009.pdf}
    % 	}
    % @article{Chamroukhi-FDA-2018,
    % 	Journal = {},
    % 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
    % 	Volume = {},
    % 	Title = {Model-Based Clustering and Classification of Functional Data},
    % 	Year = {2018},
    % 	eprint ={arXiv:1803.00276v2},
    % 	url =  {https://chamroukhi.users.lmno.cnrs.fr/papers/MBCC-FDA.pdf}
    % 	}
    "
    K <- ncol(tauijk)
    n <- nrow(phiW)
    q <- ncol(phiW)

    if (K == 1) {
      W <- matrix(nrow = (q), ncol = 0)
      piik <- ones(piik_len, 1)
      reg_irls <- 0
      LL <- 0
      loglik <- 0

    } else {
      if (is.null(Wg_init)) {
        Wg_init <- zeros(q, K - 1)
      }
      lambda <- 1e-9
      I <- diag(q * (K - 1))


      #IRLS Initialization (iter = 0)
      W_old <- Wg_init

      problik <- modele_logit(W_old, phiW, tauijk, cluster_weights)
      piik_old <- problik$probas
      loglik_old <- problik$loglik

      loglik_old <-
        loglik_old - lambda * sum(W_old ^ 2) #norm(as.vector(W_old),"2")^2
      iter <- 0
      converge <- FALSE
      max_iter <- 300
      LL <- c()
      if (verbose_IRLS) {
        message("IRLS : Iteration ", iter, "Log-likehood : ", loglik_old)
      }

      while (!converge && (iter < max_iter)) {
        #Hw_old a squared matrix of dimensions  q*(K-1) x  q*(K-1)
        hx <- q * (K - 1)
        Hw_old <- zeros(hx, hx)
        gw_old <- zeros(q, K - 1)

        #Gradient
        for (k in 1:(K - 1)) {
          if (is.null(cluster_weights)) {
            gwk <- tauijk[, k] - piik_old[, k]
          }
          else{
            gwk <- cluster_weights * (tauijk[, k] - piik_old[, k])
          }

          for (qq in 1:q) {
            vq <- phiW[, qq]
            gw_old[qq, k] <- as.numeric(t(gwk) %*% vq)
          }
        }
        gw_old <- matrix(gw_old, nrow = q * (K - 1), ncol = 1)


        #Hessienne
        for (k in 1:(K - 1)) {
          for (ell in 1:(K - 1)) {
            delta_kl <- as.numeric(k == ell)
            if (is.null(cluster_weights)) {
              gwk <- piik_old[, k] * (ones(n, 1) * delta_kl - piik_old[, ell])
            }
            else{
              gwk <-
                cluster_weights * (piik_old[, k] * (ones(n, 1) * delta_kl - piik_old[, ell]))
            }

            Hkl <- zeros(q, q)
            for (qqa in 1:q) {
              vqa <- phiW[, qqa]
              for (qqb in 1:q) {
                vqb <- phiW[, qqb]
                hwk <- t(vqb) %*% (gwk * vqa)
                Hkl[qqa, qqb] <- hwk
              }
            }
            Hw_old[(((k - 1) * q) + 1):(k * q), (((ell - 1) * q) + 1):(ell * q)] <- -Hkl
          }
        }

        # if a gaussien prior on W (lambda ~=0)
        Hw_old <- Hw_old + lambda * I
        gw_old = gw_old - lambda * as.vector(W_old)

        # Newton Raphson : W(c+1) = W(c) - H(W(c))^(-1)g(W(c))
        w <-
          as.vector(W_old) - solve(Hw_old) %*% gw_old # [(q+1)x(K-1),1]
        W <- matrix(w, q, (K - 1)) #[(q+1)*(K-1)]

        # mise a jour des probas et de la loglik
        problik <- modele_logit(W, phiW, tauijk, cluster_weights)
        piik <- problik$probas
        loglik <- problik$loglik
        loglik <-
          loglik - lambda * sum(W ^ 2) #(norm(as.vector(W_old),"2"))^2

        ##  check if Qw1(w^(t+1),w^(t))> Qw1(w^(t),w^(t))
        ##(adaptive stepsize in case of troubles with stepsize 1) Newton Raphson : W(t+1) = W(t) - stepsize*H(W)^(-1)*g(W)
        pas <- 1
        alpha <- 2

        while (loglik < loglik_old) {
          pas <- pas / alpha # pas d'adaptation de l'algo Newton raphson
          #recalcul du parametre W et de la loglik
          #Hw_old = Hw_old + lambda*I;
          w <- as.vector(W_old) - pas * solve(Hw_old) %*% gw_old
          W = matrix(w, q, K - 1)
          problik <- modele_logit(W, phiW, tauijk, cluster_weights)
          piik <- problik$probas
          loglik <- problik$loglik

          loglik <-
            loglik - lambda ** sum(W ^ 2)#(norm(as.vector(W),"2"))^2
        }

        converge1 <- abs((loglik - loglik_old) / loglik_old) <= 1e-7
        converge2 <- abs(loglik - loglik_old) <= 1e-6

        converge <- converge1 | converge2

        piik_old <- piik
        W_old <- W
        iter <- iter + 1
        LL[iter] <- loglik_old
        loglik_old <- loglik
        if (verbose_IRLS) {
          message("IRLS: Iteration ", iter, "Log-likelihood: ", loglik_old)
        }
      }

      if (converge) {
        if (verbose_IRLS) {
          message("IRLS : convergence  OK ; nbre d''iterations : ", iter)
        }
      }
      else{
        message("IRLS : pas de convergence (augmenter le nombre d''iterations > ", max_iter,")")
      }

      reg_irls <- 0
      if (lambda != 0) {
        reg_irls <- lambda * (norm(as.vector(W), "2")) ^ 2
      }

    }

    return(list(W = W, piik = piik, reg_irls = reg_irls, LL = LL, loglik = loglik))
    }
