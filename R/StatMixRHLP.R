#' A Reference Class which contains statistics of a MixRHLP model.
#'
#' StatMRHLP contains all the parameters of a [MixRHLP][ParamMixRHLP] model.
#'
#' @field h_ig Fuzzy segmentation matrix of size \eqn{(n, G)} representing the
#' post probabilities prob(curve|cluster_g)
#' @field c_ig Hard partition logical matrix of dimension \eqn{(n, G)}
#' obtained by the Maximum a posteriori (MAP) rule:
#' \eqn{c_{ig} = 1 \ \textrm{if} \ c_{ig} = \textrm{arg} \ \textrm{max}_{k} \
#' P(c_i = g | Y, W, \beta);\ 0 \ \textrm{otherwise}}{c_ig = 1 if c_ig =
#' arg max_k P(c_i = g | Y, W, \beta); 0 otherwise}, \eqn{k = 1,\dots,G}.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#' \eqn{klas(i) = k}, \eqn{k = 1,\dots,K}.
#' @field Ex_g Column vector of dimension m for each g \emph{(m,G)}. `Ex_g` is the curve expectation
#' : sum of the polynomial components \eqn{\beta_{gk} \times X_{i}}{\betak x X_i}
#' weighted by the logistic probabilities `pij_gk`:
#' \eqn{Ex_g(j) = \sum_{k = 1}^{K} pi_jik \times \beta_{gk} \times X_{i}}{Ey(i) =
#' \sum_{k=1}^K pi_jik x \betak x X_i}, \eqn{i = 1,\dots,n}.
#' @field log\_lik Numeric. Log-likelihood of the MixRHLP model.
#' @field com_loglik Numeric. Complete log-likelihood of the MixRHLP model.
#' @field stored_loglik List. Stored values of the log-likelihood at each EM
#' iteration.
#' @field stored_com_loglik List. Stored values of the complete log-likelihood at each EM iteration.
#' @field tau_ijgk Matrix of dimension \emph{(nm,K)} for each g (g=1,...,G), representing the fuzzy segmentation
#' for the cluster g prob(y_{ij}|kth_segment,cluster_g)
#' @field log_tau_ijgk the logarithmic values for tau_ijgk
#' @field BIC Numeric. Value of the BIC (Bayesian Information Criterion)
#' criterion. The formula is \eqn{BIC = log\_lik - nu \times
#' \textrm{log}(n) / 2}{BIC = log\_lik - nu x log(n) / 2} with \emph{nu} the
#' degree of freedom of the MixRHLP model.
#' @field ICL Numeric. Value of the ICL (Integrated Completed Likelihood)
#' criterion. The formula is \eqn{ICL = com\_loglik - nu \times
#' \textrm{log}(n) / 2}{ICL = com_loglik - nu x log(n) / 2} with \emph{nu} the
#' degree of freedom of the MixRHLP model.
#' @field AIC Numeric. Value of the AIC (Akaike Information Criterion)
#' criterion. The formula is \eqn{AIC = log\_lik - nu}{AIC = log\_lik - nu}.
#' @field cpu_time Numeric. Average executon time of a EM step.
#' @field log_fg_xij Loglikelihood for the n_g curves of cluster g. Matrix of dimension \eqn{(n, G)}
#' @field log_alphag_fg_xij Matrix of dimension \eqn{(n, G)}
#' @field polynomials Matrix of size \eqn{(m, K, G)} giving the values of
#' \eqn{\beta_{k} \times X_{i}}{\betak x X_i}, \eqn{i = 1,\dots,m}.
#' @field weighted_polynomials Matrix of size \eqn{(m, K, G)}
#' @seealso [ParamMixRHLP], [FData]
#' @export
StatMixRHLP <- setRefClass(
  "StatMixRHLP",
  fields = list(
    h_ig = "matrix", # h_ig = prob(curve|cluster_g) : post prob (fuzzy segmentation matrix of dim [nxG])
    c_ig = "matrix", # c_ig : Hard partition obtained by the AP rule :  c_{ig} = 1
    # if and only c_i = arg max_g h_ig (g=1,...,G)
    klas = "matrix", # klas : column vector of cluster labels
    Ex_g = "matrix",
    # Ex_g: curve expectation: sum of the polynomial components beta_gk ri weighted by
    # the logitic probabilities pij_gk: Ex_g(j) = sum_{k=1}^K pi_jgk beta_gk rj, j=1,...,m. Ex_g
    # is a column vector of dimension m for each g.
    log_lik = "numeric", # the loglikelihood of the EM or CEM algorithm
    com_loglik = "numeric", # the complete loglikelihood of the EM (computed at the convergence) or CEM algorithm
    stored_loglik = "list", # vector of stored valued of the comp-log-lik at each EM teration
    stored_com_loglik = "list",
    tau_ijgk = "array",
    # tau_ijgk prob(y_{ij}|kth_segment,cluster_g), fuzzy
    # segmentation for the cluster g. matrix of dimension
    # [nmxK] for each g  (g=1,...,G).
    log_tau_ijgk = "array",
    BIC = "numeric", # BIC value = loglik - nu*log(nm)/2.
    ICL = "numeric", # ICL value = comp-loglik_star - nu*log(nm)/2.
    AIC = "numeric", # AIC value = loglik - nu.
    cpu_time = "numeric",
    log_fg_xij = "matrix",
    log_alphag_fg_xij = "matrix",
    polynomials = "array",
    weighted_polynomials = "array"
  ),
  methods = list(
    initialize = function(paramMixRHLP = ParamMixRHLP()) {

      h_ig <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$G)
      c_ig <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$G)
      klas <<- matrix(NA, paramMixRHLP$fData$n, 1)
      Ex_g <<- matrix(NA, paramMixRHLP$fData$m, paramMixRHLP$G)
      log_lik <<- -Inf
      com_loglik <<- -Inf
      stored_loglik <<- list()
      stored_com_loglik <<- list()
      BIC <<- -Inf
      ICL <<- -Inf
      AIC <<- -Inf
      cpu_time <<- Inf
      log_fg_xij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$G)
      log_alphag_fg_xij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$G)
      polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$K, paramMixRHLP$G))
      weighted_polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$K, paramMixRHLP$G))
      tau_ijgk <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$K, paramMixRHLP$G))
      log_tau_ijgk <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$K, paramMixRHLP$G))
    },

    MAP = function() {
      "
         calculate a partition by applying the Maximum A Posteriori Bayes
         allocation rule
      "
      N <- nrow(h_ig)
      K <- ncol(h_ig)
      ikmax <- max.col(h_ig)
      ikmax <- matrix(ikmax, ncol = 1)
      c_ig <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K)
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[c_ig[, k] == 1] <<- k
      }
    },

    computeStats = function(mixParam, cpu_time_all) {

      for (g in 1:mixParam$G) {

        polynomials[, , g] <<- mixParam$phi$XBeta[1:mixParam$fData$m, ] %*% as.matrix(mixParam$betag[, , g])

        if (mixParam$K != 1 && mixParam$G != 1) {
          weighted_polynomials[, , g] <<- mixParam$pi_jgk[, , g] * polynomials[, , g]
          Ex_g[, g] <<- rowSums(weighted_polynomials[, , g])
        } else if (mixParam$K == 1 && mixParam$G != 1) {
          weighted_polynomials[, , g] <<- mixParam$pi_jgk[, g] * polynomials[, , g]
          Ex_g[, g] <<- weighted_polynomials[, , g]
        } else if (mixParam$K != 1 && mixParam$G == 1) {
          weighted_polynomials[, , g] <<- mixParam$pi_jgk * polynomials[, , g]
          Ex_g[, g] <<- matrix(rowSums(weighted_polynomials[, , g]))
        } else {
          weighted_polynomials[, , g] <<- mixParam$pi_jgk * polynomials[, , g]
          Ex_g[, g] <<- weighted_polynomials[, , g]
        }
      }

      Ex_g <<- matrix(Ex_g, nrow = mixParam$fData$m)
      cpu_time <<- mean(cpu_time_all)
      Psi <- c(as.vector(mixParam$alpha_g), as.vector(mixParam$Wg), as.vector(mixParam$betag), as.vector(mixParam$sigma2_g))
      nu <- length(Psi)
      BIC <<- log_lik - (nu * log(mixParam$fData$n) / 2)
      AIC <<- log_lik - nu

      cig_log_alphag_fg_xij <- (c_ig) * (log_alphag_fg_xij)

      com_loglik <<- sum(rowSums(cig_log_alphag_fg_xij))

      ICL <<- com_loglik - nu * log(mixParam$fData$n) / 2
    },

    CStep = function(reg_irls) {

      h_ig <<- exp(lognormalize(log_alphag_fg_xij))

      MAP() # Setting klas and c_ig

      # Compute the optimized criterion
      cig_log_alphag_fg_xij <- (c_ig) * log_alphag_fg_xij
      com_loglik <<- sum(cig_log_alphag_fg_xij) +  reg_irls
    },

    EStep = function(mixParam) {

      for (g in 1:mixParam$G) {
        alpha_g <- mixParam$alpha_g[g]
        beta_g <- mixParam$betag[, , g]
        Wg <- mixParam$Wg[, , g]
        pi_jgk <- mixParam$pi_jgk[, , g]

        if (!is.matrix(beta_g)) {
          beta_g <- matrix(beta_g)
          pi_jgk <- matrix(pi_jgk)
        }

        log_pijgk_fgk_xij <- zeros(mixParam$fData$n * mixParam$fData$m, mixParam$K)

        for (k in 1:mixParam$K) {
          beta_gk <- beta_g[, k]
          if (mixParam$variance_type == "homoskedastic") {
            sgk <- mixParam$sigma2_g[g]
          } else {
            sgk <- mixParam$sigma2_g[k, g]
          }
          z <- ((mixParam$fData$vecY - mixParam$phi$XBeta %*% beta_gk) ^ 2) / sgk
          log_pijgk_fgk_xij[, k] <- log(pi_jgk[, k]) - 0.5 * (log(2 * pi) + log(sgk)) - 0.5 * z # pdf cond c_i = g et z_i = k de xij
        }

        log_pijgk_fgk_xij <- pmin(log_pijgk_fgk_xij, log(.Machine$double.xmax))
        log_pijgk_fgk_xij <- pmax(log_pijgk_fgk_xij, log(.Machine$double.xmin))

        pijgk_fgk_xij <- exp(log_pijgk_fgk_xij)
        sumk_pijgk_fgk_xij <- rowSums(pijgk_fgk_xij) # sum over k
        log_sumk_pijgk_fgk_xij <- log(sumk_pijgk_fgk_xij) # [n*m, 1]

        log_tau_ijgk[, , g] <<- log_pijgk_fgk_xij - log_sumk_pijgk_fgk_xij %*% ones(1, mixParam$K)
        tau_ijgk[, , g] <<- exp(lognormalize(log_tau_ijgk[, , g]))

        log_fg_xij[, g] <<- rowSums(t(matrix(log_sumk_pijgk_fgk_xij, mixParam$fData$m, mixParam$fData$n))) # [n x 1]:  sum over j=1,...,m: fg_xij = prod_j sum_k pi_{jgk} N(x_{ij},mu_{gk},s_{gk))
        log_alphag_fg_xij[, g] <<- log(alpha_g) + log_fg_xij[, g] # [nxg]
      }
      log_alphag_fg_xij <<- pmin(log_alphag_fg_xij, log(.Machine$double.xmax))
      log_alphag_fg_xij <<- pmax(log_alphag_fg_xij, log(.Machine$double.xmin))

      h_ig <<- exp(lognormalize(log_alphag_fg_xij))

      log_lik <<- sum(log(rowSums(exp(log_alphag_fg_xij))))
    }

  )
)
