#' Score function for the working propensity score.
#' 
#' Function that gets the score function of the working prepensity score based
#' on the values of theta. This function uses the mathematical expression for
#' the derivative of the log likelihood. An equivalent computational approach
#' to the score function could be based on numerical gradients and could be
#' implemented as:
#' 
#' numDeriv::grad(function(theta_pars)
#'                  log_lik(theta_pars, Xmat = Xmat, Zvec = Zvec, r = r),
#'                x = theta).
#'
#' @param theta The parameters of the working propensity score model.
#' @param Xmat The covariate matrix where rows correspond to individuals and
#' columns to covariates. The first column is equal to 1 for the intercept.
#' @param Zvec The observed treatment of the units.
#' @param r Numeric. Corresponds to the ratio of the treatment probability over
#' 1 - treatment probability, where treatment probability is the probability
#' used to assign the cluster-level treatment.
#' 
m_fun <- function(theta, Xmat, Zvec, r) {
  
  # The conditional probability that Z = 1 among the recruited:
  p_Z1 <- (1 + (r * (1 + exp(Xmat %*% theta))) ^ (- 1)) ^ (-1)
  p_Z0 <- 1 - p_Z1
  
  # Creating the cluster sums of the derivative of the log-likelihood.
  if (all(Zvec == 0)) {  # All the units are controls:
    
    # Derivative of the log-likelihood for Z = 0
    der_log_Z0 <- - p_Z0 * r * exp(Xmat %*% theta)
    der_log_Z0 <- sweep(Xmat, MARGIN = 1, STATS = der_log_Z0, FUN = '*')
    s_clust <- apply(der_log_Z0, 2, sum)
    
  } else if (all(Zvec == 1)) {  # All the units are treated:
    
    # Derivative of the log-likelihood for Z = 1
    der_log_Z1 <- (p_Z1 * (r * (1 + exp(Xmat %*% theta))) ^ (- 2) *
                     r * exp(Xmat %*% theta))
    der_log_Z1 <- sweep(Xmat, MARGIN = 1, STATS = der_log_Z1, FUN = '*')
    s_clust <- apply(der_log_Z1, 2, sum)
    
  } else {
    stop('In a cluster, units should be either all controls or all treated.')
  }
  
  return(s_clust)
}