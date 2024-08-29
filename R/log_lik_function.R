#' Working propensity score log likelihood.
#' 
#' Function that calculates the log-likelihood of the treatment among the
#' recruited according to the working propensity score model.
#' 
#' @param theta The parameters of the working propensity score model.
#' @param Xmat The covariate matrix where rows correspond to individuals and
#' columns to covariates. The first column is equal to 1 for the intercept.
#' @param Zvec The observed treatment of the units.
#' @param r Numeric. Corresponds to the ratio of the treatment probability over
#' 1 - treatment probability, where treatment probability is the probability
#' used to assign the cluster-level treatment.
#' 
log_lik <- function(theta, Xmat, Zvec, r) {
  
  # The conditional probability that Z = 1 among the recruited:
  p_Z1 <- (1 + (r * (1 + exp(Xmat %*% theta))) ^ (- 1)) ^ (-1)
  
  # The log-likelihood
  log_lik <- Zvec * log(p_Z1) + (1 - Zvec) * log(1 - p_Z1)
  
  return(sum(log_lik))
}
