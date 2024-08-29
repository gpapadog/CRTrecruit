#' Working propensity score estimation.
#' 
#' Estimating the parameters of the working propensity score model which
#' expresses the conditional probability of treatment in the recruited sample.
#' 
#' @param Xobs Matrix with covariates for the recruited individuals. Rows
#' correspond to the individual, and columns to the different covariate. The
#' first column is a column of 1s for the intercept of the model.
#' @param Zobs Vector of treatment indicators for the recruited individuals.
#' @param IDobs Vector of integers with the index of the cluster to which the
#' recruited individuals belong. IDobs should include integers from 1 up to
#' the number of different clusters.
#' @param r Numeric. Corresponds to the ratio of the treatment probability over
#' 1 - treatment probability, where treatment probability is the probability
#' used to assign the cluster-level treatment.
#' @param inference Logical. Whether the function should calculate the estimate
#' of the asymptotic variance for the working propensity score parameters.
#' Defaults to TRUE.
#' 
#' @export

EstWPSpars <- function(Xobs, Zobs, IDobs, r, inference = TRUE) {
  
  # Number of clusters:
  J <- length(unique(IDobs))
  
  # Optimizing the log-likelihood to find the estimated theta:
  est_delta_coefs <- optim(par = rep(0, ncol(Xobs)), log_lik,
                           control = list(fnscale = - 1),
                           # Arguments of the log_lik function:
                           Xmat = Xobs, Zvec = Zobs, r = r)$par
  
  if (!inference) {
    return(list(est_delta_coefs = est_delta_coefs))
  }
  
  # Getting the s function per cluster:
  s_clust_est <- matrix(NA, nrow = J, ncol = ncol(Xobs))
  for (cc in 1 : J) {
    s_clust_est[cc, ] <- m_fun(est_delta_coefs,
                               Xmat = Xobs[IDobs == cc, , drop = FALSE],
                               Zvec = Zobs[IDobs == cc],
                               r = r)
    
    # Alternative:
    #   s_clust_est[cc, ] <- numDeriv::grad(function(theta_pars)
    #       log_lik(theta_pars, Xmat = Xobs[IDobs == cc, , drop = FALSE],
    #               Zvec = Zobs[IDobs == cc], r = r),
    #       x = est_delta_coefs)
    
  }
  
  # Asymptotic variance V matrix
  V_mat <- matrix(0, ncol(Xobs),  ncol(Xobs))
  for (cc in 1 : J) {
    V_mat <- V_mat + as.vector(s_clust_est[cc, ]) %*% t(as.vector(s_clust_est[cc, ]))
  }
  V_mat <- V_mat / J
  
  # For the A matrix we need the derivative of the score, which we calculate
  # numerically.
  A_mat <- matrix(0, nrow = ncol(Xobs), ncol = ncol(Xobs))
  for (cc in 1 : J) {
    deriv_clust <- numDeriv::hessian(function(theta)
      log_lik(theta, Xmat = Xobs[IDobs == cc, ], Zvec = Zobs[IDobs == cc],
              r = r),
      x = est_delta_coefs
    )
    A_mat <- A_mat - deriv_clust
  }
  A_mat <- A_mat / J
  A_mat_inv <- solve(A_mat)
  
  # Combining A and V in the asymptotic variance formula:
  asym_var <- A_mat_inv %*% V_mat %*% t(A_mat_inv) / J
  
  
  return(list(est_delta_coefs = est_delta_coefs, s_clust_est = s_clust_est,
              V_mat = V_mat, A_mat = A_mat, asym_var = asym_var))
  
}


