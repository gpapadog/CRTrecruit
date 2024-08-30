#' Estimation and inference of causal effects in cluster randomized trials
#' 
#' Estimating the causal effect on the subpopulation of the always-recruited
#' and disincentivized, the always- and the incentivized-recruited, and the
#' recruited.
#' 
#' Builds on the NewEstimateCRT function that does estimation and inference for
#' always+disincentivized and recruited only.
#' 
#' @param Zobs Vector of treatments for the recruited individuals.
#' @param Yobs Vector of outcomes for the recruited individuals.
#' @param weights Matrix of weights for each individual used for estimation.
#' Columns correspond to weights for the always-recruited, for the always+
#' incentivized-recruited, and for the recruited.
#' @param IDobs Vector of the cluster indices for the recruited individuals.
#' @param ps Whether the propensity score is known or estimated. Options are
#' 'known' or 'estimated'. Defaults to 'known'.
#' @param est_wps_pars If the PS is estimated, est_wps_pars is needed. This
#' is the output of the EstWPSpars function used to estimate the parameters.
#' This output will include quantities needed for variance estimation in the
#' A and V matrices.
#' @param inference Logical. Whether asymptotic variance and inference is
#' performed. Defaults to TRUE.
#' @param treat_prop Proportion of clusters that are treated. Defaults to 0.5.
#' 
#' @export
#'
SimultEstimateCRT <- function(Zobs, Yobs, weights, IDobs, Xobs = NULL,
                              ps = c('known', 'estimated'),
                              est_wps_pars = NULL, inference = TRUE,
                              treat_prop = 1 / 2) {
  
  # Note: If we do not need inference (inference is set to FALSE), PS set to
  # 'estimated' or 'known' does not matter.
  
  ps <- match.arg(ps)
  
  J <- max(IDobs)
  avg_cluster_size <- length(Zobs) / J  # average cluster size
  
  # Treatment randomization ratio! Need to generalize this.
  r <- treat_prop / (1 - treat_prop)
  
  if (ncol(weights) != 3) {
    stop('weights should be of 3 columns corresponding to
          always+disincentivized, always+incentivized, and recruited.')
  }

  # --------------- PART A ------------------- #
  # -------------- Estimation ---------------- #
  
  # theta1hat is the mean of w * Z * Y
  # theta2hat is the mean of w * (1 - Z) * Y
  # theta3hat is the mean of w * Z
  # theta4hat is the mean of w * (1 - Z)
  theta1hat <- colMeans(sweep(weights, MARGIN = 1, STATS = Zobs * Yobs, FUN = '*'))
  theta2hat <- colMeans(sweep(weights, MARGIN = 1, STATS = (1 - Zobs) * Yobs, FUN = '*'))
  theta3hat <- colMeans(sweep(weights, MARGIN = 1, STATS = Zobs, FUN = '*'))
  theta4hat <- colMeans(sweep(weights, MARGIN = 1, STATS = 1 - Zobs, FUN = '*'))

  WT_mu1 <- theta1hat / theta3hat
  WT_mu0 <- theta2hat / theta4hat
  estimate <- WT_mu1 - WT_mu0
  names(estimate) <- c('alw+dis', 'alw+inc', 'recr')
  
  theta_hat <- rbind(theta1hat, theta2hat, theta3hat, theta4hat)
  colnames(theta_hat) <- names(estimate)
  
  # If we do not need inference, the function will stop here and return.
  if (!inference) {
    return(list(estimate = estimate, theta_hat = theta_hat))
  }
  
  
  # The function continues only if inference is set to TRUE.
  
  # --------------- PART B ------------------- #
  # -------------- Inference ----------------- #
  
  if (ps == 'estimated' & is.null(Xobs)) {
    stop('For estimated PS, provide Xobs.')
  }
  
  # Acquiring the propensity score score functions.
  if (ps == 'estimated') {
    est_alpha <- est_wps_pars$est_delta_coefs
    num_alpha <- length(est_alpha)
    psi_alpha <- est_wps_pars$s_clust_est
  }
  
  
  # STEP 1:
  # Creating the matrix V from the asymptotic variance:
  asym_V <- matrix(0, 12, 12)  # For elements for each of 3 populations.
  if (ps == 'estimated') {
    asym_V <- matrix(0, nrow(asym_V) + num_alpha, nrow(asym_V) + num_alpha)
  }
  
  for (j in 1 : J) {
    
    wh_units <- which(IDobs == j)
    
    # Creating the estimating function as a matrix for the four quantities and
    # the three estimates. Afterwards, we convert it to a vector.
    psi <- matrix(NA, 4, 3)
    
    for (pp in 1 : 3) {  # For each estimand:
      psi[, pp] <- c(sum((Zobs * Yobs * weights[, pp])[wh_units]) - theta1hat[pp] * length(wh_units),
                     sum(((1 - Zobs) * Yobs * weights[, pp])[wh_units]) - theta2hat[pp] * length(wh_units),
                     sum((Zobs * weights[, pp])[wh_units]) - theta3hat[pp] * length(wh_units),
                     sum(((1 - Zobs) * weights[, pp])[wh_units]) - theta4hat[pp] * length(wh_units))
    }
    
    psi <- c(psi[, 1], psi[, 2], psi[, 3])
    
    # If the propensity score is estimate, I include the PS score function.
    if (ps == 'estimated') {
      psi <- c(psi, psi_alpha[j, ])
    }
    
    # The contribution of cluster j to the asymptotic matrix V:
    psi <- as.vector(psi)
    asym_V <- asym_V + psi %*% t(psi)
  }
  
  # We calculated the sum over the clusters, so we divide by the # of clusters.
  asym_V <- asym_V / J
  
  
  # STEP 2:
  # Creating the matrix A from the asymptotic variance.
  if (ps == 'known') {
    
    A_mat_inv <- diag(- 1 / avg_cluster_size, nrow = nrow(asym_V))
    
  } else {  # PS is estimated.
    
    A_mat <- matrix(0, nrow(asym_V), ncol(asym_V))
    for (qq in 1 : 12) {
      A_mat[qq, qq] <- - avg_cluster_size
    }
    
    # Bottom right of the A matrix is the Hessian from the PS estimation.
    wps_par_indices <- 13 : (12 + num_alpha)
    A_mat[wps_par_indices, wps_par_indices] <- est_wps_pars$A_mat
    
    
    # ----- Calculating the derivatives of the weights.
    
    # Useful quantity for all derivatives
    exp_Xt_alpha <- exp(Xobs %*% est_alpha)
    
    # (a) For the always-recruited:
    
    # (a1) The derivative of w1 for always-recruited:
    deriv_w1 <- sweep(t(Xobs), MARGIN = 2, FUN = '*',
                      STATS = c(- exp_Xt_alpha / (r * (1 + exp_Xt_alpha) ^ 2)))
    Z_deriv_w1 <- sweep(deriv_w1, 2, FUN = '*', STATS = Zobs)
    ZY_deriv_w1 <- sweep(Z_deriv_w1, 2, FUN = '*', STATS = Yobs)
    
    # Always-recruited for treated units - rows 1 & 3:
    A_mat[c(1, 3), wps_par_indices] <- rbind(apply(ZY_deriv_w1, 1, sum),
                                                 apply(Z_deriv_w1, 1, sum)) / J
    
    # (a0) The derivative of w0 for always-recruited is zero.
    # Always-recruited for control units remains at 0 (rows 2 & 4)
    
    # (b) For always + incentivized:
    
    # (b1) For the treated units the derivative is 0 - rows 5 & 7 remain at 0.
    # (b0) The derivative of w0 for always + incentivized:
    deriv_w0 <- r * sweep(t(Xobs), 2, FUN = '*', exp_Xt_alpha)
    Z_deriv_w0 <- sweep(deriv_w0, 2, FUN = '*', STATS = 1 - Zobs)
    ZY_deriv_w0 <- sweep(Z_deriv_w0, 2, FUN = '*', STATS = Yobs)
    
    # Always + incentivized for the control units:
    A_mat[c(6, 8), wps_par_indices] <- rbind(apply(ZY_deriv_w0, 1, sum),
                                             apply(Z_deriv_w0, 1, sum)) / J
    
    # (c) For the recruited:
    
    # (c1) The derivative of w1 is the same for always and recruited
    #      so we use what we have for the always-recruited.

    # Recruited for treated units - rows 9 & 11:
    A_mat[c(9, 11), wps_par_indices] <- A_mat[c(1, 3), wps_par_indices]
    
    # (c0) The derivative of w0 is the same for recruited and always+
    #      incentivized so we use what we have for the always+incentivized.
    
    # Recruited for control units - rows 10 & 12:
    A_mat[c(10, 12), wps_par_indices] <- A_mat[c(6, 8), wps_par_indices]
    
    
    # Inverting the matrix A:
    A_mat_inv <- solve(A_mat)
    
  }
  
  # STEP 3
  # Putting A and V together for the estimated variance:
  
  asym_var <- A_mat_inv %*% asym_V %*% t(A_mat_inv)
  asym_var <- asym_var / J
  
  
  # STEP 4
  # Using the Delta method to get the variance of the estimator:
  
  # If ps is estimated the entries for the PS are set to 0.
  anadelta <- matrix(0, nrow(asym_V), 3)
  for (pp in 1 : 3) {
    wh_rows <- 4 * (pp - 1) + 1 : 4
    anadelta[wh_rows, pp] <- c(1 / theta3hat[pp], - 1 / theta4hat[pp],
                               - theta1hat[pp] / (theta3hat[pp] ^ 2),
                               theta2hat[pp] / (theta4hat[pp] ^ 2))
  }
  
  asym_var_est <- (t(anadelta) %*% asym_var %*% anadelta)
  rownames(asym_var_est) <- names(estimate)
  colnames(asym_var_est) <- names(estimate)
  
  lb <- estimate - 1.96 * sqrt(diag(asym_var_est))
  ub <- estimate + 1.96 * sqrt(diag(asym_var_est))
  
  
  
  return(list(estimate = estimate, asym_var_est = asym_var_est,
              lb = lb, ub = ub,
              theta_hat = theta_hat, asym_var = asym_var))
  
}