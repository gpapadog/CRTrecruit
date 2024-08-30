#' Estimation and inference on the incentivized subpopulation
#' 
#' Function that receives quantities on the always+incentivized population and
#' the always population, and returns an estimate of the causal effect and a
#' confidence interval for the incentivized population. It can only be used in
#' the absence of disincentivized-recruited individuals.
#' 
#' @param prop_incen The proportion of incentivized among the always-recruited
#' and the incentivized-recruited. Numeric.
#' @param effect_estimates The effect estimates for the always-recruited and
#' the always+incentivized-recruited population. Vector of length 2.
#' @param asym_var_est 2 x 2 matrix for the asymptotic variance of the causal 
#' estimators for the effect on the always and the always+incentivized
#' populations. Defaults to NULL. If null, inference will not be performed, and
#' only the estimate will be returned.
#'    
#' @export
#' 
ComplierEffect <- function(prop_incen, effect_estimates, asym_var_est = NULL) {
  
  r <- c(estimate = NA, asym_sd = NA, asym_lb = NA, asym_ub = NA)
  
  # Weights corresponding to the always, and always+incentivized populations.
  use_weights <- c(alw = - (1 - prop_incen), alw_inc = 1) / prop_incen
  
  # Estimate: (combination of effect on always+incentivized and effect on 
  # incentivized)
  r[1] <- weighted.mean(x = effect_estimates, w = use_weights)
  
  if (is.null(asym_var_est)) {
    r <- r[1]
    return(r)
  }
  
  # Asymptotic standard deviation by using the delta method:
  use_weights <- as.vector(use_weights)
  r[2] <- sqrt(t(use_weights) %*% asym_var_est %*% use_weights)
  
  # Confidence intervals:
  r[3 : 4] <- r[1] + 1.96 * c(- 1, 1) * r[2]
  
  return(r)
}
