#' Estimation and inference on the complier subpopulation
#' 
#' Function that receives quantities on the always+complier population and the
#' always population, and returns an estimate of the causal effect and a
#' confidence interval for the complier population. It can only be used in the
#' absence of defiers.
#' 
#' @param prop_compl The proportion of compliers among the always-recruited and
#' the complier recruited. Numeric.
#' @param effect_estimates The effect estimates for the always-recruited and
#' the always+complier population. Vector of length 2.
#' @param asym_var_est 2 x 2 matrix for the asymptotic variance of the causal 
#' estimators for the effect on the always and the always+complier populations.
#' Defaults to NULL. If null, inference will not be performed, and only the
#' estimate will be returned.
#'    
#' @export
#' 
ComplierEffect <- function(prop_compl, effect_estimates, asym_var_est = NULL) {
  
  r <- c(estimate = NA, asym_sd = NA, asym_lb = NA, asym_ub = NA)
  
  # The weights corresponding to the always, and always+complier populations.
  use_weights <- c(alw = - (1 - prop_compl), alw_com = 1) / prop_compl
  
  # Estimate: (combination of effect on always+compliers and effect on compliers)
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
