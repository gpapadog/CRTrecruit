#' Bootstrap inference in CRTs with recruitment bias
#' 
#' Resamples clusters of the CRT and re-estimates quantities of interest in
#' order to acquire bootstrap-based inference.
#' 
#' @param B Number of bootstrap samples.
#' @param Xobs Matrix with covariates for the recruited individuals. Rows
#' correspond to the individual, and columns to the different covariate. The
#' first column is a column of 1s for the intercept of the model.
#' @param Zobs Vector of treatments for the recruited individuals.
#' @param Yobs Vector of outcomes for the recruited individuals.
#' @param IDobs Vector of the cluster indices for the recruited individuals.
#' IDobs should include integers from 1 up to the number of different clusters.
#' @param treat_prop Proportion of clusters that are treated.
#' @param fix_trt_num Logical. Whether we want to keep the number of treated
#' clusters fixed in the bootstrap samples. This should agree with the
#' treatment assignment mechanism. Defaults to TRUE.
#' @param prop_incen The proportion of incentivized-recruited. Bootstrap can be
#' performed without re-estimating the proportion of incentivized-recruited
#' within each bootstrap sample. Then, set prop_incen to the estimated value.
#' Defaults to NULL. If left NULL, this will not be calculated. We suggest
#' re-estimating the proportion of incentivized-recruited. based on the
#' bootstrap sample (leaving this to NULL). In simulations, we have found
#' similar inferential performance whether this is re-estimated or not.
#' @param disincentivized Logical. Whether we allow for disincentivized-
#' recruited to exist. If they exist, we do not estimate effects on the
#' incentivized-recruited. Defaults to FALSE.
#' @param verbose Logical. Whether to print bootstrap progression. Defaults to
#' TRUE.
#' 
#' @export
BootVar <- function(B, Xobs, Zobs, Yobs, IDobs, treat_prop, fix_trt_num = TRUE,
                    prop_incen = NULL, disincentivized = FALSE,
                    verbose = TRUE) {
  
  print('For estimated propensity score only.')
  
  num_alphas <- ncol(Xobs)  # Number of covariates + intercept
  J <- max(IDobs)  # Number of clusters
  
  # The cluster level treatment
  Zc <- unique(cbind(Zobs, IDobs))  
  Zc <- Zc[order(Zc[, 2]), ]  # Re-ordering according to cluster indices.
  Zc <- Zc[, 1]
  
  
  # --------- PART 1: Where bootstrap results will be saved ----------- #
  
  # Keeping track of the estimates across bootstrap samples.
  est_boot <- array(NA, dim = c(B, 5))
  dimnames(est_boot) <- list(
    boot = 1 : B,
    population = c('always+disincentivized', 'recruited',
                   'always+incentivized', 'incentivized',
                   'incentivized_propfix'))
  
  # Keeping track of the estimates of the wps across bootstraps.
  alpha_boot <- matrix(NA, nrow = B, ncol = num_alphas)
  
  # Keeping track of the bootstrap estimated proportion of incentivized-
  # recruited.
  if (!disincentivized) {
    prop_incen_boot <- rep(NA, B)
  }
  
  
  for (bb in 1 : B) {
    
    if (verbose & bb %% 100 == 0) print(bb)
    
    # --------- PART 2: Creating the bootstrap sample ----------- #
    
    if (fix_trt_num) {
      chosen_clusters <- c(sample(which(Zc == 1), sum(Zc == 1), replace = TRUE),
                           sample(which(Zc == 0), sum(Zc == 0), replace = TRUE))
    } else {
      chosen_clusters <- sample(1 : J, J, replace = TRUE)
    }
    
    boot_Zobs <- NULL
    boot_Yobs <- NULL
    boot_Xobs <- matrix(NA, nrow = 0, ncol = ncol(Xobs))
    boot_IDobs <- NULL
    
    for (j in 1 : J) {
      wh_units <- which(IDobs == chosen_clusters[j])
      boot_Zobs <- c(boot_Zobs, Zobs[wh_units])
      boot_Yobs <- c(boot_Yobs, Yobs[wh_units])
      boot_Xobs <- rbind(boot_Xobs, Xobs[wh_units, ])
      boot_IDobs <- c(boot_IDobs, rep(j, length(wh_units)))
    }
    
    
    # ------ PART 3: Effect estimates except on incentivized-recruited ------ #
    
    
    # Re-estimating the propensity score:
    boot_wps_pars <- EstWPSpars(Xobs = boot_Xobs, Zobs = boot_Zobs,
                                IDobs = boot_IDobs, treat_prop = treat_prop,
                                inference = FALSE)
    alpha_boot[bb, ] <- boot_wps_pars$est_delta_coefs
    
    
    # Acquiring the estimates of the working propensity score.
    boot_est_delta <- 1 + exp(boot_Xobs %*% alpha_boot[bb, ])
    boot_est_eind <- r * boot_est_delta / (1 + r * boot_est_delta)
    
    boot_weights <-
      cbind(alw_dis = ((1 - boot_est_eind) / boot_est_eind) * boot_Zobs + (1 - boot_Zobs),
            alw_inc = boot_Zobs + (1 - boot_Zobs) * boot_est_eind / (1 - boot_est_eind),
            recr = boot_Zobs / boot_est_eind + (1 - boot_Zobs) / (1 - boot_est_eind))
    
    boot_est <- SimultEstimateCRT(Zobs = boot_Zobs, Yobs = boot_Yobs,
                                  weights = boot_weights, IDobs = boot_IDobs,
                                  Xobs = boot_Xobs, treat_prop = treat_prop,
                                  inference = FALSE)
    
    est_boot[bb, c(1, 3, 2)] <- boot_est$estimate
    
    
    
    # --------- PART 4: Effect estimates on incentivized-recruited ----------- #
    
    if (!disincentivized) {
      
      boot_pi_trt <- mean(Zc[chosen_clusters])  # Proportion of treated clusters.
      boot_p_trt <- mean(boot_Zobs)  # Proportion of treatment in the recruited.
      boot_est_prop_incen <- 1 - ((boot_pi_trt / (1 - boot_pi_trt)) *
                                    ((1 - boot_p_trt) / boot_p_trt))
      prop_incen_boot[bb] <- boot_est_prop_incen
      
      est_boot[bb, 4] <- ComplierEffect(prop_incen = boot_est_prop_incen,
                                        effect_estimates = est_boot[bb, 1 : 2])
      
      if (!is.null(prop_incen)) {
        est_boot[bb, 5] <- ComplierEffect(prop_incen = prop_incen,
                                          effect_estimates = est_boot[bb, 1 : 2])
      }
    }
  }
  
  
  # --------- PART 5: Organizing results ----------- #
  
  # If we have disincentivized-recruited, then we do not have the results on
  # incentivized-recruited.
  if (disincentivized) {
    est_boot <- est_boot[, 1 : 3]
    return(list(alpha_boot = alpha_boot, est_boot = est_boot))
  }
  
  if (!disincentivized) {
    
    # If we do not have disincentivized-recruited, we might not fix the
    # proportion of incentivized-recruited.
    if (is.null(prop_incen)) {
      est_boot <- est_boot[, 1 : 4]
    }
    
    return(list(alpha_boot = alpha_boot, prop_incen_boot = prop_incen_boot,
                est_boot = est_boot))
  }
  
}