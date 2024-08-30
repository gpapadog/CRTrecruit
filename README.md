CRTrecruit
================
Georgia Papadogeorgou

### Inverse weighting methods to adjust for recruitment bias in cluster randomized trials

This R package performs correction for recruitment bias in cluster
randomized trials where available data are on the recruited units only.

The corresponding paper can be found here:
<https://arxiv.org/pdf/2309.07365>

## Installing CRTrecruit

Installing and using CRTrecruit in Rstudio is straightforward. You will
first need the devtools R package.

    install.packages('devtools')
    library(devtools)
    devtools::install_github("gpapadog/CRTrecruit")

### Why use it:

- Clusters are randomly assigned to treatment or controls. Once
  cluster-level assignment is perfromed, individuals are recruited in
  the study.
- Recruitment is generally not blinded to the treatment level of the
  cluster. This creates post-treatment differences in the treated and
  control recruited populations. These differences render simple
  contrasts of outcomes in treated and control recruited groups not
  causally interpretable.

### What are the situations that CRTrecruit should be used?

- In cases with post-treatment recruitment differences, CRTrecruit can
  be used to estimate causal effects on interpretable subpopulations.
- These subpopulations correspond to the combination of always- and
  defier-recruited, the always- and complier-recruited, and the
  recruited population. When defiers do not exist, such as in the case
  where treatment can be thought to only lead to increase in the
  recruitment indicator, the always- and defier-recruited population is
  the same as the always-recruited population, and we can further
  estimate effects on the complier-recruited population.
- Our approach relies on assumptions on recruitment ignorability and
  positivity.

# Using CRTrecruit

## The data set

The data set need to provide information on

- The treatment indicator of the recruited units `Zobs`
- The outcome of the recruited units `Yobs`
- The covariate matrix of the recruited units `Xobs` where the first
  column is a column of 1s corresponding to intercepts. These covariates
  are assumed to satisfy the ignorability (or the non-differential
  recruitment) assumption.
- Indices of cluster membership for the recruited units `IDobs`. These
  indices should be integers starting at 1 and going to the maximum
  number of clusters.
- The probability of treatment in the cluster assignment `treat_prop`.

A simulated data set is included in the R package. We can load it using

``` r
library(CRTrecruit)
data('dta')

# The data only for those who enrolled:
Zobs <- dta$Zobs
Yobs <- dta$Yobs
Xobs <- dta$Xobs
IDobs <- dta$IDobs
treat_prop <- dta$treat_prop
```

and investigate the data structure

``` r
dta$Zobs[1:55]  # The first 55 units
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
    ## [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

``` r
dta$Yobs[1:10]  # The first 10 units
```

    ##  [1] 0.186996 4.270205 1.561991 4.072359 1.767918 6.252805 4.531032 1.468970
    ##  [9] 3.763320 4.630092

``` r
head(dta$Xobs)
```

    ##                 X1 X2          X3         X4 X5
    ## [1,] 1  0.36594112  0  0.09725489 -0.1656856  1
    ## [2,] 1  0.06528818  1  1.57137462 -0.1656856  1
    ## [3,] 1  0.25733838  0 -0.87871270 -0.1656856  1
    ## [4,] 1 -0.64901008  0  0.03916714 -0.1656856  1
    ## [5,] 1  0.66413570  0 -0.37881680 -0.1656856  1
    ## [6,] 1 -0.79708953  1  0.62253002 -0.1656856  1

``` r
dta$IDobs[1:55]  # The first 55 units
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2
    ## [39] 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3

## The working propensity score

The estimation procedure relies on estimating the working propensity
score. We provide code to estimate it using the parametric approach in
the manuscript. Specifically:

``` r
r <- treat_prop / (1 - treat_prop)

# Estimating the propensity score:
est_wps_pars <- EstWPSpars(Xobs = Xobs, Zobs = Zobs, IDobs = IDobs,
                           treat_prop = treat_prop)

# The estimated parameters of the working propensity score:
alphas <- est_wps_pars$est_delta_coefs

# Acquiring the estimates of the working propensity score.
est_delta <- 1 + exp(Xobs %*% alphas)  # The delta function
est_eind <- r * est_delta / (1 + r * est_delta)  # The working PS.
```

## Estimating the causal effects

- Once the working propensity scores have been acquired, we can create
  the weights that are used for estimating the causal effects in the
  different populations.
- The weights are calculated as follows:

``` r
# The weights for the always+defiers, always+compliers, and recruited.
weights <- cbind(alw_def = ((1 - est_eind) / est_eind) * Zobs + (1 - Zobs),
                 alw_com = Zobs + (1 - Zobs) * est_eind / (1 - est_eind),
                 recr = Zobs / est_eind + (1 - Zobs) / (1 - est_eind))
```

### Effects on three interpretable supopulations

- These weights are fed into the `SimultEstimateCRT` function for
  acquiring the causal effect estimates on the three subpopulation, the
  always- and defier-recruited, the always- and complier-recruited, and
  the recruited populations.

``` r
simult_est <- SimultEstimateCRT(Zobs = Zobs, Yobs = Yobs, weights = weights,
                                IDobs = IDobs, Xobs = Xobs,
                                est_wps_pars = est_wps_pars,
                                ps = 'estimated', treat_prop = treat_prop)
```

- The results include the estimates, asymptotic variance, and confidence
  intervals based on the asymptotic distribution.

``` r
# The estimates:
simult_est$estimate
```

    ##  alw+def  alw+com     recr 
    ## 2.742589 2.869440 2.825215

``` r
# The asymptotic variance:
simult_est$asym_var_est
```

    ##             alw+def     alw+com        recr
    ## alw+def 0.008795245 0.007043695 0.007554217
    ## alw+com 0.007043695 0.010389458 0.009382022
    ## recr    0.007554217 0.009382022 0.008840990

``` r
# The confidence intervals' lower and upper bounds
rbind(lb = simult_est$lb, ub = simult_est$ub)
```

    ##     alw+def alw+com     recr
    ## lb 2.558774 2.66966 2.640923
    ## ub 2.926403 3.06922 3.009507

### Effects on the compliers

- It is common to assume monotonicity of the recruitment indicators,
  implying that recruitment status under treatment can only be equal or
  higher than recruitment status under control. This assumption rules
  out the presence of defiers. This means that the effect on the always-
  and defier-recruited population is simply the effect on the
  always-recruited population.
- This allows us to also estimate the effect on the complier population.
- To do so, we first need to estimate the proportion of compliers in the
  group of always- and complier-recruited units. We can do so as
  follows:

``` r
p_trt <- mean(Zobs)  # Observed proportion of treatment in the recruited.
pi_trt <- treat_prop  # The theoretical treatment proportion

# The estimated proportion of compliers among always + compliers
est_prop_compl <- 1 - (pi_trt / (1 - pi_trt)) * ((1 - p_trt) / p_trt)
```

- Then, we can use the `ComplierEffect` function to acquire the causal
  effect on the complier group.

``` r
compl <- ComplierEffect(prop_compl = est_prop_compl,
                        effect_estimates = simult_est$estimate[1:2],
                        asym_var_est = simult_est$asym_var_est[1 : 2, 1 : 2])
compl
```

    ##  estimate   asym_sd   asym_lb   asym_ub 
    ## 3.0250324 0.1620678 2.7073795 3.3426852

## Inference using the bootstrap

- Bootstrap based confidence intervals can also be acquired.
- We resample clusters, either completely randomly by setting
  `fix_trt_num = FALSE` or by keeping the number of treated clusters
  fixed `fix_trt_num = TRUE` (default).
- If we assume that there are no defiers (`defiers=FALSE`) the bootstrap
  will also calculate bootstrap estimates for the effect on the
  compliers.

``` r
boots <- BootVar(B = 500, Xobs = Xobs, Zobs = Zobs, Yobs = Yobs, IDobs = IDobs,
                 treat_prop = treat_prop, fix_trt_num = TRUE,
                 prop_compl = est_prop_compl, defiers = FALSE, verbose = TRUE)
```

    ## [1] "For estimated propensity score only."
    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500

- We create confidence intervals using the standard deviation of the
  bootstrap samples and using the normal approximation.
- Quantiles could be considered alternatively.

``` r
boot_sd <- apply(boots$est_boot, 2, sd)[1 : 4]
boot_sd <- boot_sd[c(1, 3, 2, 4)]  # Re-ordering.
effects <- c(simult_est$estimate, compliers = as.numeric(compl[1]))

effects
```

    ##   alw+def   alw+com      recr compliers 
    ##  2.742589  2.869440  2.825215  3.025032

``` r
boot_sd
```

    ##   always+defiers always+compliers        recruited        compliers 
    ##       0.07652908       0.06282128       0.06078533       0.07754785

``` r
rbind(lb = effects - 1.96 * boot_sd, ub = effects + 1.96 * boot_sd)
```

    ##     alw+def  alw+com     recr compliers
    ## lb 2.592592 2.746310 2.706076  2.873039
    ## ub 2.892586 2.992569 2.944355  3.177026
