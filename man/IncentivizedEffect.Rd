% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IncentivizedEffect_function.R
\name{IncentivizedEffect}
\alias{IncentivizedEffect}
\title{Estimation and inference on the incentivized subpopulation}
\usage{
IncentivizedEffect(prop_incen, effect_estimates, asym_var_est = NULL)
}
\arguments{
\item{prop_incen}{The proportion of incentivized among the always-recruited
and the incentivized-recruited. Numeric.}

\item{effect_estimates}{The effect estimates for the always-recruited and
the always+incentivized-recruited population. Vector of length 2.}

\item{asym_var_est}{2 x 2 matrix for the asymptotic variance of the causal 
estimators for the effect on the always and the always+incentivized
populations. Defaults to NULL. If null, inference will not be performed, and
only the estimate will be returned.}
}
\description{
Function that receives quantities on the always+incentivized population and
the always population, and returns an estimate of the causal effect and a
confidence interval for the incentivized population. It can only be used in
the absence of disincentivized-recruited individuals.
}
