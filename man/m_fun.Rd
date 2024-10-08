% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m_fun_function.R
\name{m_fun}
\alias{m_fun}
\title{Score function for the working propensity score.}
\usage{
m_fun(theta, Xmat, Zvec, r)
}
\arguments{
\item{theta}{The parameters of the working propensity score model.}

\item{Xmat}{The covariate matrix where rows correspond to individuals and
columns to covariates. The first column is equal to 1 for the intercept.}

\item{Zvec}{The observed treatment of the units.}

\item{r}{Numeric. Corresponds to the ratio of the treatment probability over
1 - treatment probability, where treatment probability is the probability
used to assign the cluster-level treatment.}
}
\description{
Function that gets the score function of the working prepensity score based
on the values of theta. This function uses the mathematical expression for
the derivative of the log likelihood. An equivalent computational approach
to the score function could be based on numerical gradients and could be
implemented as:
}
\details{
numDeriv::grad(function(theta_pars)
                 log_lik(theta_pars, Xmat = Xmat, Zvec = Zvec, r = r),
               x = theta).
}
