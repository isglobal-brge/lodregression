% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmLimits.R
\name{glmLimits}
\alias{glmLimits}
\title{glmLimits function}
\usage{
glmLimits(formula, data, LimitVariables, lod = NULL, loq = NULL,
  family = "gaussian", logLimit = FALSE, n.iter = 1000, n.chains = 1,
  alpha.corrected = 0.05, quantiles = unique(c(alpha.corrected/2, 0.025,
  0.975, 1 - (alpha.corrected/2))), quiet = FALSE, progress.bar = "text",
  dic.type = "popt", n.adapt = 1000, ...)
}
\arguments{
\item{formula}{an object of class "formula":  a symbolic description of the model to be fitted.}

\item{data}{a data.frame containing the variables in the model.}

\item{LimitVariables}{character string of names of the variables subjected to a limit.}

\item{lod}{numeric vector with the same length of LimitVariables containing the lower limit values for the variables with LOD (with the same order).}

\item{loq}{optional numeric vector with the same length of LimitVariables containing the upper limit values for the variables with limits (with the same order)}

\item{family}{character with a description of the error distribution to be used in the model. Possible values are "gaussian" (default) for continuous outcomes or "binomial" for binary outcomes.}

\item{logLimit}{a logical indicating if limit variables would be assumed distributed as log-normal.}

\item{n.iter}{number of iterations of the Markov chain to run to be used in rjags::update(), rjags::coda.samples() and rjags::dic.samples().}

\item{n.chains}{the number of parallel chains for the model to be used in rjags::jags.model().}

\item{alpha.corrected}{significance level to be used for the confident intervals.}

\item{quantiles}{numeric vector containing the percentiles in [0,1] to obtain of the estimated parameters.}

\item{progress.bar}{type of progress bar. Possible values are "text" (default), "gui", and "none" to be used in rjags::update().}

\item{dic.type}{type of penalty to use in rjags::dic.samples().}

\item{n.adapt}{the number of iterations for adaptation to be used in rjags::jags.model(). If n.adapt = 0 then no adaptation takes place.}
}
\value{
It returns an object of class "glmLimits" containing a model of class "jags".
}
\description{
This function allows you to adjust multiple Bayesian GLM models with predictors subjected to a limit of detection (LOD) specified with a BUGS-language description of the prior distribution, and a set of data.
}
\examples{
glmLimits(formula, data, LimitVariables, lod = NULL, loq = NULL, family = "gaussian", logLimit = FALSE, n.iter = 1000, n.chains = 1, alpha.corrected = 0.05, quantiles = unique(c(alpha.corrected/2, 0.025, 0.975, 1-(alpha.corrected/2))), quiet = FALSE, progress.bar = "text", dic.type = "popt", n.adapt = 1000, ...)
}
\keyword{GLM}
\keyword{LOD}
\keyword{bayesian}

