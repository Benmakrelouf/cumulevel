\name{cumulevel-package}
\alias{cumulevel-package}
\alias{cumulevel}
\docType{package}
\title{Weighted Causal Estimands for Cumulative Multilevel Exposures}
\description{
Estimation of a family of weighted causal estimands for exposures that are
cumulative and hierarchical, where attaining a level entails having attained
all preceding ones. The package makes the attribution assumption an explicit
parameter rather than an implicit consequence of dichotomising.
}
\details{
The entry point is \code{\link{weighted_level}}, which by default reports the
whole family rather than a single member, together with an automatic positivity
diagnostic (\code{\link{positivity}}). \code{\link{mixture_interval}} adds a
probabilistic sensitivity layer over the attribution weights.

Two data sets are shipped: \code{\link{nhanes_smoking}}, the real NHANES
2005--2006 pack-years and mortality data used in the reference paper, and
\code{\link{edu_synth}}, a synthetic set reproducing the structure of the
education application, whose survey file cannot be redistributed.
}
\references{
Benmakrelouf, M. R. and Tatachak, A. (2026). Weighted causal inference for
cumulative multilevel exposures: a targeted learning approach.
}
\author{Mohamed Raouf Benmakrelouf and Abdelkader Tatachak}
\keyword{package}
