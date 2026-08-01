\name{weighted_level}
\alias{weighted_level}
\title{Weighted causal estimands for a cumulative multilevel exposure}
\usage{
weighted_level(formula, data, J, pi = "family", rho = 0.5,
               learner = c("parametric", "forest"), V = 5L,
               cap = 0.025, fluctuation = c("weighted", "covariate"),
               seed = NULL)
}
\arguments{
\item{formula}{A formula \code{Y ~ E | W1 + W2 + ...}: binary outcome, ordered
cumulative exposure coded \code{1, ..., K}, and adjustment covariates.}
\item{data}{A data frame.}
\item{J}{Target level, an integer in \code{2:K}.}
\item{pi}{\code{"family"} (default) to estimate the exact, geometric and
binary members together; or a scheme name; or a numeric weight vector.}
\item{rho}{Decay rate for the geometric scheme.}
\item{learner}{\code{"parametric"} or \code{"forest"} (needs \pkg{ranger}).}
\item{V}{Number of cross-fitting folds.}
\item{cap}{Lower bound on the estimated sub-threshold propensity.}
\item{fluctuation}{\code{"weighted"} (default) or \code{"covariate"}.}
\item{seed}{Optional seed for the fold assignment.}
}
\value{
An object of class \code{"cumulevel"} with components \code{estimates},
\code{positivity}, \code{spread_pct}, \code{nuisances} and the call.
}
\description{
Estimates the family of weighted causal estimands \eqn{\Psi_J(\pi)} for an
exposure that is cumulative and hierarchical, by a level-specific weighted
targeted minimum loss-based estimator with cross-fitted nuisances.
}
\details{
\strong{Why the default returns a family.} When an exposure is cumulative,
dichotomising it at \eqn{J} attributes to level \eqn{J} whatever happens at the
levels above it. The restricted and binary estimands are not right and wrong:
they are the two boundaries of a continuum of attribution assumptions, and the
gap between them is the quantity a reader needs. Reporting a single member
hides that gap.

\strong{Positivity is checked automatically.} The distribution of
\eqn{\hat G_<(W) = \hat P(E < J \mid W)} is returned, and a warning is issued
when any observation falls below \code{cap}. This matters because the capping
term is second order when the outcome model is consistent, but becomes
\emph{first order} in the outcome-model error when the cap binds.

\strong{Why \code{J = 1} is excluded.} The sub-threshold set would be empty, the
reference regression a \eqn{0/0} ratio in every stratum, and the balancing
weight undefined.

\strong{The two fluctuation forms.} Both impose the same score equation and are
valid targeted estimators, but they solve it along different one-dimensional
paths and are therefore \emph{not} the same estimator in finite samples. They
are asymptotically equivalent; the weighted form is markedly stabler when the
balancing weight is large.
}
\references{
Benmakrelouf, M. R. and Tatachak, A. (2026). Weighted causal inference for
cumulative multilevel exposures: a targeted learning approach.
}
\seealso{\code{\link{positivity}}, \code{\link{mixture_interval}},
\code{\link{weight_scheme}}}
\examples{
data(edu_synth)
fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
fit
summary(fit)
plot(fit)

## Tobacco, real NHANES data: the decomposition changes the conclusion.
## The exact-level interval covers zero; the binary interval does not.
data(nhanes_smoking)
f2 <- weighted_level(died ~ packyears | age + sex + ethnicity + pir + bmi,
                     data = nhanes_smoking, J = 3, seed = 2005)
confint(f2)
}
