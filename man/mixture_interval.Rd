\name{mixture_interval}
\alias{mixture_interval}
\alias{print.mixture}
\title{Mixture interval over attribution assumptions}
\usage{
mixture_interval(object, prior = attribution_prior("agnostic"),
                 S = 1000L, B = 200L, level = 0.95)
}
\arguments{
\item{object}{A fitted \code{"cumulevel"} object.}
\item{prior}{An \code{"attribution_prior"} object.}
\item{S}{Number of prior draws.}
\item{B}{Number of inner-layer replicates.}
\item{level}{Confidence level.}
}
\value{
An object of class \code{"mixture"}: the mixture point estimate, the interval,
and the attribution share of the total variance.
}
\description{
Propagates attribution uncertainty jointly with sampling uncertainty, using a
two-layer Monte Carlo scheme with a common multiplier in the inner layer.
}
\details{
The variance decomposes exactly into a sampling component and an attribution
component. The first shrinks at rate \eqn{1/n}; the second is fixed by the
prior. The \emph{attribution share} therefore converges to one: at large sample
sizes what remains is the irreducible cost of not knowing \eqn{\pi}.

The inner layer uses a common multiplier because all members are computed on
the same sample from the same nuisance fits, so their estimation errors are
strongly dependent. Drawing independently across members would have the correct
variance but the wrong joint law.

\strong{This interval is prior-predictive, not a confidence interval.} It is
calibrated for a member of the family drawn from the prior, and
\emph{conservative} for the fixed mixture estimand, increasingly so as the
sample grows. Do not report it as a confidence interval for the mixture
estimand.
}
\examples{
data(edu_synth)
fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
mixture_interval(fit, attribution_prior("agnostic"), S = 100, B = 50)
}
