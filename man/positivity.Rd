\name{positivity}
\alias{positivity}
\alias{print.positivity}
\alias{plot.positivity}
\title{Positivity diagnostic for a cumulative multilevel exposure}
\usage{
positivity(object, data = NULL, J = NULL, cap = 0.025, ...)
}
\arguments{
\item{object}{A fitted \code{"cumulevel"} object, or a formula.}
\item{data}{A data frame, when \code{object} is a formula.}
\item{J}{Target level, when \code{object} is a formula.}
\item{cap}{Capping constant to compare against.}
\item{...}{Passed to \code{\link{weighted_level}}.}
}
\value{
An object of class \code{"positivity"}: quantiles of \eqn{\hat G_<}, the share
below \code{cap}, the control-side effective sample size, and a verdict.
}
\description{
Computes and grades the distribution of the estimated sub-threshold propensity
\eqn{\hat G_<(W) = \hat P(E < J \mid W)}, the single quantity on which the
identification of \eqn{\Psi_J(\pi)} rests.
}
\details{
Only \emph{one} positivity condition is required, and it concerns the
\emph{sum} \eqn{G_<}, not the individual sub-threshold probabilities: the
reference regime must be feasible in every covariate stratum. No positivity at
all is required at the levels above \eqn{J}, because they enter weighted by
their own occupation probabilities, so a stratum in which some level is
unattainable contributes nothing rather than dividing by zero.

The verdict is graded as follows.
\describe{
  \item{comfortable}{No observation below \code{cap}, first percentile at least
    four times \code{cap}. The cap never binds.}
  \item{acceptable}{No observation below \code{cap}, but the first percentile
    is close to it. Report the distribution alongside estimates.}
  \item{strained}{Some observations below \code{cap}. Condition (i) of the
    asymptotic theorem fails: the capping term is then first order in the
    outcome-model error, so consistency survives only through a correctly
    specified outcome model. Vary \code{cap} and report the sensitivity.}
}

The control-side effective sample size is reported because the treated-side
figure cannot show it: under positivity strain it is the control side that
collapses first, and schemes with high \eqn{\pi} are affected most.
}
\examples{
data(edu_synth)
p <- positivity(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
p
plot(p)
}
