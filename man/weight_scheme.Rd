\name{weight_scheme}
\alias{weight_scheme}
\title{Contribution weight schemes}
\usage{
weight_scheme(scheme = c("exact", "geometric", "binary"), J, K, rho = 0.5)
}
\arguments{
\item{scheme}{\code{"exact"}, \code{"geometric"}, \code{"binary"}, or a numeric
vector of length \code{K - J + 1}.}
\item{J}{Target level.}
\item{K}{Number of exposure levels.}
\item{rho}{Decay rate for the geometric scheme.}
}
\value{A numeric vector of length \code{K - J + 1} whose first entry is 1.}
\description{
Builds the weight vector \eqn{(\pi(J), \ldots, \pi(K))}. The anchor
\eqn{\pi(J) = 1} is imposed by every scheme: individuals who stop exactly at
the target level always contribute in full.
}
\details{
\code{"exact"} gives the restricted average treatment effect on
\eqn{\{E = J\}}; \code{"binary"} the binary average treatment effect on
\eqn{\{E \ge J\}}. These are the two boundaries of the weight space, and both
are causally valid: they answer different questions. \code{"geometric"}
interpolates with \eqn{\pi(J + k) = \rho^k}.
}
\examples{
weight_scheme("geometric", J = 3, K = 5, rho = 0.5)
weight_scheme("exact", J = 3, K = 5)
}
