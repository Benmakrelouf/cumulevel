\name{attribution_prior}
\alias{attribution_prior}
\title{Monotone prior on the attribution weights}
\usage{
attribution_prior(stance = c("agnostic", "dilution", "inheritance"),
                  a = NULL, b = NULL)
}
\arguments{
\item{stance}{\code{"agnostic"} (Beta(1,1)), \code{"dilution"} (Beta(2,4),
prior mean 1/3), or \code{"inheritance"} (Beta(4,2)).}
\item{a, b}{Beta hyperparameters, overriding \code{stance}.}
}
\value{A list of class \code{"attribution_prior"}.}
\description{
A one-parameter prior on the decay rate \eqn{\rho}, with
\eqn{\pi(J + k) = \rho^k}. Monotonicity holds for every draw by construction,
which independent priors on the components do not ensure.
}
\examples{
attribution_prior("dilution")
attribution_prior(a = 3, b = 3)
}
