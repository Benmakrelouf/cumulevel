\name{edu_synth}
\alias{edu_synth}
\docType{data}
\title{Synthetic education and far-right voting data}
\format{A data frame with 5000 rows and 5 variables:
\describe{
  \item{Y}{Binary outcome: declared sympathy for the far-right candidate.}
  \item{E}{Cumulative exposure, educational attainment on 5 ordered levels;
    level 3 is the baccalaureate.}
  \item{W1}{Opinion on globalisation, 3 categories coded 0--2.}
  \item{W2}{Age band, 5 categories.}
  \item{W3}{Sex, binary.}
}}
\source{Simulated. See \code{data-raw.R} in the package sources.}
\usage{data(edu_synth)}
\description{
A synthetic data set reproducing the \emph{structure} and the dilution
gradient of the education application, for use in examples.
}
\details{
\strong{These are not the survey data.} The application in the reference paper
uses the first wave of the \emph{Barometre Politique Francais} 2006--2007,
archived at the Centre de donnees socio-politiques (Sciences Po, Paris). That
file cannot be redistributed under its access conditions. This data set is
simulated to exhibit the same qualitative pattern: crude outcome rates
declining across levels, and a dilution gradient such that dichotomising at
level 3 over-attributes.

To reproduce the published analysis, obtain the survey file from the CDSP and
use \code{system.file("scripts", "prepare_bpf.R", package = "cumulevel")}.
}
\examples{
data(edu_synth)
table(edu_synth$E)
tapply(edu_synth$Y, edu_synth$E, mean)
}
\keyword{datasets}
