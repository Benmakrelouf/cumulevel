\name{nhanes_smoking}
\alias{nhanes_smoking}
\docType{data}
\title{NHANES 2005--2006 pack-years and 10-year mortality}
\format{A data frame with 4329 rows and 9 variables:
\describe{
  \item{seqn}{NHANES respondent sequence number.}
  \item{died}{Binary outcome: all-cause mortality within 10 years of
    examination.}
  \item{packyears}{Cumulative tobacco exposure on 5 ordered levels:
    1 = never smoker, 2 = 1--10 pack-years, 3 = 11--20, 4 = 21--40,
    5 = more than 40. Level 3 is the target in the reference analysis.}
  \item{age}{Age in years at examination, 20 to 85 (top-coded).}
  \item{sex}{Factor: Female, Male.}
  \item{ethnicity}{Factor with 5 levels.}
  \item{education}{Factor with 6 levels, from less than 9th grade to college
    graduate, plus \code{"Refused"} (2 respondents).}
  \item{pir}{Family income-to-poverty ratio, 0 to 5 (top-coded).}
  \item{bmi}{Body mass index.}
}}
\source{
National Health and Nutrition Examination Survey 2005--2006, linked to the
National Center for Health Statistics 2019 Public-Use Linked Mortality File.
Both are works of the United States government and are in the public domain.
\url{https://wwwn.cdc.gov/nchs/nhanes/} and
\url{https://www.cdc.gov/nchs/data-linkage/mortality-public.htm}
}
\usage{data(nhanes_smoking)}
\description{
Adults aged 20 and over from the NHANES 2005--2006 cycle, linked to mortality
follow-up, with tobacco exposure coded as an ordered cumulative variable.
Complete cases on the analysis variables.
}
\details{
Pack-years are cumulative by construction: reaching 21--40 pack-years entails
having passed through the lighter categories. Dichotomising at a threshold
therefore attributes to that threshold whatever happens above it, which is the
problem the package addresses.

The 2005--2006 cycle is used because the 2019 linked file gives every
participant at least about 156 months of follow-up, so no respondent is
administratively censored before the ten-year horizon and the binary endpoint
is legitimate.

Crude mortality rises monotonically across levels: 0.090, 0.101, 0.146, 0.196,
0.362. The analysis below shows that this gradient is not evenly attributable.

\strong{Two caveats.} The analysis is unweighted for the NHANES complex survey
design; a design-based extension is straightforward but not implemented here.
And two respondents have an estimated sub-threshold propensity below the
default capping constant, so \code{weighted_level} reports the positivity
verdict as strained; the estimates are essentially unchanged when the cap is
varied, but the diagnostic should be reported.
}
\examples{
data(nhanes_smoking)
tapply(nhanes_smoking$died, nhanes_smoking$packyears, mean)

## The conclusion depends on the attribution assumption
fit <- weighted_level(died ~ packyears | age + sex + ethnicity + pir + bmi,
                      data = nhanes_smoking, J = 3, seed = 2005)
confint(fit)
}
\keyword{datasets}
