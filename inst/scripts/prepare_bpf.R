## ---------------------------------------------------------------------
##  Prepare the Barometre Politique Francais 2006-2007 (wave 1) for
##  weighted_level().
##
##  The raw file CANNOT be redistributed. Obtain it from the Centre de
##  donnees socio-politiques (CDSP), Sciences Po, Paris, under its access
##  conditions:  https://cdsp.sciences-po.fr/
##
##  Place the extracted file as 'datavote.txt' in the working directory.
## ---------------------------------------------------------------------
library(cumulevel)

raw <- read.table("datavote.txt", header = TRUE)

## Exposure: Diplome, 5 ordered levels; target J = 3 (baccalaureat).
## Outcome : declared sympathy for the far-right candidate.
## Covariates follow the back-door set of Rynkiewicz et al. (2016),
## augmented with demographics.
d <- data.frame(
  Y      = raw$Lepen,
  E      = raw$Diplome,
  glob   = factor(raw$Mondialisation),
  age    = factor(raw$age),
  sex    = factor(raw$sexe),
  region = factor(raw$region),
  occ    = factor(raw$Situ.prof),
  relig  = factor(raw$religion))

fit <- weighted_level(Y ~ E | glob + age + sex + region + occ + relig,
                      data = d, J = 3, seed = 2016)
print(summary(fit))
plot(fit)

## Sensitivity over attribution assumptions
for (s in c("agnostic", "dilution", "inheritance"))
  print(mixture_interval(fit, attribution_prior(s), S = 2000))

## Bracketing the causal status of the globalisation variable:
## re-run WITHOUT it. If it is a mediator rather than a confounder,
## these are the total-effect estimates.
fit_no_glob <- weighted_level(Y ~ E | age + sex + region + occ + relig,
                              data = d, J = 3, seed = 2016)
print(fit_no_glob)
