# cumulevel

**Weighted causal estimands for cumulative multilevel exposures.**

When an exposure is cumulative — educational attainment, disease stage,
pack-years of tobacco — reaching a level entails having reached all the
preceding ones. Dichotomising it at a target level silently attributes to that
level whatever happens above it. `cumulevel` makes the attribution assumption an
explicit parameter instead of an implicit consequence of the coding.

## Install

```r
# install.packages("remotes")
remotes::install_github("benmakrelouf/cumulevel")
```

## The default reports a family, not a number

```r
library(cumulevel)
data(edu_synth)

fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 2026)
fit
```

```
Weighted causal estimands for a cumulative exposure
------------------------------------------------------------------
  target level J = 3 of K = 5   |   n = 5000   |   parametric nuisances, V = 5

    scheme estimate     se                 CI  ESS
     exact  -0.0818 0.0110 [-0.1034, -0.0603]  862
 geometric  -0.0912 0.0096 [-0.1100, -0.0724] 1665
    binary  -0.1010 0.0093 [-0.1192, -0.0827] 2034

  Attribution spread, exact to binary: +23.4%
  This is what a dichotomised analysis would attribute to level 3
  but that originates above it.

  Positivity: COMFORTABLE (min Ghat_< = 0.315, cap = 0.025)
```

The exact-level and binary members are **both causally valid**. They are the two
boundaries of a continuum of attribution assumptions, and they answer different
questions. The gap between them is the quantity a dichotomised analysis
conceals, which is why it is printed by default.

## Positivity is checked without being asked

The distribution of the sub-threshold propensity is computed on every call, and
a warning is raised when the capping constant binds. This is not decoration: the
capping term is second order when the outcome model is consistent, but becomes
**first order** in the outcome-model error when the cap binds.

```r
positivity(fit)
plot(positivity(fit))
```

## When the decomposition changes the conclusion

Real NHANES 2005–2006 data, tobacco pack-years against 10-year mortality, with
the target at 11–20 pack-years:

```r
data(nhanes_smoking)
fit <- weighted_level(died ~ packyears | age + sex + ethnicity + pir + bmi,
                      data = nhanes_smoking, J = 3, seed = 2005)
confint(fit)
```

```
                 2.5 %     97.5 %
exact     -0.010978539 0.05742758
geometric  0.008216004 0.06744590
binary     0.020299575 0.08732606
```

**The exact-level interval covers zero; the binary interval does not.** A
dichotomised analysis at 11–20 pack-years would report a significant excess
mortality — but that significance is borrowed from the heavier strata, not
earned at the target level. Only the level-specific decomposition makes this
visible, and the attribution spread is +132%.

This is the case for the framework in one table.

## What it costs not to choose

```r
mixture_interval(fit, attribution_prior("agnostic"))
```

Reports the attribution share of the total variance. The sampling component
shrinks at rate 1/n; the attribution component is fixed by the prior. At large
sample sizes what remains is the irreducible cost of not knowing the weights.

The interval is **prior-predictive, not a confidence interval**: calibrated for a
member drawn from the prior, conservative for the fixed mixture estimand.

## Real data

`nhanes_smoking` is **real**: NHANES 2005–2006 linked to the NCHS 2019
Public-Use Linked Mortality File, both public domain. Crude mortality by level
reproduces the published figures exactly (0.090, 0.101, 0.146, 0.196, 0.362).
`system.file("scripts", "prepare_nhanes.R", package = "cumulevel")` documents
how it was derived from source, so the coding of pack-years and of the endpoint
can be audited.

`edu_synth` is **synthetic**. The survey file underlying the education
application cannot be redistributed under the access conditions of the archive
that holds it; see
`system.file("scripts", "prepare_bpf.R", package = "cumulevel")` to reproduce
that analysis once the file is obtained from the CDSP (Sciences Po).

## A caution about the inflation factor

`summary()` reports `lambda`, the inflation factor built on the discrepancy
between the targeted and one-step estimators. **It is a flag, not a variance
correction.** It tracks only the portion of the second-order remainder that the
targeting step removed, and is structurally blind to first-order failures: it
can exceed one when nothing needs correcting, and equal one to three decimals
when the bias is severe. Read it alongside the positivity diagnostic, never
instead of it.

## Reference

Benmakrelouf, M. R. and Tatachak, A. (2026). *Weighted causal estimands for
cumulative multilevel exposures: identification and targeted estimation.*
Manuscript in preparation.

## License

GPL-3
