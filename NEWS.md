# cumulevel 0.1.0

* First public release.
* `weighted_level()`: level-specific weighted TMLE for cumulative multilevel exposures, doubly robust and efficient under cross-fitting. Reports the exact-level, geometric and binary members of the family by default, with the attribution spread.
* `positivity()`: positivity diagnostic computed on every call.
* `attribution_prior()` and `mixture_interval()`: probabilistic sensitivity analysis over the attribution weights.
* Data: `nhanes_smoking` (NHANES 2005-2006 linked to NCHS mortality, n = 4,329) and `edu_synth` (synthetic).
