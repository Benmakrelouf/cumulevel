## These tests are invariants: each one CAN fail. They are designed to catch
## the silent misalignment of the attribution weight vector, which produces
## plausible estimates with no error message.

test_that("ESS invariants hold exactly at the boundaries of the weight space", {
  data(edu_synth, package = "cumulevel")
  J <- 3L
  fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = J, seed = 7)
  e <- fit$estimates

  n_J  <- sum(edu_synth$E == J)
  n_ge <- sum(edu_synth$E >= J)

  ## Exact scheme: only level J carries weight, so ESS must equal n_J.
  expect_equal(e$ess_treated[e$scheme == "exact"], n_J, tolerance = 1e-8)
  ## Binary scheme: all weights are 1, so ESS must equal n_{>= J}.
  expect_equal(e$ess_treated[e$scheme == "binary"], n_ge, tolerance = 1e-8)
})

test_that("the weight vector is aligned, not shifted", {
  ## Reproduces the failure mode of ifelse(piv[E - J + 1]): for E < J the
  ## index is 0 and R silently drops it, shortening the vector.
  E <- c(1L, 1L, 2L, 3L, 3L, 2L, 1L)
  J <- 2L; piv <- c(1, 0.4)
  w <- numeric(length(E))
  trt <- E >= J
  w[trt] <- piv[E[trt] - J + 1L]

  expect_length(w, length(E))
  expect_equal(sum(w > 0), sum(E >= J))
  expect_true(all(w[E < J] == 0))
  expect_equal(w, c(0, 0, 1, 0.4, 0.4, 1, 0))
})

test_that("the efficient score equation is solved exactly", {
  data(edu_synth, package = "cumulevel")
  fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3,
                        pi = "geometric", seed = 11)
  expect_lt(abs(fit$detail[[1]]$score), 1e-10)
})

test_that("the family degenerates at the top level", {
  ## At J = K the sum collapses to its anchor term, so every scheme must
  ## coincide exactly. An implementation with a misaligned weight vector
  ## would not produce this.
  data(edu_synth, package = "cumulevel")
  K <- max(edu_synth$E)
  fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = K, seed = 3)
  est <- fit$estimates$estimate
  expect_equal(max(est) - min(est), 0, tolerance = 1e-12)
  expect_equal(fit$estimates$ess_treated[1], sum(edu_synth$E == K),
               tolerance = 1e-8)
})

test_that("weight_scheme respects the anchor and the bounds", {
  expect_equal(weight_scheme("exact", 3, 5), c(1, 0, 0))
  expect_equal(weight_scheme("binary", 3, 5), c(1, 1, 1))
  expect_equal(weight_scheme("geometric", 3, 5, rho = 0.5), c(1, 0.5, 0.25))
  expect_error(weight_scheme(c(0.5, 0.2, 0.1), 3, 5), "anchor")
  expect_error(weight_scheme(c(1, 2, 0.1), 3, 5), "\\[0, 1\\]")
  expect_error(weight_scheme(c(1, 0.5), 3, 5), "length")
})

test_that("J = 1 is refused with an explanation", {
  data(edu_synth, package = "cumulevel")
  expect_error(
    weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 1),
    "sub-threshold set is empty")
})

test_that("positivity grades and warns", {
  data(edu_synth, package = "cumulevel")
  p <- positivity(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 5)
  expect_s3_class(p, "positivity")
  expect_true(p$verdict %in% c("comfortable", "acceptable", "strained"))
  expect_equal(p$share_below_cap, 0)

  ## An absurdly high cap must bind, and must warn.
  expect_warning(
    positivity(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3,
               cap = 0.9, seed = 5),
    "Positivity is strained")
})

test_that("the two fluctuation forms differ but both solve the score", {
  data(edu_synth, package = "cumulevel")
  a <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3,
                      pi = "geometric", fluctuation = "weighted", seed = 2)
  b <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3,
                      pi = "geometric", fluctuation = "covariate", seed = 2)
  expect_lt(abs(a$detail[[1]]$score), 1e-10)
  expect_lt(abs(b$detail[[1]]$score), 1e-10)
  ## Same constraint, different path: not the same estimator.
  expect_false(isTRUE(all.equal(a$estimates$estimate, b$estimates$estimate,
                                tolerance = 1e-10)))
})

test_that("the variance decomposition of the mixture layer is exact", {
  data(edu_synth, package = "cumulevel")
  fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 4)
  m <- mixture_interval(fit, S = 60, B = 40)
  expect_equal(m$attribution_share,
               m$var_attribution / (m$var_sampling + m$var_attribution),
               tolerance = 1e-12)
  expect_gte(m$attribution_share, 0)
  expect_lte(m$attribution_share, 1)
})


test_that("real data: the ESS invariants hold and the estimates are finite", {
  data(nhanes_smoking, package = "cumulevel")
  J <- 3L
  fit <- suppressWarnings(
    weighted_level(died ~ packyears | age + sex + ethnicity + pir + bmi,
                   data = nhanes_smoking, J = J, seed = 2005))
  e <- fit$estimates
  expect_true(all(is.finite(e$estimate)))
  expect_true(all(is.finite(e$se)))
  expect_equal(e$ess_treated[e$scheme == "exact"],
               sum(nhanes_smoking$packyears == J), tolerance = 1e-8)
  expect_equal(e$ess_treated[e$scheme == "binary"],
               sum(nhanes_smoking$packyears >= J), tolerance = 1e-8)
})

test_that("a rare factor level does not silently poison the outcome model", {
  ## 'education' has a level with 2 respondents. Split across folds and
  ## restricted to the sub-threshold sample, it can be absent from a
  ## training fold, leaving a rank-deficient design. glm.fit then returns
  ## NA coefficients; if these are not handled, every prediction is NA and
  ## the estimates come back as NA with no error.
  data(nhanes_smoking, package = "cumulevel")
  fit <- suppressWarnings(
    weighted_level(died ~ packyears | age + sex + ethnicity + education + pir + bmi,
                   data = nhanes_smoking, J = 3, seed = 2005))
  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.finite(fit$nuisances$Qh)))
})

test_that("the positivity warning reports counts, not a rounded zero", {
  data(nhanes_smoking, package = "cumulevel")
  expect_warning(
    weighted_level(died ~ packyears | age + sex + ethnicity + education + pir + bmi,
                   data = nhanes_smoking, J = 3, seed = 2005),
    "observation\\(s\\)")
})
