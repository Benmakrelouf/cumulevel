#' Contribution weight schemes
#'
#' Builds the weight vector \eqn{(\pi(J), \ldots, \pi(K))} for a named scheme.
#' The anchor \eqn{\pi(J) = 1} is imposed by every scheme: individuals who stop
#' exactly at the target level always contribute in full.
#'
#' @param scheme One of \code{"exact"}, \code{"geometric"}, \code{"binary"}, or
#'   a numeric vector of length \code{K - J + 1} supplied directly.
#' @param J Target level.
#' @param K Number of exposure levels.
#' @param rho Decay rate for the geometric scheme, in \eqn{[0, 1]}.
#'
#' @return A numeric vector of length \code{K - J + 1} whose first entry is 1.
#'
#' @details
#' \code{"exact"} gives \eqn{(1, 0, \ldots, 0)}: the restricted average
#' treatment effect on \eqn{\{E = J\}}. \code{"binary"} gives
#' \eqn{(1, 1, \ldots, 1)}: the binary average treatment effect on
#' \eqn{\{E \ge J\}}. These are the two boundaries of the weight space, and
#' both are causally valid; they answer different questions.
#' \code{"geometric"} interpolates with \eqn{\pi(J + k) = \rho^k}.
#'
#' @examples
#' weight_scheme("geometric", J = 3, K = 5, rho = 0.5)
#' weight_scheme("exact", J = 3, K = 5)
#' @export
weight_scheme <- function(scheme = c("exact", "geometric", "binary"),
                          J, K, rho = 0.5) {
  if (is.numeric(scheme)) {
    if (length(scheme) != K - J + 1L)
      stop("A numeric weight vector must have length K - J + 1 = ",
           K - J + 1L, ".", call. = FALSE)
    if (!isTRUE(all.equal(scheme[1], 1)))
      stop("The anchor pi(J) must equal 1.", call. = FALSE)
    if (any(scheme < 0) || any(scheme > 1))
      stop("Weights must lie in [0, 1].", call. = FALSE)
    return(as.numeric(scheme))
  }
  scheme <- match.arg(scheme)
  k <- 0:(K - J)
  switch(scheme,
         exact     = c(1, rep(0, K - J)),
         geometric = rho^k,
         binary    = rep(1, K - J + 1L))
}

#' Weighted causal estimands for a cumulative multilevel exposure
#'
#' Estimates the family of weighted causal estimands \eqn{\Psi_J(\pi)} for an
#' exposure that is cumulative and hierarchical, by a level-specific weighted
#' targeted minimum loss-based estimator with cross-fitted nuisances.
#'
#' @param formula A formula of the form \code{Y ~ E | W1 + W2 + ...}, where
#'   \code{Y} is a binary outcome, \code{E} the ordered cumulative exposure
#'   coded \code{1, ..., K}, and the terms after \code{|} the adjustment
#'   covariates.
#' @param data A data frame.
#' @param J Target level, an integer in \code{2:K}. See Details for why
#'   \code{J = 1} is excluded.
#' @param pi Either \code{"family"} (the default), which estimates the exact,
#'   geometric and binary members together, or a scheme name accepted by
#'   \code{\link{weight_scheme}}, or a numeric weight vector.
#' @param rho Decay rate used when \code{pi} includes the geometric scheme.
#' @param learner \code{"parametric"} (logistic outcome model and multinomial
#'   logit propensities) or \code{"forest"} (cross-fitted probability random
#'   forests, requires \pkg{ranger}).
#' @param V Number of cross-fitting folds.
#' @param cap Lower bound applied to the estimated sub-threshold propensity.
#' @param fluctuation \code{"weighted"} (default) or \code{"covariate"}. Both
#'   impose the same score equation and are valid, but they are not the same
#'   estimator in finite samples; the weighted form is numerically stabler.
#' @param seed Optional integer seed for the fold assignment.
#'
#' @return An object of class \code{"cumulevel"}: a list with components
#'   \code{estimates} (a data frame with one row per weight scheme),
#'   \code{positivity} (an object of class \code{"positivity"}),
#'   \code{nuisances}, and the call.
#'
#' @details
#' \strong{Why the default returns a family.} When an exposure is cumulative,
#' dichotomising it at \eqn{J} attributes to level \eqn{J} whatever happens at
#' levels above it. The restricted and binary estimands are not right and wrong;
#' they are the two boundaries of a continuum of attribution assumptions, and
#' the gap between them is the quantity a reader needs. Reporting a single
#' member hides that gap, so \code{weighted_level()} reports the family.
#'
#' \strong{Positivity is checked automatically.} The distribution of the
#' estimated sub-threshold propensity \eqn{\hat G_<(W) = \hat P(E < J \mid W)}
#' is computed and returned, and a warning is issued when any observation falls
#' below \code{cap}. This matters because the capping term is second order when
#' the outcome model is consistent, but becomes \emph{first order} in the
#' outcome-model error when the cap binds.
#'
#' \strong{Why \code{J = 1} is excluded.} The sub-threshold set would be empty,
#' the reference regression would be a \eqn{0/0} ratio in every stratum, and
#' the balancing weight would be undefined.
#'
#' @references
#' Benmakrelouf, M. R. and Tatachak, A. (2026). Weighted causal inference for
#' cumulative multilevel exposures: a targeted learning approach.
#'
#' @seealso \code{\link{positivity}} for the diagnostic alone,
#'   \code{\link{mixture_interval}} for the sensitivity layer.
#'
#' @examples
#' data(edu_synth)
#' fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth,
#'                       J = 3, seed = 1)
#' fit
#' summary(fit)
#'
#' @export
weighted_level <- function(formula, data, J,
                           pi = "family",
                           rho = 0.5,
                           learner = c("parametric", "forest"),
                           V = 5L,
                           cap = 0.025,
                           fluctuation = c("weighted", "covariate"),
                           seed = NULL) {

  learner <- match.arg(learner)
  fluctuation <- match.arg(fluctuation)
  cl <- match.call()

  prt <- .parse_formula(formula, data)
  Y <- prt$Y; E <- prt$E; W <- prt$W
  K <- max(E)

  if (!is.numeric(J) || length(J) != 1L || J != round(J))
    stop("'J' must be a single integer.", call. = FALSE)
  J <- as.integer(J)
  if (J < 2L)
    stop("'J' must be at least 2: with J = 1 the sub-threshold set is empty, ",
         "the reference regression is undefined, and so is the balancing ",
         "weight.", call. = FALSE)
  if (J > K)
    stop("'J' cannot exceed the number of exposure levels (K = ", K, ").",
         call. = FALSE)
  if (!all(Y %in% c(0, 1)))
    stop("The outcome must be binary, coded 0/1.", call. = FALSE)
  if (sum(E < J) < 10L)
    stop("Fewer than 10 observations below the target level. ",
         "The reference arm cannot be estimated.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)
  nuis <- .fit_nuisances(W, E, Y, J, K, V = V, learner = learner)
  if (!nuis$converged)
    warning("The multinomial propensity model did not converge cleanly. ",
            "Inspect the positivity diagnostics before interpreting.",
            call. = FALSE)

  ## ---- weight schemes to report ------------------------------------
  if (identical(pi, "family")) {
    schemes <- list(exact     = weight_scheme("exact", J, K),
                    geometric = weight_scheme("geometric", J, K, rho),
                    binary    = weight_scheme("binary", J, K))
  } else if (is.character(pi) && length(pi) == 1L) {
    schemes <- stats::setNames(list(weight_scheme(pi, J, K, rho)), pi)
  } else {
    schemes <- list(custom = weight_scheme(pi, J, K))
  }

  res <- lapply(schemes, function(p)
    .target_plugin(W, E, Y, nuis, p, J, K, cap = cap,
                   fluctuation = fluctuation))

  est <- data.frame(
    scheme    = names(schemes),
    estimate  = vapply(res, `[[`, 0, "psi"),
    se        = vapply(res, `[[`, 0, "se"),
    se_robust = vapply(res, `[[`, 0, "se_robust"),
    lower     = NA_real_, upper = NA_real_,
    ess_treated = vapply(res, `[[`, 0, "ess_treated"),
    ess_control = vapply(res, `[[`, 0, "ess_control"),
    lambda    = vapply(res, `[[`, 0, "lambda"),
    row.names = NULL, stringsAsFactors = FALSE)
  est$lower <- est$estimate - stats::qnorm(0.975) * est$se_robust
  est$upper <- est$estimate + stats::qnorm(0.975) * est$se_robust

  pos <- .positivity_from(res[[1]]$G_lt, res[[1]]$rhat, cap, learner)

  ## Attribution spread: the quantity a dichotomised analysis conceals.
  spread <- NA_real_
  if (all(c("exact", "binary") %in% est$scheme)) {
    e <- est$estimate[est$scheme == "exact"]
    b <- est$estimate[est$scheme == "binary"]
    if (isTRUE(abs(e) > .Machine$double.eps)) spread <- 100 * (b / e - 1)
  }

  out <- list(call = cl, estimates = est, positivity = pos,
              spread_pct = spread, J = J, K = K, n = length(Y),
              learner = nuis$learner, V = V, cap = cap, rho = rho,
              fluctuation = fluctuation,
              weights = schemes, detail = res, nuisances = nuis,
              data = list(W = W, E = E, Y = Y))
  class(out) <- "cumulevel"
  out
}

## Parse  Y ~ E | W1 + W2  -----------------------------------------------
.parse_formula <- function(formula, data) {
  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula of the form Y ~ E | W1 + W2.",
         call. = FALSE)
  rhs <- formula[[3L]]
  if (!(is.call(rhs) && identical(as.character(rhs[[1L]]), "|")))
    stop("'formula' must separate exposure and covariates with '|', ",
         "as in  Y ~ E | W1 + W2.", call. = FALSE)

  yname <- deparse(formula[[2L]])
  ename <- deparse(rhs[[2L]])
  wform <- stats::as.formula(paste("~", deparse(rhs[[3L]])))

  if (!yname %in% names(data)) stop("Outcome '", yname, "' not found in data.",
                                    call. = FALSE)
  if (!ename %in% names(data)) stop("Exposure '", ename, "' not found in data.",
                                    call. = FALSE)

  mf <- stats::model.frame(wform, data = data, na.action = stats::na.pass)
  W  <- stats::model.matrix(wform, mf)
  W  <- W[, colnames(W) != "(Intercept)", drop = FALSE]

  Y <- data[[yname]]; E <- data[[ename]]
  ok <- stats::complete.cases(Y, E, W)
  if (any(!ok)) {
    message(sum(!ok), " incomplete observation(s) dropped.")
    Y <- Y[ok]; E <- E[ok]; W <- W[ok, , drop = FALSE]
  }
  E <- as.integer(as.numeric(E))
  if (min(E) < 1L) stop("Exposure levels must be coded 1, ..., K.", call. = FALSE)
  list(Y = as.numeric(Y), E = E, W = W)
}
