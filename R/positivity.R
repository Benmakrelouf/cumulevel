#' Positivity diagnostic for a cumulative multilevel exposure
#'
#' Computes and grades the distribution of the estimated sub-threshold
#' propensity \eqn{\hat G_<(W) = \hat P(E < J \mid W)}, which is the single
#' quantity on which the identification of \eqn{\Psi_J(\pi)} rests.
#'
#' @param object Either a fitted \code{"cumulevel"} object, or a formula, in
#'   which case \code{data} and \code{J} must be supplied.
#' @param data A data frame, when \code{object} is a formula.
#' @param J Target level, when \code{object} is a formula.
#' @param cap Capping constant to compare against.
#' @param ... Passed to \code{\link{weighted_level}} when \code{object} is a
#'   formula.
#'
#' @return An object of class \code{"positivity"} with the quantiles of
#'   \eqn{\hat G_<}, the share of observations below \code{cap}, the control
#'   side effective sample size, and a verdict.
#'
#' @details
#' Only \emph{one} positivity condition is required by the framework, and it
#' concerns the \emph{sum} \eqn{G_<}, not the individual sub-threshold
#' probabilities: the reference regime must be feasible in every covariate
#' stratum. No positivity at all is required at the levels above \eqn{J},
#' because they enter weighted by their own occupation probabilities, so a
#' stratum in which some level is unattainable contributes nothing rather than
#' dividing by zero.
#'
#' The verdict is graded as follows.
#' \describe{
#'   \item{comfortable}{no observation below \code{cap}, and the first
#'     percentile at least four times \code{cap}. The capping constant never
#'     binds and the asymptotic theory applies as stated.}
#'   \item{acceptable}{no observation below \code{cap}, but the first
#'     percentile is close to it. Report the distribution alongside estimates.}
#'   \item{strained}{some observations below \code{cap}. Condition (i) of the
#'     asymptotic theorem fails: the capping term is then \emph{first order} in
#'     the outcome-model error rather than a product of two vanishing errors,
#'     so consistency survives only through a correctly specified outcome
#'     model. Vary \code{cap} and report the sensitivity.}
#' }
#'
#' A second quantity is reported that the treated-side effective sample size
#' cannot show: the control-side effective sample size computed from the
#' balancing weights. Under positivity strain it is the control side that
#' collapses first, and schemes with high \eqn{\pi} are affected most.
#'
#' @examples
#' data(edu_synth)
#' positivity(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
#'
#' @export
positivity <- function(object, data = NULL, J = NULL, cap = 0.025, ...) {
  if (inherits(object, "cumulevel")) return(object$positivity)
  if (inherits(object, "formula")) {
    if (is.null(data) || is.null(J))
      stop("'data' and 'J' are required when 'object' is a formula.",
           call. = FALSE)
    fit <- weighted_level(object, data, J = J, pi = "exact", cap = cap, ...)
    return(fit$positivity)
  }
  stop("'object' must be a fitted model or a formula.", call. = FALSE)
}

.positivity_from <- function(G_lt, rhat, cap, learner) {
  q <- stats::quantile(G_lt, c(0, 0.01, 0.05, 0.5, 0.95, 1), names = FALSE)
  share <- mean(G_lt < cap)
  ess_c <- { w <- rhat[rhat > 0]; if (length(w)) sum(w)^2 / sum(w^2) else 0 }

  verdict <- if (share > 0) "strained" else
    if (q[2] >= 4 * cap) "comfortable" else "acceptable"

  out <- list(min = q[1], p01 = q[2], p05 = q[3], median = q[4],
              p95 = q[5], max = q[6],
              share_below_cap = share, cap = cap,
              ess_control = ess_c, n = length(G_lt),
              learner = learner, verdict = verdict, G_lt = G_lt)
  class(out) <- "positivity"

  if (share > 0) {
    nbelow <- sum(G_lt < cap)
    warning(sprintf(
      "Positivity is strained: %d observation(s) (%.2f%%) have Ghat_< below the cap (%.3f);
  the smallest is %.2e. Condition (i) of the asymptotic theorem fails, and the
  capping term becomes first order in the outcome-model error rather than a
  product of two vanishing errors. Vary 'cap' and report the sensitivity.",
      nbelow, 100 * share, cap, min(G_lt)), call. = FALSE)
  }

  out
}

#' @export
print.positivity <- function(x, ...) {
  cat("\nPositivity diagnostic  (Ghat_< = P(E < J | W),", x$learner, "fit)\n")
  cat(strrep("-", 62), "\n")
  cat(sprintf("  min %.3f | 1%% %.3f | 5%% %.3f | median %.3f | max %.3f\n",
              x$min, x$p01, x$p05, x$median, x$max))
  cat(sprintf("  below cap %.3f : %s\n", x$cap,
              if (x$share_below_cap == 0) "none"
              else sprintf("%d obs (%.2f%%), smallest %.2e",
                           round(x$share_below_cap * x$n), 100 * x$share_below_cap,
                           x$min)))
  cat(sprintf("  effective sample size, control side : %.0f  (of %d)\n",
              x$ess_control, x$n))
  cat("  verdict :", switch(x$verdict,
      comfortable = "COMFORTABLE - the cap never binds.",
      acceptable  = "ACCEPTABLE - the cap does not bind, but report the distribution.",
      strained    = "STRAINED - the cap binds; see ?positivity."), "\n\n")
  invisible(x)
}

#' @export
plot.positivity <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::hist(x$G_lt, breaks = 40, col = "grey85", border = "white",
                 main = "Sub-threshold propensity", xlab = expression(hat(G)["<"]),
                 ...)
  graphics::abline(v = x$cap, col = "firebrick", lwd = 2, lty = 2)
  graphics::mtext(sprintf("cap = %.3f", x$cap), side = 3, adj = 1,
                  cex = 0.8, col = "firebrick")
  invisible(x)
}

## ---------------------------------------------------------------------
##  Sensitivity layer
## ---------------------------------------------------------------------

#' Monotone prior on the attribution weights
#'
#' A one-parameter prior on the decay rate \eqn{\rho}, with
#' \eqn{\pi(J + k) = \rho^k}. Monotonicity is guaranteed for every draw by
#' construction, which independent priors on the components do not ensure.
#'
#' @param stance One of \code{"agnostic"} (Beta(1,1)), \code{"dilution"}
#'   (Beta(2,4), prior mean 1/3), or \code{"inheritance"} (Beta(4,2)).
#' @param a,b Beta hyperparameters, overriding \code{stance} if supplied.
#'
#' @return A list of class \code{"attribution_prior"}.
#' @examples
#' attribution_prior("dilution")
#' @export
attribution_prior <- function(stance = c("agnostic", "dilution", "inheritance"),
                              a = NULL, b = NULL) {
  if (!is.null(a) && !is.null(b)) {
    if (a <= 0 || b <= 0) stop("'a' and 'b' must be positive.", call. = FALSE)
    return(structure(list(a = a, b = b, stance = "custom"),
                     class = "attribution_prior"))
  }
  stance <- match.arg(stance)
  ab <- switch(stance, agnostic = c(1, 1), dilution = c(2, 4),
               inheritance = c(4, 2))
  structure(list(a = ab[1], b = ab[2], stance = stance),
            class = "attribution_prior")
}

#' Mixture interval over attribution assumptions
#'
#' Propagates attribution uncertainty jointly with sampling uncertainty, using
#' a two-layer Monte Carlo scheme with a common multiplier in the inner layer.
#'
#' @param object A fitted \code{"cumulevel"} object.
#' @param prior An \code{"attribution_prior"} object.
#' @param S Number of prior draws.
#' @param B Number of inner-layer replicates.
#' @param level Confidence level.
#'
#' @return An object of class \code{"mixture"}: the mixture point estimate, the
#'   interval, and the attribution share of the total variance.
#'
#' @details
#' The variance decomposes exactly into a sampling component and an attribution
#' component. The first shrinks at rate \eqn{1/n}; the second is fixed by the
#' prior. The \emph{attribution share} therefore converges to one: at large
#' sample sizes what remains is the irreducible cost of not knowing \eqn{\pi}.
#'
#' \strong{This interval is prior-predictive, not a confidence interval.} It is
#' calibrated for a member of the family drawn from the prior, and
#' \emph{conservative} for the fixed mixture estimand, increasingly so as the
#' sample grows. Do not report it as a confidence interval for the mixture
#' estimand.
#'
#' @examples
#' data(edu_synth)
#' fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
#' mixture_interval(fit, attribution_prior("agnostic"), S = 100, B = 50)
#'
#' @export
mixture_interval <- function(object, prior = attribution_prior("agnostic"),
                             S = 1000L, B = 200L, level = 0.95) {
  if (!inherits(object, "cumulevel"))
    stop("'object' must be a fitted cumulevel model.", call. = FALSE)
  if (!inherits(prior, "attribution_prior"))
    stop("'prior' must be built by attribution_prior().", call. = FALSE)

  J <- object$J; K <- object$K
  W <- object$data$W; E <- object$data$E; Y <- object$data$Y
  n <- length(Y)

  rho <- stats::rbeta(S, prior$a, prior$b)
  psi <- se <- lam <- numeric(S)
  PHI <- matrix(0, S, n)
  for (s in seq_len(S)) {
    m <- .target_plugin(W, E, Y, object$nuisances, rho[s]^(0:(K - J)),
                        J, K, cap = object$cap,
                        fluctuation = object$fluctuation)
    psi[s] <- m$psi; se[s] <- m$se; lam[s] <- m$lambda; PHI[s, ] <- m$phi
  }

  ## Common multiplier: all members come from the same sample, so their
  ## estimation errors are strongly dependent. Drawing independently would
  ## have the right variance but the wrong joint law.
  xi <- matrix(stats::rnorm(n * B), n, B)
  draws <- sweep(lam * (PHI %*% xi) / n, 1, psi, "+")

  a <- (1 - level) / 2
  ci <- unname(stats::quantile(draws, c(a, 1 - a)))
  share <- stats::var(psi) / (mean(se^2) + stats::var(psi))

  out <- list(estimate = mean(psi), lower = ci[1], upper = ci[2],
              attribution_share = share,
              var_sampling = mean(se^2), var_attribution = stats::var(psi),
              prior = prior, S = S, B = B, level = level)
  class(out) <- "mixture"
  out
}

#' @export
print.mixture <- function(x, ...) {
  cat("\nMixture interval over attribution assumptions\n")
  cat(strrep("-", 62), "\n")
  cat(sprintf("  prior              : Beta(%g, %g)  [%s]\n",
              x$prior$a, x$prior$b, x$prior$stance))
  cat(sprintf("  mixture estimate   : %+.4f\n", x$estimate))
  cat(sprintf("  %.0f%% interval       : [%+.4f, %+.4f]\n",
              100 * x$level, x$lower, x$upper))
  cat(sprintf("  attribution share  : %.1f%% of total variance\n",
              100 * x$attribution_share))
  cat("\n  Prior-predictive, not a confidence interval: calibrated for a\n")
  cat("  prior-drawn member of the family, conservative for the fixed\n")
  cat("  mixture estimand. See ?mixture_interval.\n\n")
  invisible(x)
}
