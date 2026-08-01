#' @export
print.cumulevel <- function(x, digits = 4, ...) {
  cat("\nWeighted causal estimands for a cumulative exposure\n")
  cat(strrep("-", 66), "\n")
  cat(sprintf("  target level J = %d of K = %d   |   n = %d   |   %s nuisances, V = %d\n",
              x$J, x$K, x$n, x$learner, x$V))
  cat("\n")

  e <- x$estimates
  tab <- data.frame(
    scheme   = e$scheme,
    estimate = round(e$estimate, digits),
    se       = round(e$se_robust, digits),
    CI       = sprintf("[%+.*f, %+.*f]", digits, e$lower, digits, e$upper),
    ESS      = round(e$ess_treated),
    stringsAsFactors = FALSE)
  print(tab, row.names = FALSE)

  if (!is.na(x$spread_pct)) {
    cat(sprintf("\n  Attribution spread, exact to binary: %+.1f%%\n", x$spread_pct))
    cat("  This is what a dichotomised analysis would attribute to level ",
        x$J, "\n  but that originates above it.\n", sep = "")
  }

  v <- x$positivity$verdict
  cat(sprintf("\n  Positivity: %s (min Ghat_< = %.3f, cap = %.3f)\n",
              toupper(v), x$positivity$min, x$cap))
  if (identical(v, "strained"))
    cat("  >> The cap binds. See ?positivity before interpreting.\n")
  cat("\n")
  invisible(x)
}

#' @export
summary.cumulevel <- function(object, ...) {
  structure(list(fit = object), class = "summary.cumulevel")
}

#' @export
print.summary.cumulevel <- function(x, digits = 4, ...) {
  o <- x$fit
  print.cumulevel(o, digits = digits)

  cat("Diagnostics by scheme\n")
  cat(strrep("-", 66), "\n")
  e <- o$estimates
  d <- data.frame(
    scheme      = e$scheme,
    ESS_treated = round(e$ess_treated),
    ESS_control = round(e$ess_control),
    lambda      = round(e$lambda, 3),
    se_EIF      = round(e$se, digits),
    se_robust   = round(e$se_robust, digits),
    stringsAsFactors = FALSE)
  print(d, row.names = FALSE)

  cat("\n  ESS_treated counts only the attribution weights and is blind to the\n")
  cat("  balancing weights; under positivity strain it is ESS_control that\n")
  cat("  collapses first. Report both.\n")

  cat("\n  lambda is a flag, not a variance correction. It tracks only the part\n")
  cat("  of the second-order remainder that the targeting step removed, and is\n")
  cat("  structurally blind to first-order failures. Read it alongside the\n")
  cat("  positivity diagnostic, never instead of it.\n")

  print(o$positivity)
  invisible(x)
}

#' @export
coef.cumulevel <- function(object, ...) {
  stats::setNames(object$estimates$estimate, object$estimates$scheme)
}

#' @export
confint.cumulevel <- function(object, parm, level = 0.95, ...) {
  e <- object$estimates
  z <- stats::qnorm(1 - (1 - level) / 2)
  m <- cbind(e$estimate - z * e$se_robust, e$estimate + z * e$se_robust)
  dimnames(m) <- list(e$scheme,
                      sprintf("%.1f %%", 100 * c((1 - level) / 2,
                                                 1 - (1 - level) / 2)))
  if (!missing(parm)) m <- m[parm, , drop = FALSE]
  m
}

#' Plot the estimand family
#'
#' Displays the members of the family side by side with their intervals. The
#' plot is deliberately the default output of \code{plot()} on a fitted model:
#' the gap between the exact-level and binary members is the quantity a
#' dichotomised analysis conceals.
#'
#' @param x A fitted \code{"cumulevel"} object.
#' @param main Title.
#' @param xlab Axis label.
#' @param ... Further graphical parameters.
#' @return \code{x}, invisibly.
#' @examples
#' data(edu_synth)
#' fit <- weighted_level(Y ~ E | W1 + W2 + W3, data = edu_synth, J = 3, seed = 1)
#' plot(fit)
#' @export
plot.cumulevel <- function(x, main = NULL, xlab = NULL, ...) {
  e <- x$estimates
  m <- nrow(e)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(4.5, 8, 3, 2))

  rng <- range(c(e$lower, e$upper, 0))
  pad <- diff(rng) * 0.08
  if (is.null(main)) main <- sprintf("Estimand family, target level J = %d", x$J)
  if (is.null(xlab)) xlab <- expression(hat(Psi)[J](pi))

  plot(NA, xlim = c(rng[1] - pad, rng[2] + pad), ylim = c(0.5, m + 0.5),
       yaxt = "n", ylab = "", xlab = xlab, main = main, ...)
  graphics::abline(v = 0, col = "grey70", lty = 3)
  cols <- c(exact = "black", geometric = "steelblue4",
            binary = "firebrick", custom = "darkgreen")

  for (i in seq_len(m)) {
    yy <- m - i + 1
    cc <- if (e$scheme[i] %in% names(cols)) cols[[e$scheme[i]]] else "grey30"
    graphics::segments(e$lower[i], yy, e$upper[i], yy, col = cc, lwd = 2)
    graphics::segments(e$lower[i], yy - 0.1, e$lower[i], yy + 0.1, col = cc)
    graphics::segments(e$upper[i], yy - 0.1, e$upper[i], yy + 0.1, col = cc)
    graphics::points(e$estimate[i], yy, pch = 19, col = cc, cex = 1.3)
    graphics::mtext(e$scheme[i], side = 2, at = yy, las = 1, line = 0.6,
                    cex = 0.9)
  }
  if (!is.na(x$spread_pct))
    graphics::mtext(sprintf("attribution spread: %+.1f%%", x$spread_pct),
                    side = 3, adj = 1, cex = 0.8, col = "grey30")
  invisible(x)
}
