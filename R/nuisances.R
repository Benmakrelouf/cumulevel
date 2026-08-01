## ---------------------------------------------------------------------
##  Internal: nuisance estimation with cross-fitting
## ---------------------------------------------------------------------

## Multinomial logistic regression, hand-coded (level 1 is the reference).
## Avoids a dependency on nnet; BFGS with analytic gradient.
.fit_mlogit <- function(X, E, K) {
  p <- ncol(X); n <- length(E)
  Yind <- matrix(0, n, K)
  Yind[cbind(seq_len(n), E)] <- 1
  nll <- function(par) {
    eta <- cbind(0, X %*% matrix(par, p, K - 1L))
    m <- apply(eta, 1, max)
    -sum(eta[cbind(seq_len(n), E)] - (m + log(rowSums(exp(eta - m)))))
  }
  gr <- function(par) {
    eta <- cbind(0, X %*% matrix(par, p, K - 1L))
    m <- apply(eta, 1, max)
    P <- exp(eta - m); P <- P / rowSums(P)
    as.vector(-crossprod(X, (Yind - P)[, -1, drop = FALSE]))
  }
  o <- stats::optim(rep(0, p * (K - 1L)), nll, gr, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-12))
  list(B = matrix(o$par, p, K - 1L), converged = identical(o$convergence, 0L))
}

.pred_mlogit <- function(fit, X) {
  eta <- cbind(0, X %*% fit$B)
  m <- apply(eta, 1, max)
  P <- exp(eta - m)
  P / rowSums(P)
}

## Cross-fitted nuisance estimation.
##   learner = "parametric" : logistic Q_ref + multinomial logit g
##   learner = "forest"     : ranger for both (falls back with a warning)
.fit_nuisances <- function(W, E, Y, J, K, V = 5L, learner = "parametric",
                           folds = NULL, num.trees = 300L) {
  n <- length(Y)
  X <- cbind(1, W)
  if (is.null(folds)) folds <- sample(rep(seq_len(V), length.out = n))

  use_forest <- identical(learner, "forest")
  if (use_forest && !requireNamespace("ranger", quietly = TRUE)) {
    warning("Package 'ranger' is not installed; falling back to ",
            "learner = \"parametric\".", call. = FALSE)
    use_forest <- FALSE
  }

  Qh <- rep(NA_real_, n)
  Gh <- matrix(NA_real_, n, K)
  conv <- TRUE
  aliased <- 0L

  for (v in unique(folds)) {
    te <- folds == v
    tr <- !te
    sub <- tr & (E < J)
    if (sum(sub) < 5L)
      stop("Fewer than 5 sub-threshold observations in a training fold. ",
           "Reduce 'V', or check that level J is not too low.", call. = FALSE)

    if (use_forest) {
      d1 <- data.frame(y = Y[sub], W[sub, , drop = FALSE])
      m1 <- ranger::ranger(y ~ ., data = d1, num.trees = num.trees)
      Qh[te] <- stats::predict(m1, data.frame(W[te, , drop = FALSE]))$predictions
      d2 <- data.frame(y = factor(E[tr], levels = seq_len(K)),
                       W[tr, , drop = FALSE])
      m2 <- ranger::ranger(y ~ ., data = d2, num.trees = num.trees,
                           probability = TRUE)
      Gh[te, ] <- stats::predict(m2, data.frame(W[te, , drop = FALSE]))$predictions
    } else {
      ## A factor level that is rare overall can be absent from the
      ## sub-threshold part of a training fold, leaving a constant column
      ## and a rank-deficient design. glm.fit then returns NA coefficients,
      ## which would silently propagate into every prediction. Aliased
      ## columns are dropped, as predict.glm does for aliased terms.
      f <- suppressWarnings(
        stats::glm.fit(X[sub, , drop = FALSE], Y[sub],
                       family = stats::binomial()))
      cf <- f$coefficients
      if (anyNA(cf)) {
        n_alias <- sum(is.na(cf))
        cf[is.na(cf)] <- 0
        aliased <<- aliased + n_alias
      }
      Qh[te] <- stats::plogis(drop(X[te, , drop = FALSE] %*% cf))
      mf <- .fit_mlogit(X[tr, , drop = FALSE], E[tr], K)
      conv <- conv && mf$converged
      Gh[te, ] <- .pred_mlogit(mf, X[te, , drop = FALSE])
    }
  }

  Gh <- pmax(Gh, 1e-8)
  Gh <- Gh / rowSums(Gh)

  if (aliased > 0L)
    warning(sprintf(
      "%d aliased coefficient(s) dropped from the outcome model.\n  A covariate level is too rare to be present in every training fold;\n  consider collapsing it or reducing 'V'.", aliased), call. = FALSE)
  if (!all(is.finite(Qh)))
    stop("The outcome model produced non-finite predictions. ",
         "Check for rare factor levels or collinear covariates.", call. = FALSE)

  list(Qh = Qh, Gh = Gh, folds = folds, converged = conv,
       aliased = aliased,
       learner = if (use_forest) "forest" else "parametric")
}

## Targeting step and plug-in, for one weight vector.
## This is the part that the sensitivity layer re-runs; the nuisances above
## do not depend on pi and are fitted once.
.target_plugin <- function(W, E, Y, nuis, piv, J, K, cap = 0.025,
                           fluctuation = c("weighted", "covariate")) {
  fluctuation <- match.arg(fluctuation)
  n <- length(Y)
  Qh <- nuis$Qh; Gh <- nuis$Gh

  G_lt <- rowSums(Gh[, seq_len(J - 1L), drop = FALSE])
  hW   <- drop(Gh[, J:K, drop = FALSE] %*% piv)
  rhat <- hW / pmax(G_lt, cap)
  cap_share <- mean(G_lt < cap)

  trt <- E >= J
  sub <- !trt
  ## NB: never build this with ifelse(): for E < J the index piv[E - J + 1]
  ## is zero and R silently drops it, shortening and misaligning the vector.
  wpi <- numeric(n)
  wpi[trt] <- piv[E[trt] - J + 1L]
  stopifnot(length(wpi) == n,
            sum(wpi > 0) == sum(E >= J & piv[pmax(E - J + 1L, 1L)] > 0),
            all(wpi[sub] == 0))
  Dn <- mean(wpi)

  off <- stats::qlogis(pmin(pmax(Qh, 1e-6), 1 - 1e-6))
  score <- function(eps) {
    Qs <- if (fluctuation == "covariate")
      stats::plogis(off[sub] + eps * rhat[sub])
    else
      stats::plogis(off[sub] + eps)
    sum(rhat[sub] * (Y[sub] - Qs))
  }
  eps <- tryCatch(
    stats::uniroot(score, c(-50, 50), tol = .Machine$double.eps^0.75)$root,
    error = function(e) 0)
  Qstar <- if (fluctuation == "covariate")
    stats::plogis(off + eps * rhat) else stats::plogis(off + eps)

  psi_init <- sum(wpi * (Y - Qh)) / sum(wpi)
  psi      <- sum(wpi * (Y - Qstar)) / sum(wpi)
  psi_os   <- psi_init - sum(rhat[sub] * (Y[sub] - Qh[sub])) / (n * Dn)
  Rn       <- psi - psi_os

  alpha <- wpi - rhat * sub
  phi   <- (alpha * (Y - Qstar) - psi * wpi) / Dn
  se    <- sqrt(mean(phi^2) / n)
  lambda <- sqrt(1 + n * Rn^2 / mean(phi^2))

  ess <- function(w) { w <- w[w > 0]; if (!length(w)) return(0); sum(w)^2 / sum(w^2) }

  list(psi = psi, se = se, se_robust = se * lambda,
       lambda = lambda, Rn = Rn, eps = eps,
       ess_treated = ess(wpi), ess_control = ess(rhat[sub]),
       cap_share = cap_share, score = mean(phi), phi = phi,
       G_lt = G_lt, rhat = rhat)
}
