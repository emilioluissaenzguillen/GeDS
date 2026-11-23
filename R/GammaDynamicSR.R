# Gamma-dynamic SR-pre uses the native Stage A path and a bootstrap
# under the fitted model with zero internal knots.
.geds_bootstrap_rng <- function() {
  previous <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_kind <- RNGkind()
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  function() {
    if (!identical(RNGkind(), previous_kind)) {
      do.call(RNGkind, as.list(previous_kind))
    }
    if (is.null(previous)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", previous, envir = .GlobalEnv)
    }
  }
}

.gamma_dynamic_rule <- function(stoptype, omitted, legacy, X, Y, Z, offset,
                                weights, family, initial = FALSE) {
  if (!identical(stoptype, "GD")) return(stoptype)
  valid_weights <- all(weights == 1)
  if (family %in% c("gaussian", "Gamma")) {
    valid_weights <- all(is.finite(weights) & weights >= 0)
  }
  if (identical(family, "binomial")) {
    valid_weights <- is.numeric(Y) && is.null(dim(Y)) &&
      all(is.finite(weights) & weights > 0 & weights == round(weights)) &&
      all(is.finite(Y) & Y >= 0 & Y <= 1) &&
      all(abs(Y * weights - round(Y * weights)) < 1e-10)
  }
  supported <- ncol(X) == 1L && is.null(Z) && !initial &&
    all(is.finite(offset)) && valid_weights && is.numeric(Y) &&
    length(Y) >= 4L && family %in% c("gaussian", "poisson", "Gamma", "binomial", "inverse.gaussian")
  if (supported) return(stoptype)
  if (omitted) return(legacy)
  stop("'GD' requires a univariate Gaussian or Gamma fit (optionally precision-weighted), an unweighted univariate Poisson or inverse Gaussian fit, or a univariate binomial fit with valid trial counts, with finite offsets and without parametric terms or initial knots.",
       call. = FALSE)
}

.gd_invalid_path_steps <- function(loss) {
  which(!is.finite(loss) | loss <= 0 | c(FALSE, diff(loss) >= 0))
}

.gd_draw_inverse_gaussian <- function(mu, dispersion) {
  # Michael-Schucany-Haas: use the reciprocal form of the smaller root.
  t <- mu * dispersion * stats::rchisq(length(mu), 1) / 2
  ratio <- 1 + t + sqrt(t * (t + 2))
  ifelse(stats::runif(length(mu)) <= 1 / (1 + 1 / ratio), mu / ratio, mu * ratio)
}

# The 95% type-7 quantile of stats::quantile(), computed with the same sort and
# interpolation but without its argument checks and name formatting. Missing
# values are left to stats::quantile(), which reports them as errors.
.gd_quantile95 <- function(v) {
  n <- length(v)
  if (!n || anyNA(v)) return(unname(stats::quantile(v, .95)))
  index <- 1 + (n - 1) * .95
  lo <- floor(index)
  hi <- ceiling(index)
  v <- sort(v, partial = unique(c(lo, hi)))
  q <- v[lo]
  if (index > lo && v[hi] != q) {
    h <- index - lo
    q <- (1 - h) * q + h * v[hi]
  }
  q
}

# Row-wise .gd_quantile95(), i.e. apply(W, 1L, .gd_quantile95), with one sort.
.gd_quantile95_rows <- function(W) {
  n <- ncol(W)
  if (!n || anyNA(W)) return(apply(W, 1L, function(v) unname(stats::quantile(v, .95))))
  index <- 1 + (n - 1) * .95
  lo <- floor(index)
  hi <- ceiling(index)
  # Column r of `sorted` holds row r of W in increasing order.
  sorted <- matrix(W[order(row(W), W)], nrow = n)
  q <- sorted[lo, ]
  upper <- sorted[hi, ]
  move <- index > lo & upper != q
  h <- index - lo
  q[move] <- (1 - h) * q[move] + h * upper[move]
  q
}

.gamma_dynamic_gains <- function(D, q) {
  if (length(D) - 3L < q) return(rep(NA_real_, length(D)))
  k <- seq.int(q, length(D) - 1L)
  gain <- -expm1(log(D[k + 1L] / D[k - q + 1L]) / q)
  invalid <- which(!is.finite(gain) | gain <= 0)
  if (length(invalid)) gain <- utils::head(gain, invalid[1L] - 1L)
  smoothed <- rep(NA_real_, length(D))
  if (length(gain) >= 3L) {
    y <- log(gain)
    t <- seq_along(y)
    # Endpoint of the linear regression of log gains for each prefix.
    total <- cumsum(y)
    endpoint <- exp(total / t +
      (cumsum(t * y) - (t + 1) * total / 2) / (t * (t + 1) / 6))
    smoothed[k[t] + 1L] <- endpoint
  }
  smoothed
}

.bivariate_gamma_dynamic_fit <- function(args) {
  restore_rng <- .geds_bootstrap_rng()
  on.exit(restore_rng(), add = TRUE)
  q <- args$q
  if (is.na(q) || q < 1L) stop("'q' must be a positive integer for Gamma dynamic.", call. = FALSE)
  initial <- length(args$Xintknots) + length(args$Yintknots)
  upper <- min(args$max.intknots, initial + 20L)
  if (!is.finite(upper) || upper < initial || !is.finite(args$min.intknots) ||
      args$min.intknots < 0L) {
    stop("The knot budget must accommodate the initial knots.", call. = FALSE)
  }
  lower <- max(initial, min(args$min.intknots, upper))
  prefix <- function(d) {
    bad <- which(!is.finite(d) | d <= 0 | c(FALSE, diff(d) >= 0))
    if (length(bad)) utils::head(d, bad[1L] - 1L) else d
  }
  path <- function(y, budget) {
    a <- args
    a$Z <- y
    a$max.intknots <- budget
    a$stoptype <- ".gamma_dynamic_path"
    suppressWarnings(do.call(BivariateFitter, a))
  }
  y <- args$Z
  D <- if (upper == initial) 0 else path(y, upper)
  if (!length(D) || !is.finite(D[1L]) || D[1L] < 0) {
    stop("Could not fit the initial bivariate spline basis.", call. = FALSE)
  }
  D <- if (D[1L] == 0) 0 else prefix(D)
  K <- length(D) - 1L
  lower <- min(lower - initial, K)
  selected <- K
  reason <- "path_end"
  positive <- which(args$weights > 0)
  w <- args$weights[positive]
  yp <- y[positive]
  total <- sum(w * (yp - stats::weighted.mean(yp, w))^2)
  if (upper == initial || K <= lower || D[1L] <= .Machine$double.eps * total) {
    selected <- lower
    reason <- if (upper == initial) "fixed_knots" else "no_improvement"
  } else {
    x <- cbind(args$X, args$Y)[positive, , drop = FALSE]
    x <- sweep(x, 2L, c(diff(args$Xextr), diff(args$Yextr)), "/")
    distances <- as.matrix(stats::dist(x))^2
    diag(distances) <- Inf
    p <- min(10L, nrow(x) - 1L)
    nn <- t(vapply(seq_len(nrow(x)), function(i) order(distances[i, ])[seq_len(p)], integer(p)))
    dx <- colMeans(matrix(distances[cbind(rep(seq_len(nrow(x)), p), as.vector(nn))], ncol = p))
    null <- SplineReg_biv(args$X, args$Y, y, weights = args$weights,
      InterKnotsX = args$Xintknots, InterKnotsY = args$Yintknots,
      n = 2L, Xextr = args$Xextr, Yextr = args$Yextr, fast = TRUE)
    mu <- null$predicted
    weights <- args$weights[positive]
    w <- sqrt(weights)
    varying_weights <- length(unique(weights)) > 1L
    if (varying_weights) b <- tensorProd(null$Xbasis, null$Ybasis)
    estimate <- function(response) {
      r <- response[positive]
      if (varying_weights) {
        residuals <- stats::lm.wfit(b, response, args$weights)$residuals
        r <- residuals[positive]
      }
      r <- w * r
      dy <- colMeans(matrix((r - r[nn])^2, ncol = p)) / 2
      v <- stats::.lm.fit(cbind(1, dx), dy)$coefficients[1L]
      if (is.finite(v) && v > 0) return(v)
      max(.Machine$double.eps * max(weights), stats::mad(r - r[nn[, 1L]])^2 / 2)
    }
    scale <- estimate(y)
    if (!is.finite(scale) || scale <= 0) {
      stop("Could not estimate the bivariate Gamma-dynamic noise scale.", call. = FALSE)
    }
    entry <- function(d, variance) {
      if (length(d) <= 1L) return(0)
      max((d[1L] - d[-1L]) / seq_len(length(d) - 1L)) / variance
    }
    statistic <- function(d, variance) {
      k <- seq_along(d) - 1L
      .gamma_dynamic_gains(d, q) * d[pmax(1L, k - q + 1L)] / variance
    }
    boot <- vector("list", 100L)
    accepted <- 0L
    # A shorter path contributes its available gains and zero gain thereafter.
    for (i in seq_len(3000L)) {
      response <- mu
      response[positive] <- mu[positive] + stats::rnorm(length(positive), sd = sqrt(scale) / w)
      d <- tryCatch(prefix(path(response, upper)), error = function(e) numeric())
      if (!length(d)) next
      variance <- estimate(response)
      if (!is.finite(variance) || variance <= 0) next
      accepted <- accepted + 1L
      boot[[accepted]] <- list(loss = d, scale = variance)
      if (accepted == length(boot)) break
    }
    cutoff <- function(x) unname(stats::quantile(x, .95))
    if (accepted < length(boot)) {
      warning("Bivariate Gamma-dynamic calibration failed. Using SR.", call. = FALSE)
      args$stoptype <- "SR"
      args$max.intknots <- upper
      fit <- do.call(BivariateFitter, args)
      fit$stopinfo <- list(rule = "SR", reason = "calibration_failed", bootstrap_paths = accepted)
      return(fit)
    }
    threshold <- cutoff(vapply(boot, function(b) entry(b$loss, b$scale), numeric(1)))
    if (entry(D, scale) <= threshold) {
      selected <- lower
      reason <- "entry_screen"
    } else {
      statistics <- vapply(boot, function(b) {
        z <- statistic(b$loss, b$scale)
        c(z, rep(0, max(0L, length(D) - length(z))))[seq_along(D)]
      }, numeric(length(D)))
      critical <- apply(statistics, 1L, function(z) {
        if (all(is.na(z))) NA_real_ else cutoff(z[!is.na(z)])
      })
      observed <- statistic(D, scale)
      hit <- is.finite(observed) & is.finite(critical) & observed <= critical
      confirmed <- which(utils::head(hit, -1L) & utils::tail(hit, -1L))
      if (length(confirmed)) {
        selected <- max(lower, confirmed[1L] - q)
        reason <- "dynamic_trigger"
      }
    }
  }
  args$max.intknots <- initial + selected
  args$stoptype <- ".gamma_dynamic_refit"
  fit <- do.call(BivariateFitter, args)
  fit$args$max.intknots <- upper
  fit$stopinfo <- list(rule = "GD", reason = reason, budget = upper - initial,
                       selected_additions = selected)
  fit
}

.gamma_dynamic_neighbours <- function(x) {
  # For sorted x, the p nearest neighbours of i lie within i - p, ..., i + p.
  # Rank every candidate in one pass by (observation, distance, index).
  n <- length(x)
  p <- min(10L, n - 1L)
  offsets <- c(-rev(seq_len(p)), seq_len(p))
  i <- rep(seq_len(n), each = length(offsets))
  j <- i + offsets
  valid <- j >= 1L & j <= n
  i <- i[valid]
  j <- j[valid]
  o <- order(i, abs(x[j] - x[i]), j)
  i <- i[o]
  j <- j[o]
  nn <- matrix(j[sequence(tabulate(i, n)) <= p], nrow = p)
  dx <- colMeans((x - matrix(x[nn], nrow = length(x), byrow = TRUE))^2)
  list(nn = nn, dx = dx)
}

.gamma_dynamic_select <- function(x, y, family, q, lower, upper, path,
                                  design = cbind(1, x), weights = rep(1, length(y)),
                                  allow_short_paths = FALSE, incumbent = Inf,
                                  neighbour_cache = NULL, draws = 100L,
                                  staged = FALSE, bootstrap_path = function(y, k, b) path(y, k),
                                  prefix_observed = FALSE, offset = rep(0, length(y))) {
  if (prefix_observed) {
    # Cut the observed path at its first non-improving step, as the bootstrap
    # paths are cut, so that T and W are computed on matching prefixes.
    full_path <- path
    path <- function(y, k) {
      loss <- full_path(y, k)
      bad <- .gd_invalid_path_steps(loss)
      if (length(bad) && bad[1L] > 1L) utils::head(loss, bad[1L] - 1L) else loss
    }
  }
  depth <- if (staged) min(upper, max(lower + 1L, 40L)) else upper
  D <- path(y, depth)
  K <- length(D) - 1L
  complete <- depth == upper || K < depth
  if (K < 0L || any(!is.finite(D)) || any(D < 0)) return(NULL)
  lower <- min(lower, K)
  gaussian <- family$family == "gaussian"
  identity_link <- gaussian && family$link == "identity"
  has_offset <- any(offset != 0)
  weighted <- any(weights != 1)
  positive <- weights > 0
  varying_weights <- weighted && length(unique(weights[positive])) > 1L
  numerically_exact <- FALSE
  if (gaussian && any(positive)) {
    adjusted_y <- y[positive]
    if (identity_link && has_offset) adjusted_y <- adjusted_y - offset[positive]
    w <- weights[positive]
    total <- sum(w * (adjusted_y - stats::weighted.mean(adjusted_y, w))^2)
    numerically_exact <- D[1L] <= .Machine$double.eps * total
  }
  if (K <= lower || D[1L] == 0 || numerically_exact) {
    return(list(knots = lower, reason = "no_improvement"))
  }

  poisson <- family$family == "poisson"
  binomial <- family$family == "binomial"
  unit_scale <- poisson || binomial
  inverse_gaussian <- family$family == "inverse.gaussian"
  if (sum(positive) < 4L) return(NULL)
  # A path that never reaches the incumbent's loss cannot be selected, so skip
  # the null fit, noise-scale estimate and bootstrap setup.
  if (is.finite(incumbent) &&
      min(D) > incumbent + 1e-8 * max(1, D[1L], incumbent)) {
    return(list(knots = lower, reason = "dominated"))
  }
  restore_rng <- .geds_bootstrap_rng()
  on.exit(restore_rng(), add = TRUE)
  null_fit <- function(response, bootstrap = FALSE) {
    mustart <- if (bootstrap && inverse_gaussian) rep(mean(response), length(response)) else NULL
    fit <- function(start) {
      suppressWarnings(stats::glm.fit(design, response, weights = weights, family = family,
                                      offset = offset, mustart = start))
    }
    if (inverse_gaussian) .retry_ig_zero_start(fit, mustart, response, weights) else fit(mustart)
  }
  null <- null_fit(y)
  mu <- null$fitted.values
  x <- x[positive]
  if (!unit_scale) {
    neighbours <- neighbour_cache$geometry
    if (is.null(neighbours)) {
      neighbours <- .gamma_dynamic_neighbours(x)
      if (!is.null(neighbour_cache)) neighbour_cache$geometry <- neighbours
    }
    nn <- neighbours$nn
    dx <- neighbours$dx
  }
  estimate <- function(response, fit) {
    if (unit_scale) return(1)
    if (gaussian) {
      # Varying precision weights make sqrt(w) * mean non-smooth, so remove the
      # zero-knot fit first. Constant positive weights keep the unweighted
      # residuals, so rescaling the weights or adding zero-weight rows leaves
      # the calibration unchanged.
      r <- if (varying_weights) response - fit$fitted.values else if (identity_link && has_offset) {
        response - offset
      } else if (has_offset) {
        response - fit$fitted.values
      } else response
      if (weighted) r <- r * sqrt(weights)
    } else {
      r <- (response - fit$fitted.values) / sqrt(family$variance(fit$fitted.values))
      if (weighted && family$family == "Gamma") r <- r * sqrt(weights)
    }
    r <- r[positive]
    # colMeans((r - matrix(r[nn], nrow = length(x), byrow = TRUE))^2) / 2
    dy <- gd_neighbour_semivariance(r, nn)
    v <- stats::.lm.fit(cbind(1, dx), dy)$coefficients[1L]
    if (is.finite(v) && v > 0) return(v)
    v <- stats::mad(diff(r), constant = 1.482602218505602)^2 / 2
    if (is.finite(v) && v > 0) return(v)
    if (gaussian) mean((r - mean(r))^2) else sum(r^2) / fit$df.residual
  }
  scale <- estimate(y, null)
  if (!is.finite(scale) || scale <= 0) return(NULL)
  draw <- function() {
    if (gaussian) {
      response <- mu
      response[positive] <- stats::rnorm(sum(positive), mu[positive], sqrt(scale / weights[positive]))
      return(response)
    }
    if (poisson) return(stats::rpois(length(y), mu))
    if (binomial) return(stats::rbinom(length(y), size = weights, prob = mu) / weights)
    if (inverse_gaussian) return(.gd_draw_inverse_gaussian(mu, scale))
    if (weighted) {
      response <- mu
      dispersion <- scale / weights[positive]
      response[positive] <- stats::rgamma(sum(positive), shape = 1 / dispersion,
                                          scale = mu[positive] * dispersion)
      return(response)
    }
    stats::rgamma(length(y), shape = 1 / scale, scale = mu * scale)
  }

  entry <- function(loss, variance) {
    k <- seq_len(length(loss) - 1L)
    if (!length(k)) return(0)
    max((loss[1L] - loss[-1L]) / k) / variance
  }
  observed <- entry(D, scale)
  responses <- vector("list", 150L)
  target <- min(K, 40L)
  boot <- vector("list", length(responses))
  repeat {
    accepted <- integer()
    for (b in seq_along(responses)) {
      if (is.null(boot[[b]]) || (!boot[[b]]$ended && length(boot[[b]]$D) <= target)) {
        # Generate responses in their original order and reuse them for deeper paths.
        if (is.null(responses[[b]])) responses[[b]] <- draw()
        previous_scale <- boot[[b]]$scale
        boot[b] <- list(tryCatch({
          response <- responses[[b]]
          loss <- bootstrap_path(response, target, b)
          bad <- .gd_invalid_path_steps(loss)
          if (allow_short_paths && length(bad)) loss <- utils::head(loss, bad[1L] - 1L)
          if (!length(loss) || (!allow_short_paths &&
              (length(loss) != target + 1L || length(bad)))) NULL else {
            v <- if (is.null(previous_scale)) {
              estimate(response, null_fit(response, bootstrap = TRUE))
            } else previous_scale
            if (is.finite(v) && v > 0) {
              list(D = loss, scale = v, ended = length(loss) <= target)
            }
          }
        }, error = function(e) NULL))
      }
      if (!is.null(boot[[b]])) accepted <- c(accepted, b)
      if (length(accepted) == draws) break
    }
    if (length(accepted) < draws) return(NULL)

    paths <- boot[accepted]
    lower_entry <- vapply(paths, function(b) entry(b$D, b$scale), numeric(1))
    # Beyond the fitted prefix, the reduction cannot exceed the initial loss.
    upper_entry <- if (complete && target == K) lower_entry else pmax(lower_entry,
      vapply(paths, function(b) {
        if (b$ended) 0 else b$D[1L] / ((target + 1L) * b$scale)
      }, numeric(1)))
    cutoff <- .gd_quantile95
    entry_low <- cutoff(lower_entry)
    entry_high <- cutoff(upper_entry)
    while (!complete && observed <= entry_high && D[1L] / (length(D) * scale) > entry_low) {
      depth <- min(upper, 2L * depth)
      D <- path(y, depth)
      K <- length(D) - 1L
      if (K < target || any(!is.finite(D)) || any(D < 0)) return(NULL)
      complete <- depth == upper || K < depth
      observed <- entry(D, scale)
    }
    if (complete && target == K) entry_high <- entry_low
    if (observed <= entry_low) {
      return(list(knots = lower, reason = "entry_screen"))
    }
    if (observed > entry_high && target - q >= 3L) {
      observed_gain <- .gamma_dynamic_gains(D, q)
      k <- seq.int(q + 2L, target)
      W <- vapply(paths, function(b) {
        # An exhausted bootstrap path contributes zero gain beyond its last fit.
        gain <- numeric(length(k))
        available <- k < length(b$D)
        j <- k[available]
        gain[available] <- b$D[j - q + 1L] * .gamma_dynamic_gains(b$D, q)[j + 1L] / b$scale
        gain
      }, numeric(length(k)))
      critical <- .gd_quantile95_rows(W)
      hit <- is.finite(observed_gain[k + 1L]) &
        D[k - q + 1L] * observed_gain[k + 1L] / scale <= critical
      confirmed <- which(utils::head(hit, -1L) & tail(hit, -1L))
      if (length(confirmed)) {
        return(list(knots = max(lower, k[confirmed[1L] + 1L] - q),
                    reason = "dynamic_trigger"))
      }
    }
    if (complete && target == K) return(list(knots = K, reason = "path_end"))
    while (!complete && K < 2L * target) {
      depth <- min(upper, 2L * depth)
      D <- path(y, depth)
      K <- length(D) - 1L
      if (K < target || any(!is.finite(D)) || any(D < 0)) return(NULL)
      complete <- depth == upper || K < depth
      observed <- entry(D, scale)
    }
    target <- min(K, 2L * target)
  }
}

.gamma_dynamic_fit <- function(fitter, args, family, legacy, allow_short_paths = FALSE,
                               incumbent = Inf, neighbour_cache = NULL,
                               draws = 100L, staged = FALSE, cache_initial_qr = FALSE) {
  inverse_gaussian <- family$family == "inverse.gaussian"
  # Binomial IRLS paths can lose fit at larger knot counts. Use the improving
  # prefix of both the observed and bootstrap paths instead of rejecting them.
  binomial <- family$family == "binomial"
  allow_short_paths <- allow_short_paths || inverse_gaussian || binomial
  if (!is.numeric(draws) || length(draws) != 1L || !is.finite(draws) ||
      draws < 10L || draws > 100L || draws != floor(draws)) {
    stop("'gd_draws' must be an integer between 10 and 100.", call. = FALSE)
  }
  if (!is.finite(args$q) || args$q < 1L) {
    stop("'q' must be a positive integer for Gamma dynamic.", call. = FALSE)
  }
  upper <- args$max.intknots
  inherited <- length(args$intknots_init)
  lower <- max(inherited, min(args$min.intknots, upper))
  if (!is.finite(lower) || !is.finite(upper) || lower < 0L) {
    stop("Internal-knot limits must be finite nonnegative integers.", call. = FALSE)
  }
  if (upper < inherited) {
    stop("Initial knots exceed 'max.intknots'; supply a larger limit.", call. = FALSE)
  }
  if (upper == lower) {
    choice <- list(knots = upper, reason = "fixed_knots")
  } else {
    initial_qr <- NULL
    if (cache_initial_qr && identical(fitter, UnivariateFitter) &&
        is.null(args$fit_init) && is.null(args$Z) && all(args$weights == 1) &&
        length(args$X) < 5000L && inherited < 48L) {
      # Use the same dense basis and LINPACK QR as the initial linear fit.
      initial_basis <- splines::splineDesign(
        sort(c(args$intknots_init, rep(args$extr, 2L))), args$X,
        ord = 2L, outer.ok = TRUE)
      initial_qr <- qr(initial_basis)
      if (initial_qr$rank < ncol(initial_basis)) initial_qr <- NULL
    }
    path <- function(y, cap, state = NULL) {
      path_args <- args
      path_args$Y <- y
      path_args$max.intknots <- inherited + cap
      path_args$stoptype <- if (is.null(state)) ".gamma_dynamic_path" else state
      path_args$show.iters <- FALSE
      if (!is.null(initial_qr) && is.null(state$restart)) {
        theta <- qr.coef(initial_qr, y - args$offset)
        path_args$fit_init <- list(intknots = args$intknots_init,
          pred = as.numeric(initial_basis %*% theta) + args$offset)
      }
      loss <- suppressWarnings(do.call(fitter, path_args))
      if (inherited) loss[-seq_len(inherited)] else loss
    }
    bootstrap_path <- function(y, cap, b) path(y, cap)
    if (staged || inverse_gaussian) {
      states <- new.env(parent = emptyenv())
      bootstrap_path <- function(y, cap, b) {
        key <- as.character(b)
        state <- states[[key]]
        if (is.null(state)) {
          state <- structure(new.env(parent = emptyenv()), class = "geds_gd_path")
          if (inverse_gaussian) state$inverse_gaussian_bootstrap <- TRUE
          states[[key]] <- state
        }
        path(y, cap, state)
      }
    }
    design <- if (inherited) {
      splines::splineDesign(sort(c(args$intknots_init, rep(args$extr, 2L))), args$X, ord = 2L)
    } else cbind(1, args$X)
    choice <- tryCatch(
      .gamma_dynamic_select(args$X, args$Y, family, args$q, lower - inherited,
                            upper - inherited, path, design, args$weights, allow_short_paths,
                            incumbent, neighbour_cache, draws, staged, bootstrap_path,
                            prefix_observed = binomial, offset = args$offset),
      error = identity
    )
    if (!is.null(choice) && !inherits(choice, "error")) choice$knots <- inherited + choice$knots
  }
  if (is.null(choice) || inherits(choice, "error")) {
    detail <- if (inherits(choice, "error")) paste0(" (", conditionMessage(choice), ")") else ""
    warning("Gamma-dynamic calibration failed", detail, ". Using ", legacy, ".", call. = FALSE)
    args$stoptype <- legacy
    out <- do.call(fitter, args)
    out$stopinfo$rule <- legacy
    return(out)
  }

  if (identical(choice$reason, "dominated")) {
    return(list(stopinfo = list(rule = "GD", reason = "dominated")))
  }
  if (choice$knots < lower) {
    # The valid Stage A path ended before the requested minimum.
    warning(warningCondition(
      paste0("The valid Stage A path ended with ", choice$knots,
             " internal knots, below 'min.intknots' = ", lower,
             "; the last valid fit is used."),
      class = "geds_min_intknots_warning", call = NULL))
  }
  args$max.intknots <- choice$knots
  args$stoptype <- ".gamma_dynamic_refit"
  out <- withCallingHandlers(do.call(fitter, args), warning = function(w) {
    if (inherits(w, "geds_knot_cap_warning") && choice$knots < upper) {
      invokeRestart("muffleWarning")
    }
  })
  out$args$max.intknots <- upper
  out$stopinfo$rule <- "GD"
  out$stopinfo$reason <- choice$reason
  out
}


