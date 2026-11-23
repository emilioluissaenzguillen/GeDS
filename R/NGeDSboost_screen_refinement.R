################################################################################
# NGeDSboost: shared entry screen for spline base-learners
################################################################################
.boost_screen_basis <- function(cache) {
  centered <- cache$xs - mean(cache$xs)
  norm <- sqrt(sum(centered^2))
  sites <- rle(cache$xs)
  list(ord = cache$ord, x = cache$xs, extr = cache$extr,
       linear = if (norm > 0) centered / norm else centered,
       sites = sites$values, counts = sites$lengths)
}

.boost_screen_first_knot <- function(basis, residual, scale, beta) {
  if (any(basis$counts > 1L)) {
    residual <- makeNewRes(residual, as.numeric(basis$counts))
  }
  residual[abs(residual) <= 32 * .Machine$double.eps * scale] <- 0
  if (all(residual == 0)) return(NA_real_)
  runs <- rle(sign(residual))$lengths
  ends <- cumsum(runs)
  starts <- c(1L, ends[-length(ends)] + 1L)
  cumulative <- c(0, cumsum(residual))
  means <- abs((cumulative[ends + 1L] - cumulative[starts]) / runs)
  if (max(means) == 0) return(NA_real_)
  spans <- basis$sites[ends] - basis$sites[starts]
  if (max(spans) > 0) spans <- spans / max(spans)
  knot <- Knotnew(beta * (means / max(means)) + (1 - beta) * spans,
                  residual, basis$sites, ends, rep(basis$extr, 2L), 1e-12, 2L)[1L]
  if (is.na(knot) || any(vapply(basis$extr, function(boundary) {
    isTRUE(all.equal(knot, boundary))
  }, logical(1)))) return(NA_real_)
  knot
}

# Project the candidate hinge off the line to get the one-knot RSS gain.
.boost_screen_gain <- function(basis, gradients, beta, return_line = FALSE) {
  y <- gradients[basis$ord, , drop = FALSE]
  scales <- apply(abs(y), 2L, max)
  level <- colMeans(y)
  y <- sweep(y, 2L, level)
  shift <- as.vector(crossprod(basis$linear, y))
  residuals <- y - tcrossprod(basis$linear, shift)
  rss <- colSums(residuals^2)
  gains <- numeric(ncol(y))
  for (i in which(rss > 1e-12 * colSums(y^2) & rss > 0)) {
    residual <- residuals[, i]
    knot <- .boost_screen_first_knot(basis, residual, scales[i], beta)
    if (is.na(knot)) next
    hinge <- pmax(basis$x - knot, 0)
    hinge <- hinge - mean(hinge)
    hinge <- hinge - basis$linear * sum(basis$linear * hinge)
    hinge_rss <- sum(hinge^2)
    if (hinge_rss > 0) {
      gains[i] <- min(1, sum(residual * hinge)^2 / (rss[i] * hinge_rss))
    }
  }
  if (!return_line) return(gains)
  list(gain = gains, theta = tcrossprod(basis$linear[c(1L, length(basis$x))], shift) +
         rep(level, each = 2L))
}

.boost_screen_line <- function(theta, extr, x) {
  # Keep the ordinary spline fit for degenerate or near-degenerate ranges.
  span <- diff(extr)
  if (span <= 1e-10) return(NULL)
  b1 <- diff(theta) / span
  b0 <- theta[1L] - b1 * extr[1L]
  list(Y_hat = b0 + b1 * x, knt = extr, int.knt = numeric(),
       b0 = b0, b1 = b1, theta = theta)
}

.boost_screen_calibrate <- function(bases, y, weights, family, beta) {
  restore_rng <- .geds_bootstrap_rng()
  on.exit(restore_rng(), add = TRUE)
  if (!setequal(y, c(-1, 1))) {
    stop("The shared screen requires both binomial classes.", call. = FALSE)
  }
  draws <- 199L
  probability <- family@response(family@offset(y, weights))
  if (length(probability) != 1L || !is.finite(probability) ||
      probability <= 0 || probability >= 1) {
    stop("The shared screen requires a finite intercept-only binomial fit.", call. = FALSE)
  }
  # Each bootstrap response is shared by all learners.
  gradients <- replicate(draws, {
    repeat {
      boot <- 2 * stats::rbinom(length(y), 1L, probability) - 1
      if (min(boot) != max(boot)) break
    }
    offset <- rep(family@offset(boot, weights), length.out = length(boot))
    family@ngradient(y = boot, f = offset, w = weights)
  })
  maximum <- numeric(draws)
  for (basis in bases) maximum <- pmax(maximum, .boost_screen_gain(basis, gradients, beta))
  sort(maximum, partial = 198L)[198L]
}

################################################################################
# NGeDSboost: joint coefficient refinement over fixed knots
################################################################################
.boost_svd_factor <- function(B, weights) {
  A <- B * sqrt(weights)
  s <- svd(A)
  cutoff <- max(dim(A)) * .Machine$double.eps * s$d[1L]
  keep <- s$d > cutoff
  list(u = s$u[, keep, drop = FALSE], v = s$v[, keep, drop = FALSE],
       d = s$d[keep], rank = sum(keep), cutoff = cutoff)
}

.boost_fixed_basis <- function(data, learners, knots, extr, levels = NULL, contrasts = NULL) {
  linear <- learners[vapply(learners, function(bl) bl$type == "linear", logical(1))]
  variables <- unique(unlist(lapply(linear, `[[`, "variables")))
  if (length(variables)) {
    if (!is.null(levels)) for (variable in names(levels)) {
      data[[variable]] <- factor(data[[variable]], levels = levels[[variable]])
    }
    Z <- stats::model.matrix(stats::reformulate(variables), data, contrasts.arg = contrasts)
    groups <- attr(Z, "assign")
  } else {
    Z <- matrix(1, nrow(data), 1L, dimnames = list(NULL, "(Intercept)"))
    groups <- 0L
  }
  spline <- jointSplineBasisMultivar(data, learners, knots, n = 2L,
                                    extrList = extr, prediction = TRUE)
  list(B = cbind(Z, spline$basis), parametric = variables, assign = groups,
       nparam = ncol(Z), matrices = spline$matrices,
       levels = lapply(data[variables][vapply(data[variables], is.factor, logical(1))], base::levels),
       contrasts = attr(Z, "contrasts"))
}

.boost_refinement_prediction <- function(object, model, newdata, base_learner = NULL) {
  ref <- model$refinement
  data <- newdata
  if (object$args$normalize_data) {
    variables <- intersect(names(object$args$X_mean), names(data))
    for (variable in variables) {
      data[[variable]] <- (data[[variable]] - object$args$X_mean[[variable]]) /
        object$args$X_sd[[variable]]
    }
  }
  design <- .boost_fixed_basis(data, object$args$base_learners[names(ref$anchor)],
                               ref$knots, object$args$extr, ref$levels, ref$contrasts)
  if (!identical(colnames(design$B), ref$columns)) {
    stop("New-data columns do not match the frozen boosting basis.", call. = FALSE)
  }
  delta <- ref$selected$delta
  if (!is.null(base_learner)) {
    spec <- object$args$base_learners[[base_learner]]
    if (spec$type == "linear") {
      j <- match(spec$variables, design$parametric)
      cols <- which(design$assign %in% j)
    } else {
      widths <- vapply(design$matrices, ncol, integer(1))
      j <- match(base_learner, names(widths))
      start <- design$nparam + sum(utils::head(widths, j - 1L))
      cols <- start + seq_len(widths[j])
    }
    if (is.character(ref$anchor_theta)) {
      intercept_learner <- if (length(design$parametric)) design$parametric[1L] else names(ref$anchor)[1L]
      if (base_learner == intercept_learner) cols <- c(1L, cols)
    }
    delta[-cols] <- 0
  }
  as.numeric(design$B %*% delta)
}

.boost_refine_coefficients <- function(B, y, F0, family, weights, shrinkage,
                                       phi, q, min_iterations, max_iterations) {
  factor <- .boost_svd_factor(B, weights)
  stats_family <- get_mboost_family(family@name)
  dev_y <- if (stats_family$family == "binomial") (y + 1) / 2 else y
  deviance <- function(F) {
    sum(stats_family$dev.resids(dev_y, family@response(F), weights))
  }
  F0 <- as.numeric(F0)
  delta <- numeric(ncol(B))
  coefficient_path <- list(delta)
  D <- deviance(F0)
  history <- list(list(step = 0L, deviance = D, F_hat = F0, delta = delta))
  trace <- data.frame(iteration = 0L, deviance = D, ratio = NA_real_)
  reason <- "iteration_limit"
  executed <- 0L
  F <- F0
  if (D >= 0 && D < 1e-6) {
    reason <- "near_zero_deviance"
  } else for (i in seq_len(max_iterations)) {
    gradient <- as.numeric(family@ngradient(y = y, f = F, w = weights))
    if (any(!is.finite(gradient))) {
      reason <- "nonfinite_gradient"
      break
    }
    step <- as.numeric(factor$v %*%
      (crossprod(factor$u, sqrt(weights) * gradient) / factor$d))
    next_delta <- delta + shrinkage * step
    next_F <- F0 + as.numeric(B %*% next_delta)
    next_D <- tryCatch(suppressWarnings(deviance(next_F)), error = function(e) NA_real_)
    if (!is.finite(next_D) || any(!is.finite(next_F))) {
      reason <- "nonfinite_update"
      break
    }
    executed <- i
    delta <- next_delta
    coefficient_path[[i + 1L]] <- delta
    F <- next_F
    D <- c(D, next_D)
    ratio <- NA_real_
    if (!is.null(q) && i >= q) {
      previous <- D[i + 1L - q]
      ratio <- if (previous == 0) Inf else next_D / previous
    }
    trace <- rbind(trace, data.frame(iteration = i, deviance = next_D, ratio = ratio))
    history[[length(history) + 1L]] <- list(step = i, deviance = next_D,
                                          F_hat = F, delta = delta)
    if (!is.null(q)) history <- utils::tail(history, q + 1L)
    if (next_D >= 0 && next_D < 1e-6) {
      reason <- "near_zero_deviance"
      break
    }
    if (!is.null(phi) && !is.null(q) && i >= q && i >= min_iterations && ratio >= phi) {
      reason <- "RD"
      break
    }
  }
  selected <- history[[which.min(vapply(history, `[[`, numeric(1), "deviance"))]]
  if (reason %in% c("nonfinite_gradient", "nonfinite_update")) {
    warning("Fixed-basis boosting stopped after a non-finite update; retaining the best finite model.",
            call. = FALSE)
  }
  list(selected = selected, iterations = executed, reason = reason, trace = trace,
       rank = factor$rank, cutoff = factor$cutoff, shrinkage = shrinkage,
       initial_F_hat = F0, coefficient_path = coefficient_path)
}

.boost_apply_refinement <- function(model, design, learners, delta, shrinkage) {
  # Allocate the explicit intercept to one component; no new learner is selected.
  intercept_learner <- if (length(design$parametric)) design$parametric[1L] else names(learners)[1L]
  for (variable in design$parametric) {
    j <- match(variable, design$parametric)
    increment <- delta[which(design$assign == j)]
    coefs <- model$base_learners[[variable]]$coefficients
    for (k in seq_along(increment)) coefs[[k + 1L]] <- coefs[[k + 1L]] + increment[k]
    if (variable == intercept_learner) coefs$b0 <- coefs$b0 + delta[1L]
    model$base_learners[[variable]]$coefficients <- coefs
  }
  start <- design$nparam
  for (name in names(design$matrices)) {
    width <- ncol(design$matrices[[name]])
    increment <- delta[start + seq_len(width)]
    start <- start + width
    if (name == intercept_learner) increment <- increment + delta[1L]
    learner <- model$base_learners[[name]]
    if (length(learners[[name]]$variables) == 1L) {
      k <- learner$knots
      slope <- diff(increment) / diff(k)
      slope[diff(k) == 0] <- 0
      learner$coefficients$b0 <- learner$coefficients$b0 + utils::head(increment, -1L) - slope * utils::head(k, -1L)
      learner$coefficients$b1 <- learner$coefficients$b1 + slope
    } else {
      internal <- learner$linear.int.knots
      if (length(learners[[name]]$variables) == 2L) {
        internal <- list(Xint.knt = internal$ikX, Yint.knt = internal$ikY)
      }
      # Tensor iteration contributions are multiplied by shrinkage at prediction.
      learner$iterations$phase2 <- list(int.knt = internal, coef = increment / shrinkage)
    }
    model$base_learners[[name]] <- learner
  }
  model
}
