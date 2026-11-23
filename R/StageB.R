################################################################################
# Stage B: REML-penalised refits (NGeDSboost, NGeDSgam)
################################################################################
.stageB_validate_penalty <- function(penalty, normalize_data, higher_order) {
  penalty <- match.arg(penalty, c("none", "component"))
  if (penalty != "none") {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("Install 'mgcv' to use Stage B penalisation.", call. = FALSE)
    if (isTRUE(normalize_data))
      stop("Stage B penalisation currently requires normalize_data = FALSE.", call. = FALSE)
    if (!isTRUE(higher_order))
      stop("Stage B penalisation requires higher_order = TRUE.", call. = FALSE)
  }
  penalty
}

# Four-point Gauss quadrature integrates products of cubic spline derivatives
# exactly on each nonempty knot interval. Coordinates are scaled to [0, 1].
.stageB_penalty_axis <- function(x, internal, n) {
  boundary <- range(x)
  span <- diff(boundary)
  if (!is.finite(span) || span <= 0)
    stop("Penalised spline predictors must have a nonzero range.", call. = FALSE)
  knots <- sort(c(internal, rep(boundary, n)))
  breaks <- unique(knots)
  lo <- utils::head(breaks, -1L)
  hi <- utils::tail(breaks, -1L)
  nodes <- c(-0.861136311594053, -0.339981043584856,
             0.339981043584856, 0.861136311594053)
  weights <- c(0.347854845137454, 0.652145154862546,
               0.652145154862546, 0.347854845137454)
  points <- as.vector(outer((hi - lo) / 2, nodes) + (hi + lo) / 2)
  qw <- as.vector(outer((hi - lo) / (2 * span), weights))
  roots <- lapply(0:2, function(d) {
    splines::splineDesign(knots, points, ord = n,
                         derivs = rep(d, length(points)), outer.ok = TRUE) *
      sqrt(qw) * span^d
  })
  size <- length(knots) - n
  greville <- vapply(seq_len(size), function(i) {
    mean(knots[i + seq_len(n - 1L)])
  }, numeric(1))
  list(basis = splines::splineDesign(knots, x, ord = n, outer.ok = TRUE),
    roots = roots, affine = (greville - boundary[1L]) / span)
}

# Firth replaces the curvature-penalised refit only after binomial REML fails.
# Work on the estimable design space, then map back to the original spline
# coefficients. This also defines a minimum-norm convention for aliased columns.
.stageB_firth_refit <- function(design, y, weights, B, P) {
  if (!requireNamespace("brglm2", quietly = TRUE))
    stop("Install 'brglm2' for the binomial Stage B Firth fallback.", call. = FALSE)
  decomposition <- svd(design * sqrt(weights))
  rank <- sum(decomposition$d > max(decomposition$d) * 1e-10)
  if (!rank) stop("The binomial Stage B design has no estimable directions.", call. = FALSE)
  map <- sweep(decomposition$v[, seq_len(rank), drop = FALSE], 2L,
               decomposition$d[seq_len(rank)], "/")
  x <- design %*% map
  attempts <- list()
  for (slow in c(1, 0.5, 0.1)) {
    notices <- character()
    fit <- tryCatch(withCallingHandlers(
      brglm2::brglmFit(x = x, y = y, weights = weights,
        family = stats::binomial(), intercept = TRUE, start = rep(0, rank),
        control = brglm2::brglmControl(type = "AS_mean", maxit = 1000,
                                      epsilon = 1e-10, slowit = slow)),
      warning = function(w) {
        notices <<- c(notices, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    attempts[[as.character(slow)]] <- list(
      converged = isTRUE(fit$converged), warnings = notices,
      error = if (inherits(fit, "error")) conditionMessage(fit) else NULL)
    if (!inherits(fit, "error") && isTRUE(fit$converged) &&
        all(is.finite(fit$coefficients))) {
      eta <- as.numeric(x %*% fit$coefficients)
      prob <- stats::plogis(eta)
      working_weights <- weights * prob * (1 - prob)
      q <- qr(x * sqrt(working_weights), tol = 1e-10)
      if (q$rank != rank) next
      leverage <- rowSums(qr.Q(q)^2)
      score <- crossprod(x, weights * (y - prob) + leverage * (0.5 - prob))
      if (any(!is.finite(score)) || max(abs(score)) > 1e-6) next
      theta <- setNames(as.numeric(P %*% map %*% fit$coefficients), colnames(B))
      eta <- as.numeric(B %*% theta)
      if (any(!is.finite(eta))) next
      return(list(theta = theta, eta = eta, rank = rank,
        columns = ncol(design), warnings = notices, attempts = attempts,
        score = max(abs(score))))
    }
  }
  stop("Binomial Stage B Firth fallback did not converge reliably.", call. = FALSE)
}

.stageB_penalised_fit <- function(X, Y, Z, offset, weights, learners, internal,
                                  n, family, link, penalty) {
  if (!n %in% c(3L, 4L)) stop("Penalisation is limited to orders 3 and 4.")
  if (any(offset != 0))
    stop("Stage B penalisation currently does not support formula offsets.", call. = FALSE)
  if (inherits(family, "boost_family")) {
    family <- get_mboost_family(family@name)
    if (!is.null(link)) family <- get(family$family, envir = asNamespace("stats"))(link)
    if (family$family == "binomial") Y <- (Y + 1) / 2
  }
  supported <- paste(family$family, family$link) %in%
    c("gaussian identity", "Gamma log", "poisson log", "binomial logit")
  if (!supported || !is.numeric(Y) || !is.null(dim(Y)))
    stop("Penalisation supports Gaussian identity, Gamma log, Poisson log and binary logit responses only.",
         call. = FALSE)
  if (family$family == "binomial" && any(!Y %in% c(0, 1)))
    stop("Penalised binomial fits require individual binary responses.", call. = FALSE)
  keep <- weights > 0
  w <- weights / sum(weights[keep]) * sum(keep)
  if (any(!is.finite(w)) || any(w < 0) || !any(keep))
    stop("Invalid weights for penalised Stage B fitting.", call. = FALSE)

  # Match the native coefficient order: univariate splines, tensors, linear terms.
  dims <- lengths(lapply(learners, `[[`, "variables"))
  if (any(dims > 2L))
    stop("Stage B penalisation currently supports univariate and bivariate learners only.", call. = FALSE)
  learners <- learners[order(dims)]
  blocks <- roots <- nulls <- list()
  for (name in names(learners)) {
    vars <- learners[[name]]$variables
    ik <- internal[[name]]
    if (length(vars) == 1L) ik <- list(ik)
    axes <- lapply(seq_along(vars), function(j) .stageB_penalty_axis(X[[vars[j]]], ik[[j]], n))
    if (length(vars) == 1L) {
      B <- axes[[1L]]$basis
      L <- axes[[1L]]$roots[[3L]]
      N <- cbind(1, axes[[1L]]$affine)
    } else {
      B <- tensorProd(axes[[1L]]$basis, axes[[2L]]$basis)
      a <- lapply(axes[[1L]]$roots, function(root) qr.R(qr(root)))
      b <- lapply(axes[[2L]]$roots, function(root) qr.R(qr(root)))
      L <- rbind(kronecker(a[[3L]], b[[1L]]),
                 sqrt(2) * kronecker(a[[2L]], b[[2L]]),
                 kronecker(a[[1L]], b[[3L]]))
      N <- cbind(1, rep(axes[[1L]]$affine, each = ncol(axes[[2L]]$basis)),
                 rep(axes[[2L]]$affine, times = ncol(axes[[1L]]$basis)))
    }
    colnames(B) <- paste0(name, "_", seq_len(ncol(B)))
    blocks[[name]] <- B
    roots[[name]] <- L * sqrt(sum(keep))
    nulls[[name]] <- N
  }
  if (!length(blocks))
    stop("Stage B penalisation requires at least one selected spline learner.", call. = FALSE)
  linear <- if (NCOL(Z)) stats::model.matrix(~ ., data = Z)[, -1L, drop = FALSE] else
    matrix(numeric(), nrow(X), 0L)
  B <- cbind(do.call(cbind, blocks), linear)
  sizes <- vapply(blocks, ncol, integer(1))
  p <- ncol(B)
  null <- matrix(0, p, 1L)
  null[seq_len(sizes[1L]), 1L] <- 1
  penalised <- list()
  pos <- 0L
  for (j in seq_along(blocks)) {
    ix <- pos + seq_len(sizes[j])
    pos <- pos + sizes[j]
    mu <- colSums(blocks[[j]] * w) / sum(w)
    H <- diag(sizes[j])[, -sizes[j], drop = FALSE] -
      matrix(mu[-sizes[j]], sizes[j], sizes[j] - 1L, byrow = TRUE)
    # Analytic affine null spaces avoid squaring the conditioning in L'L.
    A <- nulls[[j]][-sizes[j], -1L, drop = FALSE] -
      matrix(nulls[[j]][sizes[j], -1L], sizes[j] - 1L,
             ncol(nulls[[j]]) - 1L, byrow = TRUE)
    Q <- qr.Q(qr(A), complete = TRUE)
    nullity <- ncol(A)
    Nj <- matrix(0, p, nullity)
    Nj[ix, ] <- H %*% Q[, seq_len(nullity), drop = FALSE]
    null <- cbind(null, Nj)
    J <- Q[, -seq_len(nullity), drop = FALSE]
    decomposition <- svd(roots[[j]] %*% H %*% J, nu = 0, nv = ncol(J))
    if (any(!is.finite(decomposition$d)) || any(decomposition$d <= 0))
      stop("The curvature penalty could not be factorised reliably.", call. = FALSE)
    Pj <- matrix(0, p, ncol(J))
    Pj[ix, ] <- H %*% J %*% sweep(decomposition$v, 2L, decomposition$d, "/")
    penalised[[j]] <- Pj
  }
  if (ncol(linear)) null <- cbind(null, diag(p)[, pos + seq_len(ncol(linear)), drop = FALSE])
  # Remove redundant affine terms and unobserved isotropic penalised directions.
  design_null <- B[keep, , drop = FALSE] %*% null * sqrt(w[keep])
  qnull <- qr(design_null, tol = 1e-10)
  null <- null[, qnull$pivot[seq_len(qnull$rank)], drop = FALSE]
  reduce <- function(P) {
    D <- B[keep, , drop = FALSE] %*% P * sqrt(w[keep])
    v <- svd(D, nu = 0, nv = min(dim(D)))
    rank <- sum(v$d > max(v$d) * 1e-12)
    if (!rank) return(matrix(0, p, 0L))
    P %*% v$v[, seq_len(rank), drop = FALSE]
  }
  names(penalised) <- names(blocks)
  penalised <- lapply(penalised, reduce)
  P <- do.call(cbind, c(list(null), penalised))
  sizes_pen <- vapply(penalised, ncol, integer(1))
  S <- list()
  pos <- ncol(null)
  for (size in sizes_pen) {
    if (size) {
      Sj <- matrix(0, ncol(P), ncol(P))
      ix <- pos + seq_len(size)
      Sj[ix, ix] <- diag(size)
      S[[length(S) + 1L]] <- Sj
    }
    pos <- pos + size
  }
  if (!length(S)) stop("No estimable penalised curvature directions remain.", call. = FALSE)
  design <- B[keep, , drop = FALSE] %*% P
  y <- Y[keep]
  wt <- w[keep]
  notices <- character()
  fitted <- tryCatch(withCallingHandlers(
    mgcv::gam(y ~ design - 1, family = family, weights = wt, method = "REML",
              paraPen = list(design = c(S, list(rank = sizes_pen[sizes_pen > 0L]))),
              control = mgcv::gam.control(maxit = 100)),
    warning = function(e) {
      notices <<- c(notices, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  ), error = function(e) {
    numerical <- grepl(paste("non-finite|NA/NaN/Inf|singular|positive definite|",
      "inner loop|step size|cannot correct|infinite or missing|converg", sep = ""),
      conditionMessage(e), ignore.case = TRUE)
    if (family$family != "binomial" || !numerical) stop(e)
    e
  })
  optimisation_warnings <- notices[!grepl("non-integer.*success", notices)]
  reason <- if (inherits(fitted, "error")) conditionMessage(fitted) else NULL
  if (is.null(reason)) {
    outer <- fitted$outer.info$conv
    if (!isTRUE(fitted$converged) || length(optimisation_warnings) ||
        (!is.null(outer) && outer != "full convergence"))
      reason <- "Penalised Stage B REML did not converge reliably; inspect the data and knot configuration."
  }
  if (is.null(reason)) {
    theta <- setNames(as.numeric(P %*% stats::coef(fitted)), colnames(B))
    eta <- as.numeric(B %*% theta)
    if (any(!is.finite(eta)) || (family$family == "binomial" &&
        any(abs(eta[keep]) > -log(sqrt(.Machine$double.eps)))))
      reason <- "Penalised Stage B produced non-finite or saturated training predictions."
    else if (any(!is.finite(family$linkinv(eta))))
      reason <- "Penalised Stage B produced non-finite fitted means."
  }
  if (!is.null(reason)) {
    if (family$family != "binomial") {
      for (notice in optimisation_warnings) warning(notice, call. = FALSE)
      stop(reason, call. = FALSE)
    }
    fallback <- .stageB_firth_refit(design, y, wt, B, P)
    theta <- fallback$theta
    eta <- fallback$eta
    diagnostic <- list(structure = "none", method = "Firth", sp = numeric(),
      requested = penalty, fallback_from = "component REML", reason = reason,
      edf = NA_real_, scale = 1, rank = fallback$rank, columns = fallback$columns,
      warnings = fallback$warnings, reml_warnings = notices, converged = TRUE,
      attempts = fallback$attempts, adjusted_score = fallback$score,
      alias_convention = "Minimum norm in the REML design coordinates",
      weight_policy = "Relative weights normalised over positive-weight rows")
    warning(warningCondition(paste0("The order-", n,
      " binomial Stage B REML refit failed: ", reason,
      " Firth logistic regression is used on the same spline basis, without the curvature penalty."),
      class = "geds_stageB_firth_warning", order = n, reason = reason))
    for (notice in fallback$warnings[!grepl("non-integer.*success", fallback$warnings)])
      warning(notice, call. = FALSE)
  } else {
    sp <- setNames(fitted$sp, names(penalised)[sizes_pen > 0L])
    diagnostic <- list(structure = penalty, method = "REML", sp = sp,
      edf = sum(fitted$edf), scale = fitted$sig2, warnings = notices, converged = TRUE,
      weight_policy = "Relative weights normalised over positive-weight rows")
  }
  predicted <- family$linkinv(eta)
  residuals <- Y - predicted
  list(theta = theta, predicted = predicted, residuals = residuals,
       linear.predictors = eta,
       rss = .weighted_rss(residuals, weights),
       deviance = sum(family$dev.resids(Y, predicted, weights)),
       basis = B[, seq_len(sum(sizes)), drop = FALSE], temporary = NULL,
       nci = NULL, aci = NULL, polygon = NULL,
       penalty = diagnostic)
}

################################################################################
# Stage B: stability safeguard for unpenalised refits (NGeDSboost)
################################################################################
#' Stability of unpenalised Stage B refits
#'
#' Unpenalised quadratic and cubic Stage B fits from \code{NGeDSboost}
#' automatically simplify their averaged knots
#' when the joint spline and parametric design is excessively sensitive to the
#' training responses. Linear fits, Stage A and REML-penalised fits are unchanged.
#' \code{NGeDSgam} and standalone \code{NGeDS} and \code{GGeDS} fits do not use
#' this safeguard.
#'
#' @details
#' With the Stage A knots and prior observation weights held fixed, let
#' \eqn{h(z)} map the training responses to a weighted least-squares prediction
#' at \eqn{z}, and set \eqn{s(z) = \|h(z)\|^2}. The reference \eqn{s_0(z)} is
#' computed from the same-order, knot-free joint model, retaining parametric
#' effects. The safeguard seeks \eqn{s(z) \leq s_0(z) + 10} on training-derived
#' probe locations, so the knots may add at most 10 to the knot-free
#' sensitivity. The margin 10 is a fixed heuristic, not a
#' cross-validated tuning parameter or a guarantee of predictive accuracy.
#' With joint training design \eqn{B}, prior-weight matrix \eqn{W}, and basis
#' row \eqn{b(z)^T}, the response map is
#' \eqn{h(z) = W B (B^T W B)^+ b(z)} for identifiable predictions, where
#' \eqn{+} denotes the Moore-Penrose inverse. Known offsets do not affect this
#' response sensitivity. The reference uses the analogous map with all internal
#' knots removed, not the Stage A fit. It requires no comparison of fitted
#' responses or residual errors. The margin is additive, not
#' \eqn{10 s_0(z)}: where the knot-free model is already sensitive, for example
#' a knot-free joint learner in a small sample, knots that add little
#' sensitivity are not removed on that account.
#'
#' The probe set contains the training rows and 33 equally spaced locations per
#' interval between successive observed values or original averaged knots of
#' each spline coordinate. Other numeric coordinates are interpolated between
#' neighbouring training rows; factors retain observed levels. This set is fixed
#' before deletion and uses no test responses or test predictor locations.
#' Zero-weight rows do not enter the least-squares criterion, but their predictor
#' locations remain in the probe set and can therefore affect simplification.
#' It is a finite-grid check, not a bound over the entire predictor domain.
#' Sufficient coefficient-based bounds can certify supported designs without
#' constructing the probe grid. Inconclusive bounds leave the probe calculation
#' and knot-deletion algorithm unchanged.
#'
#' If necessary, one averaged Stage B knot at a time is deleted, choosing the
#' deletion with the smallest worst sensitivity ratio. Non-identifiable
#' prediction functionals are addressed first. Knots are not moved, and Stage A
#' knot histories are retained. Final coefficients are estimated by the usual
#' unpenalised least-squares or iteratively reweighted least-squares fitter;
#' the singular-value decomposition used by the diagnostic is not a penalty.
#' Existing coefficient limits and omission of unavailable orders still apply.
#'
#' For Gaussian identity-link fits, this sensitivity has the conditional
#' homoscedastic prediction-variance interpretation. For other families it is
#' a design-stability diagnostic based on prior weights, not on the final GLM
#' working weights, and is not a bound on response-scale prediction variance.
#' If the knot-free reference is not identifiable on the probes, the original
#' fit is retained with a \code{geds_stageB_stability_warning} warning.
#'
#' @section Diagnostics:
#' Each available unpenalised higher-order boosted fit stores a
#' \code{stability} list under \code{final_model$quadratic.fit} or
#' \code{final_model$cubic.fit}.
#' The list records the threshold, original knots, deletion path, number removed,
#' \code{target_met} and \code{max_ratio}. If \code{certified} is true,
#' \code{max_ratio} is a sufficient upper bound divided by 10, rather than the
#' sampled maximum ratio; \code{reference_max} is then unavailable.
#'
#' @name StageB_stability
#' @seealso \code{\link{NGeDSgam}}, \code{\link{NGeDSboost}}
NULL

# Geometry only: final coefficients are still estimated by the native (IR)LS
# fitter. The diagnostic uses relative prior weights, not GLM working weights.
.stageB_stability_context <- function(X, Z, weights, learners, internal, n,
                                       boundaries, expand_Z = TRUE) {
  X <- as.data.frame(X)
  spline_names <- names(learners)[vapply(learners, function(b) b$type == "GeDS", logical(1))]
  internal <- internal[intersect(names(internal), spline_names)]
  dims <- vapply(learners[names(internal)], function(b) length(b$variables), integer(1))
  internal <- internal[order(dims > 1L, dims > 2L)]
  for (bl in names(internal)) {
    vars <- learners[[bl]]$variables
    if (length(vars) > 1L) names(internal[[bl]]) <- vars
  }
  parametric <- character()
  if (NCOL(Z)) {
    Z <- as.data.frame(Z)
    parametric <- make.unique(c(names(X), paste0(".geds_param", seq_len(ncol(Z)))))
    parametric <- utils::tail(parametric, ncol(Z))
    names(Z) <- parametric
    X <- cbind(X, Z)
  }
  list(data = X, learners = learners, internal = internal, n = as.integer(n),
       boundaries = boundaries, parametric = parametric, expand_Z = expand_Z,
       weights = weights)
}

.stageB_stability_basis <- function(ctx, internal, data = ctx$data, prediction = FALSE,
                                    sparse = FALSE, block_cache = NULL) {
  blocks <- lapply(names(internal), function(bl) {
    cached <- block_cache[[bl]]
    if (!is.null(cached) && identical(cached$knots, internal[[bl]]) &&
        identical(cached$sparse, sparse)) return(cached$basis)
    vars <- ctx$learners[[bl]]$variables
    k <- if (length(vars) == 1L) list(internal[[bl]]) else internal[[bl]]
    axes <- lapply(seq_along(vars), function(j) {
      boundary <- ctx$boundaries[[vars[j]]]
      x <- data[[vars[j]]]
      if (prediction) x <- pmax(boundary[1L], pmin(boundary[2L], x))
      splines::splineDesign(sort(c(rep(boundary, ctx$n), k[[j]])),
                           x, ord = ctx$n, outer.ok = TRUE,
                           sparse = sparse && length(vars) == 1L)
    })
    basis <- if (length(axes) == 1L) axes[[1L]] else {
      block <- tensorProdND(axes)
      if (sparse) Matrix::Matrix(block, sparse = TRUE) else block
    }
    if (!is.null(block_cache)) {
      block_cache[[bl]] <- list(knots = internal[[bl]], sparse = sparse, basis = basis)
    }
    basis
  })
  if (length(ctx$parametric)) {
    Z <- data[ctx$parametric]
    if (ctx$expand_Z) {
      Z <- stats::model.matrix(~ ., Z)
      Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    } else Z <- as.matrix(Z)
    blocks <- c(blocks, list(Z))
  }
  if (length(blocks)) do.call(cbind, blocks) else matrix(numeric(), nrow(data), 0L)
}

.stageB_stability_probe_basis <- function(ctx, internal, probes, rows, sparse) {
  cache <- NULL
  if (!is.null(ctx$probe_cache)) {
    first <- as.character(rows[1L])
    cache <- ctx$probe_cache[[first]]
    if (is.null(cache)) {
      cache <- new.env(parent = emptyenv())
      ctx$probe_cache[[first]] <- cache
    }
  }
  .stageB_stability_basis(ctx, internal, probes[rows, , drop = FALSE], TRUE, sparse, cache)
}

# LAPACK's divide-and-conquer SVD (dgesdd) can fail to converge on some
# rank-deficient designs. The transposed problem, or the SVD of the pivoted QR
# factor, gives the same decomposition and usually converges.
# 'svd_fun' is an argument only so that tests can exercise the fallbacks.
.stageB_svd <- function(X, size, svd_fun = svd) {
  s <- tryCatch(svd_fun(X, nu = size, nv = size), error = function(e) NULL)
  if (!is.null(s)) return(s)
  s <- tryCatch(svd_fun(t(X), nu = size, nv = size), error = function(e) NULL)
  if (!is.null(s)) return(list(d = s$d, u = s$v, v = s$u))
  q <- qr(X)
  r <- qr.R(q)[, order(q$pivot), drop = FALSE]
  s <- svd(r, nu = min(dim(r)), nv = size)
  list(d = s$d[seq_len(size)], u = qr.Q(q) %*% s$u[, seq_len(size), drop = FALSE],
       v = s$v)
}

.stageB_stability_factor <- function(ctx, internal) {
  B <- .stageB_stability_basis(ctx, internal)
  w <- ctx$weights / mean(ctx$weights[ctx$weights > 0])
  Xw <- B * sqrt(w)
  size <- min(dim(Xw))
  s <- .stageB_svd(Xw, size)
  keep <- which(s$d > max(s$d) * max(dim(Xw)) * .Machine$double.eps)
  if (!length(keep)) stop("The Stage B design has no estimable directions.", call. = FALSE)
  U <- s$u[, keep, drop = FALSE]
  V <- s$v[, keep, drop = FALSE]
  M <- sweep(V, 2L, s$d[keep], "/")
  uniform <- all(w == 1)
  T <- if (uniform) diag(length(keep)) else crossprod(U, U * w)
  root <- if (uniform) T else chol(T)
  # Economy SVD avoids allocating p by p matrices for wide designs.
  wide <- ncol(B) > nrow(B)
  null <- if (wide) NULL else s$v[, setdiff(seq_len(ncol(B)), keep), drop = FALSE]
  list(B = B, M = M, T = T, root = root, L = if (uniform) M else M %*% t(root), V = V,
       null = null, wide = wide, rank = length(keep), columns = ncol(B),
       inverse = 1 / s$d[keep], uniform = uniform)
}

.stageB_stability_evaluate <- function(ctx, internal, probes, factor,
                                      bound = NULL, keep_map = FALSE, incumbent = NULL) {
  s <- numeric(nrow(probes))
  null_ratio <- 0
  G <- if (keep_map) matrix(0, nrow(probes), factor$rank) else NULL
  for (start in seq.int(1L, nrow(probes), by = 2048L)) {
    rows <- start:min(nrow(probes), start + 2047L)
    sparse <- factor$columns > 32L
    P <- .stageB_stability_probe_basis(ctx, internal, probes, rows, sparse)
    PV <- as.matrix(P %*% factor$V)
    map <- sweep(PV, 2L, factor$inverse, "*")
    s[rows] <- if (factor$uniform) rowSums(map^2) else rowSums((map %*% t(factor$root))^2)
    if (keep_map) G[rows, ] <- map
    if (factor$wide) {
      norm <- if (sparse) as.numeric(Matrix::rowSums(P^2)) else rowSums(P^2)
      norm <- pmax(norm, 1e-30)
      delta <- pmax(0, 1 - rowSums(PV^2) / norm)
      # Subtraction is fast for wide designs, but unreliable near zero.
      ambiguous <- which(delta < sqrt(.Machine$double.eps))
      if (length(ambiguous)) {
        remainder <- as.matrix(P[ambiguous, , drop = FALSE]) -
          PV[ambiguous, , drop = FALSE] %*% t(factor$V)
        delta[ambiguous] <- rowSums(remainder^2) / norm[ambiguous]
      }
      null_ratio <- max(null_ratio, sqrt(delta))
    } else if (NCOL(factor$null)) {
      remainder <- as.matrix(P %*% factor$null)
      norm <- if (sparse) as.numeric(Matrix::rowSums(P^2)) else rowSums(P^2)
      null_ratio <- max(null_ratio, sqrt(rowSums(remainder^2) / pmax(norm, 1e-30)))
    }
    if (!is.null(incumbent) &&
        .stageB_stability_worse(null_ratio, max(s[rows] / bound[rows]), incumbent)) {
      return(NULL)
    }
  }
  if (is.null(bound)) bound <- rep(1, length(s))
  list(s = s, ratio = max(s / bound), null_ratio = null_ratio, G = G,
       B = factor$B, M = factor$M, T = factor$T, root = factor$root,
       probes = if (keep_map) probes else NULL)
}

.stageB_stability_worse <- function(null_ratio, ratio, incumbent) {
  # Rankings compare the non-identifiable part first. Both maxima can only
  # increase as more probes are evaluated.
  null_score <- if (null_ratio > 1e-8) null_ratio else 0
  best_null <- if (incumbent$null_ratio > 1e-8) incumbent$null_ratio else 0
  null_score > best_null || (null_score == best_null &&
    (isTRUE(ratio > incumbent$ratio) || (is.na(ratio) && !is.na(incumbent$ratio))))
}

# Exhaustive deletion choice, used when the fast scores are unavailable (e.g.
# joint tensor learners). It returns the first choice with the smallest
# (non-identifiable, sensitivity) ranking. With 'state', each choice is first
# evaluated on the probes most sensitive under the current knots. A maximum
# over a subset cannot exceed the full maximum, so choices whose subset ranking
# is already worse than a fully evaluated choice are skipped. The winner is
# unchanged; only the order of full evaluations differs.
.stageB_stability_choose <- function(ctx, choices, probes, bound, state = NULL,
                                     subset_size = 512L) {
  factors <- vector("list", length(choices))
  factor_of <- function(j) {
    if (is.null(factors[[j]])) factors[[j]] <<- .stageB_stability_factor(ctx, choices[[j]]$internal)
    factors[[j]]
  }
  visit <- seq_along(choices)
  lower <- NULL
  if (!is.null(state) && length(choices) > 1L && nrow(probes) > 4L * subset_size) {
    rows <- utils::head(order(state$s / bound, decreasing = TRUE), subset_size)
    # A separate cache: the full-probe cache is keyed by chunk position.
    sub <- ctx
    sub$probe_cache <- new.env(parent = emptyenv())
    lower <- lapply(seq_along(choices), function(j) .stageB_stability_evaluate(
      sub, choices[[j]]$internal, probes[rows, , drop = FALSE], factor_of(j), bound[rows]))
    null <- vapply(lower, `[[`, numeric(1), "null_ratio")
    ratio <- vapply(lower, `[[`, numeric(1), "ratio")
    visit <- order(ifelse(null > 1e-8, null, 0), ratio, seq_along(choices))
  }
  shrink <- 1 - 1e-8  # guards the subset bound against rounding differences
  best <- NULL
  winner <- 1L
  for (j in visit) {
    if (!is.null(best) && !is.null(lower) &&
        .stageB_stability_worse(lower[[j]]$null_ratio * shrink, lower[[j]]$ratio * shrink, best)) {
      next
    }
    candidate <- .stageB_stability_evaluate(ctx, choices[[j]]$internal, probes,
      factor_of(j), bound, incumbent = best)
    if (is.null(candidate)) next
    if (is.null(best) || .stageB_stability_worse(best$null_ratio, best$ratio, candidate) ||
        (j < winner && !.stageB_stability_worse(candidate$null_ratio, candidate$ratio, best))) {
      best <- candidate
      winner <- j
    }
  }
  winner
}

.stageB_stability_probes <- function(ctx) {
  d <- ctx$data
  parts <- list(d)
  numeric_columns <- names(d)[vapply(d, is.numeric, logical(1))]
  vars <- unique(unlist(lapply(ctx$learners[names(ctx$internal)], `[[`, "variables")))
  for (v in vars) {
    k <- unlist(lapply(names(ctx$internal), function(bl) {
      bv <- ctx$learners[[bl]]$variables
      if (!v %in% bv) return(NULL)
      if (length(bv) == 1L) ctx$internal[[bl]] else ctx$internal[[bl]][[v]]
    }))
    x <- sort(unique(d[[v]]))
    if (length(x) < 2L) next
    breaks <- sort(unique(c(x, k)))
    grid <- sort(unique(unlist(lapply(seq_len(length(breaks) - 1L), function(j) {
      seq(breaks[j], breaks[j + 1L], length.out = 33L)
    }))))
    lo <- pmax(1L, pmin(length(x) - 1L, findInterval(grid, x)))
    i <- match(x[lo], d[[v]])
    j <- match(x[lo + 1L], d[[v]])
    t <- (grid - x[lo]) / (x[lo + 1L] - x[lo])
    a <- d[i, , drop = FALSE]
    b <- d[j, , drop = FALSE]
    for (z in numeric_columns) a[[z]] <- (1 - t) * a[[z]] + t * b[[z]]
    for (z in setdiff(names(d), numeric_columns)) a[[z]][t == 1] <- b[[z]][t == 1]
    a[[v]] <- grid
    parts[[length(parts) + 1L]] <- a
  }
  unique(do.call(rbind, parts))
}

.stageB_stability_candidates <- function(ctx, internal) {
  choices <- list()
  for (bl in names(internal)) {
    vars <- ctx$learners[[bl]]$variables
    axes <- if (length(vars) == 1L) list(internal[[bl]]) else internal[[bl]]
    for (a in seq_along(axes)) for (j in seq_along(axes[[a]])) {
      k <- internal
      if (length(vars) == 1L) k[bl] <- list(k[[bl]][-j]) else k[[bl]][a] <- list(k[[bl]][[a]][-j])
      choices[[length(choices) + 1L]] <- list(internal = k, learner = bl,
        axis = a, index = j, knot = axes[[a]][j])
    }
  }
  choices
}

.stageB_stability_scores <- function(ctx, internal, state, choices, bound, sparse = TRUE) {
  if (state$null_ratio > 1e-8 || any(vapply(ctx$learners[names(internal)],
      function(b) length(b$variables) != 1L, logical(1)))) return(NULL)
  if (any(vapply(internal, function(k) anyDuplicated(k) > 0L, logical(1)))) return(NULL)
  sizes <- lengths(internal) + ctx$n
  starts <- c(0L, utils::head(cumsum(sizes), -1L))
  D <- matrix(0, length(choices), ncol(state$B))
  # Removing a simple knot imposes one zero highest-derivative-jump constraint.
  for (bl in names(internal)) {
    indices <- which(vapply(choices, function(ch) ch$learner == bl, logical(1)))
    if (!length(indices)) next
    k <- internal[[bl]]
    j <- vapply(choices[indices], `[[`, integer(1), "index")
    boundary <- ctx$boundaries[[ctx$learners[[bl]]$variables]]
    left <- c(boundary[1L], k)[j]
    right <- c(k, boundary[2L])[j + 1L]
    full <- sort(c(rep(boundary, ctx$n), k))
    jump <- splines::splineDesign(full, (k[j] + right) / 2, ord = ctx$n,
                                 derivs = rep(ctx$n - 1L, length(j)), outer.ok = TRUE) -
      splines::splineDesign(full, (left + k[j]) / 2, ord = ctx$n,
                           derivs = rep(ctx$n - 1L, length(j)), outer.ok = TRUE)
    columns <- starts[match(bl, names(internal))] + seq_len(sizes[[bl]])
    for (i in seq_along(indices))
      D[indices[i], columns] <- jump[i, ] / sqrt(sum(jump[i, ]^2))
  }
  Z <- t(D %*% state$M)
  norms <- sqrt(colSums(Z^2))
  if (any(!is.finite(norms) | norms < 1e-12)) return(NULL)
  Z <- sweep(Z, 2L, norms, "/")
  uniform <- identical(state$T, diag(nrow(state$T)))
  use_sparse <- sparse && ncol(state$G) > 32L && !is.null(state$probes)
  TZ <- if (uniform) Z else state$T %*% Z
  if (use_sparse) {
    # G = P M; associating the small product first exploits local spline support.
    MZ <- state$M %*% Z
    MTZ <- if (!uniform) state$M %*% TZ
    n_probes <- nrow(state$probes)
    chunk <- if (is.null(ctx$probe_cache)) n_probes else 2048L
    A <- matrix(0, n_probes, ncol(Z))
    if (!uniform) Q <- matrix(0, n_probes, ncol(Z))
    for (start in seq.int(1L, n_probes, by = chunk)) {
      rows <- start:min(n_probes, start + chunk - 1L)
      P <- .stageB_stability_probe_basis(ctx, internal, state$probes, rows, sparse = TRUE)
      A[rows, ] <- as.matrix(P %*% MZ)
      if (!uniform) Q[rows, ] <- as.matrix(P %*% MTZ)
    }
  } else {
    A <- state$G %*% Z
    if (!uniform) Q <- state$G %*% TZ
  }
  if (uniform) Q <- A
  scores <- vapply(seq_along(choices), function(i) {
    candidate <- state$s - 2 * A[, i] * Q[, i] + A[, i]^2 * sum(Z[, i] * TZ[, i])
    danger <- which(candidate < sqrt(.Machine$double.eps) * state$s)
    if (length(danger)) {
      residual <- state$G[danger, , drop = FALSE] - tcrossprod(A[danger, i], Z[, i])
      candidate[danger] <- if (uniform) rowSums(residual^2) else
        rowSums((residual %*% t(state$root))^2)
    }
    max(pmax(0, candidate) / bound)
  }, numeric(1))
  if (!use_sparse) return(scores)
  if (!all(is.finite(scores)))
    return(.stageB_stability_scores(ctx, internal, state, choices, bound, sparse = FALSE))

  # Recheck screened winners and near ties using the original multiplication order.
  checked <- rep(FALSE, length(scores))
  repeat {
    best <- min(scores)
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(best))
    check <- which(!checked & scores <= best + tolerance)
    if (!length(check)) break
    scores[check] <- .stageB_stability_scores(ctx, internal, state, choices[check],
                                            bound, sparse = FALSE)
    checked[check] <- TRUE
  }
  scores
}

.stageB_stability_screen <- function(ctx, internal, state, choices, bound) {
  indices <- seq_along(choices)
  if (length(indices) < 8L || ncol(state$G) <= 32L || length(bound) <= 2048L) {
    return(indices)
  }

  # A maximum over a subset cannot exceed the full score. Start with the
  # probes most sensitive under the current knot configuration.
  rows <- utils::head(order(state$s / bound, decreasing = TRUE), 512L)
  probe_state <- state
  probe_state$G <- state$G[rows, , drop = FALSE]
  probe_state$s <- state$s[rows]
  lower <- .stageB_stability_scores(ctx, internal, probe_state, choices, bound[rows],
                                    sparse = FALSE)
  if (is.null(lower) || any(!is.finite(lower))) return(indices)

  candidate <- which.min(lower)
  upper <- .stageB_stability_scores(ctx, internal, state, choices[candidate], bound,
                                    sparse = FALSE)
  if (is.null(upper) || !is.finite(upper)) return(indices)
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(upper))
  indices[lower <= upper + tolerance | indices == candidate]
}

# Each B-spline block is a nonnegative partition of unity, also for tensors.
# Centre each block of a coefficient map A at its mean. Convexity and the
# triangle inequality bound ||b A|| by the norm of the summed centres plus
# the largest row deviation per block. Parametric coordinates range between
# their observed extrema (probe interpolation stays within these ranges).
.stageB_stability_certificate <- function(ctx, internal, factor) {
  if (factor$wide) return(list(sensitivity = Inf, null_ratio = Inf))
  sizes <- vapply(names(internal), function(bl) {
    k <- internal[[bl]]
    if (length(ctx$learners[[bl]]$variables) == 1L) length(k) + ctx$n else
      prod(lengths(k) + ctx$n)
  }, numeric(1))
  bound <- function(A) {
    if (!ncol(A)) return(0)
    centre <- numeric(ncol(A))
    radius <- 0
    pos <- 0L
    for (size in sizes) {
      block <- A[pos + seq_len(size), , drop = FALSE]
      mu <- colMeans(block)
      centre <- centre + mu
      radius <- radius + sqrt(max(rowSums(sweep(block, 2L, mu, "-")^2)))
      pos <- pos + size
    }
    if (pos < nrow(A)) {
      Z <- ctx$data[ctx$parametric]
      if (ctx$expand_Z) Z <- stats::model.matrix(~ ., Z)[, -1L, drop = FALSE]
      Z <- as.matrix(Z)
      lo <- apply(Z, 2L, min)
      hi <- apply(Z, 2L, max)
      block <- A[seq.int(pos + 1L, nrow(A)), , drop = FALSE]
      centre <- centre + as.numeric((lo / 2 + hi / 2) %*% block)
      radius <- radius + sum((hi / 2 - lo / 2) * sqrt(rowSums(block^2)))
    }
    sqrt(sum(centre^2)) + radius
  }
  # ||b||^2 >= sum(1 / block_size), independently of parametric effects.
  list(sensitivity = bound(factor$L)^2,
       null_ratio = bound(factor$null) / sqrt(sum(1 / sizes)))
}

.stageB_stabilise <- function(X, Z, weights, learners, internal, n,
                              boundaries = lapply(as.data.frame(X), range), expand_Z = TRUE) {
  diagnostic <- list(method = "reference", threshold = 10, probe_resolution = 33L,
    removed = 0L, original_knots = internal, knots = internal, deleted = NULL,
    target_met = NA, max_ratio = NA_real_, reference_max = NA_real_, certified = FALSE)
  if (!n %in% c(3L, 4L) || !length(internal))
    return(list(internal = internal, stability = diagnostic))
  ctx <- .stageB_stability_context(X, Z, weights, learners, internal, n, boundaries, expand_Z)
  factor <- .stageB_stability_factor(ctx, ctx$internal)
  single <- length(ctx$internal) == 1L && !length(ctx$parametric) &&
    length(ctx$learners[[names(ctx$internal)[1L]]]$variables) == 1L && factor$rank == factor$columns
  # A nonnegative partition of unity bounds sensitivity by coefficient-row norms.
  certificate <- if (single) max(rowSums(factor$L^2)) else Inf
  if (certificate <= 10) {
    diagnostic$target_met <- TRUE
    diagnostic$max_ratio <- certificate / 10
    diagnostic$certified <- TRUE
    return(list(internal = internal, stability = diagnostic))
  }
  certificate <- .stageB_stability_certificate(ctx, ctx$internal, factor)
  if (is.finite(certificate$sensitivity) &&
      certificate$sensitivity <= 10 * (1 - 1e-8) &&
      certificate$null_ratio <= 1e-10) {
    empty <- lapply(ctx$internal, function(k) if (is.list(k)) lapply(k, function(z) NULL) else NULL)
    reference_factor <- .stageB_stability_factor(ctx, empty)
    reference_certificate <- .stageB_stability_certificate(ctx, empty, reference_factor)
    # Preserve the reference-identifiability check before skipping the probes.
    if (is.finite(reference_certificate$null_ratio) &&
        reference_certificate$null_ratio <= 1e-10) {
      diagnostic$target_met <- TRUE
      diagnostic$max_ratio <- certificate$sensitivity / 10
      diagnostic$certified <- TRUE
      return(list(internal = internal, stability = diagnostic))
    }
  }
  probes <- .stageB_stability_probes(ctx)
  # One basis per learner and probe chunk, replaced only when its knots change.
  # Keep large designs on the streaming path to bound the extra memory.
  univariate <- all(vapply(ctx$learners[names(ctx$internal)], function(learner) {
    length(learner$variables) == 1L
  }, logical(1)))
  bytes_per_row <- if (factor$columns > 32L && univariate && !length(ctx$parametric)) {
    # At most n nonzero entries per row per spline block (double and row index).
    max(12 * n * length(ctx$internal), 8 * 32)
  } else 12 * factor$columns
  if (nrow(probes) * bytes_per_row <= 80 * 1024^2) {
    ctx$probe_cache <- new.env(parent = emptyenv())
  }
  empty <- lapply(ctx$internal, function(k) if (is.list(k)) lapply(k, function(z) NULL) else NULL)
  reference <- .stageB_stability_evaluate(ctx, empty, probes, .stageB_stability_factor(ctx, empty))
  diagnostic$reference_max <- max(reference$s)
  if (reference$null_ratio > 1e-8) {
    diagnostic$status <- "reference_not_identifiable"
    diagnostic$target_met <- FALSE
    warning(warningCondition(paste0("The order-", n, " Stage B stability reference is not identifiable; ",
      "the original fit and knots are retained."), class = "geds_stageB_stability_warning"))
    return(list(internal = internal, stability = diagnostic))
  }
  # Knots may add at most 10 to the knot-free reference sensitivity. A fixed
  # floor of 10 removed knots whose small increase pushed an already sensitive
  # reference (e.g. a knot-free joint learner) over the threshold.
  bound <- reference$s + 10
  k <- ctx$internal
  state <- .stageB_stability_evaluate(ctx, k, probes, factor, bound, TRUE)
  path <- list()
  repeat {
    if (state$ratio <= 1 + 1e-10 && state$null_ratio <= 1e-8) break
    choices <- .stageB_stability_candidates(ctx, k)
    if (!length(choices)) break
    eligible <- .stageB_stability_screen(ctx, k, state, choices, bound)
    j <- eligible[1L]
    if (length(eligible) > 1L || length(choices) == 1L) {
      scores <- .stageB_stability_scores(ctx, k, state, choices[eligible], bound)
      j <- if (is.null(scores)) .stageB_stability_choose(ctx, choices, probes, bound, state) else
        eligible[which.min(scores)]
    }
    ch <- choices[[j]]
    k <- ch$internal
    before <- state$ratio
    factor <- .stageB_stability_factor(ctx, k)
    state <- .stageB_stability_evaluate(ctx, k, probes, factor, bound, TRUE)
    path[[length(path) + 1L]] <- data.frame(learner = ch$learner, axis = ch$axis,
      knot = ch$knot, before = before, after = state$ratio)
  }
  diagnostic$removed <- length(path)
  diagnostic$deleted <- if (length(path)) do.call(rbind, path) else NULL
  diagnostic$max_ratio <- state$ratio
  diagnostic$target_met <- state$ratio <= 1 + 1e-10 && state$null_ratio <= 1e-8
  if (!diagnostic$target_met)
    warning(warningCondition(paste0("The order-", n, " Stage B stability criterion could not be met."),
      class = "geds_stageB_stability_warning"))
  if (length(path)) for (bl in names(k)) {
    if (is.list(k[[bl]])) names(k[[bl]]) <- names(internal[[bl]])
    internal[bl] <- k[bl]
  }
  diagnostic$knots <- internal
  list(internal = internal, stability = diagnostic)
}

