stability_fixture <- function(n = 4L) {
  x <- c(seq(0, 0.4, length.out = 20), seq(0.8, 1, length.out = 10))
  list(X = data.frame(x = x), Z = NULL, weights = rep(1, length(x)),
       learners = list(sx = list(type = "GeDS", variables = "x")),
       internal = list(sx = c(0.45, 0.45001, 0.45002, 0.45003, 0.75)), n = n)
}

test_that("the reference safeguard deletes knots and retains ordinary least squares", {
  for (n in 3:4) {
    a <- stability_fixture(n)
    set.seed(92)
    rng <- .Random.seed
    guarded <- do.call(GeDS:::.stageB_stabilise, a)
    expect_identical(.Random.seed, rng)
    expect_gt(guarded$stability$removed, 0)
    expect_true(guarded$stability$target_met)
    expect_lte(guarded$stability$max_ratio, 1 + 1e-10)
    expect_true(all(guarded$internal$sx %in% a$internal$sx))
    expect_equal(length(a$internal$sx) - length(guarded$internal$sx),
                 guarded$stability$removed)

    ctx <- do.call(GeDS:::.stageB_stability_context,
                   c(a, list(boundaries = list(x = c(0, 1)))))
    probes <- GeDS:::.stageB_stability_probes(ctx)
    evaluate <- function(k) GeDS:::.stageB_stability_evaluate(
      ctx, k, probes, GeDS:::.stageB_stability_factor(ctx, k))
    reference <- evaluate(list(sx = NULL))
    checked <- evaluate(guarded$internal)
    expect_lte(max(checked$s / (reference$s + 10)), 1 + 1e-10)
    expect_lte(checked$null_ratio, 1e-8)

    y <- sin(5 * a$X$x)
    fit <- GeDS:::SplineReg_LM(a$X$x, y, InterKnots = guarded$internal$sx,
                              n = n, extr = c(0, 1))
    design <- GeDS:::.stageB_stability_basis(ctx, guarded$internal)
    expect_equal(as.numeric(fit$predicted),
                 as.numeric(lm.fit(design, y)$fitted.values), tolerance = 1e-8)
    expect_null(fit$penalty)
  }
})

test_that("reference floor accommodates intrinsic knot-free sensitivity", {
  a <- stability_fixture()
  a$X <- data.frame(x = seq(0, 1, length.out = 60))
  a$weights <- as.numeric(a$X$x < 0.55)
  a$internal <- list(sx = NULL)
  guarded <- do.call(GeDS:::.stageB_stabilise, a)
  expect_gt(guarded$stability$reference_max, 10)
  expect_true(guarded$stability$target_met)
  expect_equal(guarded$stability$removed, 0)
  expect_identical(guarded$internal, a$internal)
})

test_that("well-supported fits and linear knots remain unchanged", {
  a <- stability_fixture()
  a$X <- data.frame(x = seq(0, 1, length.out = 80))
  a$weights <- rep(1, 80)
  a$internal <- list(sx = c(0.25, 0.5, 0.75))
  guarded <- do.call(GeDS:::.stageB_stabilise, a)
  expect_identical(guarded$internal, a$internal)
  expect_true(guarded$stability$certified)
  a$n <- 2L
  expect_identical(do.call(GeDS:::.stageB_stabilise, a)$internal, a$internal)
})

test_that("weighted sensitivity and fast deletion scores match direct calculations", {
  set.seed(93)
  X <- data.frame(x = sort(runif(50)), z = runif(50))
  Z <- data.frame(group = factor(rep(1:2, 25)))
  w <- rep(c(0, 1, 2, 3, 1), 10)
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  k <- list(sx = c(.3, .6), sz = c(.3, .6))
  ctx <- GeDS:::.stageB_stability_context(X, Z, w, learners, k, 4L, lapply(X, range))
  probes <- GeDS:::.stageB_stability_probes(ctx)
  factor <- GeDS:::.stageB_stability_factor(ctx, k)
  state <- GeDS:::.stageB_stability_evaluate(ctx, k, probes, factor, keep_map = TRUE)
  P <- GeDS:::.stageB_stability_basis(ctx, k, probes, TRUE)
  H <- P %*% factor$M %*% t(svd(factor$B * sqrt(w / mean(w[w > 0])),
                               nu = factor$rank, nv = 0)$u) *
    rep(sqrt(w / mean(w[w > 0])), each = nrow(P))
  expect_equal(state$s, unname(rowSums(H^2)), tolerance = 1e-7)
  bound <- pmax(10, seq(1, 20, length.out = nrow(probes)))
  choices <- GeDS:::.stageB_stability_candidates(ctx, k)
  scores <- GeDS:::.stageB_stability_scores(ctx, k, state, choices, bound)
  direct <- vapply(choices, function(ch) GeDS:::.stageB_stability_evaluate(
    ctx, ch$internal, probes, GeDS:::.stageB_stability_factor(ctx, ch$internal),
    bound)$ratio, numeric(1))
  expect_equal(scores, direct, tolerance = 1e-7)
  first <- GeDS:::.stageB_stabilise(X, Z, w, learners, k, 4L)
  scaled <- GeDS:::.stageB_stabilise(X, Z, 7 * w, learners, k, 4L)
  expect_identical(first$internal, scaled$internal)
  expect_equal(first$stability$max_ratio, scaled$stability$max_ratio, tolerance = 1e-9)
})

test_that("an unidentifiable reference warns and retains the original knots", {
  a <- stability_fixture()
  a$weights <- as.numeric(seq_len(30) <= 2)
  expect_warning(guarded <- do.call(GeDS:::.stageB_stabilise, a),
                 class = "geds_stageB_stability_warning")
  expect_identical(guarded$internal, a$internal)
  expect_false(guarded$stability$target_met)
  expect_identical(guarded$stability$status, "reference_not_identifiable")
})

test_that("public additive models preserve Stage A and agree on prediction", {
  set.seed(94)
  d <- data.frame(x = runif(65), z = runif(65), group = factor(rep(1:2, length.out = 65)))
  d$y <- sin(6 * d$x) + d$z + as.numeric(d$group) / 4 + rnorm(65, sd = .15)
  w <- rep(c(0, 1, 2, 1, 3), 13)
  for (fun in list(NGeDSgam, NGeDSboost)) {
    a <- list(formula = y ~ f(x) + f(z) + group, data = d, weights = w, phi = .8,
              penalty = "none")
    if (identical(fun, NGeDSboost)) a <- c(a, list(initial_learner = FALSE,
      max_iterations = 4, shrinkage = .5))
    fit <- suppressMessages(suppressWarnings(do.call(fun, a)))
    linear <- suppressMessages(suppressWarnings(do.call(fun, c(a, list(higher_order = FALSE)))))
    expect_identical(as.numeric(predict(fit, n = 2)), as.numeric(predict(linear, n = 2)))
    expect_identical(knots(fit, n = 2), knots(linear, n = 2))
    for (n in 3:4) {
      b <- fit$final_model[[if (n == 3) "quadratic.fit" else "cubic.fit"]]
      if (identical(fun, NGeDSboost)) expect_true(b$stability$target_met) else
        expect_null(b$stability)
      expect_null(b$penalty)
      expect_equal(as.numeric(predict(fit, n = n)),
                   as.numeric(predict(fit, newdata = d, n = n)), tolerance = 1e-8)
    }
  }
})

test_that("the Stage B safeguard supports non-Gaussian boosting families", {
  set.seed(95)
  d <- data.frame(x = runif(75), z = runif(75), off = rep(.05, 75))
  eta <- .3 + .4 * sin(5 * d$x) + .2 * d$z
  responses <- list(poisson = rpois(75, exp(eta)),
                    binomial = rbinom(75, 1, plogis(eta)),
                    Gamma = rgamma(75, shape = 8, scale = exp(eta) / 8))
  families <- list(poisson = mboost::Poisson(), binomial = mboost::Binomial(),
                   Gamma = mboost::GammaReg())
  for (name in names(families)) {
    d$y <- responses[[name]]
    a <- list(formula = y ~ f(x) + f(z), data = d,
      family = families[[name]], initial_learner = FALSE, shrinkage = .1,
      phi = .8, max_iterations = 3, penalty = "none")
    fit <- suppressMessages(suppressWarnings(do.call(NGeDSboost, a)))
    for (n in 3:4) {
      b <- fit$final_model[[if (n == 3) "quadratic.fit" else "cubic.fit"]]
      expect_true(b$stability$target_met)
      for (type in c("link", "response")) {
        expect_equal(as.numeric(predict(fit, n = n, type = type)),
          as.numeric(predict(fit, newdata = d, n = n, type = type)), tolerance = 1e-7)
      }
    }
  }
})

test_that("REML refits bypass the unpenalised safeguard", {
  skip_if_not_installed("mgcv")
  local_mocked_bindings(.stageB_stabilise = function(...) stop("unpenalised guard called"),
                        .package = "GeDS")
  set.seed(96)
  d <- data.frame(x = runif(65), z = runif(65))
  d$y <- sin(5 * d$x) + d$z^2 + rnorm(65, sd = .1)
  fit <- suppressMessages(suppressWarnings(NGeDSgam(y ~ f(x) + f(z), data = d,
    phi = .8, penalty = "component")))
  expect_identical(fit$final_model$cubic.fit$penalty$structure, "component")
  expect_null(fit$final_model$cubic.fit$stability)
})

test_that("NGeDSgam bypasses stability simplification for unpenalised refits", {
  local_mocked_bindings(.stageB_stabilise = function(...) stop("boost safeguard called"),
                       .package = "GeDS")
  set.seed(1)
  d <- data.frame(x = c(runif(20, 0, .4), runif(10, .8, 1)), z = runif(30))
  d$y <- sin(5 * d$x) + d$z^2 + rnorm(30, sd = .15)
  fit <- suppressMessages(suppressWarnings(NGeDSgam(y ~ f(x) + f(z), data = d,
                                                    penalty = "none")))
  for (n in 3:4) {
    key <- if (n == 3L) "quadratic.fit" else "cubic.fit"
    expect_null(fit$final_model[[key]]$stability)
    expected <- lapply(fit$internal_knots$linear.int.knots, function(k) {
      if (length(k) >= n - 1L) GeDS:::makenewknots(k, n) else NULL
    })
    expect_equal(knots(fit, n = n, options = "int"), expected)
    expect_equal(as.numeric(predict(fit, n = n)),
                 as.numeric(predict(fit, newdata = d, n = n)), tolerance = 1e-8)
  }
})

test_that("standalone fits bypass the safeguard and retain all averaged knots", {
  local_mocked_bindings(.stageB_stabilise = function(...) stop("standalone guard called"),
    .package = "GeDS")
  set.seed(97)
  d <- data.frame(x = runif(80), z = runif(80), v = runif(80))
  d <- d[order(d$x), ]
  d$off <- .2 + .1 * d$x
  d$y <- exp(.2 + .3 * sin(4 * d$x) + .2 * d$z + .1 * d$v) + d$off
  formulas <- list(y ~ f(x) + z + offset(off),
    y ~ f(x, z) + v + offset(off))
  w <- rep(c(0, 1, 2, 1), 20)
  for (v in c("x", "z", "v")) w[c(which.min(d[[v]]), which.max(d[[v]]))] <- 1
  for (dimension in seq_along(formulas)) for (generalised in c(FALSE, TRUE)) {
    formula <- formulas[[dimension]]
    weights <- w
    limit <- 4L
    fit <- suppressWarnings(if (generalised) GGeDS(formula,
      data = d, weights = weights, family = Gamma(link = "log"), beta = .5, max.intknots = limit) else
      NGeDS(formula, data = d, weights = weights, max.intknots = limit))
    linear <- fit$linear.intknots
    axes <- if (is.list(linear)) linear else list(x = linear)
    for (n in 3:4) {
      label <- if (n == 3L) "quadratic" else "cubic"
      expected <- lapply(axes, function(k) {
        if (length(k) >= n - 1L) GeDS:::makenewknots(k, n) else NULL
      })
      if (!is.list(linear)) expected <- expected[[1L]]
      expect_equal(fit[[paste0(label, ".intknots")]], expected)
      expect_null(fit[[paste0(label, ".fit")]]$stability)
      expect_equal(as.numeric(predict(fit, n = n)),
        as.numeric(predict(fit, newdata = d, n = n)), tolerance = 1e-7)
    }
  }
})

test_that("tensor simplification retains axis names and identifies prediction maps", {
  x <- stability_fixture()$X$x
  d <- expand.grid(x = x, z = seq(0, 1, length.out = 6))
  k <- list(Xk = c(.45, .45001, .45002, .45003, .75), Yk = NULL)
  a <- list(X = d, Z = NULL, weights = rep(1, nrow(d)),
    learners = list(surface = list(type = "GeDS", variables = names(d))),
    internal = list(surface = k))
  for (n in 3:4) {
    guarded <- do.call(GeDS:::.stageB_stabilise, c(a, list(n = n)))
    expect_gt(guarded$stability$removed, 0)
    expect_true(guarded$stability$target_met)
    expect_identical(names(guarded$internal$surface), names(k))
    expect_true(all(guarded$internal$surface$Xk %in% k$Xk))
    expect_null(guarded$internal$surface$Yk)
  }
})

test_that("non-Gaussian Stage B refits use the simplified knots without a penalty", {
  x <- stability_fixture()$X$x
  learners <- list(sx = list(type = "GeDS", variables = "x"))
  final <- list(base_learners = list(sx = list(
    knots = c(0, .4, .45, .45001, .45002, .45003, .75, .8, 1))))
  families <- list(mboost::Poisson(), mboost::Binomial(), mboost::GammaReg())
  responses <- list(rep(1:3, 10), rep(c(-1, 1), 15), exp(.3 * sin(5 * x)))
  for (j in seq_along(families)) {
    args <- list(predictors = data.frame(x = x), response = data.frame(y = responses[[j]]),
      weights = rep(1, 30), offset = rep(.05, 30), family = families[[j]],
      base_learners = learners)
    for (n in 3:4) {
      b <- suppressMessages(GeDS:::stageB_fit(n, args, "x", character(), "y", final,
        normalize_data = FALSE, weights = TRUE, offset = TRUE, link = TRUE))
      expect_gt(b$fit$stability$removed, 0)
      expect_true(b$fit$stability$target_met)
      expect_true(all(is.finite(b$preds)))
      expect_true(b$fit$temporary$converged)
      expect_null(b$fit$penalty)
      direct <- GeDS:::SplineReg_Multivar(args$predictors, responses[[j]],
        offset = args$offset, base_learners = learners, weights = args$weights,
        InterKnotsList = b$int.knots, n = n, family = families[[j]], link = NULL,
        only_pred = TRUE)
      expect_equal(b$preds, as.numeric(direct$predicted), tolerance = 1e-9)
    }
  }
})

test_that("repeated knots use direct candidate calculations", {
  a <- stability_fixture()
  a$internal <- list(sx = c(.3, .3, .7))
  ctx <- do.call(GeDS:::.stageB_stability_context,
                 c(a, list(boundaries = list(x = c(0, 1)))))
  k <- ctx$internal
  probes <- GeDS:::.stageB_stability_probes(ctx)
  state <- GeDS:::.stageB_stability_evaluate(ctx, k, probes,
    GeDS:::.stageB_stability_factor(ctx, k), keep_map = TRUE)
  choices <- GeDS:::.stageB_stability_candidates(ctx, k)
  expect_null(GeDS:::.stageB_stability_scores(ctx, k, state, choices, rep(10, nrow(probes))))
  guarded <- do.call(GeDS:::.stageB_stabilise, a)
  expect_true(guarded$stability$target_met)
})

test_that("sparse wide-design diagnostics agree with the direct response map", {
  set.seed(132)
  X <- data.frame(x = sort(runif(20)), z = runif(20))
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  k <- list(sx = seq(.01, .99, length.out = 40), sz = c(.3, .7))
  for (n in 3:4) for (w in list(rep(1, 20), rep(c(0, 1, 3, 1), 5))) {
    ctx <- GeDS:::.stageB_stability_context(X, NULL, w, learners, k, n, lapply(X, range))
    probes <- GeDS:::.stageB_stability_probes(ctx)
    fac <- GeDS:::.stageB_stability_factor(ctx, k)
    expect_true(fac$wide)
    state <- GeDS:::.stageB_stability_evaluate(ctx, k, probes, fac, keep_map = TRUE)
    P <- GeDS:::.stageB_stability_basis(ctx, k, probes, TRUE)
    expect_equal(state$s, unname(rowSums((P %*% fac$L)^2)), tolerance = 1e-9)
    expect_equal(state$G, unname(P %*% fac$M), tolerance = 1e-9)
    remainder <- P - (P %*% fac$V) %*% t(fac$V)
    direct <- max(sqrt(rowSums(remainder^2) / pmax(rowSums(P^2), 1e-30)))
    expect_equal(state$null_ratio, direct, tolerance = 1e-9)
  }
})

test_that("overcomplete weighted additive designs are coarsened before ordinary refitting", {
  set.seed(133)
  X <- data.frame(x = runif(18), z = runif(18))
  w <- rep(c(0, 1, 2), 6)
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  k <- lapply(X, function(x) seq(min(x), max(x), length.out = 10)[2:9])
  names(k) <- names(learners)
  y <- sin(4 * X$x) + X$z^2
  for (n in 3:4) {
    ctx <- GeDS:::.stageB_stability_context(X, NULL, w, learners, k, n, lapply(X, range))
    probes <- GeDS:::.stageB_stability_probes(ctx)
    fac <- GeDS:::.stageB_stability_factor(ctx, k)
    expect_true(fac$wide)
    initial <- GeDS:::.stageB_stability_evaluate(ctx, k, probes, fac)
    expect_gt(initial$null_ratio, 1e-8)
    guarded <- GeDS:::.stageB_stabilise(X, NULL, w, learners, k, n)
    expect_gt(guarded$stability$removed, 0)
    expect_true(guarded$stability$target_met)
    checked <- GeDS:::.stageB_stability_evaluate(ctx, guarded$internal, probes,
      GeDS:::.stageB_stability_factor(ctx, guarded$internal))
    reference <- GeDS:::.stageB_stability_evaluate(ctx, list(sx = NULL, sz = NULL), probes,
      GeDS:::.stageB_stability_factor(ctx, list(sx = NULL, sz = NULL)))
    expect_lte(checked$null_ratio, 1e-8)
    expect_lte(max(checked$s / (reference$s + 10)), 1 + 1e-10)
    b <- GeDS:::SplineReg_LM_Multivar(X, y, base_learners = learners, weights = w,
      InterKnotsList = guarded$internal, n = n, coefficients = NULL,
      linear.predictors = NULL, only_pred = TRUE)
    B <- GeDS:::.stageB_stability_basis(ctx, guarded$internal)
    expect_equal(as.numeric(b$predicted), as.numeric(B %*% b$theta), tolerance = 1e-8)
    expect_equal(as.numeric(b$predicted),
      as.numeric(lm.wfit(B, y, w)$fitted.values), tolerance = 1e-8)
    expect_null(b$penalty)
  }
})

test_that("knots may add up to 10 to an already sensitive knot-free reference", {
  skip_if_not_installed("TH.data")
  data("bodyfat", package = "TH.data", envir = environment())
  # A small Bodyfat split whose knot-free joint hip-waist cubic is already near
  # the threshold. The old floor max(10, s0) removed the stable kneebreadth knot.
  set.seed(123)
  N <- nrow(bodyfat)
  cal <- sample(seq_len(N), floor(0.4 * N))
  ev <- setdiff(seq_len(N), cal)
  set.seed(40005)
  train <- bodyfat[c(cal, sample(ev, floor(0.7 * N) - length(cal))), ]
  set.seed(123)
  invisible(capture.output(fit <- suppressWarnings(NGeDSboost(
    DEXfat ~ age + f(hipcirc, waistcirc) + f(kneebreadth), data = train,
    initial_learner = FALSE, shrinkage = 0.8, phi = 0.8, q = 2, penalty = "none"))))
  st <- fit$final_model$cubic.fit$stability
  expect_equal(st$removed, 0L)
  expect_true(st$target_met)
  expect_false(st$certified)

  X <- train[c("hipcirc", "waistcirc", "kneebreadth")]
  learners <- fit$args$base_learners[names(st$original_knots)]
  ctx <- GeDS:::.stageB_stability_context(X, train["age"], rep(1, nrow(train)), learners,
                                          st$original_knots, 4L, lapply(X, range))
  probes <- GeDS:::.stageB_stability_probes(ctx)
  evaluate <- function(k) GeDS:::.stageB_stability_evaluate(
    ctx, k, probes, GeDS:::.stageB_stability_factor(ctx, k))
  empty <- lapply(ctx$internal, function(k) if (is.list(k)) lapply(k, function(z) NULL) else NULL)
  reference <- evaluate(empty)
  original <- evaluate(ctx$internal)
  expect_gt(max(original$s / pmax(10, reference$s)), 1)
  expect_lte(max(original$s / (reference$s + 10)), 1)
  expect_equal(st$reference_max, max(reference$s), tolerance = 1e-8)
})

test_that("the stability SVD falls back when LAPACK fails to converge", {
  # dgesdd can fail on some rank-deficient designs (seen on a joint tensor with
  # rank 120 of 170 columns); the fallbacks must give the same decomposition.
  set.seed(7)
  for (dims in list(c(40, 12), c(10, 25))) {
    X <- matrix(rnorm(dims[1] * 4), dims[1], 4) %*% matrix(rnorm(4 * dims[2]), 4, dims[2])
    size <- min(dim(X))
    direct <- svd(X, nu = size, nv = size)
    fail_once <- local({ calls <- 0L; function(...) {
      calls <<- calls + 1L; if (calls == 1L) stop("dgesdd") else svd(...) } })
    always_fail <- function(...) stop("dgesdd")
    for (svd_fun in list(fail_once, always_fail)) {
      s <- GeDS:::.stageB_svd(X, size, svd_fun)
      expect_equal(s$d, direct$d, tolerance = 1e-8)
      expect_equal(dim(s$u), dim(direct$u))
      expect_equal(dim(s$v), dim(direct$v))
      expect_equal(s$u %*% (s$d * t(s$v)), X, tolerance = 1e-8)
      # Projections onto the estimable directions agree with the direct SVD.
      keep <- which(direct$d > max(direct$d) * 1e-10)
      expect_equal(tcrossprod(s$v[, keep]), tcrossprod(direct$v[, keep]), tolerance = 1e-8)
    }
  }
})
