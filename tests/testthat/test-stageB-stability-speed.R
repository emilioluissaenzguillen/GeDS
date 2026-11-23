legacy_stability_scores <- function(ctx, internal, state, choices, bound) {
  sizes <- lengths(internal) + ctx$n
  starts <- c(0L, utils::head(cumsum(sizes), -1L))
  D <- matrix(0, length(choices), ncol(state$B))
  for (i in seq_along(choices)) {
    ch <- choices[[i]]
    bl <- ch$learner
    k <- internal[[bl]]
    j <- ch$index
    boundary <- ctx$boundaries[[ctx$learners[[bl]]$variables]]
    left <- if (j == 1L) boundary[1L] else k[j - 1L]
    right <- if (j == length(k)) boundary[2L] else k[j + 1L]
    full <- sort(c(rep(boundary, ctx$n), k))
    jump <- splines::splineDesign(full, (k[j] + right) / 2, ord = ctx$n,
      derivs = ctx$n - 1L, outer.ok = TRUE) -
      splines::splineDesign(full, (left + k[j]) / 2, ord = ctx$n,
        derivs = ctx$n - 1L, outer.ok = TRUE)
    columns <- starts[match(bl, names(internal))] + seq_len(sizes[[bl]])
    D[i, columns] <- jump / sqrt(sum(jump^2))
  }
  Z <- t(D %*% state$M)
  Z <- sweep(Z, 2L, sqrt(colSums(Z^2)), "/")
  A <- state$G %*% Z
  TZ <- state$T %*% Z
  Q <- state$G %*% TZ
  vapply(seq_along(choices), function(i) {
    candidate <- state$s - 2 * A[, i] * Q[, i] + A[, i]^2 * sum(Z[, i] * TZ[, i])
    danger <- which(candidate < sqrt(.Machine$double.eps) * state$s)
    if (length(danger)) {
      residual <- state$G[danger, , drop = FALSE] - tcrossprod(A[danger, i], Z[, i])
      candidate[danger] <- rowSums((residual %*% t(state$root))^2)
    }
    max(pmax(0, candidate) / bound)
  }, numeric(1))
}

speed_stability_state <- function(n = 4L, weighted = FALSE, knots = 20L) {
  set.seed(719)
  X <- data.frame(x = seq(0, 1, length.out = 120),
                  z = sample(seq(0, 1, length.out = 120)))
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  internal <- list(sx = seq(.04, .96, length.out = knots),
                   sz = seq(.05, .95, length.out = knots))
  w <- if (weighted) rep(c(1, 0, 2, 1), 30) else rep(1, 120)
  ctx <- GeDS:::.stageB_stability_context(X, NULL, w, learners, internal, n, lapply(X, range))
  probes <- GeDS:::.stageB_stability_probes(ctx)
  state <- GeDS:::.stageB_stability_evaluate(ctx, internal, probes,
    GeDS:::.stageB_stability_factor(ctx, internal), rep(10, nrow(probes)), keep_map = TRUE)
  list(ctx = ctx, internal = internal, state = state,
    choices = GeDS:::.stageB_stability_candidates(ctx, internal), bound = rep(10, nrow(probes)))
}

test_that("sparse candidate scoring retains the legacy winner under prior weights", {
  for (n in 3:4) for (weighted in c(FALSE, TRUE)) {
    a <- speed_stability_state(n, weighted)
    expect_gt(ncol(a$state$G), 32L)
    expect_lte(a$state$null_ratio, 1e-8)
    old <- do.call(legacy_stability_scores, a)
    dense <- do.call(GeDS:::.stageB_stability_scores, c(a, list(sparse = FALSE)))
    fast <- do.call(GeDS:::.stageB_stability_scores, a)
    expect_identical(dense, old)
    expect_equal(fast, old, tolerance = 1e-10)
    expect_identical(which.min(fast), which.min(old))
  }
})

test_that("screened candidate winners are rechecked with dense arithmetic", {
  a <- speed_stability_state()
  X <- expand.grid(x = seq(0, 1, length.out = 30), z = seq(0, 1, length.out = 30))
  a$internal <- list(sx = seq(.05, .95, length.out = 15),
                     sz = seq(.05, .95, length.out = 15))
  a$ctx <- GeDS:::.stageB_stability_context(X, NULL, rep(1, nrow(X)),
    a$ctx$learners, a$internal, 4L, lapply(X, range))
  probes <- GeDS:::.stageB_stability_probes(a$ctx)
  a$bound <- rep(10, nrow(probes))
  a$state <- GeDS:::.stageB_stability_evaluate(a$ctx, a$internal, probes,
    GeDS:::.stageB_stability_factor(a$ctx, a$internal), a$bound, keep_map = TRUE)
  a$choices <- GeDS:::.stageB_stability_candidates(a$ctx, a$internal)
  score <- GeDS:::.stageB_stability_scores
  rechecks <- integer()
  local_mocked_bindings(.stageB_stability_scores = function(ctx, internal, state, choices,
    bound, sparse = TRUE) {
    if (!sparse) rechecks <<- c(rechecks, length(choices))
    score(ctx, internal, state, choices, bound, sparse)
  }, .package = "GeDS")
  fast <- do.call(score, a)
  expect_gt(length(rechecks), 0L)
  expect_true(all(rechecks > 0L))
  expect_true(any(rechecks > 1L))
  expect_identical(which.min(fast), which.min(do.call(legacy_stability_scores, a)))
})

test_that("small maps and states without probes retain the dense scoring path", {
  for (size in c(3L, 20L)) {
    a <- speed_stability_state(knots = size)
    if (size == 20L) a$state$probes <- NULL
    expected <- do.call(legacy_stability_scores, a)
    local({
      local_mocked_bindings(.stageB_stability_basis = function(...) stop("probe basis rebuilt"),
        .package = "GeDS")
      expect_identical(do.call(GeDS:::.stageB_stability_scores, a), expected)
    })
  }
})

test_that("early probe exits preserve the exhaustive deletion ranking", {
  a <- speed_stability_state(weighted = TRUE, knots = 6L)
  probes <- a$state$probes
  choices <- a$choices[c(1:5, 1L)]
  candidates <- lapply(choices, function(ch) {
    factor <- GeDS:::.stageB_stability_factor(a$ctx, ch$internal)
    GeDS:::.stageB_stability_evaluate(a$ctx, ch$internal, probes, factor, a$bound)
  })
  null <- vapply(candidates, `[[`, numeric(1), "null_ratio")
  ratios <- vapply(candidates, `[[`, numeric(1), "ratio")
  expected <- order(ifelse(null > 1e-8, null, 0), ratios)[1L]
  expect_identical(GeDS:::.stageB_stability_choose(a$ctx, choices, probes, a$bound), expected)

  # A duplicated best choice should retain the earlier position.
  tied <- rep(choices[expected], 2L)
  expect_identical(GeDS:::.stageB_stability_choose(a$ctx, tied, probes, a$bound), 1L)
})

test_that("screening sparse score products retains the exhaustive winner", {
  for (weighted in c(FALSE, TRUE)) {
    a <- speed_stability_state(weighted = weighted)
    old <- do.call(GeDS:::.stageB_stability_scores, a)
    eligible <- do.call(GeDS:::.stageB_stability_screen, a)
    a$choices <- a$choices[eligible]
    new <- do.call(GeDS:::.stageB_stability_scores, a)
    expect_identical(eligible[which.min(new)], which.min(old))
  }
})

test_that("probe locations belong to each state rather than a shared cache", {
  first <- speed_stability_state(3L)
  second <- speed_stability_state(4L, knots = 18L)
  expect_false(identical(first$state$probes, second$state$probes))
  expect_identical(which.min(do.call(GeDS:::.stageB_stability_scores, first)),
    which.min(do.call(legacy_stability_scores, first)))
  reference <- GeDS:::.stageB_stability_evaluate(second$ctx, second$internal,
    second$state$probes, GeDS:::.stageB_stability_factor(second$ctx, second$internal))
  expect_null(reference$probes)
})

test_that("coefficient certificates bound weighted additive sensitivity and null maps", {
  set.seed(721)
  X <- data.frame(x = runif(100), z = runif(100))
  Z <- data.frame(group = factor(rep(1:2, 50)), v = runif(100, -2, 3))
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  k <- list(sx = c(.3, .6), sz = c(.3, .6))
  for (n in 3:4) for (w in list(rep(1, 100), rep(c(0, 1, 2, 1), 25))) {
    ctx <- GeDS:::.stageB_stability_context(X, Z, w, learners, k, n, lapply(X, range))
    fac <- GeDS:::.stageB_stability_factor(ctx, k)
    cert <- GeDS:::.stageB_stability_certificate(ctx, k, fac)
    probes <- GeDS:::.stageB_stability_probes(ctx)
    state <- GeDS:::.stageB_stability_evaluate(ctx, k, probes, fac)
    expect_lte(max(state$s), cert$sensitivity * (1 + 1e-10))
    expect_lte(state$null_ratio, cert$null_ratio + 1e-12)
  }
})

test_that("supported additive fits skip probes without changing the retained knots", {
  set.seed(722)
  X <- data.frame(x = runif(300), z = runif(300))
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  k <- list(sx = c(.25, .5, .75), sz = c(.25, .5, .75))
  for (n in 3:4) {
    args <- list(X = X, Z = NULL, weights = rep(1, 300), learners = learners,
                 internal = k, n = n)
    old <- with_mocked_bindings(do.call(GeDS:::.stageB_stabilise, args),
      .stageB_stability_certificate = function(...) list(sensitivity = Inf, null_ratio = Inf),
      .package = "GeDS")
    fast <- with_mocked_bindings(do.call(GeDS:::.stageB_stabilise, args),
      .stageB_stability_probes = function(...) stop("probe grid constructed"), .package = "GeDS")
    expect_true(fast$stability$certified)
    expect_identical(fast$internal, old$internal)
    expect_identical(fast$stability$removed, old$stability$removed)
    expect_identical(fast$stability$target_met, old$stability$target_met)
  }
})

test_that("inconclusive certificates preserve the complete deletion path", {
  x <- c(seq(0, .4, length.out = 20), seq(.8, 1, length.out = 10))
  for (n in 3:4) {
    args <- list(X = data.frame(x = x), Z = NULL, weights = rep(1, 30),
      learners = list(sx = list(type = "GeDS", variables = "x")),
      internal = list(sx = c(.45, .45001, .45002, .45003, .75)), n = n)
    fast <- do.call(GeDS:::.stageB_stabilise, args)
    old <- with_mocked_bindings(do.call(GeDS:::.stageB_stabilise, args),
      .stageB_stability_certificate = function(...) list(sensitivity = Inf, null_ratio = Inf),
      .package = "GeDS")
    expect_gt(fast$stability$removed, 0)
    expect_identical(fast, old)
  }
})

test_that("tensor certificates bound the joint response map", {
  X <- expand.grid(x = seq(0, 1, length.out = 16), z = seq(0, 1, length.out = 16))
  learners <- list(surface = list(type = "GeDS", variables = c("x", "z")))
  k <- list(surface = list(x = .5, z = .5))
  for (n in 3:4) for (w in list(rep(1, nrow(X)), rep(c(0, 1, 2, 1), 64))) {
    ctx <- GeDS:::.stageB_stability_context(X, NULL, w, learners, k, n, lapply(X, range))
    fac <- GeDS:::.stageB_stability_factor(ctx, k)
    cert <- GeDS:::.stageB_stability_certificate(ctx, k, fac)
    grid <- expand.grid(x = seq(0, 1, length.out = 51), z = seq(0, 1, length.out = 51))
    state <- GeDS:::.stageB_stability_evaluate(ctx, k, grid, fac)
    expect_lte(max(state$s), cert$sensitivity * (1 + 1e-10))
    expect_lte(state$null_ratio, cert$null_ratio + 1e-12)
  }
})

test_that("probe-subset screening keeps the exhaustive winner for joint learners", {
  set.seed(311)
  X <- data.frame(x = runif(70), z = runif(70), u = runif(70))
  learners <- list(surface = list(type = "GeDS", variables = c("x", "z")),
                   su = list(type = "GeDS", variables = "u"))
  internal <- list(surface = list(x = c(0.2, 0.35, 0.5, 0.8), z = c(0.15, 0.3, 0.6, 0.7)),
                   su = c(0.25, 0.5, 0.75))
  for (weighted in c(FALSE, TRUE)) {
    w <- if (weighted) rep(c(0.5, 1, 2), length.out = 70) else rep(1, 70)
    ctx <- GeDS:::.stageB_stability_context(X, NULL, w, learners, internal, 4L, lapply(X, range))
    probes <- GeDS:::.stageB_stability_probes(ctx)
    expect_gt(nrow(probes), 2048L)
    empty <- list(surface = list(x = NULL, z = NULL), su = NULL)
    reference <- GeDS:::.stageB_stability_evaluate(ctx, empty, probes,
      GeDS:::.stageB_stability_factor(ctx, empty))
    bound <- reference$s + 10
    state <- GeDS:::.stageB_stability_evaluate(ctx, ctx$internal, probes,
      GeDS:::.stageB_stability_factor(ctx, ctx$internal), bound, TRUE)
    choices <- GeDS:::.stageB_stability_candidates(ctx, ctx$internal)
    exhaustive <- GeDS:::.stageB_stability_choose(ctx, choices, probes, bound)
    for (size in c(8L, 64L, 512L)) {
      expect_identical(GeDS:::.stageB_stability_choose(ctx, choices, probes, bound, state,
                                                       subset_size = size), exhaustive)
    }
    # Duplicated winners keep the earlier position, as in the exhaustive search.
    tied <- c(choices[-exhaustive], choices[exhaustive], choices[exhaustive])
    expect_identical(GeDS:::.stageB_stability_choose(ctx, tied, probes, bound, state, 64L),
                     GeDS:::.stageB_stability_choose(ctx, tied, probes, bound))
  }
})
