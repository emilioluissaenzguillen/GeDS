test_that("linear GD basis assembly matches splineDesign exactly", {
  internal <- c(0.2, 0.5, 0.500000000000001, 0.8)
  knots <- c(0, 0, internal, 1, 1)
  sites <- list(c(0, 0, 0.2, 0.2, 0.5, 0.8, 1, 1),
                seq(0, 1, length.out = 64),
                c(-0.1, 0, 0.4, 0.6, 1, 1.1))
  for (x in sites) {
    for (sparse in c(FALSE, TRUE)) {
      expected <- splines::splineDesign(knots, x, ord = 2L, derivs = 0L,
                                        outer.ok = TRUE, sparse = sparse)
      actual <- .gd_linear_basis(knots, x, ord = 2L, derivs = 0L,
                                outer.ok = TRUE, sparse = sparse)
      expect_identical(actual, expected)
    }
  }
})

test_that("standalone fits use Gamma dynamic with the family beta defaults", {
  set.seed(121)
  x <- seq(0, 1, length.out = 100)
  mu <- exp(sin(9 * x) + 2)
  cases <- list(
    list(model = NGeDS, family = NULL, y = rnorm(100, sin(9 * x), .15), beta = .5),
    list(model = GGeDS, family = gaussian(), y = rnorm(100, sin(9 * x), .15), beta = .5),
    list(model = GGeDS, family = poisson(), y = rpois(100, mu), beta = .2),
    list(model = GGeDS, family = Gamma("log"), y = rgamma(100, 25, scale = mu / 25), beta = .3)
  )
  for (case in cases) {
    args <- list(formula = y ~ f(x), data = data.frame(x, y = case$y),
                 max.intknots = 12)
    if (!is.null(case$family)) args$family <- case$family
    seed <- .Random.seed
    fit <- suppressWarnings(do.call(case$model, args))
    expect_identical(.Random.seed, seed)
    expect_identical(fit$stopinfo$rule, "GD")
    expect_equal(fit$args$beta, case$beta)
    expect_equal(fit$args$q, 2L)
    expect_equal(fit$args$max.intknots, 12L)
    for (degree in 2:4) {
      expect_true(all(is.finite(predict(fit, n = degree))))
      expect_length(predict(fit, newdata = data.frame(x = x[1:5]), n = degree), 5L)
    }

    seed <- .Random.seed
    rd <- suppressWarnings(do.call(case$model, c(args, list(stoptype = "RD", phi = .9))))
    expect_identical(.Random.seed, seed)
    expect_equal(rd$args$phi, .9)
    expect_true(all(is.finite(predict(rd, n = 2))))
  }
})

test_that("Poisson GD fixes the scale without changing its default bootstrap draws", {
  set.seed(319)
  x <- seq(-1, 1, length.out = 40)
  y <- rpois(length(x), exp(0.5 + x))
  family <- poisson()
  mu <- glm.fit(cbind(1, x), y, family = family)$fitted.values
  draws <- formals(GGeDS)$gd_draws
  expect_identical(draws, 10L)
  expected <- local({
    restore <- .geds_bootstrap_rng()
    on.exit(restore())
    replicate(draws, rpois(length(y), mu), simplify = FALSE)
  })
  observed <- list()
  local_mocked_bindings(.gamma_dynamic_neighbours = function(...) stop("Unexpected scale estimate"))
  loss <- function(y, k) 100 - seq.int(0L, k)
  bootstrap <- function(y, k, b) {
    observed[[b]] <<- y
    loss(y, k)
  }
  seed <- .Random.seed
  fit <- .gamma_dynamic_select(x, y, family, q = 2L, lower = 0L, upper = 6L,
                               path = loss, draws = draws, bootstrap_path = bootstrap)
  expect_identical(fit$reason, "entry_screen")
  expect_identical(observed, expected)
  expect_identical(.Random.seed, seed)
})

test_that("NGeDSgam uses GD for univariate smoothers unless overridden", {
  x <- seq(0, 1, length.out = 60L)
  dat <- data.frame(x, y = sin(8 * x) + 0.15 * cos(60 * x))
  calls <- 0L
  native <- .gamma_dynamic_fit
  local_mocked_bindings(.gamma_dynamic_fit = function(..., allow_short_paths) {
    calls <<- calls + 1L
    expect_true(allow_short_paths)
    native(..., allow_short_paths = allow_short_paths)
  })
  fit <- function(...) suppressWarnings(NGeDSgam(y ~ f(x), data = dat,
    internal_knots = 8L, higher_order = FALSE, ...))
  default <- fit()
  expect_gt(calls, 0L)
  expect_equal(predict(default, n = 2L), predict(fit(stoptype = "GD"), n = 2L))
  calls <- 0L
  invisible(fit(stoptype = "RD"))
  expect_identical(calls, 0L)
})

test_that("GAM GD retains valid shorter paths on the Kyphosis regression split", {
  skip_if_not_installed("rpart")
  indices <- c(34L, 15L, 59L, 11L, 5L, 69L, 16L, 23L, 38L, 33L,
    25L, 36L, 66L, 64L, 12L, 19L, 26L, 40L, 78L, 60L, 49L, 65L,
    14L, 1L, 41L, 57L, 31L, 17L, 72L, 61L, 70L, 29L, 44L, 75L,
    13L, 10L, 27L, 80L, 48L, 53L, 22L, 39L, 55L, 3L)
  dat <- rpart::kyphosis[indices, ]
  native <- .gamma_dynamic_fit
  strict <- TRUE
  rules <- character()
  local_mocked_bindings(.gamma_dynamic_fit = function(fitter, args, family, legacy,
                                                     allow_short_paths, ...) {
    expect_true(allow_short_paths)
    result <- native(fitter, args, family, legacy,
      allow_short_paths = if (strict) FALSE else allow_short_paths, ...)
    rules <<- c(rules, result$stopinfo$rule)
    result
  })
  fit <- function() {
    invisible(capture.output(result <- suppressWarnings(NGeDSgam(
      Kyphosis ~ f(Age) + f(Start), data = dat, family = binomial,
      higher_order = FALSE))))
    result
  }
  set.seed(123)
  state <- .Random.seed
  old <- fit()
  expect_true(any(rules == "RD"))
  expect_identical(.Random.seed, state)
  strict <- FALSE
  rules <- character()
  current <- fit()
  expect_gt(length(rules), 0L)
  expect_true(all(rules == "GD"))
  expect_identical(.Random.seed, state)
  expect_true(all(is.finite(predict(old, n = 2L))))
  expect_true(all(is.finite(predict(current, n = 2L))))
})

test_that("Gamma dynamic boosting is independent of the caller's seed", {
  x <- seq(0, 1, length.out = 60)
  dat <- data.frame(x, y = sin(10 * x) + cos(101 * x) / 10)
  fit <- function(seed) {
    set.seed(seed)
    state <- .Random.seed
    result <- suppressWarnings(NGeDSboost(y ~ f(x), data = dat,
      max_iterations = 3L, higher_order = FALSE))
    expect_identical(.Random.seed, state)
    result
  }
  first <- fit(1)
  second <- fit(42)
  expect_identical(knots(first, n = 2), knots(second, n = 2))
  expect_identical(predict(first, n = 2), predict(second, n = 2))
})

test_that("resumed bootstrap paths match uninterrupted losses", {
  set.seed(831)
  x <- seq(0, 1, length.out = 500)
  mu <- exp(2 + sin(10 * x))
  for (family in list(gaussian(), poisson(), Gamma("log"))) {
    fitter <- if (family$family == "gaussian") UnivariateFitter else GenUnivariateFitter
    y <- switch(family$family, gaussian = rnorm(length(x), mu),
                poisson = rpois(length(x), mu), Gamma = rgamma(length(x), 10, 10 / mu))
    args <- list(X = x, Y = y, stoptype = ".gamma_dynamic_path", max.intknots = 80L)
    if (family$family != "gaussian") args$family <- family
    full <- suppressWarnings(do.call(fitter, args))
    expect_length(full, 81L)
    args$stoptype <- structure(new.env(parent = emptyenv()), class = "geds_gd_path")
    args$max.intknots <- 40L
    prefix <- suppressWarnings(do.call(fitter, args))
    expect_identical(prefix, full[1:41])
    args$max.intknots <- 80L
    expect_identical(suppressWarnings(do.call(fitter, args)), full)
  }
})

test_that("weighted standalone Gaussian fits default to GD", {
  set.seed(123)
  x <- seq(0, 1, length.out = 80)
  w <- rep(c(0.5, 2), length.out = length(x))
  d <- data.frame(x, y = sin(9 * x) + rnorm(length(x), sd = 0.15 / sqrt(w)), w)
  padded <- rbind(data.frame(x = -1, y = 100, w = 0), d,
                  data.frame(x = 2, y = -100, w = 0))
  for (fitter in list(NGeDS, GGeDS)) {
    fit <- suppressWarnings(fitter(y ~ f(x), data = padded, weights = w,
                                   max.intknots = 8))
    gd <- suppressWarnings(fitter(y ~ f(x), data = d, weights = w,
                                  max.intknots = 8, stoptype = "GD"))
    expect_identical(fit$stopinfo$rule, "GD")
    expect_identical(gd$stopinfo$rule, "GD")
    expect_equal(knots(fit, n = 2), knots(gd, n = 2))
    expect_equal(as.numeric(predict(fit, newdata = d, n = 2)),
                 as.numeric(predict(gd, newdata = d, n = 2)), tolerance = 1e-9)
  }
})

test_that("Gamma GD supports precision weights and ignores zero-weight rows", {
  set.seed(123)
  x <- seq(0, 1, length.out = 100)
  w <- rep(c(0.5, 2), length.out = length(x))
  for (family in list(Gamma("log"), Gamma("inverse"))) {
    off <- 0.1 * cos(4 * x)
    mu <- family$linkinv(1 + 0.4 * sin(9 * x) + off)
    d <- data.frame(x, y = rgamma(length(x), shape = 20 * w, scale = mu / (20 * w)), w, off)
    padded <- rbind(data.frame(x = -1, y = 100, w = 0, off = 0), d,
                    data.frame(x = 2, y = 200, w = 0, off = 0))
    formula <- y ~ f(x) + offset(off)
    fit <- suppressWarnings(GGeDS(formula, data = padded, family = family, weights = w,
                                  max.intknots = 10L))
    gd <- suppressWarnings(GGeDS(formula, data = d, family = family, weights = w,
                                 stoptype = "GD", max.intknots = 10L))
    expect_identical(fit$stopinfo$rule, "GD")
    expect_identical(gd$stopinfo$rule, "GD")
    expect_gt(gd$Nintknots, 0)
    for (n in 2:4) {
      expect_equal(knots(fit, n = n), knots(gd, n = n))
      expect_equal(predict(fit, newdata = d, n = n), predict(gd, newdata = d, n = n),
                   tolerance = 1e-7)
    }
  }
})

test_that("Gaussian GD with a fixed offset matches subtracting it under identical draws", {
  set.seed(123)
  x <- seq(0, 1, length.out = 100)
  off <- 5 * cos(19 * x)
  y <- sin(9 * x) + rnorm(length(x), sd = 0.12)
  d <- data.frame(x, y, shifted = y + off, off)
  for (fitter in list(NGeDS, GGeDS)) {
    for (weights in list(rep(1, length(x)), rep(c(0.5, 2), length.out = length(x)))) {
      reference <- suppressWarnings(fitter(y ~ f(x), data = d, weights = weights,
                                            max.intknots = 16))
      fit <- suppressWarnings(fitter(shifted ~ f(x) + offset(off), data = d,
                                      weights = weights, max.intknots = 16))
      expect_identical(reference$stopinfo$rule, "GD")
      expect_identical(fit$stopinfo$rule, "GD")
      expect_identical(fit$stopinfo$reason, reference$stopinfo$reason)
      expect_gt(reference$Nintknots, 0)
      for (n in 2:4) {
        expect_equal(knots(fit, n = n), knots(reference, n = n), tolerance = 1e-8)
        expect_equal(as.numeric(predict(fit, newdata = d, n = n)) - off,
                     as.numeric(predict(reference, newdata = d, n = n)), tolerance = 1e-8)
      }
    }
  }
})

test_that("Poisson GD draws responses from the zero-knot model with its offset", {
  set.seed(123)
  x <- seq(0, 1, length.out = 80)
  off <- cos(15 * x)
  y <- rpois(length(x), exp(1 + x + off))
  loss <- c(30, 26, 24, 23, 22, 21, 20)
  calibration <- new.env()
  calibration$draws <- list()
  choice <- .gamma_dynamic_select(
    x, y, poisson(), q = 2L, lower = 0L, upper = 6L,
    path = function(y, k) loss, offset = off, draws = 10L,
    bootstrap_path = function(y, k, b) {
      calibration$draws[[b]] <- y
      loss
    }
  )
  reference <- glm.fit(cbind(1, x), y, offset = off, family = poisson())
  expected <- local({
    restore <- .geds_bootstrap_rng()
    on.exit(restore())
    replicate(10L, rpois(length(x), reference$fitted.values), simplify = FALSE)
  })
  expect_identical(choice$reason, "entry_screen")
  expect_identical(calibration$draws, expected)
  d <- data.frame(x, y, off)
  fit <- suppressWarnings(GGeDS(y ~ f(x) + offset(off), family = poisson(),
                                data = d, max.intknots = 6))
  expect_identical(fit$stopinfo$rule, "GD")
})

test_that("Gamma dynamic is translation invariant with matched bootstrap draws", {
  set.seed(123)
  x <- seq(0, 1, length.out = 100)
  d <- data.frame(x = x, y = sin(10 * x) + rnorm(100, sd = .1))
  fit <- function(shift) {
    set.seed(456)
    suppressWarnings(NGeDS(y ~ f(x), data = transform(d, y = y + shift),
                           beta = .6, q = 3))
  }
  reference <- fit(0)
  expect_identical(reference$stopinfo$rule, "GD")
  expect_gt(reference$Nintknots, 0)
  for (shift in c(-mean(d$y), 1e6)) {
    shifted <- fit(shift)
    expect_identical(shifted$stopinfo$rule, "GD")
    for (n in 2:4) {
      expect_equal(knots(shifted, n = n, options = "internal"),
                   knots(reference, n = n, options = "internal"),
                   tolerance = 1e-7)
      expect_lt(max(abs(predict(shifted, newdata = d, n = n) - shift -
                          predict(reference, newdata = d, n = n))), 1e-6)
    }
  }
})
test_that("bivariate GD calibrates with short paths and keeps inherited knots", {
  set.seed(41)
  d <- expand.grid(x = seq(0, 1, length.out = 8L), z = seq(0, 1, length.out = 8L))
  d$y <- with(d, x + z + x*z) + rnorm(nrow(d), sd = .1)
  path_lengths <- integer()
  native <- GeDS:::BivariateFitter
  local_mocked_bindings(BivariateFitter = function(..., stoptype) {
    result <- native(..., stoptype = stoptype)
    if (identical(stoptype, ".gamma_dynamic_path")) {
      path_lengths <<- c(path_lengths, length(result))
    }
    result
  }, .package = "GeDS")
  fit <- NGeDS(y ~ f(x, z), data = d, stoptype = "GD", higher_order = FALSE)
  expect_identical(fit$stopinfo$rule, "GD")
  expect_equal(fit$stopinfo$budget, 20L)
  expect_true(all(is.finite(predict(fit, newdata = d, n = 2L))))
  expect_gte(length(path_lengths), 101L)
  expect_true(all(path_lengths < 21L))

  before_fixed <- length(path_lengths)
  fixed <- NGeDS(y ~ f(x, z), data = d, stoptype = "GD", higher_order = FALSE,
    intknots_init = list(ikX = .4, ikY = .6), max.intknots = 2L)
  expect_identical(fixed$stopinfo$reason, "fixed_knots")
  expect_equal(unname(unlist(fixed$linear.intknots)), c(.4, .6))
  expect_length(path_lengths, before_fixed)

  invisible(capture.output(boost <- NGeDSboost(
    y ~ f(x, z) + x, data = d,
    initial_learner = FALSE, boosting_with_memory = TRUE, higher_order = FALSE,
    max_iterations = 2L, int.knots_boost = 3L
  )))
  expect_false("bivariate_stoptype" %in% names(formals(NGeDSboost)))
  expect_length(path_lengths, before_fixed)
  expect_true(all(is.finite(predict(boost, newdata = d, n = 2L))))
})

test_that("boosting calibrates univariate learners from the inherited basis", {
  set.seed(42)
  d <- data.frame(x = seq(0, 1, length.out = 60L))
  d$y <- sin(8 * d$x) + rnorm(nrow(d), sd = .1)
  selected <- list()
  native <- GeDS:::.gamma_dynamic_fit
  local_mocked_bindings(.gamma_dynamic_fit = function(fitter, args, family, legacy, allow_short_paths, ...) {
    expect_true(allow_short_paths)
    fit <- native(fitter, args, family, legacy, allow_short_paths, ...)
    selected[[length(selected) + 1L]] <<- fit$stopinfo$rule
    fit
  }, .package = "GeDS")
  invisible(capture.output(fit <- NGeDSboost(y ~ f(x), data = d,
    boosting_with_memory = TRUE, initial_learner = TRUE, int.knots_init = 1L,
    int.knots_boost = 1L, max_iterations = 2L, higher_order = FALSE)))
  expect_true("GD" %in% unlist(selected))
  expect_true(all(is.finite(predict(fit, newdata = d, n = 2L))))

  inherited <- .45
  args <- list(X = d$x, Y = d$y, weights = rep(1, nrow(d)), extr = c(0, 1),
               intknots_init = inherited, min.intknots = 0L, max.intknots = 2L,
               q = 2L, higher_order = FALSE, show.iters = FALSE)
  inherited_fit <- suppressWarnings(native(GeDS:::UnivariateFitter, args, gaussian(), "RD"))
  expect_identical(inherited_fit$stopinfo$rule, "GD")
  expect_true(inherited %in% inherited_fit$linear.intknots)
})

test_that("boosting calibration accepts shorter paths with positive decreasing RSS", {
  set.seed(43)
  x <- seq(0, 1, length.out = 30L)
  y <- rnorm(length(x))
  path <- function(response, cap) {
    if (identical(response, y)) return(c(10, 9.999, 9.998, 9.997, 9.996, 9.995))
    c(10, 5)
  }
  args <- list(x = x, y = y, family = gaussian(), q = 2L, lower = 0L, upper = 5L, path = path)
  expect_null(do.call(GeDS:::.gamma_dynamic_select, args))
  fit <- do.call(GeDS:::.gamma_dynamic_select, c(args, list(allow_short_paths = TRUE)))
  expect_identical(fit, list(knots = 0L, reason = "entry_screen"))
})

test_that("a dominated GD candidate skips bootstrap fits and preserves the random stream", {
  x <- seq(0, 1, length.out = 30L)
  y <- sin(8 * x)
  boot_calls <- 0L
  path <- function(response, cap) {
    if (!identical(response, y)) boot_calls <<- boot_calls + 1L
    c(10, 8, 7, 6)
  }
  select <- function(...) GeDS:::.gamma_dynamic_select(x, y, gaussian(), 2L,
    0L, 3L, path, allow_short_paths = TRUE, ...)
  set.seed(44)
  select()
  seed <- .Random.seed
  expect_gt(boot_calls, 0L)
  boot_calls <- 0L
  set.seed(44)
  expect_identical(select(incumbent = 5), list(knots = 0L, reason = "dominated"))
  expect_identical(boot_calls, 0L)
  expect_identical(.Random.seed, seed)
})


test_that("bivariate Gaussian GD draws are invariant to weight rescaling and zero-weight rows", {
  set.seed(321)
  d <- expand.grid(x = seq(0, 1, length.out = 8L), y = seq(0, 1, length.out = 7L))
  d$z <- with(d, sin(2 * pi * x) + cos(pi * y)) + rnorm(nrow(d), sd = 0.3)
  loss <- c(30, 20, 15, 13, 12, 11.5, 11)
  recorded <- list()
  local_mocked_bindings(BivariateFitter = function(..., Z, stoptype) {
    if (identical(stoptype, ".gamma_dynamic_path")) {
      recorded[[length(recorded) + 1L]] <<- Z
      return(loss)
    }
    list(args = list())
  })
  draws_for <- function(d, w) {
    recorded <<- list()
    .bivariate_gamma_dynamic_fit(list(
      X = d$x, Y = d$y, Z = d$z, weights = w, Xextr = c(0, 1), Yextr = c(0, 1),
      q = 2L, min.intknots = 0L, max.intknots = 6L))
    expect_length(recorded, 101L)
    recorded[-1L]
  }
  n <- nrow(d)
  reference <- draws_for(d, rep(1, n))
  expect_equal(draws_for(d, rep(3, n)), reference, tolerance = 1e-10)
  padded <- rbind(d, data.frame(x = 0.32, y = 0.41, z = 500))
  expect_equal(lapply(draws_for(padded, c(rep(1, n), 0)), head, n = n),
               reference, tolerance = 1e-10)
  w <- rep(c(0.5, 2), length.out = n)
  expect_equal(draws_for(d, 4 * w), draws_for(d, w), tolerance = 1e-10)

  # A flat response exercises the numerical floor rather than the scale estimate.
  flat <- transform(d, z = 0)
  expect_lt(max(abs(unlist(draws_for(flat, rep(3, n))) -
                    unlist(draws_for(flat, rep(1, n))))), 1e-15)
})

test_that("Gaussian GD calibration is invariant to rescaling weights and zero-weight rows", {
  set.seed(321)
  n <- 120
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.3)
  loss <- c(30, 20, 15, 13, 12, 11.5, 11)
  draws_for <- function(x, y, weights, offset = rep(0, length(y))) {
    recorded <- list()
    .gamma_dynamic_select(x, y, gaussian(), q = 2L, lower = 0L, upper = 6L,
      path = function(y, k) loss, weights = weights, offset = offset, draws = 10L,
      bootstrap_path = function(y, k, b) {
        recorded[[b]] <<- y
        loss
      })
    recorded
  }
  for (off in list(rep(0, n), 0.5 * cos(7 * x))) {
    reference <- draws_for(x, y, rep(1, n), off)
    expect_equal(draws_for(x, y, rep(3, n), off), reference, tolerance = 1e-10)
    # Zero-weight rows at the ends keep the neighbour geometry of the others.
    xz <- c(-0.5, x, 1.5)
    yz <- c(50, y, -50)
    padded <- draws_for(xz, yz, c(0, rep(1, n), 0), c(0, off, 0))
    expect_equal(lapply(padded, function(d) d[-c(1L, n + 2L)]), reference,
                 tolerance = 1e-10)
  }
  # Varying weights still use the zero-knot residuals.
  w <- rep(c(0.5, 2), length.out = n)
  expect_equal(draws_for(x, y, 4 * w), draws_for(x, y, w), tolerance = 1e-10)
})
