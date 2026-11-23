test_that("binomial GD draws proportions using the original trial counts", {
  set.seed(312)
  x <- seq(-1, 1, length.out = 36)
  trials <- rep(c(2, 3, 5, 7), length.out = length(x))
  y <- round(trials * plogis(x)) / trials
  family <- binomial()
  mu <- glm.fit(cbind(1, x), y, weights = trials, family = family)$fitted.values
  expected <- local({
    restore <- .geds_bootstrap_rng()
    on.exit(restore())
    replicate(10L, rbinom(length(y), trials, mu) / trials, simplify = FALSE)
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
                               path = loss, weights = trials, draws = 10L,
                               bootstrap_path = bootstrap)
  expect_identical(observed, expected)
  expect_identical(fit$reason, "entry_screen")
  expect_identical(.Random.seed, seed)
})

test_that("standalone binomial GD accepts Bernoulli data and grouped proportions", {
  set.seed(320)
  x <- seq(-1, 1, length.out = 90)
  p <- plogis(1.8 * cos(pi * x))
  for (trials in list(rep(1, length(x)), rep(c(2, 3, 5), length.out = length(x)))) {
    dat <- data.frame(x, trials, y = rbinom(length(x), trials, p) / trials)
    seed <- .Random.seed
    fit <- suppressWarnings(GGeDS(y ~ f(x), data = dat, weights = trials,
      family = binomial(), max.intknots = 8L))
    expect_identical(fit$stopinfo$rule, "GD")
    expect_identical(.Random.seed, seed)
    for (order in 2:4) {
      pred <- predict(fit, newdata = dat, n = order)
      expect_true(all(is.finite(pred) & pred >= 0 & pred <= 1))
    }
  }
})

test_that("quasi families and non-trial weights stay outside GD", {
  x <- matrix(seq_len(8L))
  y <- rep(c(0.5, 1), 4L)
  rule <- function(family, weights, omitted = TRUE) {
    .gamma_dynamic_rule("GD", omitted, "SR", x, y, NULL, rep(0, 8), weights, family)
  }
  expect_identical(rule("binomial", rep(2, 8)), "GD")
  expect_identical(rule("binomial", rep(3, 8)), "SR")
  expect_identical(rule("binomial", rep(1.5, 8)), "SR")
  expect_identical(rule("inverse.gaussian", rep(2, 8)), "SR")
  for (family in c("quasi", "quasibinomial", "quasipoisson")) {
    expect_identical(rule(family, rep(1, 8)), "SR")
    expect_error(rule(family, rep(1, 8), omitted = FALSE), "'GD' requires")
  }
})

test_that("inverse Gaussian GD and SR retry only a failed zero-knot start", {
  set.seed(123)
  d <- data.frame(x = seq(0, 1, length.out = 100))
  family <- inverse.gaussian()
  d$y <- .gd_draw_inverse_gaussian(family$linkinv(1 + .3 * sin(2 * pi * d$x)), .05)
  native <- SplineReg_GLM
  calls <- list()
  fail <- FALSE
  local_mocked_bindings(
    SplineReg_GLM = function(..., Y, weights, InterKnots, n, mustart = NULL, fast = FALSE) {
      if (fast && n == 2L && !length(InterKnots) && identical(as.numeric(Y), d$y)) {
        calls <<- c(calls, list(mustart))
        if (fail && is.null(mustart)) {
          stop("no valid set of coefficients has been found: please supply starting values")
        }
      }
      native(..., Y = Y, weights = weights, InterKnots = InterKnots,
             n = n, mustart = mustart, fast = fast)
    })
  for (rule in c("GD", "SR")) {
    w <- if (rule == "GD") rep(1, nrow(d)) else rep(c(.5, 2), length.out = nrow(d))
    calls <- list()
    fail <- FALSE
    invisible(suppressWarnings(GGeDS(y ~ f(x), data = d, family = family, weights = w,
                                     stoptype = rule, max.intknots = 8L)))
    expect_gt(length(calls), 0L)
    expect_true(all(vapply(calls, is.null, logical(1))))

    calls <- list()
    fail <- TRUE
    fit <- suppressWarnings(GGeDS(y ~ f(x), data = d, family = family, weights = w,
                                  stoptype = rule, max.intknots = 8L))
    expect_gt(length(calls), 1L)
    expect_null(calls[[1L]])
    expect_equal(calls[[2L]], rep(weighted.mean(d$y, w), nrow(d)))
    if (rule == "GD") expect_identical(fit$stopinfo$rule, "GD")
    expect_true(all(is.finite(predict(fit, n = 2L))))
  }
})

test_that("inverse Gaussian GD retains the prefix after a later IRLS step failure", {
  x <- seq(-1, 1, length.out = 40)
  y <- exp(0.2 + 0.1 * x + 0.3 * cos(pi * x) + 0.1 * sin(33 * x))
  family <- inverse.gaussian("log")
  native <- SplineReg_GLM
  max_zero <- GenUnivariateFitter(x, y, family = family, max.intknots = 0L,
                                 stoptype = ".gamma_dynamic_path")
  # Failures are injected through SplineReg_GLM(), which only the R loop calls.
  local_mocked_bindings(.gd_glm_fast_codes = function(family) NULL)
  local_mocked_bindings(SplineReg_GLM = function(..., InterKnots) {
    if (length(InterKnots)) stop("inner loop 1; cannot correct step size")
    native(..., InterKnots = InterKnots)
  })
  prefix <- GenUnivariateFitter(x, y, family = family, max.intknots = 6L,
                               stoptype = ".gamma_dynamic_path")
  expect_length(prefix, 1L)
  expect_identical(prefix, max_zero)
  expect_error(GenUnivariateFitter(x, y, family = family, max.intknots = 1L,
                                  stoptype = "SR"), "cannot correct step size")
})

test_that("inverse Gaussian GD calibrates on a previously failing strong-curvature case", {
  set.seed(264506)
  x <- runif(400, -1, 1)
  invisible(runif(400))
  mu <- exp(0.15 + 0.5 * x + 0.9 * cos(pi * x))
  dat <- data.frame(x, y = .gd_draw_inverse_gaussian(mu, dispersion = 0.5))
  seed <- .Random.seed
  fit <- suppressWarnings(GGeDS(y ~ f(x), family = inverse.gaussian("log"), data = dat,
                                max.intknots = 12L, gd_draws = 100L))
  expect_identical(fit$stopinfo$rule, "GD")
  expect_gt(length(knots(fit, n = 2L, options = "internal")), 0L)
  expect_identical(.Random.seed, seed)
  for (order in 2:4) {
    pred <- predict(fit, newdata = dat, n = order)
    expect_true(all(is.finite(pred) & pred > 0))
  }
})

test_that("binomial GD cuts the observed path at its first non-improving step", {
  set.seed(41)
  x <- seq(0, 1, length.out = 60)
  y <- rbinom(length(x), 1, plogis(-1 + 2 * x))
  # The observed loss improves for five knots and then breaks down.
  observed <- c(100, 60, 40, 30, 25, 22, 500, 480, 470, 460, 455)
  path <- function(y, k) utils::head(observed, k + 1L)
  # Null-like bootstrap paths with small, steady improvements.
  bootstrap <- function(y, k, b) 100 - 0.01 * seq.int(0L, k)
  select <- function(prefix) {
    .gamma_dynamic_select(x, y, binomial(), q = 2L, lower = 0L, upper = 10L, path = path,
                          allow_short_paths = TRUE, draws = 10L,
                          bootstrap_path = bootstrap, prefix_observed = prefix)
  }
  expect_identical(select(TRUE), list(knots = 5L, reason = "path_end"))
  # Without the cut, the invalid gain blocks every later check and the full path is kept.
  expect_identical(select(FALSE), list(knots = 10L, reason = "path_end"))
})

test_that("Bernoulli GD keeps valid path prefixes instead of falling back at the default cap", {
  skip_on_cran()
  set.seed(5)
  n <- 150
  dat <- data.frame(x = runif(n))
  dat$y <- rbinom(n, 1, plogis(-1 + 2 * dat$x))
  fit <- expect_no_warning(
    withCallingHandlers(GGeDS(y ~ f(x), data = dat, family = binomial()),
      warning = function(w) if (!grepl("calibration failed", conditionMessage(w))) invokeRestart("muffleWarning")))
  expect_identical(fit$stopinfo$rule, "GD")
  expect_identical(fit$stopinfo$reason, "entry_screen")
  expect_length(fit$linear.intknots, 0L)
})

irls_failure_data <- function() {
  set.seed(1009)
  x <- runif(200)
  data.frame(x = x, trials = 10, y = rbinom(200, 10, plogis(2 * sin(2 * pi * x))) / 10)
}

fail_from_third_knot <- function(message) {
  native <- SplineReg_GLM
  function(..., InterKnots, n) {
    if (n == 2 && length(InterKnots) >= 3L) stop(message)
    native(..., InterKnots = InterKnots, n = n)
  }
}

test_that("binomial GD keeps the fitted prefix when IRLS fails at a later knot", {
  dat <- irls_failure_data()
  # Failures are injected through SplineReg_GLM(), which only the R loop calls.
  local_mocked_bindings(.gd_glm_fast_codes = function(family) NULL)
  local_mocked_bindings(SplineReg_GLM = fail_from_third_knot("inner loop 1; cannot correct step size"))
  warned <- character()
  fit <- withCallingHandlers(
    GGeDS(y ~ f(x), data = dat, weights = trials, family = binomial(), max.intknots = 8L),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(fit$stopinfo$rule, "GD")
  expect_false(any(grepl("calibration failed", warned)))
  expect_lte(length(fit$linear.intknots), 2L)
})

test_that("binomial GD propagates IRLS errors it does not recognise", {
  dat <- irls_failure_data()
  local_mocked_bindings(.gd_glm_fast_codes = function(family) NULL)
  local_mocked_bindings(SplineReg_GLM = fail_from_third_knot("unrelated failure"))
  expect_error(suppressWarnings(GGeDS(y ~ f(x), data = dat, weights = trials,
                                      family = binomial(), max.intknots = 8L)),
               "unrelated failure")
})

test_that("GD warns when the valid path ends below min.intknots", {
  dat <- irls_failure_data()
  # Failures are injected through SplineReg_GLM(), which only the R loop calls.
  local_mocked_bindings(.gd_glm_fast_codes = function(family) NULL)
  local_mocked_bindings(SplineReg_GLM = fail_from_third_knot("inner loop 1; cannot correct step size"))
  fit <- NULL
  expect_warning(
    withCallingHandlers(
      fit <- GGeDS(y ~ f(x), data = dat, weights = trials, family = binomial(),
                   min.intknots = 5L, max.intknots = 8L),
      warning = function(w) if (!inherits(w, "geds_min_intknots_warning")) invokeRestart("muffleWarning")),
    "below 'min.intknots' = 5", class = "geds_min_intknots_warning")
  expect_identical(fit$stopinfo$rule, "GD")
  expect_lte(length(fit$linear.intknots), 2L)
})
