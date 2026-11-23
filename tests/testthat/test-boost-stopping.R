test_that("boosting defaults to memory and starts multiple learners at the intercept", {
  d <- data.frame(x = seq(0, 1, length.out = 25L))
  d$z <- cos(3 * d$x)
  d$y <- sin(5 * d$x) + d$z
  fit <- function(formula, ...) suppressWarnings(NGeDSboost(
    formula, data = d, int.knots_boost = 1L, max_iterations = 1L,
    higher_order = FALSE, ...
  ))
  expect_true(fit(y ~ f(x, z))$args$initial_learner)
  multi <- fit(y ~ f(x) + z)
  expect_false(multi$args$initial_learner)
  expect_named(multi$args$learner_knots_frozen_at, "f(x)")
  explicit <- fit(y ~ f(x) + z, initial_learner = TRUE, boosting_with_memory = FALSE)
  expect_true(explicit$args$initial_learner)
  expect_null(explicit$args$learner_knots_frozen_at)
})

test_that("the RD override bypasses univariate GD calibration", {
  d <- data.frame(x = seq(0, 1, length.out = 40L))
  d$y <- sin(8 * d$x)
  local_mocked_bindings(.gamma_dynamic_fit = function(...) {
    stop("Unexpected GD calibration")
  }, .package = "GeDS")
  fit <- suppressWarnings(NGeDSboost(y ~ f(x), data = d,
    stoptype = "RD", int.knots_boost = 3L, max_iterations = 2L))
  expect_identical(fit$args$stoptype, "RD")
  expect_true("f(x)" %in% names(fit$final_model$base_learners))
  expect_true(all(is.finite(predict(fit, newdata = d, n = 2L))))
})

test_that("deviance stopping is unaffected by additive constants in family risk", {
  dat <- data.frame(x = seq(-1, 1, length.out = 40))
  dat$y <- 2 + 10 * dat$x + sin(5 * dat$x)
  fam <- mboost::Gaussian()
  fam@risk <- function(y, f, w = 1) sum(w * (y - f)^2) + 1e6

  invisible(capture.output(fit <- NGeDSboost(
    y ~ x, data = dat, family = fam, initial_learner = FALSE, shrinkage = 0.1,
    q_boost = 2L, max_iterations = 4L, higher_order = FALSE
  )))
  invisible(capture.output(reference <- NGeDSboost(
    y ~ x, data = dat, initial_learner = FALSE, shrinkage = 0.1,
    q_boost = 2L, max_iterations = 4L, higher_order = FALSE
  )))

  expect_equal((length(fit$models) - 1L), 4L)
  expect_equal(N.boost.iter(fit), N.boost.iter(reference))
  expect_equal(predict(fit, n = 2), predict(reference, n = 2))
})

test_that("Poisson stopping compares deviance with the initial model deviance", {
  set.seed(123)
  x <- sort(runif(80, -1, 1))
  dat <- data.frame(x = x, y = rpois(80, exp(2 + sin(3 * x))))

  invisible(capture.output(fit <- suppressWarnings(NGeDSboost(
    y ~ f(x), data = dat, family = mboost::Poisson(),
    shrinkage = 0.005, max_iterations = 12L, higher_order = FALSE,
    boosting_with_memory = FALSE
  ))))
  dev <- vapply(fit$models, function(model) {
    sum(poisson()$dev.resids(dat$y, model$Y_hat, rep(1, nrow(dat))))
  }, numeric(1))
  ratios <- dev[3:length(dev)] / dev[1:(length(dev) - 2L)]
  if ((length(fit$models) - 1L) < 12L && is.null(fit$args$knot_cap_stop)) {
    expect_gte(tail(ratios, 1), 0.99)
  }
  if (!is.null(fit$args$knot_cap_stop)) {
    expect_equal((length(fit$models) - 1L), fit$args$knot_cap_stop$iteration)
  }
  if (length(ratios) > 1L) expect_true(all(head(ratios, -1L) < 0.99))
  expect_true(all(is.finite(dev)))
})

test_that("negative likelihood risk is not treated as zero deviance", {
  dat <- data.frame(x = seq(-1, 1, length.out = 40))
  dat$y <- 2 + dat$x + sin(5 * dat$x)
  fam <- mboost::Gaussian()
  fam@risk <- function(y, f, w = 1) sum(w * (y - f)^2) - 1000

  invisible(capture.output(fit <- NGeDSboost(
    y ~ x, data = dat, family = fam, initial_learner = FALSE,
    shrinkage = 0.1, phi_boost_exit = NULL,
    max_iterations = 3L, higher_order = FALSE
  )))

  expect_equal((length(fit$models) - 1L), 3L)
  expect_length(fit$models, 4L)
  expect_true(all(is.finite(predict(fit, n = 2))))
})

test_that("the model attaining near-zero deviance is retained", {
  dat <- data.frame(x = seq(-1, 1, length.out = 40))
  dat$y <- 2 + dat$x

  invisible(capture.output(fit <- NGeDSboost(
    y ~ x, data = dat, initial_learner = FALSE,
    shrinkage = 1, phi_boost_exit = NULL,
    max_iterations = 3L, higher_order = FALSE
  )))

  expect_equal((length(fit$models) - 1L), 1L)
  expect_length(fit$models, 2L)
  expect_equal(as.numeric(predict(fit, n = 2)), dat$y, tolerance = 1e-10)
})

test_that("Gamma boosting can continue beyond a negative likelihood risk", {
  set.seed(123)
  x <- sort(runif(500, -2, 2))
  mu <- exp(40 * x / (1 + 100 * x^2) + 4)
  dat <- data.frame(x = x, y = rgamma(500, shape = mu, rate = 0.1))

  invisible(capture.output(fit <- suppressWarnings(NGeDSboost(
    y ~ f(x), data = dat, family = mboost::GammaReg(),
    beta = 0.5, phi = 0.99, q = 2L, shrinkage = 0.1,
    phi_boost_exit = NULL, max_iterations = 3L, higher_order = FALSE
  ))))

  expect_equal((length(fit$models) - 1L), 3L)
  expect_length(fit$models, 4L)
  expect_true(all(is.finite(predict(fit, n = 2))))
})
