test_that("the shared first-knot score is the relative RSS gain", {
  x <- rep(seq(-1, 1, length.out = 30), each = 2)
  y <- x^2 + sin(2 * x)
  basis <- GeDS:::.boost_screen_basis(list(ord = seq_along(x), xs = x, extr = range(x)))
  line <- lm(y ~ x)
  knot <- GeDS:::.boost_screen_first_knot(basis, residuals(line), max(abs(y)), 0.5)
  hinge <- lm(y ~ x + I(pmax(x - knot, 0)))
  expect_equal(GeDS:::.boost_screen_gain(basis, cbind(y, x), 0.5),
               c(1 - deviance(hinge) / deviance(line), 0))
  screened <- GeDS:::.boost_screen_gain(basis, cbind(y, x), 0.5, return_line = TRUE)
  expect_identical(screened$gain, GeDS:::.boost_screen_gain(basis, cbind(y, x), 0.5))
  expect_equal(as.numeric(screened$theta[, 1]), unname(predict(line, data.frame(x = range(x)))),
               tolerance = 1e-12)
})

test_that("reusing a screened line preserves the boosting path and predictions", {
  x <- rep(seq(-1, 1, length.out = 20), each = 2)
  ord <- c(seq(40, 2, by = -2), seq(1, 39, by = 2))
  dat <- data.frame(y = factor(ifelse(x > 0 & seq_along(x) %% 3 != 0, "yes", "no")),
                    x = x, twin = x)[ord, ]
  fit <- function() {
    NGeDSboost(y ~ f(x) + f(twin), data = dat, family = mboost::Binomial(),
               shrinkage = 0.1, phi_boost_exit = 0.95, max_iterations = 6,
               high_dim_fit = TRUE)
  }
  native_fit <- GeDS:::.gamma_dynamic_fit
  native_calls <- 0L
  local_mocked_bindings(
    .boost_screen_calibrate = function(...) 1,
    .gamma_dynamic_fit = function(...) {
      native_calls <<- native_calls + 1L
      native_fit(...)
    },
    .package = "GeDS"
  )
  reused <- fit()
  expect_identical(native_calls, 0L)
  refit <- with_mocked_bindings(fit(), .boost_screen_line = function(...) NULL,
                                .package = "GeDS")
  expect_gt(native_calls, 0L)

  selected <- function(object) {
    vapply(object$models[-1L], function(m) m$best_bl$name, character(1))
  }
  expect_true(all(selected(reused) == "f(x)"))
  expect_identical(selected(reused), selected(refit))
  expect_identical(reused$iters, refit$iters)
  expect_identical(reused$final_model$model_name, refit$final_model$model_name)
  expect_equal(lengths(knots(reused, n = 2, options = "internal")),
               lengths(knots(refit, n = 2, options = "internal")))
  expect_equal(predict(reused, n = 2), predict(refit, n = 2), tolerance = 1e-12)
  newdata <- data.frame(x = c(-1.4, 0.15, 1.7), twin = c(-1.4, 0.15, 1.7))
  expect_equal(predict(reused, newdata = newdata, n = 2),
               predict(refit, newdata = newdata, n = 2), tolerance = 1e-12)
})

test_that("near-degenerate screened predictors keep the standard fit path", {
  expect_null(GeDS:::.boost_screen_line(c(0, 1), c(2, 2 + 1e-12), 2))
})

test_that("high-dimensional mode uses a 0.95 stopping default unless overridden", {
  threshold <- formals(NGeDSboost)$phi_boost_exit
  expect_identical(eval(threshold, list(high_dim_fit = FALSE)), 0.99)
  expect_identical(eval(threshold, list(high_dim_fit = TRUE)), 0.95)

  x <- seq(-1, 1, length.out = 30)
  dat <- data.frame(y = factor(ifelse(x > 0 & seq_along(x) %% 5 != 0, "yes", "no")), x = x)
  local_mocked_bindings(.boost_screen_calibrate = function(...) 1, .package = "GeDS")
  fit <- function(...) {
    NGeDSboost(y ~ f(x), data = dat, family = mboost::Binomial(),
               initial_learner = FALSE, high_dim_fit = TRUE, shrinkage = 0.1,
               max_iterations = 3, ...)
  }
  default <- fit()
  explicit <- fit(phi_boost_exit = 0.95)
  early <- fit(phi_boost_exit = 0.01)
  expect_identical(default$iters, explicit$iters)
  expect_identical(predict(default, n = 2), predict(explicit, n = 2))
  expect_equal(early$iters, 2L)
  expect_gt(default$iters, early$iters)
})

test_that("the screen is opt-in and rejected learners can still win as lines", {
  x <- seq(-1, 1, length.out = 30)
  dat <- data.frame(y = factor(ifelse(x > 0 & seq_along(x) %% 5 != 0, "yes", "no")), x = x)
  family <- mboost::Binomial()
  local_mocked_bindings(.boost_screen_calibrate = function(...) 1, .package = "GeDS")
  fit <- NGeDSboost(y ~ f(x), data = dat, family = family, initial_learner = FALSE,
                    max_iterations = 3, shrinkage = 0.1,
                    high_dim_fit = TRUE)
  expect_equal(fit$args$shared_screen_cutoff, 1)
  expect_length(fit$models, 4)
  expect_true(all(lengths(knots(fit, n = 2, options = "internal")) == 0))
  expect_null(fit$final_model$refinement)
  expect_null(fit$final_model$quadratic.fit)
  expect_null(fit$final_model$cubic.fit)
  expect_equal(predict(fit, n = 2), predict(fit, newdata = dat, n = 2))

  default <- suppressWarnings(NGeDSboost(y ~ f(x), data = dat, family = family,
                                        initial_learner = FALSE, max_iterations = 1,
                                        int.knots_boost = 0))
  expect_null(default$args$shared_screen_cutoff)
  expect_false(is.null(default$final_model$refinement))
  expect_false(is.null(default$final_model$quadratic.fit))
  expect_false(is.null(default$final_model$cubic.fit))
  override <- suppressWarnings(NGeDSboost(
    y ~ f(x), data = dat, family = family, initial_learner = FALSE,
    max_iterations = 3, shrinkage = 0.1, high_dim_fit = TRUE, higher_order = TRUE
  ))
  expect_false(is.null(override$final_model$quadratic.fit))
  expect_false(is.null(override$final_model$cubic.fit))
  expect_null(override$final_model$refinement)
  expect_equal(predict(override, n = 2), predict(fit, n = 2))
  expect_error(NGeDSboost(y ~ f(x), data = transform(dat, y = x), high_dim_fit = TRUE),
               "requires unweighted binomial logit")
  expect_error(NGeDSboost(y ~ f(x), data = dat, high_dim_fit = TRUE, stoptype = "RD"),
               "requires stoptype")
})

test_that("the shared cutoff is reproducible for fixed bootstrap draws", {
  x <- seq_len(20)
  y <- rep(c(-1, 1), 10)
  bases <- lapply(list(x, sin(x)), function(x) {
    ord <- order(x)
    GeDS:::.boost_screen_basis(list(ord = ord, xs = x[ord], extr = range(x)))
  })
  get_cutoff <- function(bases) {
    set.seed(123)
    GeDS:::.boost_screen_calibrate(bases, y, rep(1, 20), mboost::Binomial(), 0.5)
  }
  cutoff <- get_cutoff(bases)
  expect_equal(cutoff, get_cutoff(rev(bases)))
  expect_gt(cutoff, 0)
  expect_lt(cutoff, 1)
})

test_that("a first-knot gain above the cutoff is passed on to GD", {
  x <- seq(-1, 1, length.out = 30)
  dat <- data.frame(y = factor(ifelse(abs(x) > 0.5, "yes", "no")), x = x)
  limits <- numeric()
  native_fit <- GeDS:::.gamma_dynamic_fit
  local_mocked_bindings(
    .boost_screen_calibrate = function(...) 0,
    .gamma_dynamic_fit = function(fitter, args, ...) {
      limits <<- c(limits, args$max.intknots)
      args$max.intknots <- 1L
      native_fit(fitter, args, ...)
    },
    .package = "GeDS"
  )
  suppressWarnings(NGeDSboost(y ~ f(x), data = dat, family = mboost::Binomial(),
                              initial_learner = FALSE, max_iterations = 1,
                              higher_order = FALSE, high_dim_fit = TRUE))
  expect_length(limits, 1)
  expect_gt(limits, 1)
})
