library("testthat")
library("GeDS")

test_that("f() rejects more than two covariates", {
  expect_error(
    f(1:5, 1:5, 1:5)
  )
})

test_that("NGeDS rejects invalid beta values", {
  set.seed(123)
  N <- 50
  X <- sort(runif(N))
  Y <- sin(2 * pi * X) + rnorm(N, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  expect_error(
    NGeDS(Y ~ f(X), data = dat, beta = -0.1)
  )

  expect_error(
    NGeDS(Y ~ f(X), data = dat, beta = 1.1)
  )
})

test_that("NGeDS rejects weights of incorrect length", {
  set.seed(123)
  N <- 50
  X <- sort(runif(N))
  Y <- sin(2 * pi * X) + rnorm(N, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  expect_error(
    NGeDS(Y ~ f(X), data = dat, weights = rep(1, N - 1))
  )
})

test_that("methods warn and default to n = 3 for unavailable spline orders", {
  set.seed(123)
  N <- 50
  X <- sort(runif(N))
  Y <- sin(2 * pi * X) + rnorm(N, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  Gmod <- NGeDS(Y ~ f(X), data = dat, phi = 0.9)

  pred_n3 <- predict(Gmod, newdata = dat, n = 3)
  coef_n3 <- coef(Gmod, n = 3)
  dev_n3  <- deviance(Gmod, n = 3)
  knot_n3 <- knots(Gmod, n = 3)

  expect_warning(
    pred_n5 <- predict(Gmod, newdata = dat, n = 5),
    "'n' incorrectly specified"
  )
  expect_equal(pred_n5, pred_n3)

  expect_warning(
    coef_n5 <- coef(Gmod, n = 5),
    "'n' incorrectly specified"
  )
  expect_equal(coef_n5, coef_n3)

  expect_warning(
    dev_n5 <- deviance(Gmod, n = 5),
    "'n' incorrectly specified"
  )
  expect_equal(dev_n5, dev_n3)

  expect_warning(
    knot_n5 <- knots(Gmod, n = 5),
    "'n' incorrectly specified"
  )
  expect_equal(knot_n5, knot_n3)
})

test_that("GeDSgam methods warn and default to n = 3 for unavailable spline orders", {
  data(mtcars)
  mtcars$cyl <- factor(mtcars$cyl)

  Gmodgam <- NULL
  invisible(capture.output({
    Gmodgam <- suppressWarnings(
      NGeDSgam(
        mpg ~ cyl + f(wt),
        data = mtcars,
        family = gaussian,
        phi_gam_exit = 0.8,
        phi = 0.8
      )
    )
  }))

  pred_n3 <- predict(Gmodgam, newdata = mtcars, n = 3)
  coef_n3 <- coef(Gmodgam, n = 3)
  dev_n3  <- deviance(Gmodgam, n = 3)
  knot_n3 <- knots(Gmodgam, n = 3)

  expect_warning(
    pred_n5 <- predict(Gmodgam, newdata = mtcars, n = 5),
    "incorrectly specified"
  )
  expect_equal(pred_n5, pred_n3)

  expect_warning(
    coef_n5 <- coef(Gmodgam, n = 5),
    "incorrectly specified"
  )
  expect_equal(coef_n5, coef_n3)

  expect_warning(
    dev_n5 <- deviance(Gmodgam, n = 5),
    "incorrectly specified"
  )
  expect_equal(dev_n5, dev_n3)

  expect_warning(
    knot_n5 <- knots(Gmodgam, n = 5),
    "incorrectly specified"
  )
  expect_equal(knot_n5, knot_n3)
})

test_that("GeDSboost methods warn and default to n = 3 for unavailable spline orders", {
  data(mtcars)
  mtcars$cyl <- factor(mtcars$cyl)

  Gmodboost <- NULL
  invisible(capture.output({
    Gmodboost <- suppressWarnings(
      NGeDSboost(
        mpg ~ cyl + f(wt),
        data = mtcars,
        family = mboost::Gaussian(),
        phi_boost_exit = 0.8,
        phi = 0.8
      )
    )
  }))

  pred_n3 <- predict(Gmodboost, newdata = mtcars, n = 3)
  coef_n3 <- coef(Gmodboost, n = 3)
  dev_n3  <- deviance(Gmodboost, n = 3)
  knot_n3 <- knots(Gmodboost, n = 3)

  expect_warning(
    pred_n5 <- predict(Gmodboost, newdata = mtcars, n = 5),
    "incorrectly specified"
  )
  expect_equal(pred_n5, pred_n3)

  expect_warning(
    coef_n5 <- coef(Gmodboost, n = 5),
    "incorrectly specified"
  )
  expect_equal(coef_n5, coef_n3)

  expect_warning(
    dev_n5 <- deviance(Gmodboost, n = 5),
    "incorrectly specified"
  )
  expect_equal(dev_n5, dev_n3)

  expect_warning(
    knot_n5 <- knots(Gmodboost, n = 5),
    "incorrectly specified"
  )
  expect_equal(knot_n5, knot_n3)
})

test_that("Derive rejects derivative order not lower than spline order", {
  set.seed(123)
  N <- 50
  X <- sort(runif(N))
  Y <- sin(2 * pi * X) + rnorm(N, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  Gmod <- NGeDS(Y ~ f(X), data = dat, phi = 0.9)

  expect_error(
    Derive(Gmod, x = c(0.25, 0.5, 0.75), order = 3, n = 3)
  )
})

test_that("GGeDS rejects invalid family specifications", {
  set.seed(123)
  N <- 50
  X <- sort(runif(N))
  Y <- rpois(N, lambda = exp(sin(2 * pi * X)))
  dat <- data.frame(X = X, Y = Y)

  expect_error(
    GGeDS(Y ~ f(X), data = dat, family = "not_a_family")
  )
})

test_that("NGeDSgam prediction fails for unsupported base learner names", {
  data(mtcars)
  mtcars$cyl <- factor(mtcars$cyl)

  Gmodgam <- suppressWarnings(
    NGeDSgam(
      mpg ~ cyl + f(wt),
      data = mtcars,
      family = gaussian,
      phi_gam_exit = 0.8,
      phi = 0.8
    )
  )

  expect_error(
    predict(Gmodgam, newdata = mtcars, n = 2, base_learner = "f(not_in_model)")
  )
})

test_that("NGeDSboost prediction fails for unsupported base learner names", {
  data(mtcars)
  mtcars$cyl <- factor(mtcars$cyl)

  Gmodboost <- suppressWarnings(
    NGeDSboost(
      mpg ~ cyl + f(wt),
      data = mtcars,
      family = mboost::Gaussian(),
      phi_boost_exit = 0.8,
      phi = 0.8
    )
  )

  expect_error(
    predict(Gmodboost, newdata = mtcars, n = 2, base_learner = "f(not_in_model)")
  )
})

test_that("GeDSgam and GeDSboost predictions support n = 'all'", {
  data(mtcars)
  mtcars$cyl <- factor(mtcars$cyl)

  Gmodgam <- NULL
  invisible(capture.output({
    Gmodgam <- suppressWarnings(
      NGeDSgam(
        mpg ~ cyl + f(wt),
        data = mtcars,
        family = gaussian,
        phi_gam_exit = 0.8,
        phi = 0.8
      )
    )
  }))

  Gmodboost <- NULL
  invisible(capture.output({
    Gmodboost <- suppressWarnings(
      NGeDSboost(
        mpg ~ cyl + f(wt),
        data = mtcars,
        family = mboost::Gaussian(),
        phi_boost_exit = 0.8,
        phi = 0.8
      )
    )
  }))

  pred_gam <- predict(Gmodgam, newdata = mtcars, n = "all")
  pred_boost <- predict(Gmodboost, newdata = mtcars, n = "all")

  expect_named(pred_gam, c("pred_linear", "pred_quadratic", "pred_cubic"))
  expect_named(pred_boost, c("pred_linear", "pred_quadratic", "pred_cubic"))

  expect_true(all(vapply(pred_gam, length, integer(1)) == nrow(mtcars)))
  expect_true(all(vapply(pred_boost, length, integer(1)) == nrow(mtcars)))

  expect_true(all(vapply(pred_gam, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(vapply(pred_boost, function(x) all(is.finite(x)), logical(1))))
})


