library("testthat")
library("GeDS")

test_that("f() combines more than two covariates", {
  result <- f(1:5, 6:10, 11:15)

  expect_equal(dim(result), c(5L, 3L))
  expect_equal(unname(result), unname(cbind(1:5, 6:10, 11:15)))
})

test_that("Knotnew handles empty residual clusters", {
  result <- GeDS:::Knotnew(
    weights = numeric(),
    residuals = numeric(),
    x = numeric(),
    dcum = numeric(),
    oldknots = numeric(),
    tol = 1e-12,
    support_order = 2L
  )

  expect_equal(result, c(NA_real_, NA_real_))
})

test_that("Knotnew uses the CRAN-era strict support intervals", {
  x <- 1:5
  residuals <- rep(1, length(x))
  oldknots <- c(1, 5, 1, 5, 3)

  left <- GeDS:::Knotnew(
    weights = c(0.1, 1, 0.2, 0.3),
    residuals = residuals,
    x = x,
    dcum = c(1, 2, 4, 5),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )
  left_R <- GeDS:::Knotnew_R_clean(
    wht = c(0.1, 1, 0.2, 0.3),
    restmp = residuals,
    x = x,
    dcm = c(1, 2, 4, 5),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )

  right <- GeDS:::Knotnew(
    weights = c(0.1, 0.2, 0.3, 1, 0.4),
    residuals = residuals,
    x = x,
    dcum = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )
  right_R <- GeDS:::Knotnew_R_clean(
    wht = c(0.1, 0.2, 0.3, 1, 0.4),
    restmp = residuals,
    x = x,
    dcm = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )

  expect_equal(left, c(NA_real_, NA_real_))
  expect_equal(c(left_R$newknot, left_R$cluster_index), left)
  expect_equal(right, c(NA_real_, NA_real_))
  expect_equal(c(right_R$newknot, right_R$cluster_index), right)
})

test_that("Knotnew checks every CRAN-era support interval independently", {
  x <- c(0, 0.25, 0.5, 3, 4)
  weights <- c(0.1, 0.7, 0.8, 1, 0.9)
  residuals <- rep(1, length(x))
  oldknots <- c(0, 4, 0, 4, 1, 2)

  result <- GeDS:::Knotnew(
    weights = weights,
    residuals = residuals,
    x = x,
    dcum = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )
  result_R <- GeDS:::Knotnew_R_clean(
    wht = weights,
    restmp = residuals,
    x = x,
    dcm = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )

  # No proposal puts an observation strictly inside every support interval.
  expect_equal(result, c(NA_real_, NA_real_))
  expect_equal(c(result_R$newknot, result_R$cluster_index), result)
})

test_that("Knotnew may reuse observations across overlapping support intervals", {
  x <- c(0.5, 2, 3.5)
  weights <- c(0.1, 1, 0.2)
  residuals <- rep(1, length(x))
  oldknots <- c(0, 4, 0, 4, 1, 3)

  result <- GeDS:::Knotnew(
    weights = weights,
    residuals = residuals,
    x = x,
    dcum = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )
  result_R <- GeDS:::Knotnew_R_clean(
    wht = weights,
    restmp = residuals,
    x = x,
    dcm = seq_along(x),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )

  # For sorted knots 0, 0, 1, 2, 3, 4, 4, x = 0.5 supports the
  # first two overlapping intervals and x = 3.5 supports the last two.
  expect_equal(result, c(2, 2))
  expect_equal(c(result_R$newknot, result_R$cluster_index), result)
})

test_that("Knotnew reports exhaustion when every proposal lacks support", {
  x <- c(0, 1.5, 3)
  weights <- c(0.1, 1)
  original_weights <- weights
  residuals <- c(1, 2, 1)
  oldknots <- c(0, 3, 0, 3, 1)

  result <- GeDS:::Knotnew(
    weights = weights,
    residuals = residuals,
    x = x,
    dcum = c(1, 3),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )
  result_R <- GeDS:::Knotnew_R_clean(
    wht = weights,
    restmp = residuals,
    x = x,
    dcm = c(1, 3),
    oldknots = oldknots,
    tol = 1e-12,
    support_order = 2L
  )

  # The first proposal leaves (0, 1) empty; the boundary proposal also fails.
  expect_equal(result, c(NA_real_, NA_real_))
  expect_equal(result_R, list(newknot = NA_real_, cluster_index = NA_integer_))
  expect_identical(weights, original_weights)
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


