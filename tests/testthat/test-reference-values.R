library("testthat")
library("GeDS")

test_that("NGeDS returns stable reference values", {
  set.seed(123)
  N <- 80
  X <- sort(runif(N, min = -2, max = 2))
  Y <- (10 * X / (1 + 100 * X^2)) * 4 + 4 + rnorm(N, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  Gmod <- suppressWarnings(
    NGeDS(Y ~ f(X), data = dat, beta = 0.6, phi = 0.95, stoptype = "RD")
  )

  ref <- list(
    deviance = round(c(
      linear = deviance(Gmod, n = 2),
      quadratic = deviance(Gmod, n = 3),
      cubic = deviance(Gmod, n = 4)
    ), 6),
    pred_first_5 = round(predict(Gmod, newdata = dat[1:5, ], n = 3), 6),
    n_internal_knots = c(
      linear = length(knots(Gmod, n = 2, options = "internal")),
      quadratic = length(knots(Gmod, n = 3, options = "internal")),
      cubic = length(knots(Gmod, n = 4, options = "internal"))
    )
  )

  expect_equal(ref$deviance,
               c(linear = 0.532337, quadratic = 0.570885, cubic = 0.622715),
               tolerance = 1e-6)
  expect_equal(ref$pred_first_5,
               c(3.769455, 3.780928, 3.785291, 3.785762, 3.785793),
               tolerance = 1e-6)
  expect_equal(ref$n_internal_knots,
               c(linear = 10L, quadratic = 9L, cubic = 8L))
})

test_that("GGeDS Poisson fit returns stable reference values", {
  set.seed(123)
  N <- 80
  X <- sort(runif(N, min = -2, max = 2))
  eta <- sin(X) + 1
  Y <- rpois(N, lambda = exp(eta))
  dat <- data.frame(X = X, Y = Y)

  Gmod <- suppressWarnings(
    GGeDS(Y ~ f(X), data = dat, family = poisson(), beta = 0.2, phi = 0.95, stoptype = "SR")
  )

  ref <- list(
    deviance = round(c(
      linear = deviance(Gmod, n = 2),
      quadratic = deviance(Gmod, n = 3),
      cubic = deviance(Gmod, n = 4)
    ), 6),
    response_pred_first_5 = round(predict(Gmod, newdata = dat[1:5, ], n = 3), 6),
    link_pred_first_5 = round(predict(Gmod, newdata = dat[1:5, ], n = 3, type = "link"), 6)
  )

  expect_equal(ref$deviance,
               c(linear = 79.842135, quadratic = 85.539446, cubic = 79.997988),
               tolerance = 1e-6)
  expect_equal(ref$response_pred_first_5,
               c(0.469763, 0.529133, 0.576, 0.585782, 0.586556),
               tolerance = 1e-6)
  expect_equal(ref$link_pred_first_5,
               c(-0.755527, -0.636516, -0.551648, -0.534808, -0.533487),
               tolerance = 1e-6)
})

test_that("NGeDSgam returns stable reference values", {
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
        phi = 0.8,
        stoptype = "RD", penalty = "none"
      )
    )
  }))

  ref <- list(
    deviance = round(c(
      linear = deviance(Gmodgam, n = 2),
      quadratic = deviance(Gmodgam, n = 3),
      cubic = deviance(Gmodgam, n = 4)
    ), 6),
    pred_first_5 = round(predict(Gmodgam, newdata = mtcars[1:5, ], n = 3), 6),
    local_scoring_iterations = Gmodgam$iters$local_scoring
  )

  expect_equal(ref$deviance,
               c(linear = 144.647406, quadratic = 159.528346, cubic = 147.858943),
               tolerance = 1e-6)
  expect_equal(ref$pred_first_5,
               c(21.912536, 20.710742, 26.256086, 19.261310, 16.508585),
               tolerance = 1e-6)
  expect_equal(ref$local_scoring_iterations, 1L)
  # With no Stage B internal knots, the spline is a global polynomial.
  for (n in c(3, 4)) {
    polynomial <- lm(mpg ~ cyl + poly(wt, n - 1, raw = TRUE), data = mtcars)
    expect_equal(predict(Gmodgam, newdata = mtcars, n = n),
                 unname(predict(polynomial, mtcars)), tolerance = 1e-10)
  }
})

test_that("NGeDSboost returns stable reference values", {
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
        phi = 0.8,
        stoptype = "RD", penalty = "none"
      )
    )
  }))

  ref <- list(
    deviance = round(c(
      linear = deviance(Gmodboost, n = 2),
      quadratic = deviance(Gmodboost, n = 3),
      cubic = deviance(Gmodboost, n = 4)
    ), 6),
    pred_first_5 = round(predict(Gmodboost, newdata = mtcars[1:5, ], n = 3), 6),
    boosting_iterations = N.boost.iter(Gmodboost)
  )

  expect_equal(ref$deviance,
               c(linear = 172.616114, quadratic = 159.528346, cubic = 147.858943),
               tolerance = 1e-6)
  expect_equal(ref$pred_first_5,
               c(21.912536, 20.710742, 26.256086, 19.261310, 16.508585),
               tolerance = 1e-6)
  expect_equal(ref$boosting_iterations, 5L)
  # With unit shrinkage, Gaussian Phase 2 reaches LS on the selected basis.
  internal <- knots(Gmodboost, n = 2, options = "internal")[["f(wt)"]]
  linear_reference <- lm(mpg ~ cyl + splines::bs(wt, degree = 1, knots = internal), data = mtcars)
  expect_equal(deviance(Gmodboost, n = 2), deviance(linear_reference), tolerance = 1e-10)
  for (n in 3:4) {
    polynomial <- lm(mpg ~ cyl + poly(wt, n - 1, raw = TRUE), data = mtcars)
    expect_equal(predict(Gmodboost, newdata = mtcars, n = n),
                 unname(predict(polynomial, mtcars)), tolerance = 1e-10)
  }
})
