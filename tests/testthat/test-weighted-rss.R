library("testthat")
library("GeDS")

test_that("Normal spline fitters report the weighted residual sum of squares", {
  x <- seq(0, 1, length.out = 40)
  y <- sin(5 * x) + 0.2 * x^2
  weights <- seq(0.25, 2.5, length.out = length(x))

  univariate <- lapply(c(FALSE, TRUE), function(fast) {
    SplineReg_LM(
      X = x, Y = y, weights = weights,
      InterKnots = c(0.3, 0.65), n = 2L,
      only_pred = TRUE, fast = fast
    )
  })
  for (fit in univariate) {
    expect_equal(fit$rss, sum(weights * fit$residuals^2), tolerance = 1e-13)
  }

  grid <- expand.grid(
    x = seq(0, 1, length.out = 8),
    y = seq(0, 1, length.out = 7)
  )
  response <- with(grid, sin(4 * x) + cos(3 * y) + x * y^2)
  biv_weights <- seq(0.3, 2.2, length.out = NROW(grid))

  bivariate <- lapply(c(FALSE, TRUE), function(fast) {
    SplineReg_biv(
      X = grid$x, Y = grid$y, Z = response, weights = biv_weights,
      InterKnotsX = 0.5, InterKnotsY = 0.45, n = 2L,
      fast = fast
    )
  })
  for (fit in bivariate) {
    expect_equal(fit$rss,
                 sum(biv_weights * fit$residuals^2),
                 tolerance = 1e-13)
  }
})

test_that("unit weights retain the ordinary residual sum of squares", {
  x <- seq(0, 1, length.out = 40)
  y <- sin(5 * x) + 0.2 * x^2

  fit <- SplineReg_LM(
    X = x, Y = y, InterKnots = c(0.3, 0.65), n = 2L,
    only_pred = TRUE
  )

  expect_equal(fit$rss, sum(fit$residuals^2), tolerance = 1e-13)
})

test_that("repeated-site residuals retain legacy GeDS averaging", {
  contributions <- c(2, -1, 3, 4, 8, -2)
  recurrence <- c(3, 2, 1)

  expect_equal(
    GeDS:::makeNewRes(contributions, recurrence),
    c(mean(contributions[1:3]), mean(contributions[4:5]), contributions[6]),
    tolerance = 1e-15
  )
})

test_that("an initial univariate fit enters Stage A with weighted RSS", {
  x <- seq(0, 1, length.out = 50)
  y <- sin(5 * x) + 0.2 * x^2
  weights <- seq(0.25, 2.5, length.out = length(x))
  intknots <- c(0.3, 0.65)
  initial <- SplineReg_LM(
    X = x, Y = y, weights = weights,
    InterKnots = intknots, n = 2L, only_pred = TRUE
  )

  fit <- suppressWarnings(
    UnivariateFitter(
      X = x, Y = y, weights = weights,
      max.intknots = 3L, q = 1L, show.iters = FALSE,
      stoptype = "RD", higher_order = FALSE,
      fit_init = list(
        pred = as.numeric(initial$predicted),
        intknots = intknots,
        coef = initial$theta
      ),
      only_pred = TRUE
    )
  )

  initial_iteration <- length(intknots) + 1L
  expect_equal(fit$rss[initial_iteration],
               sum(weights * initial$residuals^2),
               tolerance = 1e-13)
})

test_that("multidimensional Normal fits report weighted RSS", {
  grid <- expand.grid(
    x0 = seq(0, 1, length.out = 6),
    x1 = seq(0, 1, length.out = 5),
    x2 = seq(0, 1, length.out = 4)
  )
  coordinates <- as.matrix(grid)
  response <- with(grid, sin(3 * x0) + x1^2 - x2 + x0 * x1 * x2)
  weights <- seq(0.4, 1.8, length.out = NROW(coordinates))

  stage_a <- GeDS:::stageAOneStepND(
    coordinates = coordinates,
    response = response,
    upper.bounds = rep(list(c(0.5, 1) + 1e-15), 3L),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    weights = weights
  )
  expect_equal(stage_a$rss,
               sum(weights * stage_a$residuals^2),
               tolerance = 1e-13)

  stage_b <- GeDS:::fitTensorSplineND(
    coordinates = coordinates,
    response = response,
    intknots = rep(list(0.5), 3L),
    spline.order = 2L,
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    weights = weights
  )
  expect_equal(stage_b$rss,
               sum(weights * stage_b$residuals^2),
               tolerance = 1e-13)
})

test_that("the additive Normal spline fitter reports weighted RSS", {
  x <- seq(0, 1, length.out = 30)
  response <- sin(4 * x) + x
  weights <- seq(0.3, 2, length.out = length(x))

  fit <- GeDS:::SplineReg_LM_Multivar(
    X = data.frame(x = x), Y = response, Z = NULL,
    base_learners = list(s = list(type = "GeDS", variables = "x")),
    weights = weights, InterKnotsList = list(s = c(0.35, 0.7)),
    n = 2L, coefficients = NULL, linear.predictors = NULL,
    only_pred = TRUE
  )

  expect_equal(fit$rss,
               sum(weights * fit$residuals^2),
               tolerance = 1e-13)
})

test_that("legacy bivariate Stage B forwards fitting weights", {
  set.seed(42)
  x <- runif(90)
  y <- runif(90)
  W <- cbind(z = rnorm(90))
  response <- sin(4 * x) + cos(3 * y) + 0.8 * W[, 1] +
    rnorm(90, sd = 0.15)
  weights <- seq(0.25, 2.5, length.out = length(x))

  fit <- suppressWarnings(
    BivariateFitter(
      X = x, Y = y, Z = response, W = W, weights = weights,
      indicator = NULL, beta = 0.5, phi = 0.99,
      min.intknots = 0L, max.intknots = 5L, q = 2L,
      Xextr = c(0, 1), Yextr = c(0, 1),
      show.iters = FALSE, stoptype = "RD"
    )
  )
  direct <- SplineReg_biv(
    X = x, Y = y, Z = response, W = W, weights = weights,
    InterKnotsX = fit$linear.intknots$Xk,
    InterKnotsY = fit$linear.intknots$Yk,
    Xextr = c(0, 1), Yextr = c(0, 1), n = 2L
  )

  expect_equal(fit$linear.fit$predicted,
               direct$predicted,
               tolerance = 1e-11)
  expect_equal(fit$linear.fit$rss,
               sum(weights * fit$linear.fit$residuals^2),
               tolerance = 1e-12)
})

test_that("additive fitting criteria use their fitting weights", {
  dat <- data.frame(
    y = c(0, 1, 2, 6, 8),
    x = 0:4
  )
  weights <- c(0.2, 0.5, 1, 2, 3)
  base_learners <- list(x = list(type = "linear", variables = "x"))

  component <- GeDS:::componentwise_fit(
    bl_name = "x", response = "y", data = dat,
    model_formula_template = "y ~ ", family = gaussian(),
    weights = weights, base_learners = base_learners, m = 1L,
    internal_knots = 0L, beta = 0.5, phi = 0.99, q = 2L
  )
  expect_equal(component$ssr,
               sum(weights * component$resid^2),
               tolerance = 1e-13)

  initial <- list(x = list(coefficients = list(b0 = 0, b1 = 0)))
  backfit <- GeDS:::backfitting(
    z = dat$y, base_learners = base_learners,
    base_learners_list = initial, data = dat, wz = weights,
    phi_gam_exit = 0.99, q_gam = 1L, iter = 1L,
    internal_knots = 2L, beta = 0.5, phi = 0.99, q = 2L
  )
  reference <- lm(y ~ x, data = dat, weights = weights)
  expect_equal(backfit$z_hat,
               unname(fitted(reference)),
               tolerance = 1e-10)
})
