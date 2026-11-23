test_that("standalone Gaussian splines retain interiors and extend fixed boundaries", {
  data <- data.frame(x = seq(0, 1, length.out = 80))
  data$z <- sin(seq_len(nrow(data)))
  data$off <- .1 + .2 * data$x
  data$y <- 2 + data$x + .3 * data$z + data$off
  fit <- suppressWarnings(NGeDS(y ~ f(x) + z + offset(off), data = data,
                                min.intknots = 4, max.intknots = 6, q = 1))
  newdata <- data.frame(x = c(-.2, .4, 1.2), z = c(.2, .3, .4), off = .5)
  clipped <- newdata
  clipped$x <- pmax(0, pmin(1, clipped$x))
  for (n in 2:4) {
    expect_equal(predict(fit, newdata = data, n = n), predict(fit, n = n))
    pred <- suppressWarnings(predict(fit, newdata = newdata, n = n))
    at_boundary <- predict(fit, newdata = clipped, n = n)
    if (n == 2L) {
      expect_equal(pred - at_boundary, newdata$x - clipped$x, tolerance = 1e-7)
    } else {
      expect_equal(pred, at_boundary)
    }
    shifted <- newdata
    shifted$off <- shifted$off + .25
    expect_equal(suppressWarnings(predict(fit, shifted, n = n)) - pred, rep(.25, 3))
    expect_equal(pred, vapply(seq_len(nrow(newdata)), function(i) {
      suppressWarnings(predict(fit, newdata[i, , drop = FALSE], n = n))
    }, numeric(1)))
  }
})

test_that("standalone bivariate splines continue bilinearly or clamp by order", {
  data <- expand.grid(x = seq(0, 1, length.out = 8), z = seq(0, 1, length.out = 8))
  data$y <- 2 + data$x + data$z + data$x * data$z
  fit <- suppressWarnings(NGeDS(y ~ f(x, z), data = data,
                                min.intknots = 2, max.intknots = 4, q = 1))
  newdata <- data.frame(x = c(-.2, .4, 1.2), z = c(1.2, .3, -.2))
  clipped <- as.data.frame(lapply(newdata, function(x) pmax(0, pmin(1, x))))
  for (n in 2:4) {
    pred <- suppressWarnings(predict(fit, newdata, n = n))
    expect_equal(predict(fit, data, n = n), predict(fit, n = n))
    if (n == 2L) {
      expect_equal(pred, 2 + newdata$x + newdata$z + newdata$x * newdata$z, tolerance = 1e-7)
    } else {
      expect_equal(pred, predict(fit, clipped, n = n))
    }
  }
})

test_that("Gamma inverse-link standalone predictions clamp all orders before inversion", {
  data <- data.frame(x = seq(0, 1, length.out = 90))
  data$z <- sin(seq_len(nrow(data)))
  data$off <- .1 + .1 * data$x
  data$y <- 1 / (2 + .4 * data$x + .1 * data$z + data$off)
  fit <- suppressWarnings(GGeDS(y ~ f(x) + z + offset(off), data = data,
    family = Gamma(), min.intknots = 4, max.intknots = 6, q = 1))
  newdata <- data.frame(x = c(-.5, .5, 1.5), z = .2, off = .3)
  clipped <- newdata
  clipped$x <- pmax(0, pmin(1, clipped$x))
  for (n in 2:4) {
    eta <- suppressWarnings(predict(fit, newdata, n = n, type = "link"))
    expect_equal(eta, predict(fit, clipped, n = n, type = "link"))
    expect_equal(suppressWarnings(predict(fit, newdata, n = n)), 1 / eta)
    expect_true(all(eta > 0))
    shifted <- newdata
    shifted$off <- shifted$off + .25
    expect_equal(suppressWarnings(predict(fit, shifted, n = n, type = "link")) - eta, rep(.25, 3))
    shifted$off <- -100
    expect_error(suppressWarnings(predict(fit, shifted, n = n)), "positive finite linear predictor")
  }
})

test_that("Gamma inverse-link bivariate splines clamp every coordinate at every order", {
  data <- expand.grid(x = seq(0, 1, length.out = 8), z = seq(0, 1, length.out = 8))
  data$y <- 1 / (2 + .3 * data$x + .2 * data$z)
  fit <- suppressWarnings(GGeDS(y ~ f(x, z), data = data, family = Gamma(),
    min.intknots = 4, max.intknots = 6, q = 1))
  newdata <- data.frame(x = c(-.5, .4, 1.5), z = c(1.5, .2, -.5))
  clipped <- as.data.frame(lapply(newdata, function(x) pmax(0, pmin(1, x))))
  for (n in 2:4) {
    eta <- suppressWarnings(predict(fit, newdata, n = n, type = "link"))
    expect_equal(eta, predict(fit, clipped, n = n, type = "link"))
    expect_equal(suppressWarnings(predict(fit, newdata, n = n)), 1 / eta)
  }
})

test_that("Gamma inverse-link additive full and component predictions clamp linear smoothers", {
  data <- data.frame(x = seq(0, 1, length.out = 70))
  data$y <- 1 / (2 + .4 * data$x + .1 * sin(3 * data$x))
  fit <- NULL
  invisible(capture.output(fit <- suppressWarnings(NGeDSgam(
    y ~ f(x), data = data, family = Gamma(), phi = .7, q = 1, penalty = "none"))))
  newdata <- data.frame(x = c(-.5, .5, 1.5))
  clipped <- data.frame(x = c(0, .5, 1))
  for (n in 2:4) {
    expect_equal(predict(fit, newdata, n = n), predict(fit, clipped, n = n))
    expect_equal(predict(fit, newdata, n = n, base_learner = "f(x)"),
                 predict(fit, clipped, n = n, base_learner = "f(x)"))
    eta <- predict(fit, newdata, n = n, type = "link")
    expect_equal(predict(fit, newdata, n = n), 1 / eta)
  }
})

test_that("boundary exceptions depend on both family and link", {
  expect_true(GeDS:::.prediction_clamp_linear(Gamma()))
  expect_false(GeDS:::.prediction_clamp_linear(Gamma(link = "log")))
  expect_false(GeDS:::.prediction_clamp_linear(gaussian(link = "inverse")))
  expect_false(GeDS:::.prediction_clamp_linear(NULL))
  boost <- structure(list(args = list(family = mboost::GammaReg())), class = "GeDSboost")
  expect_false(GeDS:::.prediction_clamp_linear(GeDS:::.prediction_family(boost)))
  knots <- c(0, 0, .5, 1, 1)
  x <- matrix(c(-1, .5, 2), ncol = 1)
  expect_equal(GeDS:::.prediction_spline_basis(knots, x, 2, clamp_linear = TRUE),
               splines::splineDesign(knots, c(0, .5, 1), ord = 2))
})

test_that("Gamma inverse-link GAMs clamp bivariate smoothers without clamping linear terms", {
  data <- expand.grid(x = seq(0, 1, length.out = 8), z = seq(0, 1, length.out = 8))
  data$u <- sin(seq_len(nrow(data)))
  data$y <- 1 / (2 + .3 * data$x + .2 * data$z + .05 * data$u)
  fit <- NULL
  invisible(capture.output(fit <- suppressWarnings(NGeDSgam(
    y ~ f(x, z) + u, data = data, family = Gamma(), phi = .7, q = 1, penalty = "none"))))
  newdata <- data.frame(x = c(-.5, .4, 1.5), z = c(1.5, .2, -.5), u = .3)
  clipped <- newdata
  clipped$x <- pmax(0, pmin(1, clipped$x))
  clipped$z <- pmax(0, pmin(1, clipped$z))
  for (n in 2:4) {
    expect_equal(predict(fit, newdata, n = n), predict(fit, clipped, n = n))
    expect_equal(predict(fit, newdata, n = n, base_learner = "f(x, z)"),
                 predict(fit, clipped, n = n, base_learner = "f(x, z)"))
    shifted <- newdata
    shifted$u <- shifted$u + .2
    expect_equal(predict(fit, shifted, n = n, base_learner = "f(x, z)"),
                 predict(fit, newdata, n = n, base_learner = "f(x, z)"))
    expect_equal(predict(fit, shifted, n = n, type = "link") -
                   predict(fit, newdata, n = n, type = "link"),
                 predict(fit, shifted, n = n, base_learner = "u") -
                   predict(fit, newdata, n = n, base_learner = "u"), tolerance = 1e-7)
  }
})
