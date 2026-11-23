test_that("prediction bases preserve interiors, clamp higher orders and continue linear segments", {
  for (n in 2:4) {
    knots <- c(rep(0, n), 0.2, 0.6, rep(1, n))
    x <- c(-0.5, 0, 0.3, 1, 1.5)
    basis <- GeDS:::.prediction_spline_basis(knots, x, n)
    inside <- x >= 0 & x <= 1
    expect_equal(basis[inside, ], splines::splineDesign(knots, x[inside], ord = n))
    expect_equal(rowSums(basis), rep(1, length(x)), tolerance = 1e-12)
    greville <- vapply(seq_len(length(knots) - n), function(j) {
      mean(knots[(j + 1):(j + n - 1)])
    }, numeric(1))
    target <- if (n == 2) x else pmax(0, pmin(1, x))
    expect_equal(as.numeric(basis %*% greville), target, tolerance = 1e-12)
    individual <- do.call(rbind, lapply(x, function(z) GeDS:::.prediction_spline_basis(knots, z, n)))
    expect_equal(basis, individual)
    for (boundary in c(0, 1)) {
      probe <- boundary + c(-1e-8, 0, 1e-8)
      values <- GeDS:::.prediction_spline_basis(knots, probe, n) %*% greville
      expect_lt(max(abs(diff(as.numeric(values)))), 1e-6)
    }
  }
})

test_that("boundary rules also preserve normalised predictions", {
  set.seed(123)
  data <- expand.grid(x = seq(2, 5, length.out = 6), z = seq(10, 20, length.out = 6))
  data$y <- data$x + .2 * data$z + rnorm(nrow(data), sd = .1)
  newdata <- data.frame(x = c(1, 3, 6), z = c(8, 15, 22))
  clipped <- data.frame(x = c(2, 3, 5), z = c(10, 15, 20))
  for (fun in list(NGeDSgam, NGeDSboost)) {
    fit <- NULL
    invisible(capture.output(fit <- suppressWarnings(fun(y ~ f(x, z), data = data,
      phi = .7, q = 1, normalize_data = TRUE, penalty = "none"))))
    for (n in 2:4) {
      expect_equal(predict(fit, newdata = data, n = n), predict(fit, n = n), tolerance = 1e-6)
      expect_equal(predict(fit, newdata = newdata, n = n), vapply(1:3, function(i) {
        predict(fit, newdata = newdata[i, , drop = FALSE], n = n)
      }, numeric(1)), tolerance = 1e-8)
      if (n > 2) expect_equal(predict(fit, newdata = newdata, n = n),
                              predict(fit, newdata = clipped, n = n), tolerance = 1e-8)
    }
  }
})

test_that("joint linear prediction retains the fitted tensor and its mixed term", {
  learners <- list(joint = list(type = "GeDS", variables = c("x", "z")))
  internal <- list(joint = list(ikX = NULL, ikY = NULL))
  newdata <- data.frame(x = c(-1, 0.4, 2), z = c(2, 0.2, -1))
  for (n in 2:4) {
    basis <- GeDS:::jointSplineBasisMultivar(newdata, learners, internal, n,
      list(x = c(0, 1), z = c(0, 1)), prediction = TRUE)$basis
    # Coefficients of x*z in a tensor product of Bernstein bases.
    theta <- as.numeric(GeDS:::tensorProd(matrix(seq(0, 1, length.out = n), nrow = 1),
                                         matrix(seq(0, 1, length.out = n), nrow = 1)))
    x <- if (n == 2) newdata$x else pmax(0, pmin(1, newdata$x))
    z <- if (n == 2) newdata$z else pmax(0, pmin(1, newdata$z))
    expect_equal(as.numeric(basis %*% theta), x * z, tolerance = 1e-12)
  }
})

test_that("higher-dimensional component predictions follow the same boundary rule", {
  set.seed(123)
  data <- expand.grid(x = seq(0, 1, length.out = 5), z = seq(0, 1, length.out = 5),
                      w = seq(0, 1, length.out = 5))
  data$y <- data$x + data$z + data$w + rnorm(nrow(data), sd = .05)
  fit <- NULL
  invisible(capture.output(fit <- NGeDSgam(y ~ f(x, z, w), data = data, phi = .7,
                                         q = 1, penalty = "none")))
  newdata <- data.frame(x = c(-.2, 1.2), z = c(.5, 1.3), w = c(.2, -.1))
  clipped <- newdata
  clipped[] <- lapply(clipped, function(x) pmax(0, pmin(1, x)))
  for (n in 2:4) {
    pred <- function(d) predict(fit, newdata = d, n = n, base_learner = "f(x, z, w)")
    expect_equal(pred(newdata), vapply(1:2, function(i) pred(newdata[i, , drop = FALSE]), numeric(1)))
    if (n > 2) {
      expect_equal(pred(newdata), pred(clipped))
      expect_equal(pred(newdata), predict(fit, newdata = newdata, n = n))
    }
  }
})

test_that("GAM and boost full and component predictions use fixed boundary rules", {
  set.seed(123)
  data <- expand.grid(x = seq(0, 1, length.out = 7), z = seq(0, 1, length.out = 7))
  data$u <- seq_len(nrow(data)) / nrow(data)
  data$y <- 2 + data$x + data$z + data$x * data$z + sin(3 * data$u) + rnorm(nrow(data), sd = .05)
  newdata <- data.frame(x = c(-.2, .4, 1.2), z = c(.3, .4, 1.3), u = c(-.1, .5, 1.1))
  clamped <- newdata
  clamped[] <- lapply(clamped, function(x) pmax(0, pmin(1, x)))
  clamped$u <- pmax(min(data$u), clamped$u)
  for (fun in list(NGeDSgam, NGeDSboost)) {
    for (penalty in c("none", "component")) {
      if (penalty != "none") skip_if_not_installed("mgcv")
      fit <- NULL
      invisible(capture.output(fit <- suppressWarnings(fun(
        y ~ f(x, z) + f(u), data = data, phi = .7, q = 1, penalty = penalty
      ))))
      for (n in 2:4) {
        expect_equal(predict(fit, newdata = data, n = n), predict(fit, n = n), tolerance = 1e-6)
        for (bl in c(NA_character_, "f(x, z)", "f(u)")) {
          pred <- function(d) {
            if (is.na(bl)) predict(fit, newdata = d, n = n)
            else predict(fit, newdata = d, n = n, base_learner = bl)
          }
          expect_equal(pred(newdata), vapply(seq_len(nrow(newdata)), function(i) {
            pred(newdata[i, , drop = FALSE])
          }, numeric(1)), tolerance = 1e-8)
          if (n > 2) expect_equal(pred(newdata), pred(clamped), tolerance = 1e-8)
        }
      }
    }
  }
})
