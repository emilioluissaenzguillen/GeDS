bivariate_offset_data <- function() {
  set.seed(123)
  data <- expand.grid(x = seq(0, 1, length.out = 10), v = seq(0, 1, length.out = 10))
  data$off <- runif(nrow(data), -.4, .4)
  data$w <- rep(c(.5, 1, 2, 0), length.out = nrow(data))
  data$signal <- 2 + sin(3 * data$x) + .4 * data$v
  data
}

test_that("weighted bivariate Gaussian offsets affect knot selection and every spline order", {
  data <- bivariate_offset_data()
  data$y <- data$signal + data$off + rnorm(nrow(data), sd = .05)
  data$adjusted <- data$y - data$off
  fit <- suppressWarnings(NGeDS(y ~ f(x, v) + offset(off), data = data,
    weights = w, phi = .9, q = 1, max.intknots = 6, stoptype = "RD"))
  adjusted <- suppressWarnings(NGeDS(adjusted ~ f(x, v), data = data,
    weights = w, phi = .9, q = 1, max.intknots = 6, stoptype = "RD"))
  expect_equal(fit$args$offset, data$off)
  expect_equal(fit$stored, adjusted$stored, tolerance = 1e-8)
  expect_equal(fit$rss, adjusted$rss, tolerance = 1e-8)
  for (n in 2:4) {
    pred <- predict(fit, data, n = n)
    expect_equal(pred, predict(fit, n = n))
    expect_equal(pred - data$off, predict(adjusted, data, n = n), tolerance = 1e-8)
    expect_equal(knots(fit, n = n), knots(adjusted, n = n), tolerance = 1e-8)
    kn <- knots(fit, n = n, options = "all")
    B <- tensorProd(splines::splineDesign(kn$Xk, data$x, ord = n),
                    splines::splineDesign(kn$Yk, data$v, ord = n))
    direct <- lm.wfit(B, data$adjusted, data$w)
    theta <- coef(direct)
    theta[is.na(theta)] <- 0
    expect_equal(pred, as.numeric(B %*% theta) + data$off, tolerance = 1e-8)
    shifted <- transform(data, off = off + .25)
    expect_equal(predict(fit, shifted, n = n) - pred, rep(.25, nrow(data)), tolerance = 1e-10)
    stage <- fit[[c("linear.fit", "quadratic.fit", "cubic.fit")[n - 1L]]]
    expect_equal(as.numeric(stage$predicted), pred)
    expect_equal(stage$rss, sum(data$w * (data$y - pred)^2), tolerance = 1e-8)
  }
})

test_that("bivariate GLM offsets and weights reach Stage A and every Stage B fit", {
  data <- bivariate_offset_data()
  data$y <- rpois(nrow(data), exp(.4 + .3 * data$x + .2 * data$v + data$off))
  seen <- list()
  original <- SplineReg_biv_GLM
  testthat::local_mocked_bindings(SplineReg_biv_GLM = function(...) {
    args <- list(...)
    seen[[length(seen) + 1L]] <<- list(order = args$n, offset = args$offset,
                                       weights = args$weights)
    original(...)
  }, .package = "GeDS")
  fit <- suppressWarnings(GGeDS(y ~ f(x, v) + offset(off), data = data,
    family = poisson(), weights = w, phi = .9, q = 1, max.intknots = 6))
  expect_equal(fit$args$offset, data$off)
  expect_true(all(2:4 %in% vapply(seen, `[[`, numeric(1), "order")))
  expect_true(all(vapply(seen, function(args) identical(args$offset, data$off), logical(1))))
  expect_true(all(vapply(seen, function(args) identical(args$weights, data$w), logical(1))))
  for (n in 2:4) {
    eta <- predict(fit, data, n = n, type = "link")
    mu <- predict(fit, data, n = n)
    expect_equal(mu, predict(fit, n = n))
    expect_equal(mu, exp(eta))
    kn <- knots(fit, n = n, options = "all")
    B <- tensorProd(splines::splineDesign(kn$Xk, data$x, ord = n),
                    splines::splineDesign(kn$Yk, data$v, ord = n))
    direct <- glm.fit(B, data$y, family = poisson(), weights = data$w,
                      offset = data$off, intercept = FALSE)
    expect_equal(mu, as.numeric(direct$fitted.values), tolerance = 1e-6)
    shifted <- transform(data, off = off + log(2))
    expect_equal(predict(fit, shifted, n = n, type = "link") - eta,
                 rep(log(2), nrow(data)), tolerance = 1e-10)
    expect_equal(predict(fit, shifted, n = n), 2 * mu, tolerance = 1e-8)
    stage <- fit[[c("linear.fit", "quadratic.fit", "cubic.fit")[n - 1L]]]
    expect_equal(as.numeric(stage$predicted), mu)
    expect_equal(stage$rss, sum(poisson()$dev.resids(data$y, mu, data$w)), tolerance = 1e-8)
  }
})

test_that("bivariate fits without offsets agree with explicitly zero offsets", {
  data <- bivariate_offset_data()
  data$zero <- 0
  data$y <- data$signal + rnorm(nrow(data), sd = .1)
  for (fun in list(NGeDS, GGeDS)) {
    for (weighted in c(FALSE, TRUE)) {
      args <- list(data = data, phi = .9, q = 1, max.intknots = 5)
      if (weighted) args$weights <- data$w
      if (identical(fun, GGeDS)) args$family <- gaussian()
      fit <- suppressWarnings(do.call(fun, c(list(formula = y ~ f(x, v)), args)))
      zero <- suppressWarnings(do.call(fun, c(list(formula = y ~ f(x, v) + offset(zero)), args)))
      expect_equal(fit$rss, zero$rss)
      expect_equal(fit$stored, zero$stored)
      for (n in 2:4) {
        expect_equal(predict(fit, data, n = n), predict(zero, data, n = n))
        expect_equal(coef(fit, n = n), coef(zero, n = n))
      }
      legacy <- fit
      legacy$args$offset <- NULL
      expect_equal(predict(legacy, n = 2), predict(fit, n = 2))
    }
  }
})

test_that("bivariate prediction sums formula offsets once and respects new-data order", {
  data <- bivariate_offset_data()
  data$other <- .1 * data$x
  data$y <- data$signal + data$off + data$other + rnorm(nrow(data), sd = .1)
  fit <- suppressWarnings(NGeDS(y ~ f(x, v) + offset(off) + offset(other), data = data,
    phi = .9, q = 1, max.intknots = 5))
  expect_equal(fit$args$offset, data$off + data$other)
  index <- c(91, 13, 48)
  newdata <- data[index, ]
  shifted <- transform(newdata, off = off + .2, other = other - .05)
  for (n in 2:4) {
    pred <- predict(fit, newdata, n = n)
    expect_equal(predict(fit, shifted, n = n) - pred, rep(.15, 3), tolerance = 1e-10)
    expect_equal(pred, predict(fit, data, n = n)[index])
    expect_equal(pred, vapply(seq_len(nrow(newdata)), function(i) {
      predict(fit, newdata[i, , drop = FALSE], n = n)
    }, numeric(1)))
  }
})

test_that("bivariate fitting engines reject malformed offsets", {
  data <- bivariate_offset_data()
  for (fun in list(BivariateFitter, GenBivariateFitter)) {
    for (bad in list(0, rep(NA_real_, nrow(data)), rep(Inf, nrow(data)))) {
      expect_error(fun(X = data$x, Y = data$v, Z = data$signal, W = NULL,
                       offset = bad), "one finite value per observation")
    }
  }
})

test_that("bivariate inverse-Gamma offsets are added after spline clamping", {
  data <- bivariate_offset_data()
  data$y <- 1 / (2 + .2 * data$x + .3 * data$v + data$off) + rnorm(nrow(data), sd = .003)
  fit <- suppressWarnings(GGeDS(y ~ f(x, v) + offset(off), data = data,
    family = Gamma(), weights = w, phi = .9, q = 1, max.intknots = 6))
  newdata <- data.frame(x = c(-.2, .4, 1.2), v = c(1.2, .3, -.2), off = .2)
  clipped <- newdata
  clipped$x <- pmax(0, pmin(1, clipped$x))
  clipped$v <- pmax(0, pmin(1, clipped$v))
  for (n in 2:4) {
    eta <- suppressWarnings(predict(fit, newdata, n = n, type = "link"))
    expect_equal(eta, predict(fit, clipped, n = n, type = "link"))
    expect_equal(suppressWarnings(predict(fit, newdata, n = n)), 1 / eta)
    shifted <- transform(newdata, off = off + .1)
    expect_equal(suppressWarnings(predict(fit, shifted, n = n)), 1 / (eta + .1))
    shifted$off <- -100
    expect_error(suppressWarnings(predict(fit, shifted, n = n)), "positive finite linear predictor")
  }
})
