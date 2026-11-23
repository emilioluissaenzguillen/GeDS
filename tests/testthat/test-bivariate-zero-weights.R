biv_zero_weight_data <- function(seed, n = 300) {
  set.seed(seed)
  x <- runif(n); y <- runif(n)
  z <- sin(2 * pi * x) * cos(pi * y) + rnorm(n, sd = 0.2)
  data.frame(x = x, y = y, z = z)
}

biv_zero_weight_summary <- function(fit, newdata) {
  list(knots = knots(fit, n = 2L, options = "internal"),
       pred = as.numeric(predict(fit, newdata = newdata, n = 2L)))
}

expect_same_biv_fit <- function(a, b) {
  expect_equal(a$knots, b$knots)
  expect_equal(a$pred, b$pred, tolerance = 1e-6)
}

test_that("zero-weight rows do not change bivariate NGeDS fits", {
  d <- biv_zero_weight_data(4)
  one <- rep(1, nrow(d))
  pads <- list(
    inside = data.frame(x = 0.5, y = 0.5, z = 100),
    outside = data.frame(x = 1.3, y = -0.2, z = 100),
    duplicate = data.frame(x = d$x[10], y = d$y[10], z = 100)
  )
  for (stoptype in c("SR", "GD")) {
    fit <- function(data, w) suppressWarnings(
      NGeDS(z ~ f(x, y), data = data, weights = w, stoptype = stoptype))
    ref <- biv_zero_weight_summary(fit(d, one), d)
    for (pad in pads) {
      padded <- fit(rbind(d, pad), c(one, 0))
      expect_same_biv_fit(biv_zero_weight_summary(padded, d), ref)
    }
  }
})

test_that("bivariate GD is translation invariant with matched bootstrap draws", {
  d <- biv_zero_weight_data(4, n = 120)
  for (w in list(rep(1, nrow(d)), rep(c(.5, 2), length.out = nrow(d)))) {
    fit_model <- function(data, weights) suppressWarnings(
      NGeDS(z ~ f(x, y), data = data, weights = weights, stoptype = "GD",
            max.intknots = 8L))
    set.seed(123)
    seed <- .Random.seed
    reference <- fit_model(d, w)
    expect_identical(.Random.seed, seed)
    expect_gt(length(unlist(knots(reference, n = 2L, options = "internal"))), 0L)
    shifted_data <- transform(d, z = z + 1e8)
    shifted <- fit_model(shifted_data, w)
    padded_data <- rbind(shifted_data, data.frame(x = .5, y = .5, z = 1e100))
    padded <- fit_model(padded_data, c(w, 0))
    for (fit in list(shifted, padded)) {
      expect_identical(fit$stopinfo$reason, reference$stopinfo$reason)
      for (n in 2:4) {
        expect_equal(knots(fit, n = n), knots(reference, n = n), tolerance = 1e-6)
        difference <- as.numeric(predict(fit, newdata = d, n = n)) - 1e8 -
          as.numeric(predict(reference, newdata = d, n = n))
        expect_lt(max(abs(difference)), 1e-5)
      }
    }
    expect_identical(.Random.seed, seed)
  }
})

test_that("zero-weight rows do not change bivariate GGeDS fits", {
  set.seed(11)
  n <- 300
  d <- data.frame(x = runif(n), y = runif(n))
  d$z <- rpois(n, exp(1 + sin(2 * pi * d$x) * d$y))
  one <- rep(1, n)
  fit <- function(data, w) suppressWarnings(
    GGeDS(z ~ f(x, y), data = data, family = poisson(), weights = w))
  ref <- biv_zero_weight_summary(fit(d, one), d)
  for (pad in list(data.frame(x = 0.5, y = 0.5, z = 60),
                   data.frame(x = 1.3, y = -0.2, z = 60))) {
    padded <- fit(rbind(d, pad), c(one, 0))
    expect_same_biv_fit(biv_zero_weight_summary(padded, d), ref)
  }
})
