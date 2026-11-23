test_that("GGeDS ignores zero-weight boundary rows when selecting knots", {
  set.seed(610)
  d <- data.frame(x = seq(0, 1, length.out = 60), z = rnorm(60), y = rnorm(60))
  for (side in c("left", "right")) {
    w <- if (side == "right") c(rep(1, 40), rep(0, 20)) else c(rep(0, 20), rep(1, 40))
    fit <- suppressWarnings(GGeDS(y ~ f(x), data = d, weights = w, phi = .8, stoptype = "SR"))
    reference <- suppressWarnings(GGeDS(y ~ f(x), data = d[w > 0, ], phi = .8, stoptype = "SR"))
    expect_equal(fit$args$extr, range(d$x[w > 0]))
    expect_equal(knots(fit, n = 2), knots(reference, n = 2))
    for (n in 2:4) {
      expect_equal(predict(fit, newdata = d[w > 0, ], n = n), predict(reference, n = n),
                   tolerance = 1e-8)
      expect_length(suppressWarnings(predict(fit, n = n)), nrow(d))
      expect_true(all(is.finite(suppressWarnings(predict(fit, n = n)))))
    }
  }
})

test_that("zero-weight rows do not alter generalized fits with offsets and repeated sites", {
  families <- list(gaussian(), poisson(), binomial(), Gamma(link = "inverse"), Gamma(link = "log"))
  fit_names <- c("linear.fit", "quadratic.fit", "cubic.fit")
  for (family in families) {
    set.seed(620)
    d <- data.frame(x = rep(seq(0, 1, length.out = 45), each = 2), z = rnorm(90))
    d$o <- .02 * d$z
    eta <- 1.5 + .2 * sin(5 * d$x) + .04 * d$z + d$o
    mu <- family$linkinv(eta)
    d$y <- switch(family$family,
      gaussian = rnorm(90, mu, .1), poisson = rpois(90, mu),
      binomial = rbinom(90, 1, mu), Gamma = rgamma(90, shape = 12, scale = mu / 12))
    w <- rep(c(1, 2, 3), length.out = 90)
    w[c(1:10, 81:90)] <- 0
    fit <- suppressWarnings(GGeDS(y ~ f(x) + z + offset(o), family = family,
                                 data = d, weights = w, phi = .8, stoptype = "SR"))
    reference <- suppressWarnings(GGeDS(y ~ f(x) + z + offset(o), family = family,
                                       data = d[w > 0, ], weights = w[w > 0], phi = .8, stoptype = "SR"))
    for (n in 2:4) {
      expect_equal(knots(fit, n = n), knots(reference, n = n))
      expect_equal(coef(fit, n = n, onlySpline = FALSE),
                   coef(reference, n = n, onlySpline = FALSE), tolerance = 1e-7)
      expect_equal(predict(fit, newdata = d[w > 0, ], n = n),
                   predict(reference, n = n), tolerance = 1e-7)
      expect_equal(as.numeric(fit[[fit_names[n - 1]]]$predicted),
                   as.numeric(suppressWarnings(predict(fit, newdata = d, n = n))), tolerance = 1e-7)
      expect_true(all(is.finite(suppressWarnings(predict(fit, n = n)))))
    }
  }
})

test_that("explicit GGeDS boundaries are retained without counting inactive support", {
  set.seed(621)
  d <- data.frame(x = seq(0, 1, length.out = 80))
  d$y <- sin(5 * d$x) + rnorm(80, sd = .1)
  w <- c(rep(1, 60), rep(0, 20))
  fit <- suppressWarnings(GGeDS(y ~ f(x), data = d, weights = w,
                               Xextr = c(-.1, 1.1), phi = .8, stoptype = "SR"))
  reference <- suppressWarnings(GGeDS(y ~ f(x), data = d[w > 0, ],
                                     Xextr = c(-.1, 1.1), phi = .8, stoptype = "SR"))
  expect_equal(fit$args$extr, c(-.1, 1.1))
  for (n in 2:4) {
    expect_equal(knots(fit, n = n), knots(reference, n = n))
    expect_equal(predict(fit, newdata = d[w > 0, ], n = n),
                 predict(reference, n = n), tolerance = 1e-8)
    expect_equal(as.numeric(predict(fit, n = n)),
      as.numeric(predict(fit, newdata = d, n = n)), tolerance = 1e-8)
  }
  d$y[w == 0] <- d$y[w == 0] + 100
  altered <- suppressWarnings(GGeDS(y ~ f(x), data = d, weights = w,
                                     Xextr = c(-.1, 1.1), phi = .8, stoptype = "SR"))
  for (n in 2:4) {
    expect_equal(knots(altered, n = n), knots(fit, n = n))
    expect_equal(as.numeric(predict(altered, n = n)),
                 as.numeric(predict(fit, n = n)), tolerance = 1e-8)
  }
})

test_that("GGeDS diagnoses weighted parametric aliasing", {
  set.seed(622)
  d <- data.frame(x = seq(0, 1, length.out = 60), y = rnorm(60))
  d$z <- d$x
  expect_error(suppressWarnings(GGeDS(y ~ f(x) + z, data = d)),
               "weighted design matrix is rank deficient")
})
