test_that("knot placement is stable under response translations", {
  set.seed(123)
  x <- seq(0, 1, length.out = 100)
  repeated_x <- rep(seq(0, 1, length.out = 50), each = 2)
  cases <- list(
    smooth = data.frame(x = x, y = sin(10 * x) + rnorm(100, sd = .1)),
    piecewise = data.frame(x = x, y = 1 + x + 2 * pmax(x - .4, 0) - pmax(x - .7, 0)),
    repeated = data.frame(x = repeated_x, y = sin(10 * repeated_x) + rnorm(100, sd = .1)),
    weighted = data.frame(x = x, y = sin(10 * x) + rnorm(100, sd = .1),
                          w = sample(c(0, .5, 1, 3), 100, replace = TRUE)),
    parametric_offset = data.frame(x = x, z = rnorm(100), off = .3 * sin(2 * x))
  )
  cases$parametric_offset$y <- with(cases$parametric_offset,
    sin(10 * x) + 2 * z + off + rnorm(100, sd = .1))

  for (case in names(cases)) {
    dat <- cases[[case]]
    formula <- if (case == "parametric_offset") y ~ f(x) + z + offset(off) else y ~ f(x)
    fit_model <- function(dat) {
      args <- list(formula = formula, data = dat, beta = .6, phi = .995, q = 3,
                   stoptype = "RD")
      if (case == "weighted") args$weights <- dat$w
      suppressWarnings(do.call(NGeDS, args))
    }
    reference <- fit_model(dat)
    for (shift in c(-mean(dat$y), 1e6)) {
      shifted <- dat
      shifted$y <- dat$y + shift
      fit <- fit_model(shifted)
      for (n in 2:4) {
        k1 <- as.numeric(unlist(knots(reference, n = n, options = "internal")))
        k2 <- as.numeric(unlist(knots(fit, n = n, options = "internal")))
        expect_length(k2, length(k1))
        expect_equal(k2, k1, tolerance = 1e-7, info = case)
        p1 <- as.numeric(predict(reference, newdata = dat, n = n))
        p2 <- as.numeric(predict(fit, newdata = dat, n = n)) - shift
        expect_lt(max(abs(p1 - p2)), 1e-6)
      }
    }
  }
})

test_that("zero-weight responses do not inflate the geometric tolerance", {
  set.seed(123)
  dat <- data.frame(x = seq(0, 1, length.out = 80))
  dat$y <- sin(10 * dat$x) + rnorm(80, sd = .1)
  w <- rep(c(0, .5, 1, 3), length.out = nrow(dat))
  reference <- suppressWarnings(NGeDS(y ~ f(x), data = dat, weights = w,
                                      beta = .6, phi = .9, q = 3))
  dat$y[w == 0] <- 1e15
  fit <- suppressWarnings(NGeDS(y ~ f(x), data = dat, weights = w,
                                beta = .6, phi = .9, q = 3))
  for (n in 2:4) {
    expect_equal(knots(fit, n = n, options = "internal"),
                 knots(reference, n = n, options = "internal"), tolerance = 1e-10)
    expect_equal(predict(fit, newdata = dat[w > 0, ], n = n),
                 predict(reference, newdata = dat[w > 0, ], n = n), tolerance = 1e-10)
  }
})

test_that("numerically exact fits stop without normalizing zero cluster means", {
  x <- seq(0, 1, length.out = 100)
  for (y in list(rep(0, 100), rep(5, 100), 2 + 3 * x)) {
    for (shift in c(0, -mean(y), 1e6)) {
      dat <- data.frame(x = x, y = y + shift)
      fit <- suppressWarnings(NGeDS(y ~ f(x), data = dat, beta = .6, phi = .995, q = 3))
      expect_equal(fit$iters, 1)
      expect_equal(fit$Nintknots, 0)
      expect_true(all(is.finite(fit$rss)))
      for (n in 2:4) {
        expect_length(unlist(knots(fit, n = n, options = "internal")), 0)
        expect_lt(max(abs(as.numeric(predict(fit, n = n)) - dat$y)), 1e-7)
      }
      expect_equal(fit$dev.linear,
                   sum((dat$y - as.numeric(predict(fit, n = 2)))^2), tolerance = 1e-12)
    }
  }
})

test_that("a unit Gaussian boosting update agrees with a single GeDS learner", {
  set.seed(123)
  dat <- data.frame(x = seq(0, 1, length.out = 100))
  dat$y <- 3 + sin(10 * dat$x) + rnorm(100, sd = .1)
  # Stop before saturation to isolate residual geometry from the separate
  # knot-history row selection used when no further insertion is admissible.
  plain <- suppressWarnings(NGeDS(y ~ f(x), data = dat, beta = .6, phi = .9, q = 3,
                                  stoptype = "RD"))
  invisible(capture.output(boost <- suppressWarnings(NGeDSboost(
    y ~ f(x), data = dat, beta = .6, phi = .9, q = 3,
    initial_learner = FALSE, shrinkage = 1, max_iterations = 1, stoptype = "RD"
  ))))
  for (n in 2:4) {
    expect_equal(as.numeric(unlist(knots(boost, n = n, options = "internal"))),
                 as.numeric(unlist(knots(plain, n = n, options = "internal"))),
                 tolerance = 1e-7)
    expect_lt(max(abs(as.numeric(predict(boost, newdata = dat, n = n)) -
                       as.numeric(predict(plain, newdata = dat, n = n)))), 1e-6)
  }
})
