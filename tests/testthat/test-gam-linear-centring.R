# Linear (order-2) GeDSgam newdata predictions must reproduce the fitted
# local-scoring predictor. Backfitting centres with IRLS weights, so with
# non-constant weights the old intercept, mean(z), shifted newdata predictions
# by a constant on the link scale.

# Only the known "fitted probabilities 0 or 1" warning from the GD binomial
# fit's internal glm object is muffled; any other warning still surfaces.
fit_quiet <- function(...) {
  out <- NULL
  invisible(utils::capture.output(out <- withCallingHandlers(
    NGeDSgam(..., higher_order = FALSE),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })))
  out
}

test_that("weighted Gaussian linear fit matches an independent weighted lm", {
  set.seed(1)
  n <- 100
  d <- data.frame(x = runif(n), w = rep(c(1, 5), each = n / 2))
  d$y <- 1 + 2 * d$x + c(rnorm(n / 2, 2, .3), rnorm(n / 2, 0, .3))
  ref <- stats::lm(y ~ x, data = d, weights = w)
  # The fixture must exercise the old defect: mean(y) differs from the
  # mean of the weighted fit, so the old intercept would have been wrong.
  expect_gt(abs(mean(d$y) - mean(fitted(ref))), 0.1)

  fit <- fit_quiet(y ~ f(x), data = d, weights = d$w, internal_knots = 0L,
                   stoptype = "RD")
  new <- data.frame(x = seq(0.05, 0.95, length.out = 7))
  expect_equal(as.numeric(predict(fit, newdata = d, n = 2)),
               as.numeric(fitted(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(predict(fit, newdata = new, n = 2)),
               as.numeric(predict(ref, newdata = new)), tolerance = 1e-8)
})

expect_newdata_matches_fit <- function(fit, data) {
  expect_equal(as.numeric(predict(fit, newdata = data, n = 2, type = "link")),
               as.numeric(fit$final_model$Y_hat$eta), tolerance = 1e-8)
  expect_equal(as.numeric(predict(fit, newdata = data, n = 2, type = "response")),
               as.numeric(predict(fit, n = 2, type = "response")), tolerance = 1e-8)
}

test_that("binomial GeDSgam linear newdata predictions match the fit", {
  skip_if_not_installed("rpart")
  d <- subset(rpart::kyphosis, Number > 2)
  for (stoptype in c("GD", "RD")) {
    set.seed(123)
    fit <- fit_quiet(Kyphosis ~ f(Age) + f(Start), data = d, family = binomial(),
                     stoptype = stoptype)
    expect_newdata_matches_fit(fit, d)
  }
})

test_that("Poisson GeDSgam linear newdata predictions match the fit", {
  set.seed(1)
  n <- 150
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- rpois(n, exp(1 + sin(4 * d$x1) + d$x2))
  expect_newdata_matches_fit(fit_quiet(y ~ f(x1) + f(x2), data = d, family = poisson()), d)
  expect_newdata_matches_fit(fit_quiet(y ~ f(x1) + x2, data = d, family = poisson()), d)
})

# Training-data link predictions (no newdata) must be the full fitted predictor,
# intercept and formula offset included. Previously the intercept was omitted.
test_that("linear GeDSgam training link predictions include intercept and offset", {
  set.seed(2)
  n <- 200
  d <- data.frame(x = runif(n), off = log(runif(n, 1, 5)))
  d$y <- rpois(n, exp(0.5 + 1.5 * d$x + d$off))
  fit <- fit_quiet(y ~ f(x) + offset(off), data = d, family = poisson(),
                   internal_knots = 0L, stoptype = "RD")
  ref <- stats::glm(y ~ x + offset(off), data = d, family = poisson())

  link <- as.numeric(predict(fit, n = 2, type = "link"))
  # Independent reference: the Poisson GLM. Local scoring stops slightly
  # before full convergence, hence the absolute tolerance.
  expect_lt(max(abs(link - stats::predict(ref))), 0.01)
  expect_equal(exp(link), as.numeric(predict(fit, n = 2, type = "response")),
               tolerance = 1e-10)
})

test_that("linear centring also works with normalization, zero weights and joint smoothers", {
  set.seed(123)
  d <- data.frame(x = runif(90), v = runif(90), w = rep(c(0, 1, 3), 30))
  d$y <- 2 + d$x + .5 * d$v + rnorm(90, sd = .3)
  for (normalize in c(FALSE, TRUE)) {
    for (form in list(y ~ f(x) + v, y ~ f(x, v))) {
      fit <- fit_quiet(form, data = d, weights = d$w, normalize_data = normalize,
                       stoptype = "RD", internal_knots = 2)
      expect_newdata_matches_fit(fit, d)
      index <- c(80, 12, 41)
      expect_equal(predict(fit, d[index, ], n = 2), predict(fit, n = 2)[index],
                   tolerance = 1e-8)
    }
  }
})

test_that("formula offsets enter GeDSgam predictions exactly once at every order", {
  set.seed(123)
  d <- data.frame(x = runif(100), off = runif(100, -.4, .4))
  d$y <- rpois(100, exp(.4 + d$x + d$off))
  fit <- suppressMessages(NGeDSgam(y ~ f(x) + offset(off), data = d,
    family = poisson(), internal_knots = 0, stoptype = "RD", penalty = "none"))
  index <- c(83, 17, 46)
  new <- d[index, c("x", "off")]
  shifted <- transform(new, off = off + log(2))
  for (n in 2:4) {
    eta <- predict(fit, n = n, type = "link")
    mu <- predict(fit, n = n)
    expect_equal(exp(eta), mu, tolerance = 1e-8)
    expect_equal(predict(fit, new, n = n, type = "link"), eta[index], tolerance = 1e-8)
    expect_equal(predict(fit, new, n = n), mu[index], tolerance = 1e-8)
    expect_equal(predict(fit, shifted, n = n), 2 * mu[index], tolerance = 1e-8)
    if (n > 2L) {
      reference <- glm(y ~ poly(x, n - 1L) + offset(off), data = d, family = poisson())
      expect_equal(mu, as.numeric(fitted(reference)), tolerance = 1e-7)
    }
  }
  expect_error(predict(fit, transform(new, off = Inf), n = 2), "finite numeric")
})

test_that("old linear GAM objects retain the documented intercept fallback", {
  expect_equal(GeDS:::.gam_linear_intercept(list(z = c(1, 2, 6))), 3)
  expect_equal(GeDS:::.gam_linear_intercept(list(z = c(1, 2, 6), z_hat = c(2, 3, 7))), 4)
})
