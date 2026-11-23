test_that("default spline caps balance positive-weight rows and distinct support", {
  x <- rep(seq_len(20), each = 2)
  w <- as.numeric(x <= 12)
  cap <- GeDS:::.geds_spline_cap(x, w)
  expect_equal(cap$max.coef, 12)
  expect_equal(cap$max.intknots, 10L)
  expect_equal(GeDS:::.geds_spline_cap(x, w, 1L)$max.intknots, 9L)
  expect_equal(GeDS:::.geds_spline_cap(x, 3 * w), cap)
  expect_equal(GeDS:::.geds_spline_cap(x[w > 0], w[w > 0]), cap)
  expect_equal(GeDS:::.geds_spline_cap(1:3, rep(1, 3))$max.intknots, 0L)
  expect_equal(GeDS:::.geds_spline_cap(1:40, rep(1, 40))$max.intknots, 18L)
  repeated <- rep(1:5, each = 10)
  expect_equal(GeDS:::.geds_spline_cap(repeated, rep(1, 50))$max.intknots, 3L)
  expect_equal(GeDS:::.geds_spline_cap(repeated, rep(1, 50), 1L)$max.intknots, 3L)
  expect_equal(GeDS:::.geds_spline_cap(1:3, rep(1, 3), 2L)$max.coef, 4L)
  iris_x <- iris$Petal.Width[iris$Species != "virginica"]
  expect_equal(GeDS:::.geds_spline_cap(iris_x, rep(1, 100))$max.intknots, 13L)

  grid <- expand.grid(x = 1:10, z = 1:10)
  expect_error(GeDS:::.geds_spline_cap(grid, rep(1, 100)), "univariate only")
})

test_that("standalone fits use the revised cap with repeated covariates", {
  set.seed(614)
  d <- data.frame(x = rep(seq(0, 1, length.out = 15), each = 6))
  d$y <- sin(d$x) + rnorm(nrow(d), sd = .1)
  for (fun in list(NGeDS, GGeDS)) {
    fit <- suppressWarnings(fun(y ~ f(x), data = d, phi = .8))
    expect_equal(fit$args$max.intknots, 13L)
    for (n in 2:4) expect_true(all(is.finite(predict(fit, n = n))))
  }
})

test_that("bivariate cap checks the prospective tensor product on each axis", {
  flags <- GeDS:::.geds_bivariate_cap_flags(FALSE, FALSE, 3L, 1L, 0L, 18)
  expect_false(flags$flagX)
  expect_true(flags$flagY)
  expect_true(flags$capped)
  flags <- GeDS:::.geds_bivariate_cap_flags(FALSE, FALSE, 5L, 5L, 0L, 50)
  expect_true(flags$flagX)
  expect_true(flags$flagY)
  flags <- GeDS:::.geds_bivariate_cap_flags(FALSE, FALSE, 5L, 5L, 0L, Inf)
  expect_false(flags$capped)
})

test_that("univariate defaults include parametric columns and allow overrides", {
  set.seed(610)
  d <- data.frame(x = seq(0, 1, length.out = 60), z = rnorm(60), y = rnorm(60))
  for (fun in list(NGeDS, GGeDS)) {
    capped <- suppressWarnings(fun(y ~ f(x), data = d, phi = 1, q = 1))
    expect_equal(capped$args$max.intknots, 28L)
    expect_lte(length(coef(capped, n = 2, onlySpline = FALSE)), 30)
    parametric <- suppressWarnings(fun(y ~ f(x) + z, data = d, phi = 1, q = 1))
    expect_equal(parametric$args$max.intknots, 27L)
    expect_lte(length(coef(parametric, n = 2, onlySpline = FALSE)), 30)
    explicit <- suppressWarnings(fun(y ~ f(x), data = d, max.intknots = 40, phi = 1, q = 1))
    expect_equal(explicit$args$max.intknots, 40L)
    w <- rep(1, 60)
    w[seq(2, 59, by = 3)] <- 0
    weighted <- suppressWarnings(fun(y ~ f(x), data = d, weights = w, phi = .8))
    expect_equal(weighted$args$max.intknots, 18L)
  }
})

test_that("bivariate defaults retain the previous knot limit without a coefficient cap", {
  set.seed(611)
  d <- data.frame(x = runif(80), z = runif(80), a = rnorm(80), y = rnorm(80))
  for (fun in list(NGeDS, GGeDS)) {
    fit <- suppressMessages(suppressWarnings(fun(
      y ~ f(x, z) + a, data = d, phi = .8, q = 1, stoptype = "RD"
    )))
    expect_equal(fit$args$max.coef, Inf)
    expect_equal(fit$args$max.intknots, length(unique(as.matrix(d[c("x", "z")]))) - 3L)
    expect_true(all(is.finite(predict(fit, n = 2))))
    explicit <- suppressMessages(suppressWarnings(fun(
      y ~ f(x, z), data = d, max.intknots = 8, phi = .8
    )))
    expect_equal(explicit$args$max.coef, Inf)
    expect_equal(explicit$args$max.intknots, 8L)
  }
})

test_that("nonbinding caps retain original predictions for every spline order", {
  set.seed(612)
  d <- data.frame(x = seq(-2, 2, length.out = 100))
  d$y <- sin(d$x) + rnorm(100, sd = .1)
  for (fun in list(NGeDS, GGeDS)) {
    set.seed(613)
    fit <- suppressWarnings(fun(y ~ f(x), data = d, phi = .8))
    set.seed(613)
    old_limit <- suppressWarnings(fun(y ~ f(x), data = d, phi = .8, max.intknots = 98))
    for (n in 2:4) expect_equal(predict(fit, n = n), predict(old_limit, n = n))
  }
})

test_that("GAM and boosting cached paths use per-smoother defaults", {
  set.seed(613)
  d <- data.frame(x = seq(0, 1, length.out = 40), z = runif(40), y = rnorm(40))
  capture_limits <- function(expr) {
    limits <- integer()
    original <- GeDS:::UnivariateFitter
    local_mocked_bindings(UnivariateFitter = function(...) {
      limits <<- c(limits, list(...)$max.intknots)
      original(...)
    }, .package = "GeDS")
    suppressMessages(suppressWarnings(force(expr)))
    limits
  }
  gam_limits <- capture_limits(NGeDSgam(y ~ f(x) + f(z), data = d, phi = .8, stoptype = "RD"))
  expect_gt(length(gam_limits), 0)
  expect_true(all(gam_limits == 18L))
  explicit <- capture_limits(NGeDSgam(y ~ f(x), data = d, internal_knots = 25, phi = .8, stoptype = "RD"))
  expect_true(all(explicit == 25L))
  gam_gd_limits <- capture_limits(NGeDSgam(y ~ f(x) + f(z), data = d, higher_order = FALSE))
  expect_gt(length(gam_gd_limits), 0)
  expect_true(all(gam_gd_limits <= 18L))
  boost_limits <- capture_limits(NGeDSboost(
    y ~ f(x) + f(z), data = d, initial_learner = FALSE, max_iterations = 2,
    phi = .8, higher_order = FALSE
  ))
  expect_gt(length(boost_limits), 0)
  expect_true(all(boost_limits <= 18L))
  explicit <- capture_limits(NGeDSboost(
    y ~ f(x), data = d, initial_learner = FALSE, max_iterations = 2,
    int.knots_boost = 25, phi = .8, higher_order = FALSE
  ))
  expect_true(all(explicit <= 18L))
  repeated <- data.frame(x = rep(seq(0, 1, length.out = 15), each = 6))
  repeated$y <- sin(repeated$x) + rnorm(nrow(repeated), sd = .1)
  gam_limits <- capture_limits(NGeDSgam(y ~ f(x), data = repeated, phi = .8, stoptype = "RD"))
  expect_gt(length(gam_limits), 0)
  expect_true(all(gam_limits == 13L))
  boost_limits <- capture_limits(NGeDSboost(
    y ~ f(x), data = repeated, initial_learner = FALSE, max_iterations = 2,
    phi = .8, higher_order = FALSE
  ))
  expect_gt(length(boost_limits), 0)
  expect_true(all(boost_limits <= 13L))
})
