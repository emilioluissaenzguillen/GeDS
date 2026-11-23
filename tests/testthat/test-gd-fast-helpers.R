test_that("fast 95% quantiles match stats::quantile exactly", {
  set.seed(5)
  for (r in seq_len(300)) {
    n <- sample(c(1:5, 10, 11, 20, 100, 101), 1)
    v <- if (r %% 3 == 0) round(rnorm(n), 1) else if (r %% 3 == 1) rexp(n) else rep(2, n)
    expect_identical(GeDS:::.gd_quantile95(v), unname(stats::quantile(v, .95)), info = r)
    W <- matrix(if (r %% 2) round(rnorm(4 * n), 1) else rnorm(4 * n), nrow = 4)
    expect_identical(GeDS:::.gd_quantile95_rows(W),
                     apply(W, 1L, function(z) unname(stats::quantile(z, .95))), info = r)
  }
  expect_error(GeDS:::.gd_quantile95(c(1, NA)), "missing values")
  expect_error(GeDS:::.gd_quantile95_rows(matrix(c(1, NA, 2, 3), 2)), "missing values")
})

test_that("neighbour semivariance matches the R expression exactly", {
  set.seed(6)
  for (r in seq_len(100)) {
    x <- sort(if (r %% 2) round(runif(sample(c(4, 12, 300), 1)), 2) else runif(sample(c(4, 12, 300), 1)))
    geometry <- GeDS:::.gamma_dynamic_neighbours(x)
    res <- rnorm(length(x)) * 10^sample(-8:3, 1)
    nn <- geometry$nn
    expected <- colMeans((res - matrix(res[nn], nrow = length(x), byrow = TRUE))^2) / 2
    expect_identical(GeDS:::gd_neighbour_semivariance(res, nn), expected, info = r)
  }
})

test_that("the R references match linear_spline_ls and gd_neighbour_semivariance", {
  set.seed(8)
  for (r in seq_len(30)) {
    n <- sample(c(20, 200), 1)
    x <- sort(runif(n))
    y <- sin(6 * x) + rnorm(n, sd = 0.2)
    w <- if (r %% 2) rep(1, n) else rexp(n)
    t <- c(0, sort(runif(sample(1:6, 1), 0.1, 0.9)), 1)
    fast <- GeDS:::linear_spline_ls(x, y, w, t)
    ref <- GeDS:::.linear_spline_ls_R(x, y, w, t)
    if (is.null(fast)) next
    expect_equal(fast$theta, ref$theta, tolerance = 1e-8, info = r)
    expect_equal(fast$fitted, ref$fitted, tolerance = 1e-8, info = r)
    nn <- GeDS:::.gamma_dynamic_neighbours(x)$nn
    expect_identical(GeDS:::gd_neighbour_semivariance(y, nn),
                     GeDS:::.gd_neighbour_semivariance_R(y, nn), info = r)
  }
  expect_null(GeDS:::.linear_spline_ls_R(1:5 / 5, 1:5, rep(1, 5), c(0, 0.5, 0.5, 1)))
})
