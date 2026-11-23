test_that("fast spline fits preserve coefficient order after QR pivoting", {
  for (large in c(FALSE, TRUE)) {
    x <- seq(0, 1, length.out = if (large) 201 else 21)
    internal <- if (large) c(seq(.03, .97, length.out = 50), .2101, .2102, .2103) else c(.21, .22, .23)
    offset <- .1 * x
    y <- 1 + 2 * x + sin(4 * x) + offset
    for (weighted in c(FALSE, TRUE)) {
      weights <- if (weighted) rep(c(1, 2, 0), length.out = length(x)) else rep(1, length(x))
      fit <- suppressWarnings(SplineReg_LM(
        X = x, Y = y, offset = offset, weights = weights,
        InterKnots = internal, n = 2L, fast = TRUE
      ))
      basis <- as.matrix(fit$basis)
      reference <- stats::lm.wfit(basis, y - offset, weights)
      theta <- reference$coefficients
      theta[is.na(theta)] <- 0
      expected <- as.numeric(basis %*% theta) + offset

      expect_lt(reference$rank, ncol(basis))
      expect_equal(as.numeric(fit$theta), unname(theta), tolerance = 1e-10)
      expect_equal(fit$predicted, expected, tolerance = 1e-10)
      expect_equal(fit$residuals, y - expected, tolerance = 1e-10)
      expect_equal(fit$rss, sum(weights * (y - expected)^2), tolerance = 1e-10)
    }
  }
})

test_that("light weighted regression restores pivoted coefficients and residuals", {
  x <- seq(0, 1, length.out = 30)
  basis <- cbind(1, x, x, x^2)
  y <- 1 + x + sin(3 * x)
  weights <- rep(c(0, .5, 2), length.out = length(x))
  fit <- GeDS:::lm.wfit.light(basis, y, weights)
  reference <- stats::lm.wfit(basis, y, weights)
  theta <- reference$coefficients
  theta[is.na(theta)] <- 0

  expect_true(fit$pivoted)
  expect_equal(unname(fit$coefficients), unname(theta), tolerance = 1e-12)
  expect_equal(as.numeric(fit$residuals), as.numeric(y - basis %*% theta), tolerance = 1e-12)
})

test_that("full-rank QR coefficients remain unchanged", {
  x <- seq(0, 1, length.out = 30)
  basis <- cbind(1, x, x^2)
  y <- sin(3 * x)
  raw <- stats::.lm.fit(basis, y)
  fit <- GeDS:::lm.wfit.light(basis, y, rep(1, length(x)))

  expect_false(raw$pivoted)
  expect_identical(fit$coefficients, raw$coefficients)
})

test_that("GD path regression solves the tridiagonal linear-spline system directly", {
  x <- seq(0, 1, length.out = 201)
  internal <- seq(.02, .98, length.out = 50)
  y <- sin(9 * x) + .1 * cos(80 * x)
  w <- rep(c(0, 1, 2), length.out = length(x))
  ref <- SplineReg_LM(x, y, weights = w, InterKnots = internal, n = 2L)

  local_mocked_bindings(lm.wfit.light = function(...) stop("Unexpected dense fallback"),
                        .package = "GeDS")
  fit <- GeDS:::.spline_reg_lm_fast(x, y, NULL, rep(0, length(x)), w,
                                    internal, 2L, range(x), gd_path = TRUE)
  expect_null(fit$basis)
  expect_equal(fit$predicted, as.numeric(ref$predicted), tolerance = 1e-10)
  expect_equal(fit$rss, ref$rss, tolerance = 1e-10)
})

test_that("GD path regression falls back to QR for degenerate knots", {
  x <- seq(0, 1, length.out = 101)
  y <- sin(6 * x)
  w <- rep(1, length(x))
  repeated <- c(.3, .3, .7)
  expect_null(GeDS:::linear_spline_ls(x, y, w, c(0, repeated, 1)))
  fit <- GeDS:::.spline_reg_lm_fast(x, y, NULL, rep(0, length(x)), w,
                                    repeated, 2L, range(x), gd_path = TRUE)
  ref <- SplineReg_LM(x, y, weights = w, InterKnots = repeated, n = 2L)
  expect_false(is.null(fit$basis))
  expect_equal(fit$rss, ref$rss, tolerance = 1e-10)

  # Three knots inside one gap between data sites leave the middle hat
  # function without observations, so the system is singular.
  empty <- c(.501, .502, .503)
  expect_null(GeDS:::linear_spline_ls(x, y, w, c(0, empty, 1)))
})
