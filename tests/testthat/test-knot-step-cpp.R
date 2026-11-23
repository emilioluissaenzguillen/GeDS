test_that("C++ Stage A knot step matches the R reference implementation", {
  set.seed(42)
  for (r in seq_len(400)) {
    n <- sample(c(6, 15, 60, 400), 1)
    X <- sort(if (r %% 3 == 0) round(runif(n), 1) else runif(n))
    weights <- if (r %% 4 == 0) replace(rexp(n), sample(n, n %/% 3), 0) else rep(1, n)
    if (!any(weights > 0)) weights[1] <- 1
    positive <- weights > 0
    xr <- rle(X[positive])
    repeated <- any(xr$lengths > 1L)
    residuals <- rnorm(n) + if (r %% 2) sin(9 * X) else 0
    if (r %% 5 == 0) residuals[sample(n, n %/% 4)] <- 1e-18   # numerical zeros
    if (r %% 17 == 0) residuals[] <- 0                        # nothing left to place
    scale <- if (repeated) GeDS:::makeNewRes(weights[positive], as.numeric(xr$lengths)) else
      weights[positive]
    geometry_tol <- 32 * .Machine$double.eps * max(abs(residuals[positive]))
    extr <- range(X[positive])
    intknots <- if (r %% 6 == 0) NULL else runif(sample(0:8, 1), extr[1], extr[2])
    beta <- c(.1, .5, .9)[r %% 3 + 1]
    args <- list(residuals, positive, weights[positive], xr$lengths, repeated, scale,
                 geometry_tol, xr$values, beta, extr, as.numeric(intknots), 1e-12, 2L)
    expected <- do.call(GeDS:::.geds_knot_step_R, args)
    expect_identical(do.call(GeDS:::geds_knot_step, args), expected,
                     info = paste("case", r))
  }
})
