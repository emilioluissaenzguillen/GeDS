# library("GeDS")
#
# ## Rebuild the basis on the stored knots and check that theta_f reproduces f,
# ## i.e. that the Greville-based coefficient mapping is correct.
# expect_theta_f_reproduces_f <- function(knots, n, x, theta_f, f, tol = 1e-6) {
#   B <- splines::splineDesign(knots, x = x, ord = n,
#                              derivs = rep(0, length(x)), outer.ok = TRUE)
#   testthat::expect_equal(as.numeric(B %*% theta_f), f, tolerance = tol)
# }
#
# ## --- univariate Normal (NGeDS) ----------------------------------------------
#
# test_that("separateLinear.GeDS gives an exact, identifiable decomposition", {
#   set.seed(1)
#   N <- 300
#   X <- sort(runif(N, -2, 2))
#   Y <- 2 * X + cos(2 * X) + rnorm(N, sd = 0.2)          # cos(2x) is ~orthogonal to x
#   Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995)
#   d <- separateLinear(Gmod, n = 3)
#
#   expect_s3_class(d, "separateLinear")
#   expect_identical(d$scale, "response")
#
#   # (a) fit unchanged
#   expect_equal(as.numeric(predict(Gmod, n = 3)), d$fitted_spline, tolerance = 1e-6)
#   # (b) reconstruction
#   expect_equal(d$beta0 + d$beta1 * X + d$f, d$fitted_spline, tolerance = 1e-8)
#   # (c) f orthogonal to {1, X}
#   expect_equal(sum(d$f), 0, tolerance = 1e-6)
#   expect_equal(sum(X * d$f), 0, tolerance = 1e-6)
#   # (d) theta_f reproduces f  (the Greville mapping)
#   expect_theta_f_reproduces_f(d$knots, d$n, X, d$theta_f, d$f)
#   # (e) beta1 is exactly the OLS slope of the fitted curve on X ...
#   expect_equal(d$beta1, unname(coef(lm(d$fitted_spline ~ X))[2]), tolerance = 1e-6)
#   # ... and recovers the true slope (looser; tune if flaky)
#   expect_equal(d$beta1, 2, tolerance = 0.2)
# })
#
# test_that("separateLinear.GeDS works for each order n", {
#   set.seed(2)
#   N <- 250
#   X <- sort(runif(N, -2, 2))
#   Y <- 1.5 * X + sin(2 * X) + rnorm(N, sd = 0.2)
#   Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.99)
#   for (nn in 2:4) {
#     d <- separateLinear(Gmod, n = nn)
#     expect_theta_f_reproduces_f(d$knots, d$n, X, d$theta_f, d$f)
#     expect_equal(sum(X * d$f), 0, tolerance = 1e-6)
#   }
# })
#
# ## --- univariate GLM (GGeDS) -------------------------------------------------
#
# test_that("separateLinear.GeDS handles GLM fits on the link scale", {
#   set.seed(3)
#   N <- 400
#   X <- sort(runif(N, -2, 2))
#   eta <- 2 + 0.5 * X + 0.6 * cos(2 * X)                 # linear + curvature on the link scale
#   Y <- rpois(N, exp(eta))
#   Gmod <- GGeDS(Y ~ f(X), family = poisson(), beta = 0.6, phi = 0.99)
#   d <- separateLinear(Gmod, n = 3)
#
#   expect_identical(d$scale, "link")
#   # fitted_spline is the link-scale curve
#   expect_equal(as.numeric(predict(Gmod, n = 3, type = "link")),
#                d$fitted_spline, tolerance = 1e-6)
#   # algebraic identities hold regardless of fit quality
#   expect_equal(sum(d$f), 0, tolerance = 1e-6)
#   expect_equal(sum(X * d$f), 0, tolerance = 1e-6)
#   expect_theta_f_reproduces_f(d$knots, d$n, X, d$theta_f, d$f)
#   expect_gt(sqrt(sum(d$f^2)), 0)                        # genuine curvature recovered
#   # log-scale slope is sensible (loose)
#   expect_equal(d$beta1, 0.5, tolerance = 0.2)
# })
#
# ## --- additive: boosting & GAM -----------------------------------------------
#
# test_that("separateLinear.GeDSboost decomposes each univariate learner", {
#   set.seed(4)
#   N <- 300
#   x1 <- runif(N, -2, 2); x2 <- runif(N, -2, 2)
#   y  <- 1.5 * x1 + (0.5 * x2 + cos(2 * x2)) + rnorm(N, sd = 0.2)
#   dat <- data.frame(y, x1, x2)
#
#   b <- NGeDSboost(y ~ f(x1) + f(x2), data = dat)
#   d <- separateLinear(b, n = 3)
#
#   expect_identical(d$scale, "response")
#   expect_setequal(names(d$learners), c("f(x1)", "f(x2)"))
#   for (nm in names(d$learners)) {
#     l <- d$learners[[nm]]; x <- dat[[l$variable]]
#     expect_equal(sum(x * l$f), 0, tolerance = 1e-6)
#     expect_theta_f_reproduces_f(l$knots, l$n, x, l$theta_f, l$f)
#   }
#   # x1 enters linearly (small curvature), x2 nonlinearly (large curvature)
#   expect_lt(sqrt(sum(d$learners[["f(x1)"]]$f^2)),
#             sqrt(sum(d$learners[["f(x2)"]]$f^2)))
#   expect_equal(d$learners[["f(x1)"]]$beta1, 1.5, tolerance = 0.3)
# })
#
# test_that("separateLinear.GeDSboost works with a single base-learner", {
#   # Guards the single-learner coefficient-name-prefix path.
#   set.seed(10)
#   N <- 250
#   X <- sort(runif(N, -2, 2))
#   Y <- 1.2 * X + cos(2 * X) + rnorm(N, sd = 0.2)
#   b <- NGeDSboost(Y ~ f(X), data = data.frame(X, Y))
#   d <- separateLinear(b, n = 3)
#   expect_length(d$learners, 1)
#   l <- d$learners[[1]]
#   expect_theta_f_reproduces_f(l$knots, l$n, X, l$theta_f, l$f)
# })
#
# test_that("separateLinear.GeDSgam decomposes each univariate smoother", {
#   set.seed(5)
#   N <- 300
#   x1 <- runif(N, -2, 2); x2 <- runif(N, -2, 2)
#   y  <- 1.5 * x1 + (0.5 * x2 + cos(2 * x2)) + rnorm(N, sd = 0.2)
#   dat <- data.frame(y, x1, x2)
#
#   g <- NGeDSgam(y ~ f(x1) + f(x2), data = dat, phi = 0.9)
#   d <- separateLinear(g, n = 3)
#
#   expect_identical(d$scale, "response")
#   for (nm in names(d$learners)) {
#     l <- d$learners[[nm]]; x <- dat[[l$variable]]
#     expect_equal(sum(x * l$f), 0, tolerance = 1e-6)
#     expect_theta_f_reproduces_f(l$knots, l$n, x, l$theta_f, l$f)
#   }
#   expect_lt(sqrt(sum(d$learners[["f(x1)"]]$f^2)),
#             sqrt(sum(d$learners[["f(x2)"]]$f^2)))
# })
#
# ## --- guards / errors --------------------------------------------------------
#
# test_that("X + f(X) is rejected as non-identifiable", {
#   set.seed(6)
#   N <- 200
#   X <- sort(runif(N, -2, 2))
#   Y <- 2 * X + cos(2 * X) + rnorm(N, sd = 0.2)
#   expect_error(NGeDS(Y ~ X + f(X), data = data.frame(X, Y), beta = 0.6, phi = 0.99),
#                "identifiable")
# })
#
# test_that("separateLinear rejects bivariate GeDS fits", {
#   mod <- structure(list(type = "LM - Biv"), class = "GeDS")
#   expect_error(separateLinear(mod, n = 3), "univariate")
# })
#
# test_that("separateLinear rejects normalize_data = TRUE for boosting", {
#   set.seed(8)
#   N <- 200
#   x1 <- runif(N, -2, 2); x2 <- runif(N, -2, 2)
#   y  <- 1.5 * x1 + cos(2 * x2) + rnorm(N, sd = 0.2)
#   # fit exists only to test rejection; its internal CI warnings are incidental
#   b <- NGeDSboost(y ~ f(x1) + f(x2), data = data.frame(y, x1, x2),
#                   normalize_data = TRUE)
#
#   expect_error(separateLinear(b, n = 3), "normalize_data")
# })
#
# ## --- print ------------------------------------------------------------------
#
# test_that("print.separateLinear works", {
#   set.seed(9)
#   N <- 200
#   X <- sort(runif(N, -2, 2))
#   Y <- 2 * X + cos(2 * X) + rnorm(N, sd = 0.2)
#   d <- separateLinear(NGeDS(Y ~ f(X), beta = 0.6, phi = 0.99), n = 3)
#   expect_output(print(d), "decomposition")
#   expect_invisible(print(d))
# })
#
