separated_stageB <- function(quasi = FALSE, n = 4L) {
  x <- if (quasi) rep(c(-1, 0, 1), each = 30) else seq(-1, 1, length.out = 90)
  y <- as.numeric(x > 0)
  if (quasi) y[x == 0] <- rep(0:1, 15)
  d <- data.frame(x = x)
  GeDS:::.stageB_penalised_fit(d, y, d[FALSE], rep(0, 90), rep(1, 90),
    list("f(x)" = list(type = "GeDS", variables = "x")),
    list("f(x)" = numeric()), n, binomial(), NULL, "component")
}

test_that("separated binomial REML refits use an explicit Firth fallback", {
  skip_if_not_installed("mgcv")
  for (quasi in c(FALSE, TRUE)) for (n in 3:4) {
    expect_warning(fit <- withCallingHandlers(separated_stageB(quasi, n),
      warning = function(w) {
        if (grepl("brglmFit: fitted probabilities numerically 0 or 1", conditionMessage(w)))
          invokeRestart("muffleWarning")
      }), class = "geds_stageB_firth_warning")
    expect_identical(fit$penalty$method, "Firth")
    expect_identical(fit$penalty$structure, "none")
    expect_identical(fit$penalty$requested, "component")
    expect_length(fit$penalty$sp, 0L)
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(fit$predicted > 0 & fit$predicted < 1))
    expect_lt(fit$penalty$adjusted_score, 1e-6)
  }
})

test_that("failed REML is replaced rather than retried with another curvature penalty", {
  skip_if_not_installed("mgcv")
  calls <- 0L
  testthat::local_mocked_bindings(gam = function(...) {
    calls <<- calls + 1L
    list(converged = FALSE)
  }, .package = "mgcv")
  expect_warning(fit <- separated_stageB(), class = "geds_stageB_firth_warning")
  expect_equal(calls, 1L)
  expect_identical(fit$penalty$method, "Firth")
  expect_match(fit$penalty$reason, "did not converge reliably")
})

test_that("successful REML fits do not invoke Firth", {
  d <- data.frame(x = seq(0, 1, length.out = 100))
  y <- rep(c(0, 1, 0, 1, 1), 20)
  args <- list(d, y, d[FALSE], rep(0, 100), rep(1, 100),
    list("f(x)" = list(type = "GeDS", variables = "x")),
    list("f(x)" = numeric()), 3L, binomial(), NULL, "component")
  ref <- do.call(.stageB_penalised_fit, args)
  expect_identical(ref$penalty$method, "REML")
  local_mocked_bindings(.stageB_firth_refit = function(...) stop("Firth must not run"),
                       .package = "GeDS")
  fit <- do.call(.stageB_penalised_fit, args)
  expect_identical(fit$theta, ref$theta)
  expect_identical(fit$predicted, ref$predicted)
  expect_identical(fit$penalty, ref$penalty)
})

test_that("unrelated REML errors are not hidden by Firth", {
  local_mocked_bindings(gam = function(...) stop("Injected programming error"), .package = "mgcv")
  expect_error(separated_stageB(), "Injected programming error")
})

test_that("Firth is not used for non-binomial REML failures", {
  local_mocked_bindings(gam = function(...) list(converged = FALSE), .package = "mgcv")
  local_mocked_bindings(.stageB_firth_refit = function(...) stop("Firth must not run"), .package = "GeDS")
  d <- data.frame(x = seq(0, 1, length.out = 50))
  expect_error(.stageB_penalised_fit(d, sin(d$x), d[FALSE], rep(0, 50), rep(1, 50),
    list("f(x)" = list(type = "GeDS", variables = "x")),
    list("f(x)" = numeric()), 3L, gaussian(), NULL, "component"), "did not converge reliably")
})

test_that("Firth handles aliases and positive weights on the fixed design", {
  x <- seq(-1, 1, length.out = 70)
  B <- cbind(1, x, 1)
  y <- as.numeric(x > 0)
  w <- rep(c(0, .5, 1, 2, 1), 14)
  keep <- w > 0
  fit <- .stageB_firth_refit(B[keep, ], y[keep], w[keep], B, diag(3))
  expect_warning(ref <- brglm2::brglmFit(x = cbind(1, x[keep]), y = y[keep], weights = w[keep],
    family = binomial(), intercept = TRUE,
    control = brglm2::brglmControl(type = "AS_mean", epsilon = 1e-10, maxit = 1000)),
    "non-integer.*success")
  expect_identical(fit$rank, 2L)
  expect_equal(fit$eta[keep], as.numeric(ref$linear.predictors), tolerance = 1e-7)
  expect_equal(fit$eta, as.numeric(cbind(1, x) %*% ref$coefficients), tolerance = 1e-7)
  expect_lt(fit$score, 1e-6)
})
