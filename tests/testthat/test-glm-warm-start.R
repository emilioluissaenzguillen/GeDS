warm_start_data <- function() {
  set.seed(1009)
  x <- runif(200)
  y <- rbinom(200, 10, plogis(2 * sin(2 * pi * x))) / 10
  o <- order(x)
  list(x = x[o], y = y[o], w = rep(10, 200))
}

test_that("a diverging IRLS warm start is replaced by the default start", {
  d <- warm_start_data()
  extr <- range(d$x)
  one <- SplineReg_GLM(X = d$x, Y = d$y, Z = NULL, weights = d$w, InterKnots = 0.6371134,
                       n = 2, extr = extr, family = binomial(), fast = TRUE)
  # Warm start: the one-knot fit carried into the two-knot basis.
  inits <- newknot.guess(0.6371134, extr, one$theta, 0.03006273)
  knots <- c(0.03006273, 0.6371134)
  basis <- splines::splineDesign(sort(c(knots, rep(extr, 2))), d$x, ord = 2)
  plain <- suppressWarnings(glm.fit(basis, d$y, weights = d$w, family = binomial(),
                                    mustart = binomial()$linkinv(basis %*% inits),
                                    intercept = FALSE))
  expect_gt(plain$deviance, 100 * one$temporary$deviance)  # glm.fit diverges from it
  cold <- glm.fit(basis, d$y, weights = d$w, family = binomial(), intercept = FALSE)
  for (fast in c(TRUE, FALSE)) {
    fit <- suppressWarnings(SplineReg_GLM(X = d$x, Y = d$y, Z = NULL, weights = d$w,
      InterKnots = knots, n = 2, extr = extr, family = binomial(), inits = inits, fast = fast))
    deviance <- if (fast) fit$temporary$deviance else fit$rss
    expect_equal(deviance, cold$deviance, tolerance = 1e-8)
    expect_lt(deviance, one$temporary$deviance)
  }
})

test_that("the GLM Stage A deviance path does not jump after a warm-started knot", {
  d <- warm_start_data()
  path <- suppressWarnings(GenUnivariateFitter(X = d$x, Y = d$y, weights = d$w,
    family = binomial(), beta = 0.1, max.intknots = 12L, q = 2L,
    stoptype = ".gamma_dynamic_path"))
  expect_length(path, 13L)
  expect_true(all(diff(path) <= 1e-8 * path[-length(path)]))
})

test_that("a converged warm start is kept unchanged", {
  set.seed(7)
  x <- sort(runif(150))
  y <- rpois(150, exp(1 + sin(4 * x)))
  first <- SplineReg_GLM(X = x, Y = y, Z = NULL, InterKnots = 0.5, n = 2,
                         family = poisson(), fast = TRUE)
  inits <- newknot.guess(0.5, range(x), first$theta, 0.25)
  basis <- splines::splineDesign(sort(c(0.25, 0.5, rep(range(x), 2))), x, ord = 2)
  plain <- glm.fit(basis, y, family = poisson(), intercept = FALSE,
                   mustart = poisson()$linkinv(basis %*% inits))
  fit <- SplineReg_GLM(X = x, Y = y, Z = NULL, InterKnots = c(0.25, 0.5), n = 2,
                       family = poisson(), inits = inits, fast = TRUE)
  expect_identical(fit$theta, plain$coefficients)
  expect_identical(fit$temporary$iter, plain$iter)
})

test_that("only the returned fit's warnings are reported after a cold retry", {
  d <- warm_start_data()
  extr <- range(d$x)
  one <- SplineReg_GLM(X = d$x, Y = d$y, Z = NULL, weights = d$w, InterKnots = 0.6371134,
                       n = 2, extr = extr, family = binomial(), fast = TRUE)
  inits <- newknot.guess(0.6371134, extr, one$theta, 0.03006273)
  knots <- c(0.03006273, 0.6371134)
  basis <- splines::splineDesign(sort(c(knots, rep(extr, 2))), d$x, ord = 2)
  messages <- function(expr) {
    out <- character()
    withCallingHandlers(expr, warning = function(w) {
      out <<- c(out, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    out
  }
  warm <- messages(glm.fit(basis, d$y, weights = d$w, family = binomial(), intercept = FALSE,
                           mustart = binomial()$linkinv(basis %*% inits)))
  expect_true(any(grepl("did not converge", warm)))
  cold <- messages(glm.fit(basis, d$y, weights = d$w, family = binomial(), intercept = FALSE))
  for (fast in c(TRUE, FALSE)) {
    reported <- messages(SplineReg_GLM(X = d$x, Y = d$y, Z = NULL, weights = d$w,
      InterKnots = knots, n = 2, extr = extr, family = binomial(), inits = inits, fast = fast))
    expect_identical(reported, cold)
  }
})
