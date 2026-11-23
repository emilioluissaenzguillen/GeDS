test_that("rank-deficient weighted Stage B retains the QR coefficient solution", {
  x <- seq(0, 1, length.out = 60)
  X <- data.frame(x = x, z = x^2)
  Z <- data.frame(group = factor(rep(1:3, 20)))
  weights <- rep(c(1, 2, 3), 20)
  weights[x > 0.55] <- 0
  offset <- 0.2 + x / 3
  response <- sin(5 * x) + x^2 + as.numeric(Z$group) / 4 + offset
  learners <- list(sx = list(type = "GeDS", variables = "x"),
                   sz = list(type = "GeDS", variables = "z"))
  internal <- list(sx = c(0.25, 0.5, 0.75), sz = c(0.2, 0.5, 0.8))

  for (n in 3:4) {
    fit <- GeDS:::SplineReg_LM_Multivar(
      X = X, Y = response, Z = Z, offset = offset,
      base_learners = learners, weights = weights,
      InterKnotsList = internal, n = n,
      coefficients = NULL, linear.predictors = NULL, only_pred = TRUE
    )
    design <- cbind(fit$basis, model.matrix(~ ., Z)[, -1, drop = FALSE])
    reference <- lm(I(response - offset) ~ design - 1, weights = weights)
    theta <- unname(coef(reference))
    expect_true(anyNA(theta))
    theta[is.na(theta)] <- 0
    expect_equal(unname(fit$theta), theta, tolerance = 1e-12)
    reconstructed <- as.numeric(design %*% fit$theta) + offset
    for (rows in list(weights > 0, weights == 0)) {
      expect_equal(as.numeric(fit$predicted)[rows], reconstructed[rows],
                   tolerance = 1e-10)
    }
    expect_equal(as.numeric(fit$predicted),
                 unname(fitted(reference)) + offset, tolerance = 1e-10)
    expect_equal(fit$rss, sum(weights * (response - reconstructed)^2),
                 tolerance = 1e-10)
    expect_equal(as.numeric(fit$residuals), response - reconstructed,
                 tolerance = 1e-10)
  }
})

test_that("weighted additive higher-order training and newdata predictions agree", {
  set.seed(1729)
  d <- data.frame(x = rnorm(100), z = rnorm(100),
                  group = factor(rep(1:4, 25)))
  d$y <- 3 * sin(d$x) + d$z^2 + as.numeric(d$group) / 4 + rnorm(100)
  weights <- drop(rmultinom(1, 100, rep(1, 100)))
  expect_true(any(weights == 0))

  for (fun in list(NGeDSgam, NGeDSboost)) {
    args <- list(formula = y ~ f(x) + f(z) + group,
                 data = d, weights = weights, phi = 0.9, q = 2,
                 penalty = "none")
    if (identical(fun, NGeDSboost)) {
      args <- c(args, list(initial_learner = FALSE, shrinkage = 0.1))
    }
    fit <- suppressWarnings(suppressMessages(do.call(fun, args)))
    for (n in 3:4) {
      stage_b <- fit$final_model[[if (n == 3) "quadratic.fit" else "cubic.fit"]]
      expect_true(anyNA(coef(stage_b$temporary)))
      for (type in c("response", "link")) {
        stored <- as.numeric(predict(fit, n = n, type = type))
        explicit <- as.numeric(predict(fit, newdata = d, n = n, type = type))
        for (rows in list(weights > 0, weights == 0)) {
          expect_equal(stored[rows], explicit[rows], tolerance = 1e-8)
        }
      }
    }
  }
})
