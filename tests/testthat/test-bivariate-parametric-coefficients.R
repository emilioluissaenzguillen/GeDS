bivariate_parametric_data <- function() {
  set.seed(123)
  data <- expand.grid(x = seq(0, 1, length.out = 10), v = seq(0, 1, length.out = 12))
  data$z <- rnorm(nrow(data))
  data$group <- factor(rep(c("A", "B", "C", "B"), length.out = nrow(data)))
  data$off <- runif(nrow(data), -.15, .15)
  data$w <- rep(c(.5, 1, 2, 0), length.out = nrow(data))
  data
}

test_that("bivariate Stage B retains numeric and factor coefficients for all orders", {
  data <- bivariate_parametric_data()
  eta <- 1 + .3 * sin(2 * data$x) + .2 * data$v + .4 * data$z +
    .2 * (data$group == "B") - .1 * (data$group == "C") + data$off
  families <- list(NGeDS = NULL, Gaussian = gaussian(), Poisson = poisson(),
                   Binomial = binomial(), Gamma_log = Gamma(link = "log"),
                   Gamma_inverse = Gamma())
  for (kind in names(families)) {
    data$w <- rep(c(.5, 1, 2, 0), length.out = nrow(data))
    data$y <- switch(kind,
      NGeDS = eta + rnorm(nrow(data), sd = .1),
      Gaussian = eta + rnorm(nrow(data), sd = .1),
      Poisson = rpois(nrow(data), exp(eta)),
      Binomial = rbinom(nrow(data), 20, plogis(eta)) / 20,
      Gamma_log = rgamma(nrow(data), shape = 20, scale = exp(eta) / 20),
      Gamma_inverse = rgamma(nrow(data), shape = 20, scale = 1 / (2 + .1 * eta) / 20))
    if (kind == "Binomial") data$w <- rep(c(20, 40, 60, 0), length.out = nrow(data))
    args <- list(formula = y ~ f(x, v) + z + group + offset(off), data = data,
      weights = data$w, phi = .9, q = 1, max.intknots = 6)
    fun <- if (kind == "NGeDS") NGeDS else GGeDS
    if (kind != "NGeDS") args$family <- families[[kind]]
    fit <- suppressWarnings(do.call(fun, args))
    expect_equal(colnames(fit$args$W), c("z", "groupB", "groupC"))
    for (n in 2:4) {
      stage <- fit[[c("linear.fit", "quadratic.fit", "cubic.fit")[n - 1L]]]
      B <- tensorProd(stage$Xbasis, stage$Ybasis)
      design <- cbind(B, fit$args$W)
      coefficients <- coef(fit, n = n, onlySpline = FALSE)
      expect_length(coefficients, ncol(design))
      expect_length(coef(fit, n = n, onlySpline = TRUE), ncol(B))
      expect_equal(tail(names(coefficients), 3), colnames(fit$args$W))
      if (kind == "NGeDS") {
        direct <- lm.wfit(design, data$y - data$off, data$w)
      } else {
        # Inverse-Gamma IRLS may propose invalid means before step halving.
        direct <- withCallingHandlers(
          glm.fit(design, data$y, family = families[[kind]],
                  weights = data$w, offset = data$off, intercept = FALSE,
                  mustart = as.numeric(stage$predicted)),
          warning = function(w) {
            if (kind == "Gamma_inverse" && identical(conditionMessage(w), "NaNs produced")) {
              invokeRestart("muffleWarning")
            }
          }
        )
        expect_true(isTRUE(direct$converged))
        expect_true(all(is.finite(direct$fitted.values)))
        expect_true(isTRUE(families[[kind]]$validmu(direct$fitted.values)))
      }
      theta <- coef(direct)
      theta[is.na(theta)] <- 0
      expect_equal(unname(coefficients), unname(theta), tolerance = 1e-6)
      link <- as.numeric(design %*% theta) + data$off
      expected <- if (kind == "NGeDS") link else families[[kind]]$linkinv(link)
      expect_equal(predict(fit, data, n = n), expected, tolerance = 1e-6)
      expect_equal(predict(fit, data, n = n), predict(fit, n = n))
      expect_equal(as.numeric(stage$predicted), predict(fit, data, n = n))
      shifted <- data
      shifted$z <- shifted$z + .2
      expect_equal(predict(fit, shifted, n = n, type = "link") -
                     predict(fit, data, n = n, type = "link"),
                   rep(.2 * unname(coefficients["z"]), nrow(data)), tolerance = 1e-8)
    }
  }
})

test_that("both bivariate fitters retain every Stage A parametric coefficient", {
  for (generalized in c(FALSE, TRUE)) {
    data <- bivariate_parametric_data()
    eta <- 1 + sin(3 * data$x) + .4 * data$v + .3 * data$z + .2 * (data$group == "B") + data$off
    data$y <- if (generalized) rpois(nrow(data), exp(eta)) else eta + rnorm(nrow(data), sd = .1)
    recorded <- list()
    original <- if (generalized) SplineReg_biv_GLM else SplineReg_biv
    wrapper <- function(...) {
      args <- list(...)
      result <- original(...)
      if (isTRUE(args$fast)) recorded[[length(recorded) + 1L]] <<- result$theta
      result
    }
    # Scope each mock separately so the two fitting engines remain independent.
    fit <- local({
      if (generalized) {
        testthat::local_mocked_bindings(SplineReg_biv_GLM = wrapper, .package = "GeDS")
        suppressWarnings(GGeDS(y ~ f(x, v) + z + group + offset(off), data = data,
          family = poisson(), weights = w, phi = .9, q = 1, max.intknots = 6))
      } else {
        testthat::local_mocked_bindings(SplineReg_biv = wrapper, .package = "GeDS")
        suppressWarnings(NGeDS(y ~ f(x, v) + z + group + offset(off), data = data,
          weights = w, phi = .9, q = 1, max.intknots = 6))
      }
    })
    expected_columns <- (ncol(fit$stored$previousX) - 2) *
      (ncol(fit$stored$previousY) - 2) + ncol(fit$args$W)
    expect_equal(ncol(fit$coefficients), expected_columns)
    expect_equal(nrow(fit$coefficients), length(recorded))
    for (i in seq_len(nrow(fit$coefficients))) {
      expected <- as.numeric(recorded[[i]])
      expect_gte(ncol(fit$coefficients), length(expected))
      expect_equal(as.numeric(fit$coefficients[i, seq_along(expected)]), expected)
      if (length(expected) < ncol(fit$coefficients)) {
        expect_true(all(is.na(fit$coefficients[i, (length(expected) + 1):ncol(fit$coefficients)])))
      }
    }
  }
})

test_that("an empty bivariate parametric matrix leaves fits unchanged", {
  data <- bivariate_parametric_data()
  data$y <- rpois(nrow(data), exp(1 + .3 * data$x + .2 * data$v))
  base <- list(X = data$x, Y = data$v, Z = data$y, weights = data$w,
    indicator = NULL, family = poisson(), phi = .9, q = 1,
    max.intknots = 6, show.iters = FALSE)
  without <- suppressWarnings(do.call(GenBivariateFitter, c(base, list(W = NULL))))
  empty <- suppressWarnings(do.call(GenBivariateFitter,
    c(base, list(W = matrix(numeric(0), nrow(data), 0)))))
  expect_equal(empty$stored, without$stored)
  expect_equal(empty$coefficients, without$coefficients)
  for (name in c("linear.fit", "quadratic.fit", "cubic.fit")) {
    expect_equal(empty[[name]]$theta, without[[name]]$theta)
    expect_equal(empty[[name]]$predicted, without[[name]]$predicted)
  }
})
