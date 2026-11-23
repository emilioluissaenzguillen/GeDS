test_that("standalone and additive fits reject unsupported joint dimensions", {
  set.seed(512)
  d <- data.frame(x1 = runif(40), x2 = runif(40), x3 = runif(40), x4 = runif(40))
  d$y <- 1 + sin(d$x1) + d$x2 + rnorm(40, sd = .1)
  for (formula in list(y ~ f(x1, x2, x3), y ~ f(x1, x2, x3, x4))) {
    for (rule in c("GD", "RD", "SR")) {
      expect_error(NGeDS(formula, data = d, stoptype = rule), "at most two covariates")
      expect_error(GGeDS(formula, data = d, stoptype = rule), "at most two covariates")
      expect_error(NGeDSgam(formula, data = d, stoptype = rule), "at most two covariates")
      if (rule %in% c("GD", "RD"))
        expect_error(NGeDSboost(formula, data = d, stoptype = rule), "at most two covariates")
    }
  }
  expect_equal(unname(GeDS:::f(d$x1)), unname(as.matrix(d["x1"])))
  expect_equal(unname(GeDS:::f(d$x1, d$x2)), unname(as.matrix(d[c("x1", "x2")])))
  expect_error(GeDS:::f(d$x1, d$x2, d$x3), "at most two covariates")
})

test_that("multiple separate spline components remain supported", {
  set.seed(312)
  d <- data.frame(x1 = runif(60), x2 = runif(60), x3 = runif(60), x4 = runif(60))
  d$y <- 1 + d$x1 + d$x2 * d$x3 + sin(3 * d$x4) + rnorm(60, sd = .1)
  formula <- y ~ f(x1) + f(x2, x3) + f(x4)
  parsed <- GeDS:::read.formula.gam(formula, d, type = "gam")
  expect_equal(lengths(lapply(parsed$base_learners, `[[`, "variables")), c(1L, 2L, 1L),
               ignore_attr = TRUE)
  for (fun in list(NGeDSgam, NGeDSboost)) {
    invisible(capture.output(fit <- suppressWarnings(suppressMessages(fun(
      formula, data = d, max_iterations = 2L, higher_order = FALSE
    )))))
    expect_true(all(is.finite(predict(fit, n = 2))))
    expect_equal(predict(fit, newdata = d, n = 2), predict(fit, n = 2), tolerance = 1e-8)
  }
})

test_that("the deferred fitter and class are absent from the submission namespace", {
  ns <- asNamespace("GeDS")
  for (name in c("MultivariateFitter", "GenMultivariateFitter", "predict.GeDSfitND",
                 "coef.GeDSfitND", "print.GeDSfitND", "summary.GeDSfitND"))
    expect_false(exists(name, envir = ns, inherits = FALSE))
  expect_false(any(getNamespaceInfo(ns, "S3methods")[, 2L] == "GeDSfitND"))
  methods <- getNamespaceInfo(ns, "S3methods")
  for (class in c("GeDS", "GeDSgam", "GeDSboost")) for (generic in c("print", "summary"))
    expect_true(any(methods[, 1L] == generic & methods[, 2L] == class))
})
