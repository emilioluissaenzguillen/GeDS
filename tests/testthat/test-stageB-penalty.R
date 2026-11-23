penalty_fixture <- function() {
  set.seed(813)
  d <- data.frame(x = runif(100), z = runif(100), a = rnorm(100))
  d$y <- sin(6 * d$x) + d$z^2 + 0.1 * d$a + rnorm(100, sd = 0.25)
  d
}

test_that("component REML is the additive default and leaves Stage A unchanged", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  for (fun in list(NGeDSgam, NGeDSboost)) {
    args <- list(formula = y ~ f(x) + f(z) + a, data = d,
                 phi = 0.8, max_iterations = 4L)
    original <- suppressMessages(do.call(fun, args))
    none <- suppressMessages(do.call(fun, c(args, list(penalty = "none"))))
    expect_identical(predict(original, n = 2), predict(none, n = 2))
    expect_identical(original$internal_knots, none$internal_knots)
    model <- suppressMessages(do.call(fun, c(args, list(penalty = "component"))))
    for (n in 2:4) expect_identical(predict(original, n = n), predict(model, n = n))
    expect_identical(original$internal_knots, model$internal_knots)
    for (n in 3:4) {
      fit <- model$final_model[[if (n == 3) "quadratic.fit" else "cubic.fit"]]
      expect_equal(fit$penalty$structure, "component")
      expect_equal(fit$penalty$method, "REML")
      expect_true(fit$penalty$converged)
      expect_true(all(is.finite(fit$theta)))
      expect_true(all(is.finite(fit$penalty$sp)))
      expect_length(fit$penalty$sp, 2L)
      expect_equal(predict(model, newdata = d, n = n), predict(model, n = n),
                   tolerance = 1e-8)
      expect_error(confint(model, n = n), "not implemented")
    }
  }
})

test_that("linear-only and single-learner boosting defaults remain unpenalised", {
  d <- penalty_fixture()
  for (fun in list(NGeDSgam, NGeDSboost)) {
    args <- list(formula = y ~ f(x) + f(z), data = d,
                 higher_order = FALSE, max_iterations = 3L)
    default <- suppressMessages(do.call(fun, args))
    none <- suppressMessages(do.call(fun, c(args, list(penalty = "none"))))
    expect_identical(predict(default, n = 2), predict(none, n = 2))
  }
  args <- list(formula = y ~ f(x), data = d, max_iterations = 3L)
  default <- suppressMessages(do.call(NGeDSboost, args))
  none <- suppressMessages(do.call(NGeDSboost, c(args, list(penalty = "none"))))
  for (n in 2:4) expect_identical(predict(default, n = n), predict(none, n = n))
  explicit <- suppressMessages(do.call(NGeDSboost, c(args, list(penalty = "component"))))
  expect_identical(predict(default, n = 2), predict(explicit, n = 2))
  expect_equal(explicit$final_model$cubic.fit$penalty$structure, "component")
})

test_that("curvature roots retain the affine null space and tensor coefficient order", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  for (n in 3:4) {
    ax <- GeDS:::.stageB_penalty_axis(d$x, c(0.2, 0.5, 0.8), n)
    expect_equal(as.numeric(ax$roots[[3L]] %*% rep(1, ncol(ax$basis))),
                 rep(0, nrow(ax$roots[[3L]])), tolerance = 1e-9)
    expect_equal(as.numeric(ax$roots[[3L]] %*% ax$affine),
                 rep(0, nrow(ax$roots[[3L]])), tolerance = 1e-9)
    model <- suppressMessages(NGeDSgam(y ~ f(x, z) + a, data = d,
                                      phi = 0.7, max_iterations = 3,
                                      penalty = "component"))
    expect_equal(predict(model, newdata = d, n = n), predict(model, n = n),
                 tolerance = 1e-8)
  }
})

test_that("relative weights, zero weights and supported families are handled", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  learners <- list("f(x)" = list(type = "GeDS", variables = "x"),
                   "f(z)" = list(type = "GeDS", variables = "z"))
  internal <- list("f(x)" = c(0.3, 0.6), "f(z)" = c(0.3, 0.6))
  w <- rep(c(0, 1, 2, 3), length.out = nrow(d))
  refit <- function(y, family, weights = w, n = 4, penalty = "component") {
    GeDS:::.stageB_penalised_fit(d[c("x", "z")], y, d["a"],
                               rep(0, nrow(d)), weights, learners, internal,
                               n, family, NULL, penalty)
  }
  a <- refit(d$y, gaussian())
  b <- refit(d$y, gaussian(), weights = w * 5)
  expect_equal(a$predicted, b$predicted, tolerance = 1e-8)
  altered <- d$y
  altered[w == 0] <- altered[w == 0] + 1000
  expect_equal(refit(altered, gaussian())$theta, a$theta, tolerance = 1e-8)
  set.seed(17)
  mu <- exp(0.1 + d$x + d$z)
  cases <- list(list(y = rgamma(100, shape = 4, scale = mu / 4), family = Gamma("log")),
                list(y = rpois(100, mu), family = poisson()),
                list(y = rbinom(100, 1, plogis(d$x - d$z)), family = binomial()),
                list(y = rbinom(100, 1, plogis(d$x - d$z)) * 2 - 1,
                     family = mboost::Binomial()))
  for (case in cases) for (n in 3:4) {
    fit <- refit(case$y, case$family, n = n)
    expect_true(all(is.finite(fit$predicted)))
    expect_true(fit$penalty$converged)
  }
  expect_error(refit(mu, Gamma("inverse")), "supports Gaussian")
  expect_error(refit(rep(0.5, 100), binomial()), "individual binary")
})

test_that("invalid opt-in requests fail explicitly", {
  skip_if_not_installed("mgcv")
  expect_error(NGeDSgam(y ~ f(x), penalty = "invalid"), "arg")
  expect_error(NGeDSgam(y ~ f(x), penalty = "component", normalize_data = TRUE),
               "normalize_data")
  expect_error(NGeDSboost(y ~ f(x), penalty = "component", higher_order = FALSE),
               "higher_order")
  for (fun in list(NGeDSgam, NGeDSboost)) {
    expect_error(fun(y ~ f(x), penalty = "common"), "arg.*component")
  }
})

test_that("wide weighted bases are regularised without penalty Gram matrices", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  w <- rep(c(0, 1), 50)
  learners <- list("f(x)" = list(type = "GeDS", variables = "x"),
                   "f(z)" = list(type = "GeDS", variables = "z"))
  internal <- lapply(learners, function(bl) seq(0.01, 0.99, length.out = 75L))
  for (n in 3:4) {
    fit <- GeDS:::.stageB_penalised_fit(
      d[c("x", "z")], d$y, d[FALSE], rep(0, 100), w,
      learners, internal, n, gaussian(), NULL, "component"
    )
    expect_true(fit$penalty$converged)
    expect_true(all(is.finite(fit$predicted)))
    expect_length(fit$theta, 2L * (75L + n))
  }
})

test_that("non-Gaussian GAM and boosting predictions use the correct response scale", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  set.seed(27)
  mu <- exp(0.1 + d$x)
  cases <- list(
    list(y = rgamma(100, 4, rate = 4 / mu), gam = Gamma("log"), boost = mboost::GammaReg()),
    list(y = rpois(100, mu), gam = poisson(), boost = mboost::Poisson()),
    list(y = rbinom(100, 1, plogis(d$x)), gam = binomial(), boost = mboost::Binomial())
  )
  for (case in cases) {
    d$y <- case$y
    for (type in c("gam", "boost")) {
      fun <- if (type == "gam") NGeDSgam else NGeDSboost
      # The native Stage A binary IRLS can warn before the guarded Stage B refit.
      fit <- withCallingHandlers(
        suppressMessages(fun(y ~ f(x) + f(z), data = d,
                             family = case[[type]], phi = 0.7,
                             max_iterations = 3L, penalty = "component")),
        warning = function(w) {
          if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
      for (n in 3:4) {
        expect_equal(predict(fit, newdata = d, n = n), predict(fit, n = n),
                     tolerance = 1e-8)
        eta <- predict(fit, newdata = d, n = n, type = "link")
        expect_equal(case$gam$linkinv(eta), predict(fit, n = n), tolerance = 1e-8)
      }
    }
  }
})

test_that("the transformed REML fit agrees with a direct curvature penalty", {
  skip_if_not_installed("mgcv")
  d <- penalty_fixture()
  for (n in 3:4) {
    axis <- GeDS:::.stageB_penalty_axis(d$x, c(0.2, 0.4, 0.7), n)
    design <- axis$basis
    y <- d$y
    S <- crossprod(axis$roots[[3L]]) * nrow(d)
    direct <- mgcv::gam(y ~ design - 1, method = "REML",
                         paraPen = list(design = list(S, rank = ncol(design) - 2L)))
    fit <- GeDS:::.stageB_penalised_fit(
      d["x"], y, d[FALSE], rep(0, nrow(d)), rep(1, nrow(d)),
      list("f(x)" = list(type = "GeDS", variables = "x")),
      list("f(x)" = c(0.2, 0.4, 0.7)), n, gaussian(), NULL, "component"
    )
    expect_equal(fit$predicted, as.numeric(fitted(direct)), tolerance = 1e-6)
    expect_equal(fit$penalty$edf, sum(direct$edf), tolerance = 1e-6)
  }
})
