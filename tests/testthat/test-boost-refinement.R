test_that("fixed-basis SVD uses minimum-norm updates in rank-deficient weighted fits", {
  B <- cbind(1, c(0, 0, 1, 1), c(0, 1, 0, 0))
  w <- c(1, 0, 1, 1)
  y <- c(1, 3, 2, 2)
  ref <- GeDS:::.boost_refine_coefficients(
    B, y, rep(0, 4), mboost::Gaussian(), w, .5, .99, 2L, 0L, 1L
  )
  expect_equal(ref$rank, 2L)
  expect_equal(ref$selected$delta, c(.5, .5, 0), tolerance = 1e-12)
  expect_equal(ref$selected$F_hat, c(.5, .5, 1, 1), tolerance = 1e-12)
  perm <- c(3, 1, 2)
  other <- GeDS:::.boost_refine_coefficients(
    B[, perm], y, rep(0, 4), mboost::Gaussian(), w, .5, .99, 2L, 0L, 1L
  )
  expect_equal(ref$selected$F_hat, other$selected$F_hat, tolerance = 1e-12)
})

test_that("Phase 2 resets RD and keeps the original shrinkage and final-window selection", {
  B <- cbind(1, seq(-1, 1, length.out = 20))
  y <- 2 + B[, 2] + sin(4 * B[, 2])
  F0 <- rep(2, 20)
  ref <- GeDS:::.boost_refine_coefficients(
    B, y, F0, mboost::Gaussian(), rep(1, 20), .4, .99, 2L, 0L, 20L
  )
  delta <- .4 * qr.coef(qr(B), y - F0)
  expect_equal(ref$trace$deviance[2], sum((y - F0 - B %*% delta)^2))
  expect_true(all(is.na(ref$trace$ratio[1:2])))
  expect_equal(ref$trace$ratio[-(1:2)],
               tail(ref$trace$deviance, -2) / head(ref$trace$deviance, -2))
  expect_equal(ref$selected$deviance, min(tail(ref$trace$deviance, 3)))
  expect_identical(ref$reason, "RD")
  expect_equal(ref$shrinkage, .4)
})

test_that("memory freezes each learner after two knotless wins and keeps it eligible", {
  x <- seq(0, 1, length.out = 60)
  d <- data.frame(x = x, z = x^2, y = 2 * x + sin(6 * x))
  native_fit <- GeDS:::componentwise_fit
  frozen_winner <- FALSE
  local_mocked_bindings(componentwise_fit = function(bl_name, m, internal_knots,
                                                    default_spline_cap, knots_frozen, ...) {
    winner <- if (m %in% c(4L, 7L)) "f(z)" else "f(x)"
    force_knot <- bl_name == "f(x)" && m == 2L
    fit <- native_fit(bl_name = bl_name, m = m,
                      internal_knots = if (force_knot || knots_frozen) 1L else 0L,
                      default_spline_cap = FALSE, knots_frozen = knots_frozen, ...)
    if (bl_name == "f(x)" && m == 6L) frozen_winner <<- knots_frozen
    fit$ssr <- if (bl_name == winner) 1 else 2
    fit
  }, .package = "GeDS")

  invisible(capture.output(fit <- NGeDSboost(
    y ~ f(x) + f(z), data = d, initial_learner = FALSE,
    boosting_with_memory = TRUE, shrinkage = .2, int.knots_boost = 3L,
    phi_boost_exit = NULL, max_iterations = 9L, higher_order = FALSE
  )))
  expect_identical(fit$args$learner_knots_frozen_at, c("f(x)" = 5L, "f(z)" = 7L))
  expect_identical(fit$args$knot_stability_stop$iteration, 7L)
  expect_true(frozen_winner)
  expect_equal(fit$models$model2$best_bl$n_intknots, 1L)
  expect_equal(fit$models$model6$best_bl$n_intknots, 0L)
})

test_that("family gradients and actual deviances use their native response scales", {
  for (fam in list(mboost::Poisson(), mboost::GammaReg(), mboost::Binomial())) {
    fs <- GeDS:::get_mboost_family(fam@name)
    y <- if (fs$family == "binomial") rep(c(-1, 1), 10) else seq(1, 3, length.out = 20)
    B <- cbind(1, seq(-1, 1, length.out = 20))
    F0 <- rep(if (fs$family == "binomial") 0 else log(2), 20)
    # Gamma's nuisance estimate is updated by the native gradient itself.
    g <- fam@ngradient(y = y, f = F0, w = rep(1, 20))
    expected <- F0 + .1 * as.numeric(B %*% qr.coef(qr(B), g))
    ref <- GeDS:::.boost_refine_coefficients(
      B, y, F0, fam, rep(1, 20), .1, .99, 2L, 0L, 1L
    )
    expect_equal(ref$selected$F_hat, expected, tolerance = 1e-10)
    dev_y <- if (fs$family == "binomial") (y + 1) / 2 else y
    expect_equal(ref$selected$deviance,
                 sum(fs$dev.resids(dev_y, fam@response(expected), rep(1, 20))))
  }
})

test_that("refinement preserves knots and higher orders and updates stored coefficients", {
  set.seed(247)
  d <- data.frame(x = sort(runif(80, -1, 1)))
  d$y <- 2 + sin(4 * d$x) + rnorm(80, sd = .1)
  fit <- suppressMessages(suppressWarnings(NGeDSboost(
    y ~ f(x), data = d, initial_learner = FALSE, shrinkage = .5,
    max_iterations = 3L
  )))
  phase1 <- fit$models[[fit$final_model$model_name]]
  expect_equal(knots(fit, n = 2, options = "internal")[[1]],
               GeDS:::get_internal_knots(phase1$base_learners[["f(x)"]]$knots))
  expect_equal(predict(fit, n = 2), predict(fit, newdata = d, n = 2), tolerance = 1e-8)
  expect_equal(predict(fit, type = "link", n = 2), fit$final_model$F_hat)
  expect_equal(fit$iters, length(fit$models) - 1L + fit$final_model$refinement$iterations)
  expect_equal(N.boost.iter(fit), fit$iters)
  local_mocked_bindings(.boost_refine_coefficients = function(B, y, F0, family, weights,
                                                            shrinkage, phi, q, ...) {
    list(selected = list(step = 0L, delta = numeric(ncol(B)), F_hat = F0,
                         deviance = sum((y - F0)^2)), iterations = 0L)
  }, .package = "GeDS")
  reference <- suppressMessages(suppressWarnings(NGeDSboost(
    y ~ f(x), data = d, initial_learner = FALSE, shrinkage = .5,
    max_iterations = 3L
  )))
  expect_equal(fit$models, reference$models)
  for (n in 3:4) expect_equal(predict(fit, n = n), predict(reference, n = n))
  expect_equal(fit$internal_knots, reference$internal_knots)
})

test_that("additive, tensor, factor and normalized predictions include refinement", {
  set.seed(248)
  d <- data.frame(x = runif(70), z = runif(70), a = factor(rep(1:2, 35)))
  d$y <- 2 + 3 * d$x * d$z + as.numeric(d$a) + rnorm(70, sd = .1)
  for (form in list(y ~ f(x) + a, y ~ f(x, z) + a)) {
    for (normalize in c(FALSE, TRUE)) {
      fit <- suppressMessages(suppressWarnings(NGeDSboost(
        form, data = d, shrinkage = .5, initial_learner = FALSE,
        normalize_data = normalize, max_iterations = 2L, higher_order = FALSE
      )))
      expect_equal(predict(fit, n = 2), predict(fit, newdata = d, n = 2), tolerance = 1e-7)
      expect_true(all(is.finite(predict(fit, n = 2))))
      for (learner in names(fit$final_model$base_learners)) {
        whole <- predict(fit, newdata = d[1:3, ], n = 2, base_learner = learner)
        single <- vapply(1:3, function(i) {
          predict(fit, newdata = d[i, , drop = FALSE], n = 2, base_learner = learner)
        }, numeric(1))
        expect_equal(whole, single, tolerance = 1e-8)
      }
      if (normalize) {
        expect_equal(deviance(fit, n = 2), sum(fit$args$weights * (d$y - predict(fit, n = 2))^2))
      }
    }
  }
})

test_that("repeated knots preserve anchored binomial predictions", {
  d <- droplevels(subset(iris, Species != "virginica"))
  for (normalize in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(NGeDSboost(
      Species ~ f(Sepal.Length, Sepal.Width) + f(Petal.Length) + f(Petal.Width),
      data = d, family = mboost::Binomial(link = "probit"),
      initial_learner = FALSE, normalize_data = normalize,
      phi_boost_exit = .8, higher_order = FALSE
    )))
    expect_equal(predict(fit, n = 2), predict(fit, newdata = d, n = 2), tolerance = 1e-10)
    expect_equal(predict(fit, n = 2, type = "link"),
                 predict(fit, newdata = d, n = 2, type = "link"), tolerance = 1e-10)
  }
})

test_that("tensor component contributions retain the refinement intercept", {
  set.seed(248)
  grid <- seq(-1, 1, length.out = 8L)
  d <- expand.grid(x = grid, z = grid)
  d$a <- factor(rep(1:2, each = 32L))
  d$y <- with(d, 2 + 10 * x * z + as.numeric(a)) + rnorm(nrow(d), sd = .1)
  for (normalize in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(NGeDSboost(
      y ~ f(x, z) + a,
      data = d, initial_learner = TRUE, normalize_data = normalize,
      max_iterations = 3L, higher_order = FALSE, int.knots_boost = 2L
    )))
    expect_true("f(x, z)" %in% names(fit$final_model$base_learners))
    parts <- lapply(names(fit$final_model$base_learners), function(bl) {
      predict(fit, newdata = d, n = 2L, base_learner = bl)
    })
    total <- Reduce(`+`, parts)
    if (normalize) total <- total * fit$args$Y_sd + fit$args$Y_mean
    expect_equal(total, predict(fit, newdata = d, n = 2L), tolerance = 1e-8)
  }
})

test_that("nonfinite fixed-basis updates are rejected without changing the anchor", {
  fam <- mboost::Gaussian()
  fam@ngradient <- function(y, f, w = 1) rep(Inf, length(y))
  expect_warning(ref <- GeDS:::.boost_refine_coefficients(
    matrix(1, 5, 1), 1:5, rep(3, 5), fam, rep(1, 5), .1, .99, 2L, 0L, 3L
  ), "non-finite")
  expect_equal(ref$selected$F_hat, rep(3, 5))
  expect_equal(ref$iterations, 0L)
})

test_that("weighted rank-deficient additive fits predict consistently in both phases", {
  set.seed(250)
  d <- data.frame(x = rep(seq(-1, 1, length.out = 12), each = 3))
  d$z <- d$x
  d$y <- 2 + sin(4 * d$x) + rnorm(36, sd = .1)
  w <- rep(c(0, 1, 2), 12)
  fit <- suppressMessages(suppressWarnings(NGeDSboost(
    y ~ f(x) + z, data = d, weights = w, initial_learner = FALSE,
    shrinkage = .5, max_iterations = 4L, higher_order = FALSE
  )))
  expect_equal(predict(fit, n = 2), predict(fit, newdata = d, n = 2), tolerance = 1e-8)
  expect_true(all(is.finite(predict(fit, n = 2))))
  expect_equal(fit$final_model$dev, sum(w * (d$y - predict(fit, n = 2))^2))
})

test_that("refined family fits preserve native link predictions and offsets", {
  for (fam in list(mboost::Poisson(), mboost::GammaReg(), mboost::Binomial())) {
    set.seed(249)
    d <- data.frame(x = sort(runif(80, -.5, .5)))
    fs <- GeDS:::get_mboost_family(fam@name)
    d$y <- switch(fs$family,
      poisson = rpois(80, exp(2 + d$x)),
      Gamma = rgamma(80, shape = 10, scale = exp(d$x)),
      binomial = factor(rbinom(80, 1, plogis(3 * d$x))))
    for (initial in c(FALSE, TRUE)) {
      fit <- suppressMessages(suppressWarnings(NGeDSboost(
        y ~ f(x), data = d, family = fam, initial_learner = initial,
        shrinkage = .1, max_iterations = 3L, higher_order = FALSE
      )))
      for (type in c("response", "link")) {
        expect_equal(predict(fit, n = 2, type = type),
                     predict(fit, newdata = d, n = 2, type = type), tolerance = 1e-7)
      }
      native_y <- fit$args$response[[1]]
      dev_y <- if (fs$family == "binomial") (native_y + 1) / 2 else native_y
      expect_equal(fit$final_model$dev,
                   sum(fs$dev.resids(dev_y, predict(fit, n = 2), fit$args$weights)))
    }
  }
})
