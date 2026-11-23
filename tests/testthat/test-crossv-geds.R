library("testthat")
library("GeDS")

test_that("crossv_GeDS runs on a small NGeDS example", {
  set.seed(123)
  N <- 40
  X <- sort(runif(N, min = -1, max = 1))
  Y <- sin(pi * X) + rnorm(N, sd = 0.05)
  dat <- data.frame(X = X, Y = Y)

  param <- list(
    beta_grid = c(0.5),
    phi_grid = c(0.9),
    q_grid = c(2)
  )

  cv <- suppressWarnings(
    crossv_GeDS(
      Y ~ f(X),
      data = dat,
      model_fun = NGeDS,
      n = 3,
      n_folds = 2,
      parameters = param,
      n_cores = 1L
    )
  )

  expect_true(is.list(cv))
  expect_true(all(c("best_params", "results") %in% names(cv)))
  expect_s3_class(cv$best_params, "data.frame")
  expect_s3_class(cv$results, "data.frame")
  expect_true(nrow(cv$best_params) >= 1)
  expect_true(nrow(cv$results) >= 1)
  expect_true("cv_mse" %in% names(cv$results))
  expect_true(all(is.finite(cv$results$cv_mse)))
})

test_that("crossv_GeDS scores binomial NGeDSgam fits by held-out deviance", {
  skip_if_not_installed("rpart")
  data("kyphosis", package = "rpart")
  set.seed(123)
  dat <- kyphosis[sample(seq_len(nrow(kyphosis)), size = 52L), ]
  param <- list(beta_grid = 0.5, phi_grid = 0.7, q_grid = 2)

  set.seed(123)
  cv <- suppressWarnings(crossv_GeDS(
    Kyphosis ~ Number + f(Age), data = dat, model_fun = NGeDSgam,
    n = 2L, n_folds = 5L, parameters = param, n_cores = 1L,
    family = binomial, q_gam = 1L
  ))

  set.seed(123)
  shuffled <- unlist(lapply(split(seq_len(nrow(dat)), dat$Kyphosis),
                            function(idx) idx[sample.int(length(idx))]),
                     use.names = FALSE)
  fold_id <- integer(nrow(dat))
  fold_id[shuffled] <- rep(sample.int(5L), length.out = nrow(dat))
  losses <- numeric(nrow(dat))
  for (fold in 1:5) {
    suppressWarnings(invisible(capture.output({
      model <- NGeDSgam(Kyphosis ~ Number + f(Age),
                        data = dat[fold_id != fold, ], family = binomial(),
                        beta = 0.5, phi = 0.7, q = 2, q_gam = 1L, stoptype = "RD",
                        higher_order = FALSE)
    })))
    pred <- predict(model, newdata = dat[fold_id == fold, ], n = 2)
    losses[fold_id == fold] <- binomial()$dev.resids(
      as.numeric(dat$Kyphosis[fold_id == fold]) - 1L, pred,
      rep(1, sum(fold_id == fold)))
  }

  expect_true("cv_deviance" %in% names(cv$results))
  expect_false("cv_mse" %in% names(cv$results))
  expect_equal(cv$results$cv_deviance, mean(losses), tolerance = 1e-8)
  expect_equal(cv$best_params$phi, 0.7)
  expect_equal(cv$best_params$q, 2)
})

test_that("crossv_GeDS scores Poisson GGeDS fits by deviance", {
  set.seed(17)
  x <- seq(0, 1, length.out = 45)
  dat <- data.frame(x = x, y = rpois(length(x), lambda = exp(0.5 + x)))
  param <- list(beta_grid = 0.2, phi_grid = 0.8, q_grid = 2L)

  set.seed(18)
  cv <- suppressWarnings(crossv_GeDS(
    y ~ f(x), data = dat, model_fun = GGeDS, n = 2L, n_folds = 2L,
    parameters = param, n_cores = 1L, family = poisson()
  ))

  expect_true("cv_deviance" %in% names(cv$results))
  expect_true(is.finite(cv$results$cv_deviance))
})

test_that("boosting CV does not reward failed folds", {
  dat <- data.frame(x = seq_len(6), y = seq_len(6))
  parameters <- list(int.knots_init_grid = 2L, shrinkage_grid = 0.6,
                     beta_grid = 0.5, phi_grid = 0.8, q_grid = 1L)

  expect_error(
    GeDS:::cross_validate.GeDSboost(
      y ~ f(x), data = dat, model_fun = function(...) stop("forced fit failure"),
      n = 2L, n_folds = 2L, parameters = parameters, n_cores = 1L
    ),
    "All cross-validation candidates failed"
  )
})

test_that("omitted CV grids use fitting defaults", {
  grids <- GeDS:::prepare_cv_grids(
    list(phi_grid = c(0.8, 0.9)),
    c("beta_grid", "phi_grid", "q_grid")
  )
  expect_identical(grids$tuned, "phi")
  expect_true(is.na(grids$values$beta_grid))
  expect_true(is.na(grids$values$q_grid))
  expect_identical(grids$values$phi_grid, c(0.8, 0.9))
  expect_error(GeDS:::prepare_cv_grids(list(phi_grid = numeric()), "phi_grid"),
               "finite numeric values")
  expect_error(GeDS:::prepare_cv_grids(list(phii_grid = 0.8), "phi_grid"),
               "Unknown parameter grid")
})

test_that("CV results show only supplied parameter grids", {
  set.seed(71)
  dat <- data.frame(x = seq(0, 1, length.out = 40))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  cv <- suppressWarnings(crossv_GeDS(
    y ~ f(x), data = dat, model_fun = NGeDS,
    n = 2L, n_folds = 2L,
    parameters = list(phi_grid = 0.9), n_cores = 1L
  ))

  expect_true(is.finite(cv$results$cv_mse))
  expect_true("phi" %in% names(cv$results))
  expect_false(any(c("beta", "q") %in% names(cv$results)))
  expect_false(any(c("beta", "q") %in% names(cv$best_params)))

  cv_default <- suppressWarnings(crossv_GeDS(
    y ~ f(x), data = dat, model_fun = NGeDS,
    n = 2L, n_folds = 2L, n_cores = 1L
  ))
  expect_true(is.finite(cv_default$results$cv_mse))
  expect_false(any(c("beta", "phi", "q") %in% names(cv_default$results)))
})

test_that("Poisson CV retains the GGeDS family-specific beta default", {
  set.seed(17)
  dat <- data.frame(x = seq(0, 1, length.out = 45))
  dat$y <- rpois(nrow(dat), lambda = exp(0.5 + dat$x))

  cv <- suppressWarnings(crossv_GeDS(
    y ~ f(x), data = dat, model_fun = GGeDS,
    n = 2L, n_folds = 2L,
    parameters = list(phi_grid = 0.8), n_cores = 1L,
    family = poisson()
  ))

  expect_true(is.finite(cv$results$cv_deviance))
  expect_false(any(c("beta", "q") %in% names(cv$results)))
})

test_that("boosting phi tuning uses RD unless a rule is supplied", {
  set.seed(123)
  dat <- data.frame(x = seq(0, 1, length.out = 40))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.1)
  run_cv <- function(...) {
    set.seed(123)
    crossv_GeDS(y ~ f(x), data = dat, model_fun = NGeDSboost,
                n = 2L, n_folds = 2L, n_cores = 1L,
                parameters = list(phi_grid = c(0.3, 0.99)),
                initial_learner = FALSE, max_iterations = 3L, ...)
  }
  automatic <- run_cv()
  explicit <- run_cv(stoptype = "RD")
  # Compare fitted scores and model sizes, excluding wall-clock timings.
  cols <- setdiff(names(automatic$results), "avg_fitting_time")
  expect_equal(automatic$results[cols], explicit$results[cols], tolerance = 1e-12)
  expect_true(all(is.finite(automatic$results$cv_mse)))
  expect_gt(diff(range(automatic$results$avg_intknots_univbl)), 0)
})

test_that("boosting CV rejects phi tuning under explicit GD before starting workers", {
  local_mocked_bindings(
    makeCluster = function(...) stop("Workers must not start"),
    .package = "GeDS"
  )
  dat <- data.frame(x = 1:10, y = 1:10)
  for (rule in c("GD", "G")) {
    expect_error(crossv_GeDS(
      y ~ f(x), data = dat, model_fun = NGeDSboost,
      parameters = list(phi_grid = 0.9), n_cores = 1L, stoptype = rule
    ), "Choose a fixed-threshold stopping rule to tune 'phi'", fixed = TRUE)
  }
})

test_that("boosting CV omits default parameter columns", {
  set.seed(81)
  dat <- data.frame(x = seq(0, 1, length.out = 40))
  dat$y <- sin(pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  cv <- crossv_GeDS(y ~ f(x), data = dat, model_fun = NGeDSboost,
                    n = 2L, n_folds = 2L,
                    parameters = list(shrinkage_grid = 0.6), n_cores = 1L,
                    initial_learner = FALSE, max_iterations = 10L)

  expect_true(is.finite(cv$results$cv_mse))
  expect_true("shrinkage" %in% names(cv$results))
  expect_false(any(c("int.knots_init", "beta", "phi", "q") %in%
                   names(cv$results)))

  cv_default <- crossv_GeDS(y ~ f(x), data = dat,
                            model_fun = NGeDSboost, n = 2L, n_folds = 2L,
                            n_cores = 1L, initial_learner = FALSE,
                            max_iterations = 10L)
  expect_true(is.finite(cv_default$results$cv_mse))
  expect_false(any(c("int.knots_init", "shrinkage", "beta", "phi", "q") %in%
                   names(cv_default$results)))
})

