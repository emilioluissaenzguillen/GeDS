library("testthat")
library("GeDS")

test_that("shapeConstrain.GeDS imposes increasing fitted values", {
  skip_if_not_installed("quadprog")

  set.seed(123)
  N <- 120
  X <- sort(runif(N, 0, 1))
  Y <- exp(2 * X) + rnorm(N, sd = 0.15)
  dat <- data.frame(X = X, Y = Y)

  Gmod <- suppressWarnings(
    NGeDS(Y ~ f(X), data = dat, beta = 0.6, phi = 0.99)
  )

  Gmod_sc <- shapeConstrain(Gmod, n = 3, shape_constraint = "increasing")
  pred_sc <- predict(Gmod_sc, newdata = dat, n = 3)

  expect_s3_class(Gmod_sc, "shapeConstrainedGeDS")
  expect_equal(Gmod_sc$shape_constraint$constraint, "increasing")
  expect_true(all(diff(pred_sc) >= -1e-8))
  expect_true(all(is.finite(pred_sc)))
})

test_that("shapeConstrain.GeDS rejects incompatible constraints", {
  skip_if_not_installed("quadprog")

  set.seed(123)
  X <- sort(runif(60, 0, 1))
  Y <- exp(X) + rnorm(60, sd = 0.1)
  dat <- data.frame(X = X, Y = Y)

  Gmod <- suppressWarnings(
    NGeDS(Y ~ f(X), data = dat, beta = 0.6, phi = 0.95)
  )

  expect_error(
    shapeConstrain(
      Gmod,
      n = 3,
      shape_constraint = c("increasing", "decreasing")
    ),
    "Cannot impose both increasing and decreasing constraints"
  )
})

test_that("shapeConstrain.GeDSgam constrains a selected smoother", {
  skip_if_not_installed("quadprog")

  set.seed(123)
  N <- 80
  X <- sort(runif(N, 0, 1))
  Y <- exp(2 * X) + rnorm(N, sd = 0.15)
  dat <- data.frame(X = X, Y = Y)

  Gmodgam <- suppressWarnings(
    NGeDSgam(Y ~ f(X), data = dat, beta = 0.6, phi = 0.95,
             higher_order = FALSE)
  )

  bl_name <- names(Gmodgam$final_model$base_learners)[1L]
  Gmodgam_sc <- shapeConstrain(
    Gmodgam,
    n = 2,
    shape_constraint = "increasing",
    base_learner = bl_name
  )

  pred_sc <- predict(Gmodgam_sc, newdata = dat, n = 2,
                     base_learner = bl_name)

  expect_s3_class(Gmodgam_sc, "shapeConstrainedGeDSgam")
  expect_equal(Gmodgam_sc$shape_constraint$base_learner, bl_name)
  expect_true(all(diff(pred_sc) >= -1e-8))
  expect_true(all(is.finite(pred_sc)))
})

test_that("shapeConstrain.GeDSboost constrains a selected base learner", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("mboost")

  set.seed(123)
  N <- 80
  X <- sort(runif(N, 0, 1))
  Y <- exp(2 * X) + rnorm(N, sd = 0.15)
  dat <- data.frame(X = X, Y = Y)

  Gmodboost <- suppressWarnings(
    NGeDSboost(Y ~ f(X), data = dat, beta = 0.6, phi = 0.95,
               max_iterations = 2L, higher_order = FALSE)
  )

  bl_name <- names(Gmodboost$final_model$base_learners)[1L]
  Gmodboost_sc <- shapeConstrain(
    Gmodboost,
    n = 2,
    shape_constraint = "increasing",
    base_learner = bl_name
  )

  pred_sc <- predict(Gmodboost_sc, newdata = dat, n = 2,
                     base_learner = bl_name)

  expect_s3_class(Gmodboost_sc, "shapeConstrainedGeDSboost")
  expect_equal(Gmodboost_sc$shape_constraint$base_learner, bl_name)
  expect_true(all(diff(pred_sc) >= -1e-8))
  expect_true(all(is.finite(pred_sc)))
})
