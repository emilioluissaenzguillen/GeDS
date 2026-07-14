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
})

