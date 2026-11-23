library("testthat")
library("GeDS")

test_that("NGeDSboost uses serial fallback when workers are unavailable", {
  worker_count <- GeDS:::.ngedsboost_worker_count

  expect_identical(worker_count(999L, cores = 8L), 1L)
  expect_identical(worker_count(1000L, cores = NA_integer_), 1L)
  expect_identical(worker_count(1000L, cores = 1L), 1L)
  expect_identical(worker_count(1000L, cores = 2L), 1L)
  expect_identical(worker_count(1000L, cores = 3L), 2L)
  expect_identical(worker_count(1000L, cores = 4L), 3L)
  expect_identical(worker_count(1000L, cores = 16L), 3L)
})

test_that("parallel NGeDSboost preserves the first learner under ties", {
  skip_on_cran()
  skip_if(GeDS:::.ngedsboost_worker_count(1000L) == 1L)

  set.seed(921)
  x <- seq(-2, 2, length.out = 34)
  y <- sin(x) + rnorm(length(x), sd = 0.1)
  predictors <- cbind(matrix(rep(x, 900), nrow = length(x)),
                      matrix(0, nrow = length(x), ncol = 100))
  colnames(predictors) <- paste0("x", seq_len(ncol(predictors)))
  dat <- data.frame(y = y, predictors, check.names = FALSE)

  fit <- suppressWarnings(suppressMessages(NGeDSboost(
    y ~ ., data = dat, initial_learner = FALSE,
    min_iterations = 1L, max_iterations = 3L,
    higher_order = FALSE, shrinkage = 0.1
  )))

  selected <- vapply(fit$models[-1L], function(model) model$best_bl$name,
                     character(1))
  expect_identical(unname(selected), rep("x1", 3L))
  expect_true(all(is.finite(predict(fit, n = 2))))
})
