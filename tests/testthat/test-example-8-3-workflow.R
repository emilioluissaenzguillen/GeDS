library("testthat")
library("GeDS")

test_that("reduced Example 8.3-style NGeDSboost classification workflow runs", {
  set.seed(123)

  # Reduced high-dimensional binary classification setting.
  # This is intentionally much smaller than the manuscript example,
  # so that it remains suitable for routine package testing.
  N <- 80
  P <- 30

  X <- matrix(rnorm(N * P), nrow = N, ncol = P)
  colnames(X) <- paste0("x", seq_len(P))

  eta <- 0.8 * X[, 1] - 0.6 * X[, 2] + 0.4 * X[, 3]
  prob <- plogis(eta)
  Y <- rbinom(N, size = 1, prob = prob)

  dat <- data.frame(
    Y = factor(Y, levels = c(0, 1)),
    X
  )

  set.seed(321)
  train_id <- sample(seq_len(N), size = floor(0.7 * N))
  train <- dat[train_id, ]
  test  <- dat[-train_id, ]

  form <- as.formula(
    paste("Y ~", paste(colnames(X), collapse = " + "))
  )

  Gmodboost <- NULL
  invisible(capture.output({
    Gmodboost <- suppressWarnings(
      NGeDSboost(
        form,
        data = train,
        family = mboost::Binomial(link = "logit"),
        phi_boost_exit = 0.8,
        phi = 0.8,
        initial_learner = FALSE
      )
    )
  }))

  pred <- predict(Gmodboost, newdata = test, n = 2, type = "response")

  y_test <- as.integer(as.character(test$Y))
  eps <- 1e-8
  pred_clip <- pmin(pmax(pred, eps), 1 - eps)

  mean_binomial_deviance <- -2 * mean(
    y_test * log(pred_clip) + (1 - y_test) * log(1 - pred_clip)
  )

  expect_s3_class(Gmodboost, "GeDSboost")
  expect_equal(length(pred), nrow(test))
  expect_true(all(is.finite(pred)))
  expect_true(all(pred >= 0 & pred <= 1))
  expect_true(is.finite(mean_binomial_deviance))
  expect_true(N.boost.iter(Gmodboost) >= 0)
})

