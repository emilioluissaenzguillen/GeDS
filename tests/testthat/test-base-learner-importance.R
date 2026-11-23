# Check the public helper against independently computed losses, not its internals.
importance_reference <- function(fit, loss, boosting_only) {
  first <- if (fit$args$initial_learner && !boosting_only) 1L else 2L
  contributions <- c(loss(NULL) - loss(fit$models[[1L]]),
                     -diff(vapply(fit$models, loss, numeric(1))))
  selected <- vapply(fit$models, function(model) {
    if (is.null(model$best_bl$name)) NA_character_ else model$best_bl$name
  }, character(1))
  keep <- seq_along(fit$models) >= first
  learners <- names(fit$args$base_learners)
  expected <- vapply(learners, function(learner) {
    sum(contributions[keep & !is.na(selected) & selected == learner])
  }, numeric(1))
  list(values = expected, selected = selected[keep])
}

test_that("Gaussian importance attributes the actual mean-risk decrease", {
  dat <- data.frame(x = seq(-1, 1, length.out = 40))
  dat$z <- sin(seq_len(nrow(dat)))
  dat$y <- 3 + 2 * dat$x - 1.5 * dat$z + sin(5 * dat$x)

  for (initial in c(FALSE, TRUE)) {
    invisible(capture.output(fit <- NGeDSboost(
      y ~ x + z, data = dat, family = mboost::Gaussian(),
      initial_learner = initial, shrinkage = 0.2, phi_boost_exit = NULL,
      max_iterations = 4L, higher_order = FALSE
    )))
    loss <- function(model) {
      prediction <- if (is.null(model)) 0 else model$F_hat
      mean((dat$y - prediction)^2)
    }
    for (only in c(FALSE, TRUE)) {
      reference <- importance_reference(fit, loss, only)
      actual <- bl_imp(fit, boosting_iter_only = only)
      expect_equal(as.numeric(actual), unname(reference$values), tolerance = 1e-12)
      expect_equal(sum(actual), sum(reference$values), tolerance = 1e-12)
      expect_identical(names(actual), names(fit$args$base_learners))
      expect_s3_class(actual, "varimp")
      expect_equal(attr(actual, "selfreqs"), vapply(names(actual), function(v) {
        mean(reference$selected == v)
      }, numeric(1)), ignore_attr = TRUE)
      expect_length(attr(actual, "variable_names"), length(actual))
    }
  }
})

test_that("binomial importance agrees with stable Bernoulli log loss", {
  set.seed(927)
  dat <- data.frame(x = seq(-1, 1, length.out = 60))
  dat$z <- sin(seq_len(nrow(dat)))
  dat$y <- factor(rbinom(nrow(dat), 1, plogis(0.4 + dat$x - dat$z)))
  signed_y <- 2 * as.integer(as.character(dat$y)) - 1

  for (initial in c(FALSE, TRUE)) {
    invisible(capture.output(fit <- NGeDSboost(
      y ~ x + z, data = dat, family = mboost::Binomial(),
      initial_learner = initial, shrinkage = 0.2, phi_boost_exit = NULL,
      max_iterations = 4L, higher_order = FALSE
    )))
    loss <- function(model) {
      margin <- if (is.null(model)) 0 else model$F_hat
      a <- -2 * signed_y * margin
      mean(pmax(a, 0) + log1p(exp(-abs(a))))
    }
    for (only in c(FALSE, TRUE)) {
      reference <- importance_reference(fit, loss, only)
      actual <- bl_imp(fit, boosting_iter_only = only)
      expect_equal(as.numeric(actual), unname(reference$values), tolerance = 1e-12)
      expect_equal(sum(actual), sum(reference$values), tolerance = 1e-12)
      expect_true(all(is.finite(actual)))
    }
  }
})

test_that("a single boosting step has no artificial terminal risk reduction", {
  dat <- data.frame(x = seq(-1, 1, length.out = 30))
  dat$y <- 4 + 2 * dat$x
  invisible(capture.output(fit <- NGeDSboost(
    y ~ x, data = dat, initial_learner = FALSE, shrinkage = 0.2,
    phi_boost_exit = NULL, max_iterations = 1L, higher_order = FALSE
  )))
  expect_length(fit$models, 2L)
  expected <- mean((dat$y - fit$models$model0$F_hat)^2) -
    mean((dat$y - fit$models$model1$F_hat)^2)
  expect_equal(as.numeric(bl_imp(fit)), expected, tolerance = 1e-12)
  expect_equal(bl_imp(fit), bl_imp(fit, boosting_iter_only = TRUE))
})
