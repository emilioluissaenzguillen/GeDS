test_that("GAM selection uses the final local-scoring window, not the global minimum", {
  dat <- data.frame(x = seq_len(20L), y = rep(0:1, 10L))
  eta <- c(0, 1, 2, 3)
  local_mocked_bindings(backfitting = function(z, base_learners_list, iter, ...) {
    base_learners_list$x$coefficients <- list(b0 = eta[iter], b1 = 0)
    list(z_hat = rep(eta[iter], length(z)),
         base_learners_list = base_learners_list, n_iters = 1L)
  })
  fit <- function(q_gam) {
    invisible(capture.output(result <- NGeDSgam(y ~ x, data = dat,
      family = binomial(), q_gam = q_gam, min_iterations = 4L,
      max_iterations = 4L, higher_order = FALSE)))
    result
  }
  expected_deviance <- vapply(eta, function(e)
    sum(binomial()$dev.resids(dat$y, rep(plogis(e), nrow(dat)), rep(1, nrow(dat)))),
    numeric(1))
  expect_identical(which.min(expected_deviance), 1L)
  for (q_gam in c(1L, 2L, 5L)) {
    result <- fit(q_gam)
    candidates <- tail(seq_along(eta), q_gam + 1L)
    selected <- candidates[which.min(expected_deviance[candidates])]
    expect_equal(result$iters$local_scoring, 4L)
    expect_identical(result$final_model$model_name, paste0("iter", selected))
    expect_equal(result$final_model$dev, expected_deviance[selected])
    expect_equal(predict(result, n = 2L), rep(plogis(eta[selected]), nrow(dat)))
  }
  unrestricted <- fit(NULL)
  expect_identical(unrestricted$final_model$model_name, "iter1")
  expect_equal(unrestricted$final_model$dev, expected_deviance[1L])
})

test_that("a single Gaussian local-scoring fit remains eligible", {
  dat <- data.frame(x = seq(-1, 1, length.out = 20L))
  dat$y <- dat$x + sin(3 * dat$x)
  result <- suppressWarnings(NGeDSgam(y ~ f(x), data = dat,
    q_gam = 2L, higher_order = FALSE, stoptype = "RD"))
  expect_identical(result$final_model$model_name, "iter1")
  expect_equal(result$iters$local_scoring, 1L)
  expect_true(all(is.finite(predict(result, n = 2L))))
})
