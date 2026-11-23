test_that("bivariate learner conversion preserves the selected fit, not its history", {
  set.seed(823)
  dat <- data.frame(x = runif(100), v = runif(100))
  dat$y <- sin(5 * dat$x) + dat$v^2
  for (cap in c(0L, 8L)) {
    fits <- list(
      suppressWarnings(NGeDS(y ~ f(x, v), data = dat, stoptype = "SR",
                             max.intknots = cap)),
      suppressWarnings(GGeDS(y ~ f(x, v), data = transform(dat, y = exp(y)),
                             family = Gamma(link = "log"), stoptype = "SR",
                             max.intknots = cap))
    )
    for (fit in fits) {
      expected <- as.numeric(predict(fit, newdata = dat, n = 2))
      theta <- as.numeric(fit$linear.fit$theta)
      # A rejected later fit must not be resurrected during learner conversion.
      fit$stored$previousX[] <- NA_real_
      fit$stored$previousY[] <- NA_real_
      fit$coefficients[] <- 1e15
      response <- if (fit$type == "GLM - Biv") exp(dat$y) else dat$y
      converted <- GeDS:::predict_GeDS_linear(fit, dat$x, dat$v, response)
      expect_equal(converted$theta, theta)
      expect_equal(converted$int.knt$Xint.knt, fit$linear.intknots$Xk)
      expect_equal(converted$int.knt$Yint.knt, fit$linear.intknots$Yk)
      expect_equal(converted$knt$Xk, sort(c(fit$args$Xextr, fit$linear.intknots$Xk)))
      expect_equal(converted$knt$Yk, sort(c(fit$args$Yextr, fit$linear.intknots$Yk)))
      expect_equal(as.numeric(converted$Y_hat), expected, tolerance = 1e-8)
    }
  }
})

test_that("airquality backfitting keeps the bivariate fitter's rank-safe selection", {
  dat <- na.omit(datasets::airquality)
  dat$Ozone <- dat$Ozone^(1/3)
  dat <- dat[order(dat$Month, dat$Day), ]
  cutoff <- round(.71 * nrow(dat))
  convert <- GeDS:::predict_GeDS_linear
  checked <- 0L
  local_mocked_bindings(predict_GeDS_linear = function(Gmod, X, Y, Z) {
    result <- convert(Gmod, X, Y, Z)
    if (identical(Gmod$type, "LM - Biv")) {
      checked <<- checked + 1L
      expect_equal(as.numeric(result$theta), as.numeric(Gmod$linear.fit$theta))
      expect_equal(result$int.knt$Xint.knt, Gmod$linear.intknots$Xk)
      expect_equal(result$int.knt$Yint.knt, Gmod$linear.intknots$Yk)
    }
    result
  }, .package = "GeDS")
  set.seed(135)
  invisible(capture.output(fit <- suppressWarnings(suppressMessages(
    NGeDSgam(Ozone ~ f(Solar.R) + f(Wind, Temp), data = dat[seq_len(cutoff), ],
              higher_order = FALSE)
  ))))
  expect_gt(checked, 0L)
  pred <- as.numeric(predict(fit, dat[-seq_len(cutoff), ], n = 2))
  expect_true(all(is.finite(pred)))
  expect_lt(mean((dat$Ozone[-seq_len(cutoff)] - pred)^2), 1)
  expect_equal(as.numeric(predict(fit, dat[seq_len(cutoff), ], n = 2)),
               as.numeric(predict(fit, n = 2)), tolerance = 1e-8)
})
