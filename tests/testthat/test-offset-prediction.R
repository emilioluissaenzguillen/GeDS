test_that("formula offsets are counted once in univariate predictions", {
  x <- seq(-1, 1, length.out = 60)
  off <- seq(0.1, 0.5, length.out = 60)
  data <- data.frame(x = x, off = off, y = 2 + sin(2 * x) + off)

  model <- suppressWarnings(NGeDS(
    y ~ f(x) + offset(off), data = data, higher_order = FALSE,
    phi = 0.9, max.intknots = 4
  ))

  baseline <- predict(model, newdata = data, n = 2)
  shifted <- transform(data, off = off + 0.5)
  expect_equal(predict(model, newdata = shifted, n = 2) - baseline,
               rep(0.5, nrow(data)), tolerance = 1e-10)
  expect_equal(predict(model, n = 2), baseline, tolerance = 1e-10)
})
