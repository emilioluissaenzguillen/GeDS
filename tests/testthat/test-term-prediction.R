test_that("term predictions retain their documented matrix shape", {
  x <- seq(-1, 1, length.out = 60)
  z <- rep(c(-1, 0, 1), length.out = length(x))
  data <- data.frame(x = x, z = z, y = 2 + sin(2 * x) + 0.3 * z)

  model <- suppressWarnings(NGeDS(
    y ~ f(x) + z, data = data, higher_order = FALSE,
    phi = 0.9, max.intknots = 4
  ))

  terms <- predict(model, newdata = data, type = "terms", n = 2)
  expect_true(is.matrix(terms))
  expect_equal(dim(terms), c(nrow(data), 2L))
  expect_equal(colnames(terms), c("Spline", "z"))
  expect_equal(unname(rowSums(terms)), predict(model, newdata = data, type = "link", n = 2),
               tolerance = 1e-10)
})
