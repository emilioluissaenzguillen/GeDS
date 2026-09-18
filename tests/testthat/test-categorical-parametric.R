library("testthat")
library("GeDS")

test_that("NGeDS supports categorical parametric covariates", {
  data <- mtcars
  data$cyl <- factor(data$cyl)

  model <- suppressWarnings(
    NGeDS(mpg ~ f(wt) + cyl, data = data, phi = 0.9)
  )

  expect_identical(model$znames, c("cyl6", "cyl8"))
  expect_true(is.numeric(model$args$Z))
  expect_identical(colnames(model$args$Z), model$znames)

  coefficients <- coef(model, n = 2, onlySpline = FALSE)
  expect_true(all(model$znames %in% names(coefficients)))

  fitted <- predict(model, n = 2)
  predicted <- predict(model, newdata = data, n = 2)
  expect_equal(predicted[order(data$wt)], fitted, tolerance = 1e-10)

  terms <- predict(model, newdata = data, n = 2, type = "terms")
  expect_length(terms, nrow(data) * (length(model$znames) + 1L))
})

test_that("GGeDS uses the same factor encoding for new data", {
  data <- mtcars
  data$cyl <- factor(data$cyl)
  data$count <- round(data$mpg)

  model <- suppressWarnings(
    GGeDS(count ~ f(wt) + cyl, data = data, family = poisson(),
          beta = 0.2, phi = 0.9)
  )

  expect_identical(colnames(model$args$Z), model$znames)
  expect_equal(
    predict(model, newdata = data, n = 2)[order(data$wt)],
    predict(model, n = 2),
    tolerance = 1e-8
  )
})

test_that("parametric interactions use their encoded design columns", {
  data <- mtcars
  data$cyl <- factor(data$cyl)

  parsed <- GeDS:::read.formula(mpg ~ f(wt) + cyl:am, data)

  expect_identical(colnames(parsed$Z), c("cyl4:am", "cyl6:am", "cyl8:am"))
  expect_true(is.numeric(parsed$Z))
})
