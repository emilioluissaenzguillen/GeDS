library("testthat")
library("GeDS")

test_that("GAM and boost Stage B drops knots only where averaging is impossible", {
  internal <- list(NULL, 0.3, c(0.3, 0.6), c(0.2, 0.5, 0.8))
  full <- lapply(internal, function(k) c(0, k, 1))
  learners <- list(
    univ = list(type = "GeDS", variables = "x1"),
    biv = list(type = "GeDS", variables = c("x1", "x2")),
    joint = list(type = "GeDS", variables = paste0("x", 1:4))
  )
  model <- list(base_learners = list(
    univ = list(knots = full[[2]]),
    biv = list(knots = list(Xk = full[[2]], Yk = full[[4]])),
    joint = list(knots = setNames(full, paste0("x", 1:4)))
  ))
  compute <- function(n) {
    result <- NULL
    invisible(capture.output(result <- GeDS:::compute_avg_int.knots(
      model, learners, X_sd = NULL, X_mean = NULL,
      normalize_data = FALSE, n = n
    )))
    result
  }
  linear <- compute(2)
  expect_equal(linear$univ, internal[[2]])
  expect_equal(linear$biv, list(ikX = internal[[2]], ikY = internal[[4]]))
  expect_equal(linear$joint, setNames(internal, paste0("x", 1:4)))
  for (n in c(3, 4)) {
    higher <- compute(n)
    expect_null(higher$univ)
    expect_null(higher$biv$ikX)
    expect_equal(higher$biv$ikY, GeDS:::makenewknots(internal[[4]], n))
    expect_equal(higher$joint, GeDS:::stageBKnotsND(linear$joint, n))
  }
})

test_that("IRIS - NGeDSgam predictions consistency", {
  # Prepare the iris dataset
  iris_subset <- subset(iris, Species %in% c("setosa", "versicolor"))
  iris_subset$Species <- factor(iris_subset$Species)

  # Compute ranges from iris_subset
  ranges <- lapply(iris_subset[, 1:4], range)

  # Filter iris for observations within these ranges
  within_ranges <- with(iris,
                        Sepal.Length >= ranges$Sepal.Length[1] & Sepal.Length <= ranges$Sepal.Length[2] &
                          Sepal.Width  >= ranges$Sepal.Width[1]  & Sepal.Width  <= ranges$Sepal.Width[2] &
                          Petal.Length >= ranges$Petal.Length[1] & Petal.Length <= ranges$Petal.Length[2] &
                          Petal.Width  >= ranges$Petal.Width[1]  & Petal.Width  <= ranges$Petal.Width[2]
  )
  iris_in_range <- iris[within_ranges, ]
  combined <- rbind(iris_subset, iris_in_range)
  new_rows <- !duplicated(combined)[(nrow(iris_subset) + 1):nrow(combined)]
  iris_in_range_new <- iris_in_range[new_rows, ]

  for (normalize in c(TRUE, FALSE)) {
    # Run NGeDSgam silently
    invisible(capture.output({
      Gmodgam <- suppressWarnings(
        NGeDSgam(Species ~ f(Sepal.Length, Sepal.Width) + f(Petal.Length) + Petal.Width,
                 data = iris_subset, family = binomial(link = "cauchit"), normalize_data = normalize,
                 phi_gam_exit = 0.8, penalty = "none")
        )
      }))

    # Check that prediction differences are essentially zero
    expect_equal(
      predict(Gmodgam, newdata = iris_subset, type = "response", n = 2),
      Gmodgam$predictions$pred_linear,
      tolerance = 1e-6
    )

    expect_equal(
      predict(Gmodgam, newdata = iris_subset, type = "response", n = 3),
      Gmodgam$predictions$pred_quadratic,
      tolerance = 1e-6
    )

    expect_equal(
      predict(Gmodgam, newdata = iris_subset, type = "response", n = 4),
      Gmodgam$predictions$pred_cubic,
      tolerance = 1e-6
    )

    for (ord in 2:4) {
      expect_equal(
        predict(Gmodgam, newdata = iris_subset, type = "response", n = ord),
        predict(Gmodgam, newdata = rbind(iris_subset, iris_in_range_new), type = "response", n = ord)[1:nrow(iris_subset)],
        tolerance = 1e-6
      )
    }
  }
})

test_that("joint NGeDSgam predicts after an SR-based learner fit", {
  set.seed(123)
  data <- data.frame(
    x1 = runif(400), x2 = runif(400),
    x3 = runif(400), x4 = runif(400)
  )
  data$y <- with(data,
                 2 * sin(pi * x1) * sin(pi * x2) *
                   sin(pi * x3) * sin(pi * x4)) +
    rnorm(nrow(data), sd = 0.15)

  invisible(capture.output(
    fit <- suppressWarnings(NGeDSgam(
      y ~ f(x1, x2, x3, x4), data = data, penalty = "none"
    ))
  ))

  expect_gt(sum(lengths(fit$internal_knots$linear.int.knots[[1]])), 0L)
  expect_equal(predict(fit, newdata = data, n = 2),
               fit$predictions$pred_linear)
  expect_length(predict(fit, newdata = data, n = 4), nrow(data))
})

test_that("joint NGeDSboost prediction retains learners above two dimensions", {
  set.seed(713)
  data <- data.frame(
    x1 = runif(120), x2 = runif(120),
    x3 = runif(120), x4 = runif(120)
  )
  data$y <- with(data, sin(2 * pi * x1) + x2 * x3 - x4)

  invisible(capture.output(
    fit <- suppressWarnings(NGeDSboost(
      y ~ f(x1, x2, x3, x4), data = data,
      higher_order = FALSE
    ))
  ))

  expect_equal(
    predict(fit, newdata = data, n = 2),
    fit$predictions$pred_linear,
    tolerance = 1e-7
  )
  expect_true(all(is.finite(predict(fit, newdata = data, n = 2))))
})

test_that("omitted joint Stage B orders fail clearly at prediction time", {
  set.seed(123)
  data <- data.frame(
    x1 = runif(80), x2 = runif(80),
    x3 = runif(80), x4 = runif(80)
  )
  data$y <- with(data,
                 sin(2 * pi * x1) * sin(2 * pi * x2) * sin(2 * pi * x3))

  invisible(capture.output(
    fit <- suppressWarnings(NGeDSboost(
      y ~ f(x1, x2, x3, x4), data = data,
      initial_learner = FALSE
    ))
  ))

  expect_null(fit$final_model$cubic.fit)
  expect_error(predict(fit, newdata = data, n = 4),
               "joint tensor basis requires")
  expect_error(predict(fit, n = 4),
               "joint tensor basis requires")
})

test_that("IRIS - NGeDSboost predictions consistency", {
  # Prepare the iris dataset
  iris_subset <- subset(iris, Species %in% c("setosa", "versicolor"))
  iris_subset$Species <- factor(iris_subset$Species)

  # Compute ranges from iris_subset
  ranges <- lapply(iris_subset[, 1:4], range)

  # Filter iris for observations within these ranges
  within_ranges <- with(iris,
                        Sepal.Length >= ranges$Sepal.Length[1] & Sepal.Length <= ranges$Sepal.Length[2] &
                          Sepal.Width  >= ranges$Sepal.Width[1]  & Sepal.Width  <= ranges$Sepal.Width[2] &
                          Petal.Length >= ranges$Petal.Length[1] & Petal.Length <= ranges$Petal.Length[2] &
                          Petal.Width  >= ranges$Petal.Width[1]  & Petal.Width  <= ranges$Petal.Width[2]
  )
  iris_in_range <- iris[within_ranges, ]
  combined <- rbind(iris_subset, iris_in_range)
  new_rows <- !duplicated(combined)[(nrow(iris_subset) + 1):nrow(combined)]
  iris_in_range_new <- iris_in_range[new_rows, ]

  for (normalize in c(TRUE, FALSE)) {
    for (init_learner in c(TRUE, FALSE)) {

      # Run NGeDSboost silently
      invisible(capture.output({
        Gmodboost <- suppressWarnings(
          NGeDSboost(Species ~ f(Sepal.Length, Sepal.Width) + f(Petal.Length) + f(Petal.Width),
                     data = iris_subset, family = mboost::Binomial(link = "probit"),
                     initial_learner = init_learner, normalize_data = normalize,
                     phi_boost_exit = 0.8, penalty = "none")
          )
        }))

      # Check that prediction differences are essentially zero
      expect_equal(
        predict(Gmodboost, newdata = iris_subset, type = "response", n = 2),
        Gmodboost$predictions$pred_linear,
        tolerance = 1e-6
      )

      expect_equal(
        predict(Gmodboost, newdata = iris_subset, type = "response", n = 3),
        Gmodboost$predictions$pred_quadratic,
        tolerance = 1e-6
      )

      expect_equal(
        predict(Gmodboost, newdata = iris_subset, type = "response", n = 4),
        Gmodboost$predictions$pred_cubic,
        tolerance = 1e-6
      )

      for (ord in 2:4) {
        expect_equal(
          predict(Gmodboost, newdata = iris_subset, type = "response", n = ord),
          predict(Gmodboost, newdata = rbind(iris_subset, iris_in_range_new), type = "response", n = ord)[1:nrow(iris_subset)],
          tolerance = 1e-6
        )
      }

    }
  }
})

### MTCARS - NGeDSgam Checks
test_that("MTCARS - NGeDSgam predictions consistency", {
  data(mtcars)
  # Convert specified variables to factors
  categorical_vars <- c("cyl", "vs", "am", "gear", "carb")
  mtcars[categorical_vars] <- lapply(mtcars[categorical_vars], factor)

  for (normalize in c(TRUE, FALSE)) {
    invisible(capture.output({
      Gmodgam <- suppressWarnings(
        NGeDSgam(mpg ~ cyl + f(disp, hp) + f(drat) + f(wt) + f(qsec) + vs + am + gear + carb,
                 data = mtcars, family = gaussian, normalize_data = normalize,
                 penalty = "none")
        )
      }))

    # Check prediction differences for orders 2, 3, and 4
    expect_equal(
      predict(Gmodgam, newdata = mtcars, type = "response", n = 2),
      Gmodgam$predictions$pred_linear,
      tolerance = 1e-6
    )
    expect_equal(
      predict(Gmodgam, newdata = mtcars, type = "response", n = 3),
      Gmodgam$predictions$pred_quadratic,
      tolerance = 1e-6
    )
    expect_equal(
      predict(Gmodgam, newdata = mtcars, type = "response", n = 4),
      Gmodgam$predictions$pred_cubic,
      tolerance = 1e-6
    )

    # Check that the sum of base learner contributions equals the overall prediction
    for (ord in 2:4) {
      pred1 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "cyl")
      pred2 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "f(disp, hp)")
      pred3 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "f(drat)")
      pred4 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "f(wt)")
      pred5 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "f(qsec)")
      pred6 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "vs")
      pred7 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "am")
      pred8 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "gear")
      pred9 = predict(Gmodgam, n = ord, newdata = mtcars,  base_learner = "carb")

      # The combined B-spline coefficients are not uniquely identifiable for
      # this rank-deficient model. Use the offset from the fitted backfitting
      # representation, which is also used by the full prediction method.
      pred0 <- if (ord == 2) mean(Gmodgam$final_model$Y_hat$z) else 0
      sum <-  pred0 + pred1 + pred2 + pred3 + pred4 + pred5 + pred6 + pred7 + pred8 + pred9

      if (Gmodgam$args$normalize_data && ord == 2) {
        sum <- sum * Gmodgam$args$Y_sd + Gmodgam$args$Y_mean
      }

      expect_equal(
        sum,
        predict(Gmodgam, newdata = mtcars, type = "response", n = ord),
        tolerance = 1e-6
      )
    }
  }
})

### MTCARS - NGeDSboost Checks
test_that("MTCARS - NGeDSboost predictions consistency", {
  data(mtcars)
  # Convert specified variables to factors
  categorical_vars <- c("cyl", "vs", "am", "gear", "carb")
  mtcars[categorical_vars] <- lapply(mtcars[categorical_vars], factor)

  for (normalize in c(TRUE, FALSE)) {
    for (init_learner in c(TRUE, FALSE)) {
      invisible(capture.output({
        Gmodboost <- suppressWarnings(
          NGeDSboost(mpg ~ cyl + f(disp, hp) + f(drat) + f(wt) + f(qsec) + vs + am + gear + carb,
                     data = mtcars, family = mboost::Gaussian(), initial_learner = init_learner,
                     normalize_data = normalize, max_iterations = 3L, penalty = "none")
          )
        }))

      # The tensor-size safeguard can omit higher orders in this 32-row model.
      fit_keys <- c("linear.fit", "quadratic.fit", "cubic.fit")
      pred_keys <- c("pred_linear", "pred_quadratic", "pred_cubic")
      available <- !vapply(Gmodboost$final_model[fit_keys], is.null, logical(1))
      expect_true(available[1])
      for (ord in (2:4)[!available]) {
        expect_error(
          predict(Gmodboost, newdata = mtcars, type = "response", n = ord),
          "joint tensor basis requires.*positive weight"
        )
      }
      for (ord in (2:4)[available]) {
        expect_equal(
          predict(Gmodboost, newdata = mtcars, type = "response", n = ord),
          Gmodboost$predictions[[pred_keys[ord - 1]]],
          tolerance = 1e-6
        )
      }

      # Check that the sum of base learner contributions equals the overall prediction
      for (ord in (2:4)[available]) {
        pred1 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "cyl")
        pred2 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "f(disp, hp)")
        pred3 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "f(drat)")
        pred4 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "f(wt)")
        pred5 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "f(qsec)")
        pred6 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "vs")
        pred7 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "am")
        pred8 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "gear")
        pred9 = predict(Gmodboost, n = ord, newdata = mtcars,  base_learner = "carb")

        if(!Gmodboost$args$initial_learner && !Gmodboost$args$normalize_data && ord == 2) {
          pred0 <- mean(mtcars$mpg)
        } else {
          pred0 <- 0
        }

        sum <- pred1+pred2+pred3+pred4+pred5+pred6+pred7+pred8+pred9

        if (ord == 2) {
          b0 <- unname(Gmodboost$final_model$linear.fit$theta["b0"])
          if (is.na(b0)) b0 <- 0
        } else {
          b0 <- 0
        }

        if (!Gmodboost$args$initial_learner && ord == 2) {
          pred0 <- mean(mtcars$mpg)
          if (Gmodboost$args$normalize_data) {
            pred0 <- (pred0 - Gmodboost$args$Y_mean) / Gmodboost$args$Y_sd
          }
        } else {
          pred0 <- 0
        }

        sum <- pred0 + b0 + sum

        if (Gmodboost$args$normalize_data && ord == 2) {
          sum <- sum * Gmodboost$args$Y_sd + Gmodboost$args$Y_mean
        }

        expect_equal(
          sum,
          predict(Gmodboost, newdata = mtcars, type = "response", n = ord),
          tolerance = 1e-6
        )
      }
    }
  }
})

test_that("joint Stage B keeps an early prefix of fitted knots", {
  internal <- list(x1 = c(0.2, 0.5), x2 = c(0.3, 0.7),
                   x3 = numeric(), x4 = numeric())
  order <- list(
    list(dimension = "x2", knot = 0.7),
    list(dimension = "x1", knot = 0.5),
    list(dimension = "x1", knot = 0.2),
    list(dimension = "x2", knot = 0.3)
  )
  chosen <- GeDS:::.stageB_joint_prefix(internal, order, 4L, 1000L)
  expect_equal(chosen, list(x1 = c(0.2, 0.5), x2 = 0.7,
                            x3 = NULL, x4 = NULL))
  expect_equal(internal$x1, c(0.2, 0.5))
})

test_that("joint boosting retains its linear knots and fits a smaller cubic tensor", {
  set.seed(123)
  data <- data.frame(
    x1 = runif(1000), x2 = runif(1000),
    x3 = runif(1000), x4 = runif(1000)
  )
  data$y <- with(data,
                 2 * sin(pi * x1) * sin(pi * x2) *
                   sin(pi * x3) * sin(pi * x4) *
                   (1 + 0.15 * sin(4 * pi * x1) * sin(4 * pi * x2))) +
    rnorm(nrow(data), sd = 0.15)
  set.seed(99)
  invisible(capture.output(
    fit <- suppressWarnings(NGeDSboost(
      y ~ f(x1, x2, x3, x4), data = data,
      initial_learner = FALSE, max_iterations = 3
    ))
  ))
  set.seed(99)
  invisible(capture.output(
    linear_only <- suppressWarnings(NGeDSboost(
      y ~ f(x1, x2, x3, x4), data = data,
      initial_learner = FALSE, max_iterations = 3,
      higher_order = FALSE
    ))
  ))
  expect_gt(sum(lengths(fit$internal_knots$linear.int.knots[[1]])), 1L)
  expect_equal(lengths(fit$internal_knots$linear.int.knots[[1]]),
               c(x1 = 1L, x2 = 1L, x3 = 1L, x4 = 1L))
  expect_identical(fit$internal_knots$linear.int.knots,
                   linear_only$internal_knots$linear.int.knots)
  expect_equal(predict(fit, newdata = data, n = 2),
               predict(linear_only, newdata = data, n = 2))
  expect_lt(prod(lengths(fit$internal_knots$cubic.int.knots[[1]]) + 4L),
            nrow(data))
  expect_false(is.null(fit$final_model$cubic.fit))
  expect_length(predict(fit, newdata = data, n = 4), nrow(data))
})
