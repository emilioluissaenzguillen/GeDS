library(testthat)
library(GeDS) 

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
                 phi_gam_exit = 0.8)
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
                     phi_boost_exit = 0.8)
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
                 data = mtcars, family = gaussian, normalize_data = normalize)
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

      if (ord == 2) {
        b0 = Gmodgam$final_model$Linear.Fit$Theta["b0"]
        alpha = if(Gmodgam$args$normalize_data) 0 else mean(mtcars$mpg)
        pred0 = alpha + b0
      } else {
        pred0 = 0
      }
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
                     normalize_data = normalize)
          )
        }))

      # Check prediction differences for orders 2, 3, and 4
      expect_equal(
        predict(Gmodboost, newdata = mtcars, type = "response", n = 2),
        Gmodboost$predictions$pred_linear,
        tolerance = 1e-6
      )
      expect_equal(
        predict(Gmodboost, newdata = mtcars, type = "response", n = 3),
        Gmodboost$predictions$pred_quadratic,
        tolerance = 1e-6
      )
      expect_equal(
        predict(Gmodboost, newdata = mtcars, type = "response", n = 4),
        Gmodboost$predictions$pred_cubic,
        tolerance = 1e-6
      )

      # Check that the sum of base learner contributions equals the overall prediction
      for (ord in 2:4) {
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

        sum <- pred0 + pred1+pred2+pred3+pred4+pred5+pred6+pred7+pred8+pred9

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

