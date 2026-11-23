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

test_that("multiple parametric contributions reconstruct univariate and bivariate links", {
  set.seed(123)
  data <- expand.grid(x = seq(0, 1, length.out = 10), v = seq(0, 1, length.out = 10))
  data$z1 <- rnorm(nrow(data))
  data$z2 <- rnorm(nrow(data))
  data$group <- factor(rep(c("A", "B", "C", "B"), length.out = nrow(data)))
  data$off <- runif(nrow(data), -.2, .2)
  data$w <- rep(c(.5, 1, 2, 0), length.out = nrow(data))
  eta <- 1 + sin(2 * data$x) + .3 * data$v + .3 * data$z1 - .2 * data$z2 +
    .2 * (data$group == "B") - .1 * (data$group == "C") + data$off
  for (dimension in 1:2) for (generalized in c(FALSE, TRUE)) {
    data$y <- if (generalized) rpois(nrow(data), exp(eta)) else eta + rnorm(nrow(data), sd = .1)
    formula <- if (dimension == 1L) y ~ f(x) + z1 + z2 + group + offset(off) else
      y ~ f(x, v) + z1 + z2 + group + offset(off)
    args <- list(formula = formula, data = data, weights = data$w,
                 phi = .9, q = 1, max.intknots = 4)
    if (generalized) args$family <- poisson()
    fit <- suppressWarnings(do.call(if (generalized) GGeDS else NGeDS, args))
    for (n in 2:4) {
      terms <- predict(fit, data, n = n, type = "terms")
      expect_true(is.matrix(terms))
      expect_equal(dim(terms), c(nrow(data), 5L))
      expect_equal(colnames(terms), c("Spline", "z1", "z2", "groupB", "groupC"))
      expect_equal(unname(rowSums(terms)) + data$off,
                   predict(fit, data, n = n, type = "link"), tolerance = 1e-10)
      beta <- coef(fit, n = n, onlySpline = FALSE)[fit$znames]
      design <- model.matrix(~ z1 + z2 + group, data)[, fit$znames, drop = FALSE]
      expect_equal(terms[, fit$znames, drop = FALSE], sweep(design, 2L, beta, "*"),
                   ignore_attr = TRUE, tolerance = 1e-10)
      training_data <- if (dimension == 1L) data[order(data$x), ] else data
      training_terms <- predict(fit, n = n, type = "terms")
      expect_equal(unname(rowSums(training_terms)) + fit$args$offset,
                   predict(fit, n = n, type = "link"), tolerance = 1e-10)
      expect_equal(unname(training_terms),
                   unname(predict(fit, training_data, n = n, type = "terms")))
      for (index in list(c(19, 3), c(19, 3, 45))) {
        batch <- data[index, , drop = FALSE]
        batch_terms <- predict(fit, batch, n = n, type = "terms")
        individual <- do.call(rbind, lapply(seq_len(nrow(batch)), function(i) {
          predict(fit, batch[i, , drop = FALSE], n = n, type = "terms")
        }))
        expect_equal(unname(batch_terms), unname(individual), tolerance = 1e-10)
        expect_equal(unname(batch_terms), unname(terms[index, , drop = FALSE]))
        reverse <- predict(fit, batch[nrow(batch):1, , drop = FALSE], n = n, type = "terms")
        expect_equal(unname(batch_terms), unname(reverse[nrow(batch):1, , drop = FALSE]))
      }
      shifted <- transform(data, off = off + .25)
      expect_equal(predict(fit, shifted, n = n, type = "terms"), terms)
      expect_equal(predict(fit, shifted, n = n, type = "link") -
                     predict(fit, data, n = n, type = "link"), rep(.25, nrow(data)))
    }
  }
})

test_that("spline-only term matrices reconstruct links with offsets", {
  data <- expand.grid(x = seq(0, 1, length.out = 8), v = seq(0, 1, length.out = 8))
  data$off <- .1 + .2 * data$v
  data$y <- 2 + sin(2 * data$x) + .3 * data$v + data$off
  for (formula in list(y ~ f(x) + offset(off), y ~ f(x, v) + offset(off))) {
    fit <- suppressWarnings(NGeDS(formula, data = data, phi = .9, q = 1, max.intknots = 4))
    for (n in 2:4) {
      terms <- predict(fit, data, n = n, type = "terms")
      expect_equal(dim(terms), c(nrow(data), 1L))
      expect_identical(colnames(terms), "Spline")
      expect_equal(unname(rowSums(terms)) + data$off,
                   predict(fit, data, n = n, type = "link"), tolerance = 1e-10)
    }
  }
})
