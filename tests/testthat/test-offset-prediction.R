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

test_that("generalized spline fitting passes offsets on the link scale", {
  x <- seq(-1.5, 1.5, length.out = 90)
  exposure <- seq(1.1, 2, length.out = 90)
  off <- log(exposure)
  z <- rep(c(-1, 0, 1), length.out = length(x))
  set.seed(321)
  y <- rpois(length(x), exposure * exp(0.4 + sin(x) + 0.2 * z))
  expect_true(any(y == 0))

  model <- suppressWarnings(GGeDS(
    y ~ f(x) + z + offset(off),
    data = data.frame(x, z, off, y), family = poisson(),
    beta = 0.2, phi = 0.9, max.intknots = 4, higher_order = FALSE
  ))
  expect_true(all(is.finite(predict(model, n = 2))))

  spline <- splines::splineDesign(
    knots = sort(c(0, rep(range(x), each = 2))), x = x,
    ord = 2, outer.ok = TRUE
  )
  design <- cbind(spline, z)
  direct <- stats::glm.fit(design, y, family = poisson(), offset = off,
                           intercept = FALSE)
  fast <- SplineReg_GLM(
    X = x, Y = y, Z = z, offset = off, InterKnots = 0,
    n = 2, family = poisson(), fast = TRUE
  )
  full <- SplineReg_GLM(
    X = x, Y = y, Z = z, offset = off, InterKnots = 0,
    n = 2, family = poisson()
  )
  expect_equal(unname(fast$theta), unname(direct$coefficients),
               tolerance = 1e-8)
  expect_equal(unname(full$theta), unname(direct$coefficients),
               tolerance = 1e-8)
  expect_equal(as.numeric(full$predicted), as.numeric(direct$fitted.values),
               tolerance = 1e-8)
})

test_that("mortality exposure offsets fit with generalized starting values", {
  model <- GGeDS(
    Deaths ~ f(Age) + offset(log(Exposure)), data = EWmortality,
    family = quasipoisson(), phi = 0.99, beta = 0.1, q = 3,
    stoptype = "LR"
  )
  expect_true(all(is.finite(predict(model, n = 3))))
})
