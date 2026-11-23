library("testthat")
library("GeDS")

expect_same_knot_result <- function(cpp, ref, expected_knot,
                                    expected_weight, expected_bounds) {
  expect_equal(cpp$Dim.newknot, expected_knot, tolerance = 1e-14)
  expect_equal(ref$Dim.newknot, expected_knot, tolerance = 1e-14)
  expect_equal(cpp$weightDim, expected_weight, tolerance = 1e-14)
  expect_equal(ref$weightDim, expected_weight, tolerance = 1e-14)
  expect_false(cpp$flagDim)
  expect_false(ref$flagDim)
  expect_equal(c(cpp$dcumInf, cpp$dcumSup), expected_bounds)
  expect_equal(c(ref$dcumInf, ref$dcumSup), expected_bounds)
}

test_that("bivariate knot helpers agree on admissible residual clusters", {
  coordinate <- c(0.1, 0.2, 0.7, 0.8)
  residuals <- c(1, 2, -1, -2)

  fixtures <- list(
    weighted_cluster = list(
      cluster_ends = c(2L, 4L),
      cluster_weights = c(0.8, 0.3),
      old_knots = c(0, 1),
      coordinate = coordinate,
      residuals = residuals,
      expected_knot = 1 / 6,
      expected_weight = 0.8,
      expected_bounds = c(1, 2)
    ),
    blocked_best_cluster = list(
      cluster_ends = c(2L, 4L),
      cluster_weights = c(0.8, 0.3),
      old_knots = c(0, 0.15, 1),
      coordinate = coordinate,
      residuals = residuals,
      expected_knot = 23 / 30,
      expected_weight = 0.3,
      expected_bounds = c(3, 4)
    ),
    valid_singleton = list(
      cluster_ends = c(1L, 3L),
      cluster_weights = c(0.9, 0.2),
      old_knots = c(0, 1),
      coordinate = c(0.25, 0.6, 0.8),
      residuals = c(2, -1, -1),
      expected_knot = 0.25,
      expected_weight = 0.9,
      expected_bounds = c(1, 1)
    )
  )

  for (fixture in fixtures) {
    original_weights <- fixture$cluster_weights + 0
    cpp <- GeDS:::findNewDimKnot(
      fixture$cluster_ends,
      fixture$cluster_weights,
      fixture$old_knots,
      fixture$coordinate,
      fixture$residuals
    )
    ref <- GeDS:::findNewDimKnot_R(
      fixture$cluster_ends,
      fixture$cluster_weights,
      fixture$old_knots,
      fixture$coordinate,
      fixture$residuals
    )

    expect_equal(fixture$cluster_weights, original_weights)

    expect_same_knot_result(
      cpp, ref,
      fixture$expected_knot,
      fixture$expected_weight,
      fixture$expected_bounds
    )
  }
})

test_that("knot helpers flag the absence of an admissible cluster", {
  coordinate <- c(0.1, 0.2, 0.7, 0.8)
  residuals <- c(1, 2, -1, -2)

  cpp <- GeDS:::findNewDimKnot(
    c(2L, 4L),
    c(0.8, 0.3),
    c(0, 0.15, 0.75, 1),
    coordinate,
    residuals
  )
  ref <- GeDS:::findNewDimKnot_R(
    c(2L, 4L),
    c(0.8, 0.3),
    c(0, 0.15, 0.75, 1),
    coordinate,
    residuals
  )

  for (result in list(cpp, ref)) {
    expect_true(result$flagDim)
    expect_true(is.na(result$weightDim))
    expect_equal(result$Dim.newknot, 23 / 30, tolerance = 1e-14)
    expect_equal(c(result$dcumInf, result$dcumSup), c(3, 4))
  }

  expect_equal(ref, cpp, tolerance = 1e-14)
})

test_that("dimension-agnostic knot helpers validate vector lengths", {
  cpp <- GeDS:::findNewDimKnot
  ref <- GeDS:::findNewDimKnot_R

  for (helper in list(cpp, ref)) {
    expect_error(
      helper(c(2L, 4L), 0.8, c(0, 1), 1:4, c(1, 2, -1, -2)),
      "weights.*same length"
    )
    expect_error(
      helper(c(2L, 4L), c(0.8, 0.3), c(0, 1), 1:4, c(1, 2, -1)),
      "values.*same length"
    )
  }
})

test_that("dimension-indexed knot placement matches the corrected legacy path", {
  set.seed(2468)
  X <- runif(100, 0, 3)
  Y <- runif(100, 0, 3)
  residuals <- sin(2 * X) * cos(2 * Y) + rnorm(100, sd = 0.05)
  matr <- cbind(X = X, Y = Y, residual = residuals)
  nint <- as.integer(sqrt(NROW(matr)))

  strip_metadata <- function(values) {
    upper <- seq(0, 3, length.out = nint + 1L)[-1L] + 1e-15
    list(upper = upper)
  }

  metadata <- list(X = strip_metadata(X), Y = strip_metadata(Y))
  dimensions <- list(
    X = list(fixed = Y, order = order(Y, X), strips = metadata$Y),
    Y = list(fixed = X, order = order(X, Y), strips = metadata$X)
  )

  for (Dim in names(dimensions)) {
    args <- dimensions[[Dim]]
    common <- list(
      Dim = Dim,
      Dim.intknots = NULL,
      matr = matr,
      indicator = NULL,
      FixedDim = args$fixed,
      ordFixedDim = args$order,
      nintFixedDim = nint,
      upperFixedDim = args$strips$upper,
      beta = 0.5
    )

    current <- do.call(GeDS:::placeKnot, common)
    legacy <- do.call(GeDS:::placeKnot_biv_legacy, common)
    expect_equal(current, legacy, tolerance = 1e-14)
  }
})

test_that("repeated fixed-coordinate values retain one strip ID per row", {
  fixed.values <- c(0.1, 0.1, 0.8, 0.8)
  upper.bounds <- c(0.45, 0.8) + 1e-15

  expect_equal(
    GeDS:::fixedDimStripId(fixed.values, upper.bounds),
    c(1L, 1L, 2L, 2L)
  )

  X <- c(0.1, 0.2, 0.3, 0.4)
  residuals <- c(1, 2, -1, -2)
  matr <- cbind(X = X, Y = fixed.values, residual = residuals)
  common <- list(
    Dim = "X", Dim.intknots = NULL, matr = matr, indicator = NULL,
    FixedDim = fixed.values, ordFixedDim = order(fixed.values, X),
    nintFixedDim = 2L, upperFixedDim = upper.bounds, beta = 0.5
  )

  expect_equal(
    do.call(GeDS:::placeKnot, common),
    do.call(GeDS:::placeKnot_biv_legacy, common),
    tolerance = 1e-14
  )
})

test_that("fixed-coordinate cells extend strip membership beyond two dimensions", {
  fixed.coordinates <- rbind(
    c(0.1, 1),
    c(0.8, 1),
    c(0.1, 9),
    c(0.8, 9),
    c(0.8, 9)
  )
  upper.bounds <- list(
    c(0.5, 1.0) + 1e-15,
    c(5.0, 10.0) + 1e-15
  )

  expect_equal(
    GeDS:::fixedDimCellId(fixed.coordinates, upper.bounds),
    c(1L, 2L, 3L, 4L, 4L)
  )
  expect_equal(
    GeDS:::fixedDimCellId(fixed.coordinates[, 1, drop = FALSE], upper.bounds[1]),
    GeDS:::fixedDimStripId(fixed.coordinates[, 1], upper.bounds[[1]])
  )
})

test_that("placeKnotND places knots along every dimension of a 3D design", {
  coordinates <- as.matrix(expand.grid(
    X1 = c(0.1, 0.2, 0.7, 0.8),
    X2 = c(0.25, 0.75),
    X3 = c(2.5, 7.5)
  ))
  residuals <- rep(c(1, 2, -1, -2), 4L)
  bounds <- list(
    c(0.5, 1.0) + 1e-15,
    c(0.5, 1.0) + 1e-15,
    c(5.0, 10.0) + 1e-15
  )
  ranges <- list(c(0, 1), c(0, 1), c(0, 10))
  expected.knots <- c(23 / 30, 0.5, 5.0)

  for (target.index in seq_len(NCOL(coordinates))) {
    fixed.index <- setdiff(seq_len(NCOL(coordinates)), target.index)
    cell.id <- GeDS:::fixedDimCellId(
      coordinates[, fixed.index, drop = FALSE],
      bounds[fixed.index]
    )
    expected <- GeDS:::placeDimKnot(
      dim.index = target.index,
      intknots = NULL,
      coordinates = coordinates,
      residuals = residuals,
      strip.id = cell.id,
      beta = 0.5,
      dim.range = ranges[[target.index]]
    )
    actual <- GeDS:::placeKnotND(
      target.index = target.index,
      intknots = NULL,
      coordinates = coordinates,
      residuals = residuals,
      fixed.bounds = bounds[fixed.index],
      beta = 0.5,
      dim.range = ranges[[target.index]]
    )

    expect_equal(actual, expected, tolerance = 1e-14)
    expect_equal(actual$Dim.newknot, expected.knots[target.index],
                 tolerance = 1e-14)
    expect_false(actual$flagDim)
  }
})

test_that("placeKnotND validates its target and fixed-coordinate bounds", {
  coordinates <- cbind(X1 = 1:4, X2 = 5:8, X3 = 9:12)
  residuals <- c(1, 1, -1, -1)

  expect_error(
    GeDS:::placeKnotND(4L, NULL, coordinates, residuals,
                       list(8, 12), 0.5),
    "target.index"
  )
  expect_error(
    GeDS:::placeKnotND(1L, NULL, coordinates, residuals,
                       list(8), 0.5),
    "one vector per non-target"
  )
})

test_that("NGeDSgam passes explicit stopping rules to both Airquality smoothers", {
  dat <- utils::head(stats::na.omit(datasets::airquality), 40L)
  dat$Ozone <- dat$Ozone^(1 / 3)
  univariate <- UnivariateFitter
  bivariate <- NGeDS
  rules <- character()
  local_mocked_bindings(
    UnivariateFitter = function(..., stoptype) {
      rules <<- c(rules, paste0("univariate:", stoptype))
      univariate(..., stoptype = stoptype)
    },
    NGeDS = function(..., stoptype) {
      rules <<- c(rules, paste0("bivariate:", stoptype))
      bivariate(..., stoptype = stoptype)
    }
  )

  for (rule in c("RD", "SR")) {
    rules <- character()
    fit <- suppressWarnings(NGeDSgam(
      Ozone ~ f(Solar.R) + f(Wind, Temp), data = dat,
      stoptype = rule, internal_knots = 4L, max_iterations = 1L,
      higher_order = FALSE
    ))
    expect_setequal(rules, paste(c("univariate", "bivariate"), rule, sep = ":"))
    expect_true(all(is.finite(predict(fit, n = 2))))
  }

  expect_error(
    NGeDSgam(Ozone ~ f(Solar.R, Wind, Temp), data = dat, stoptype = "GD"),
    "at most two covariates"
  )
})

test_that("unsupported joint formulas do not resolve f in the formula environment", {
  set.seed(207)
  dat <- data.frame(x1 = runif(70), x2 = runif(70), x3 = runif(70))
  dat$y <- with(dat, x1 * x2 * x3 + rnorm(nrow(dat), sd = 0.05))
  isolated <- new.env(parent = baseenv())
  formula <- stats::as.formula("y ~ f(x1, x2, x3)", env = isolated)
  expect_error(NGeDS(formula, data = dat, max.intknots = 2L,
                     higher_order = FALSE), "at most two covariates")
  expect_false(exists("f", envir = isolated, inherits = TRUE))
})

test_that("Normal bivariate fitter uses all repeated-coordinate rows", {
  set.seed(123)
  X <- round(runif(400, min = 0, max = 3), 2)
  Y <- round(runif(400, min = 0, max = 3), 2)
  truth <- sin(2 * X) * sin(2 * Y)
  Z <- truth + rnorm(400, mean = 0, sd = 0.1)

  fit <- BivariateFitter(
    X = X, Y = Y, Z = Z, W = NULL, indicator = NULL,
    beta = 0.5, phi = 0.99, q = 2L,
    Xextr = c(0, 3), Yextr = c(0, 3),
    show.iters = FALSE
  )

  expect_equal(fit$iters, 14L)
  expect_equal(unlist(fit$Nintknots), c(X = 6L, Y = 5L))
  expect_equal(
    fit$linear.intknots$Xk,
    c(
      0.577359731594400, 0.690000000000000, 1.345157479213650,
      2.151389961254150, 2.450000000000000, 2.689592873309130
    ),
    tolerance = 1e-10
  )
  expect_equal(
    fit$linear.intknots$Yk,
    c(
      0.710000000000000, 0.919455489138973, 1.346201589686330,
      2.042629618542620, 2.517974268802460
    ),
    tolerance = 1e-10
  )
  expect_equal(
    c(fit$dev.linear, fit$dev.quadratic, fit$dev.cubic),
    c(4.10763189876677, 3.43691063876918, 3.39792646258526),
    tolerance = 1e-8
  )
  expect_equal(
    mean((truth - fit$quadratic.fit$predicted)^2),
    0.00170307701088141,
    tolerance = 1e-10
  )
})

test_that("Poisson bivariate fitter uses all repeated-coordinate rows", {
  set.seed(123)
  X <- round(runif(400, min = 0, max = 3), 2)
  Y <- round(runif(400, min = 0, max = 3), 2)
  truth <- exp(sin(2 * X) + sin(2 * Y))
  Z <- rpois(400, truth)

  fit <- GenBivariateFitter(
    X = X, Y = Y, Z = Z, W = NULL, family = poisson(),
    indicator = NULL, beta = 0.2, phi = 0.99, q = 2L,
    Xextr = c(0, 3), Yextr = c(0, 3),
    show.iters = FALSE
  )

  expect_equal(fit$iters, 9L)
  expect_equal(unlist(fit$Nintknots), c(X = 3L, Y = 3L))
  expect_equal(
    fit$linear.intknots$Xk,
    c(0.486723147249971, 1.696303216665780, 2.299802393514910),
    tolerance = 1e-10
  )
  expect_equal(
    fit$linear.intknots$Yk,
    c(0.194524675852381, 1.552392887826920, 1.870281149003880),
    tolerance = 1e-10
  )
  expect_equal(
    c(fit$dev.linear, fit$dev.quadratic, fit$dev.cubic),
    c(401.735849159871, 355.745631469214, 354.908531255439),
    tolerance = 1e-8
  )
  expect_equal(
    mean((truth - fit$quadratic.fit$predicted)^2),
    0.133785742941616,
    tolerance = 1e-10
  )
})

test_that("bivariate SR is the default and rank-deficient RD selections are avoided", {
  set.seed(501)
  train <- data.frame(x = runif(300), z = runif(300))
  train$y <- with(train, sin(2 * pi * x) * sin(2 * pi * z)) +
    rnorm(nrow(train), sd = 0.2)

  default <- NGeDS(y ~ f(x, z), data = train, Xextr = c(0, 1),
                   Yextr = c(0, 1), higher_order = FALSE)
  explicit.sr <- NGeDS(y ~ f(x, z), data = train, Xextr = c(0, 1),
                       Yextr = c(0, 1), stoptype = "SR", higher_order = FALSE)
  expect_equal(coef(default, n = 2L), coef(explicit.sr, n = 2L))
  expect_identical(default$linear.intknots, explicit.sr$linear.intknots)

  warnings <- character()
  rd <- withCallingHandlers(
    NGeDS(y ~ f(x, z), data = train, Xextr = c(0, 1),
          Yextr = c(0, 1), stoptype = "RD", higher_order = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  design <- GeDS:::tensorProd(rd$linear.fit$Xbasis,
                             rd$linear.fit$Ybasis)
  expect_equal(qr(design)$rank, NCOL(design))
  expect_true(any(grepl("rank deficient", warnings)))
  expect_lt(max(abs(coef(rd, n = 2L))), 10)
})

test_that("bivariate GAM and boosting learners select SR with the Stage A defaults", {
  d <- expand.grid(x = seq(0, 1, length.out = 8L), z = seq(0, 1, length.out = 8L))
  d$y <- with(d, sin(2 * pi * x) * sin(2 * pi * z) + x)
  native <- GeDS:::BivariateFitter
  settings <- list()
  local_mocked_bindings(BivariateFitter = function(..., stoptype, phi, q) {
    settings[[length(settings) + 1L]] <<- list(stoptype = stoptype, phi = phi, q = q)
    native(..., stoptype = stoptype, phi = phi, q = q)
  }, .package = "GeDS")
  expected <- list(stoptype = "SR", phi = 0.99, q = 2L)

  invisible(capture.output(gam <- suppressWarnings(NGeDSgam(
    y ~ f(x, z) + x, data = d, higher_order = FALSE, internal_knots = 6L
  ))))
  expect_true(all(is.finite(predict(gam, newdata = d, n = 2L))))
  expect_gt(length(settings), 0L)
  expect_true(all(vapply(settings, identical, logical(1), y = expected)))

  settings <- list()
  invisible(capture.output(boost <- suppressWarnings(NGeDSboost(
    y ~ f(x, z) + x, data = d, initial_learner = FALSE, higher_order = FALSE,
    max_iterations = 2L, int.knots_boost = 6L
  ))))
  expect_true(all(is.finite(predict(boost, newdata = d, n = 2L))))
  expect_gt(length(settings), 0L)
  expect_true(all(vapply(settings, identical, logical(1), y = expected)))
})
