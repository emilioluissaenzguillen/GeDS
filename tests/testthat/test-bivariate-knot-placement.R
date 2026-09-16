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

test_that("a two-dimensional Stage A step matches the bivariate components", {
  set.seed(987)
  X <- runif(100, 0, 3)
  Y <- runif(100, 0, 3)
  Z <- sin(2 * X) * cos(2 * Y) + rnorm(100, sd = 0.05)
  coordinates <- cbind(X = X, Y = Y)
  nint <- as.integer(sqrt(NROW(coordinates)))
  upper <- seq(0, 3, length.out = nint + 1L)[-1L] + 1e-15

  generic <- GeDS:::stageAOneStepND(
    coordinates = coordinates,
    response = Z,
    upper.bounds = list(upper, upper),
    intknots = list(NULL, NULL),
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    beta = 0.5
  )
  bivariate <- SplineReg_biv(
    X = X, Y = Y, Z = Z, W = NULL,
    InterKnotsX = NULL, InterKnotsY = NULL,
    Xextr = c(0, 3), Yextr = c(0, 3), n = 2L, fast = TRUE
  )

  expect_equal(generic$design,
               GeDS:::tensorProd(bivariate$Xbasis, bivariate$Ybasis),
               tolerance = 1e-14)
  expect_equal(generic$predicted, bivariate$predicted, tolerance = 1e-12)
  expect_equal(generic$residuals, bivariate$residuals, tolerance = 1e-12)

  matr <- cbind(X = X, Y = Y, residual = bivariate$residuals)
  expected <- list(
    X = GeDS:::placeKnot(
      "X", NULL, matr, NULL, Y, order(Y, X), nint, upper, 0.5
    ),
    Y = GeDS:::placeKnot(
      "Y", NULL, matr, NULL, X, order(X, Y), nint, upper, 0.5
    )
  )
  expect_equal(generic$candidates, expected, tolerance = 1e-12)
})

test_that("one three-dimensional Stage A step selects the best candidate", {
  coordinates <- as.matrix(expand.grid(
    X1 = seq(0.05, 0.95, length.out = 5L),
    X2 = seq(0.05, 0.95, length.out = 5L),
    X3 = seq(0.05, 0.95, length.out = 5L)
  ))
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * X1) + 0.7 * cos(2 * pi * X2) +
      0.4 * sin(4 * pi * X3) + X1 * X2
  )
  upper <- rep(list(c(0.5, 1) + 1e-15), 3L)

  step <- GeDS:::stageAOneStepND(
    coordinates = coordinates,
    response = response,
    upper.bounds = upper,
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    placement.ranges = rep(list(c(0, 1)), 3L),
    beta = 0.5
  )

  expect_identical(step$solver, "tensor-mesh-qr")
  expect_null(step$design)
  expect_equal(step$design.dim, c(125L, 8L))
  expect_length(step$candidates, 3L)
  expect_true(all(!vapply(step$candidates, `[[`, logical(1), "flagDim")))
  expect_equal(step$selected$index, 3L)
  expect_equal(step$selected$dimension, "X3")
  expect_equal(step$selected$newknot,
               step$candidates[[3L]]$Dim.newknot,
               tolerance = 1e-14)
  expect_equal(step$selected$weight,
               max(vapply(step$candidates, `[[`, numeric(1), "weightDim")),
               tolerance = 1e-14)
})

test_that("the early two-dimensional Stage A loop matches BivariateFitter", {
  set.seed(123)
  X <- round(runif(400, min = 0, max = 3), 2)
  Y <- round(runif(400, min = 0, max = 3), 2)
  truth <- sin(2 * X) * sin(2 * Y)
  Z <- truth + rnorm(400, mean = 0, sd = 0.1)
  nint <- as.integer(sqrt(length(Z)))
  upper <- seq(0, 3, length.out = nint + 1L)[-1L] + 1e-15

  loop <- GeDS:::stageALoopND(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    upper.bounds = list(upper, upper),
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    beta = 0.5,
    max.steps = 4L
  )
  fit <- BivariateFitter(
    X = X, Y = Y, Z = Z, W = NULL, indicator = NULL,
    beta = 0.5, phi = 0.99, q = 2L,
    Xextr = c(0, 3), Yextr = c(0, 3),
    show.iters = FALSE
  )

  extract_internal <- function(row) {
    values <- row[!is.na(row)]
    if (length(values) <= 4L) return(NULL)
    values[3L:(length(values) - 2L)]
  }
  fitted.iteration <- loop$completed + 1L
  expected <- list(
    X = extract_internal(fit$stored$previousX[fitted.iteration, ]),
    Y = extract_internal(fit$stored$previousY[fitted.iteration, ])
  )

  expect_equal(loop$completed, 4L)
  expect_equal(vapply(loop$insertions, `[[`, character(1), "dimension"),
               c("X", "X", "Y", "Y"))
  expect_equal(unname(loop$intknots$X), unname(expected$X), tolerance = 1e-12)
  expect_equal(unname(loop$intknots$Y), unname(expected$Y), tolerance = 1e-12)

  sr <- GeDS:::stageALoopND(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    upper.bounds = list(upper, upper),
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    beta = 0.5,
    max.steps = 300L,
    stop.rule = "SR",
    phi = 0.99,
    q = 2L
  )
  expect_equal(length(sr$history), fit$iters)
  expect_equal(sr$rss, fit$rss, tolerance = 1e-10)
  expect_equal(sr$stop.reason, "SR threshold reached")
  expect_equal(sr$selected.iteration, fit$iters - 2L)
  expect_equal(sr$selected.intknots$X, fit$linear.intknots$Xk,
               tolerance = 1e-12)
  expect_equal(sr$selected.intknots$Y, fit$linear.intknots$Yk,
               tolerance = 1e-12)
  expect_gte(tail(sr$phis.star, 1L), 0.99)

  stage_b <- GeDS:::stageBND(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    intknots = sr$selected.intknots,
    coordinate.ranges = list(c(0, 3), c(0, 3))
  )
  expected_knots <- list(
    linear = list(X = fit$linear.intknots$Xk, Y = fit$linear.intknots$Yk),
    quadratic = list(X = fit$quadratic.intknots$Xk,
                     Y = fit$quadratic.intknots$Yk),
    cubic = list(X = fit$cubic.intknots$Xk, Y = fit$cubic.intknots$Yk)
  )
  expected_fits <- list(
    linear = fit$linear.fit,
    quadratic = fit$quadratic.fit,
    cubic = fit$cubic.fit
  )
  expect_equal(stage_b$intknots, expected_knots, tolerance = 1e-12)
  for (fit_name in names(expected_fits)) {
    expect_equal(stage_b$fits[[fit_name]]$rss,
                 expected_fits[[fit_name]]$rss,
                 tolerance = 1e-10)
    expect_equal(stage_b$fits[[fit_name]]$predicted,
                 as.numeric(expected_fits[[fit_name]]$predicted),
                 tolerance = 1e-10)
  }

  multivariate <- GeDS:::MultivariateFitter(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    beta = 0.5,
    phi = 0.99,
    q = 2L,
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    nint = 20L
  )
  expect_s3_class(multivariate, "GeDSfitND")
  expect_equal(multivariate$iters, fit$iters)
  expect_equal(multivariate$selected.iteration, fit$iters - 2L)
  expect_equal(multivariate$Nintknots,
               c(X = length(fit$linear.intknots$Xk),
                 Y = length(fit$linear.intknots$Yk)))
  expect_equal(
    c(multivariate$dev.linear,
      multivariate$dev.quadratic,
      multivariate$dev.cubic),
    c(fit$dev.linear, fit$dev.quadratic, fit$dev.cubic),
    tolerance = 1e-10
  )
  expect_equal(multivariate$cubic.fit$predicted,
               as.numeric(fit$cubic.fit$predicted),
               tolerance = 1e-10)
  expected_predictions <- list(
    `2` = fit$linear.fit$predicted,
    `3` = fit$quadratic.fit$predicted,
    `4` = fit$cubic.fit$predicted
  )
  for (spline_order in 2:4) {
    fit_name <- c("linear", "quadratic", "cubic")[spline_order - 1L]
    expect_equal(
      predict(multivariate, n = spline_order),
      as.numeric(expected_predictions[[as.character(spline_order)]]),
      tolerance = 1e-10
    )
    expect_equal(
      knots(multivariate, n = spline_order, options = "internal"),
      stage_b$fits[[fit_name]]$intknots,
      tolerance = 1e-12
    )
    expect_equal(
      knots(multivariate, n = spline_order),
      stage_b$fits[[fit_name]]$full.knots,
      tolerance = 1e-12
    )
    expect_equal(
      deviance(multivariate, n = spline_order),
      stage_b$fits[[fit_name]]$rss,
      tolerance = 1e-10
    )
    coefficients <- coef(multivariate, n = spline_order)
    basis_index <- attr(coefficients, "basis.index")
    expect_equal(
      as.numeric(coefficients),
      stage_b$fits[[fit_name]]$coefficients,
      tolerance = 1e-12
    )
    expect_named(basis_index, c("X", "Y"))
    expect_equal(NROW(basis_index), length(coefficients))
    expect_equal(
      as.numeric(stage_b$fits[[fit_name]]$design %*% coefficients),
      as.numeric(expected_predictions[[as.character(spline_order)]]),
      tolerance = 1e-10
    )
  }
  reordered <- data.frame(Y = Y[1:10], X = X[1:10])
  expect_equal(
    predict(multivariate, newdata = reordered, n = 4L),
    as.numeric(fit$cubic.fit$predicted)[1:10],
    tolerance = 1e-10
  )

  fit_rd <- BivariateFitter(
    X = X, Y = Y, Z = Z, W = NULL, indicator = NULL,
    beta = 0.5, phi = 0.9, q = 2L,
    Xextr = c(0, 3), Yextr = c(0, 3),
    show.iters = FALSE, stoptype = "RD"
  )
  multivariate_rd <- GeDS:::MultivariateFitter(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    beta = 0.5,
    phi = 0.9,
    q = 2L,
    stoptype = "RD",
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    nint = 20L
  )
  expect_equal(multivariate_rd$stageA$stop.reason, "RD threshold reached")
  expect_equal(multivariate_rd$Nintknots,
               c(X = length(fit_rd$linear.intknots$Xk),
                 Y = length(fit_rd$linear.intknots$Yk)))
  expect_equal(
    predict(multivariate_rd, n = 3L),
    as.numeric(fit_rd$quadratic.fit$predicted),
    tolerance = 1e-10
  )
})

test_that("weighted two-dimensional Stage A and B match BivariateFitter", {
  set.seed(123)
  X <- round(runif(400, min = 0, max = 3), 2)
  Y <- round(runif(400, min = 0, max = 3), 2)
  Z <- sin(2 * X) * sin(2 * Y) + rnorm(400, sd = 0.1)
  weights <- seq(0.25, 2.25, length.out = length(Z))
  upper <- seq(0, 3, length.out = 21L)[-1L] + 1e-15

  generic <- GeDS:::stageALoopND(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    upper.bounds = list(upper, upper),
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    weights = weights,
    beta = 0.5,
    max.steps = 300L,
    stop.rule = "SR",
    phi = 0.99,
    q = 2L
  )
  legacy <- BivariateFitter(
    X = X, Y = Y, Z = Z, W = NULL, weights = weights,
    indicator = NULL, beta = 0.5, phi = 0.99, q = 2L,
    Xextr = c(0, 3), Yextr = c(0, 3),
    show.iters = FALSE, stoptype = "SR"
  )

  extract_internal <- function(row) {
    values <- row[!is.na(row)]
    if (length(values) <= 4L) return(NULL)
    values[3L:(length(values) - 2L)]
  }

  expect_equal(length(generic$history), legacy$iters)
  expect_equal(generic$rss, legacy$rss, tolerance = 1e-10)
  expect_equal(generic$selected.iteration, legacy$iters - 2L)
  expect_equal(generic$selected.intknots$X,
               legacy$linear.intknots$Xk,
               tolerance = 1e-12)
  expect_equal(generic$selected.intknots$Y,
               legacy$linear.intknots$Yk,
               tolerance = 1e-12)

  for (iteration in seq_along(generic$history)) {
    direct <- SplineReg_biv(
      X = X, Y = Y, Z = Z, W = NULL, weights = weights,
      InterKnotsX = extract_internal(legacy$stored$previousX[iteration, ]),
      InterKnotsY = extract_internal(legacy$stored$previousY[iteration, ]),
      Xextr = c(0, 3), Yextr = c(0, 3), n = 2L,
      fast = TRUE
    )
    expect_equal(generic$history[[iteration]]$coefficients,
                 as.numeric(direct$theta),
                 tolerance = 1e-10)
    expect_equal(generic$history[[iteration]]$predicted,
                 as.numeric(direct$predicted),
                 tolerance = 1e-10)
    expect_equal(generic$history[[iteration]]$rss,
                 sum(weights * generic$history[[iteration]]$residuals^2),
                 tolerance = 1e-12)
  }

  generic_stage_b <- GeDS:::stageBND(
    coordinates = cbind(X = X, Y = Y),
    response = Z,
    intknots = generic$selected.intknots,
    coordinate.ranges = list(c(0, 3), c(0, 3)),
    weights = weights
  )
  legacy_stage_b <- list(
    linear = legacy$linear.fit,
    quadratic = legacy$quadratic.fit,
    cubic = legacy$cubic.fit
  )

  for (fit_name in names(legacy_stage_b)) {
    expect_equal(generic_stage_b$fits[[fit_name]]$coefficients,
                 as.numeric(legacy_stage_b[[fit_name]]$theta),
                 tolerance = 1e-10)
    expect_equal(generic_stage_b$fits[[fit_name]]$predicted,
                 as.numeric(legacy_stage_b[[fit_name]]$predicted),
                 tolerance = 1e-10)
    expect_equal(generic_stage_b$fits[[fit_name]]$rss,
                 legacy_stage_b[[fit_name]]$rss,
                 tolerance = 1e-10)
  }
})

test_that("the three-dimensional Stage A loop inserts successive knots", {
  coordinates <- as.matrix(expand.grid(
    X1 = seq(0.05, 0.95, length.out = 5L),
    X2 = seq(0.05, 0.95, length.out = 5L),
    X3 = seq(0.05, 0.95, length.out = 5L)
  ))
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * X1) + 0.7 * cos(2 * pi * X2) +
      0.4 * sin(4 * pi * X3) + X1 * X2
  )

  loop <- GeDS:::stageALoopND(
    coordinates = coordinates,
    response = response,
    upper.bounds = rep(list(c(0.5, 1) + 1e-15), 3L),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    placement.ranges = rep(list(c(0, 1)), 3L),
    beta = 0.5,
    max.steps = 5L
  )

  expect_equal(loop$completed, 5L)
  expect_equal(loop$stop.reason, "maximum steps reached")
  expect_equal(vapply(loop$insertions, `[[`, integer(1), "index"),
               c(3L, 2L, 1L, 1L, 3L))
  expect_equal(lengths(loop$intknots), c(X1 = 2L, X2 = 1L, X3 = 2L))
  expect_true(all(diff(loop$rss) <= 1e-10))

  sr <- GeDS:::stageALoopND(
    coordinates = coordinates,
    response = response,
    upper.bounds = rep(list(c(0.5, 1) + 1e-15), 3L),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    placement.ranges = rep(list(c(0, 1)), 3L),
    beta = 0.5,
    max.steps = 50L,
    stop.rule = "SR",
    phi = 0.98,
    q = 2L
  )
  expect_equal(sr$stop.reason, "SR threshold reached")
  expect_equal(sr$selected.iteration, 8L)
  expect_equal(sr$completed, 9L)
  expect_gte(tail(sr$phis.star, 1L), 0.98)
  expect_equal(sum(lengths(sr$selected.intknots)), 7L)

  stage_b <- GeDS:::stageBND(
    coordinates = coordinates,
    response = response,
    intknots = sr$selected.intknots,
    coordinate.ranges = rep(list(c(0, 1)), 3L)
  )
  expect_equal(lengths(stage_b$intknots$linear),
               c(X1 = 4L, X2 = 1L, X3 = 2L))
  expect_equal(lengths(stage_b$intknots$quadratic),
               c(X1 = 3L, X2 = 0L, X3 = 1L))
  expect_equal(lengths(stage_b$intknots$cubic),
               c(X1 = 2L, X2 = 0L, X3 = 0L))
  expect_equal(
    lapply(stage_b$fits, `[[`, "design.dim"),
    list(linear = c(125L, 72L), quadratic = c(125L, 72L),
         cubic = c(125L, 96L))
  )
  # The selected X1 bases are wider than the five mesh levels, so the
  # rank-checked mesh solver must fall back to dense least squares here.
  expect_identical(
    vapply(stage_b$fits, `[[`, character(1), "solver"),
    c(linear = "dense", quadratic = "dense", cubic = "dense")
  )
  expect_equal(
    vapply(stage_b$fits, `[[`, numeric(1), "rss"),
    c(linear = 0.189164136790582,
      quadratic = 1.42582737244196,
      cubic = 1.41785997282093),
    tolerance = 1e-10
  )

  default_grid <- GeDS:::makeGridBoundsND(coordinates)
  expect_equal(default_grid$nint, c(X1 = 5L, X2 = 5L, X3 = 5L))

  fit <- GeDS:::MultivariateFitter(
    coordinates = coordinates,
    response = response,
    beta = 0.5,
    phi = 0.98,
    q = 2L,
    max.steps = 50L,
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    placement.ranges = rep(list(c(0, 1)), 3L),
    nint = 2L
  )
  expect_s3_class(fit, "GeDSfitND")
  expect_equal(fit$dimensions, c("X1", "X2", "X3"))
  expect_equal(fit$iters, 10L)
  expect_equal(fit$selected.iteration, 8L)
  expect_equal(fit$stageA$stop.reason, "SR threshold reached")
  expect_true(any(vapply(
    fit$stageA$history,
    function(step) identical(step$solver, "tensor-mesh-qr"),
    logical(1)
  )))
  expect_true(any(vapply(
    fit$stageA$history,
    function(step) identical(step$solver, "dense"),
    logical(1)
  )))
  expect_equal(fit$Nintknots, c(X1 = 4L, X2 = 1L, X3 = 2L))
  expect_equal(
    c(fit$dev.linear, fit$dev.quadratic, fit$dev.cubic),
    c(0.189164136790582, 1.42582737244196, 1.41785997282093),
    tolerance = 1e-10
  )
  print_output <- capture.output(print_result <- print(fit))
  expect_identical(print_result, fit)
  expect_true(any(grepl("Experimental multivariate Normal GeDS fit", print_output)))
  expect_true(any(grepl("Dimensions \\(3\\): X1, X2, X3", print_output)))
  expect_true(any(grepl("X1=4, X2=1, X3=2", print_output)))
  expect_true(any(grepl("Least-squares solver", print_output)))
  summary_output <- capture.output(summary_result <- summary(fit))
  expect_identical(summary_result, fit)
  expect_true(any(grepl("Stopping rule: SR", summary_output)))
  expect_true(any(grepl("Termination: SR threshold reached", summary_output)))
  for (spline_order in 2:4) {
    fit_name <- c("linear", "quadratic", "cubic")[spline_order - 1L]
    expect_equal(
      predict(fit, n = spline_order),
      fit$stageB$fits[[fit_name]]$predicted,
      tolerance = 1e-10
    )
    coefficients <- coef(fit, n = spline_order)
    basis_index <- attr(coefficients, "basis.index")
    basis_counts <- vapply(
      fit$stageB$fits[[fit_name]]$basis.matrices,
      NCOL,
      integer(1)
    )
    expect_named(basis_index, c("X1", "X2", "X3"))
    expect_equal(as.integer(basis_index[1L, ]), c(1L, 1L, 1L))
    expect_equal(
      as.integer(basis_index[NROW(basis_index), ]),
      unname(basis_counts)
    )
    expect_equal(
      as.numeric(GeDS:::tensorProdND(
        fit$stageB$fits[[fit_name]]$basis.matrices
      ) %*% coefficients),
      fit$stageB$fits[[fit_name]]$predicted,
      tolerance = 1e-10
    )
  }
  newdata <- as.data.frame(coordinates[1:8, c("X3", "X1", "X2")])
  expect_equal(
    predict(fit, newdata = newdata, n = 3L),
    fit$quadratic.fit$predicted[1:8],
    tolerance = 1e-10
  )
  expect_error(
    predict(fit, newdata = newdata[, c("X1", "X2")], n = 3L),
    "missing fitted dimensions"
  )
})

test_that("tensor-mesh QR agrees with dense least squares", {
  coordinates <- as.matrix(expand.grid(
    x0 = seq(0, 1, length.out = 4L),
    x1 = seq(0, 1, length.out = 5L),
    x2 = seq(0, 1, length.out = 3L)
  ))
  set.seed(912)
  coordinates <- coordinates[sample(NROW(coordinates)), , drop = FALSE]
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * x0) + x1 * x2 + 0.2 * x0 * x1 * x2
  )
  basis_matrices <- lapply(seq_len(NCOL(coordinates)), function(j) {
    splines::splineDesign(
      knots = rep(range(coordinates[, j]), 2L),
      x = coordinates[, j], ord = 2L, outer.ok = TRUE
    )
  })

  mesh_fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates, response, basis_matrices, rep(1, NROW(coordinates))
  )
  dense_design <- GeDS:::tensorProdND(basis_matrices)
  dense_coefficients <- as.numeric(stats::coef(.lm.fit(dense_design, response)))
  dense_coefficients[is.na(dense_coefficients)] <- 0

  expect_identical(mesh_fit$solver, "tensor-mesh-qr")
  expect_null(mesh_fit$design)
  expect_equal(mesh_fit$design.dim, dim(dense_design))
  expect_equal(mesh_fit$rank, NCOL(dense_design))
  expect_false(mesh_fit$rank.deficient)
  expect_equal(mesh_fit$effective.nobs, NROW(coordinates))
  expect_equal(mesh_fit$coefficients, dense_coefficients, tolerance = 1e-11)
  expect_equal(
    mesh_fit$predicted,
    as.numeric(dense_design %*% dense_coefficients),
    tolerance = 1e-11
  )

  nonuniform_weights <- seq(1, 2, length.out = NROW(coordinates))
  nonuniform_weights[c(2L, 19L, 41L)] <- 0
  weighted_fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates, response, basis_matrices, nonuniform_weights
  )
  weighted_reference <- stats::lm.wfit(
    dense_design, response, nonuniform_weights
  )
  weighted_coefficients <- as.numeric(stats::coef(weighted_reference))
  weighted_coefficients[is.na(weighted_coefficients)] <- 0
  expect_identical(weighted_fit$solver, "dense")
  expect_equal(weighted_fit$design, dense_design, tolerance = 0)
  expect_equal(weighted_fit$effective.nobs, sum(nonuniform_weights > 0))
  expect_equal(weighted_fit$rank, weighted_reference$rank)
  expect_identical(
    weighted_fit$rank.deficient,
    weighted_reference$rank < NCOL(dense_design)
  )
  expect_equal(
    weighted_fit$coefficients, weighted_coefficients, tolerance = 1e-11
  )
  expect_equal(
    weighted_fit$predicted,
    as.numeric(dense_design %*% weighted_coefficients),
    tolerance = 1e-11
  )

  incomplete_design <- dense_design[-1L, , drop = FALSE]
  incomplete_fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates[-1, ], response[-1],
    lapply(basis_matrices, function(basis) basis[-1, , drop = FALSE]),
    rep(1, NROW(coordinates) - 1L)
  )
  incomplete_reference <- stats::lm.fit(incomplete_design, response[-1L])
  incomplete_coefficients <- as.numeric(stats::coef(incomplete_reference))
  incomplete_coefficients[is.na(incomplete_coefficients)] <- 0
  expect_identical(incomplete_fit$solver, "dense")
  expect_false(incomplete_fit$mesh$complete)
  expect_equal(incomplete_fit$rank, incomplete_reference$rank)
  expect_equal(
    incomplete_fit$predicted,
    as.numeric(incomplete_design %*% incomplete_coefficients),
    tolerance = 1e-11
  )

  aliased_bases <- basis_matrices
  aliased_bases[[1L]] <- cbind(
    aliased_bases[[1L]], aliased_bases[[1L]][, 1L]
  )
  aliased_design <- GeDS:::tensorProdND(aliased_bases)
  aliased_fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates, response, aliased_bases, rep(1, NROW(coordinates))
  )
  aliased_reference <- stats::lm.fit(aliased_design, response)
  aliased_coefficients <- as.numeric(stats::coef(aliased_reference))
  aliased_coefficients[is.na(aliased_coefficients)] <- 0
  expect_identical(aliased_fit$solver, "dense")
  expect_true(aliased_fit$mesh$complete)
  expect_true(aliased_fit$rank.deficient)
  expect_equal(aliased_fit$rank, aliased_reference$rank)
  expect_true(all(is.finite(aliased_fit$coefficients)))
  expect_equal(
    aliased_fit$predicted,
    as.numeric(aliased_design %*% aliased_coefficients),
    tolerance = 1e-11
  )

  expect_error(
    GeDS:::fitTensorLeastSquaresND(
      coordinates, response, basis_matrices, rep(0, NROW(coordinates))
    ),
    "at least one positive value"
  )
})

test_that("zero-weight rows do not affect multivariate knot placement", {
  coordinates <- as.matrix(expand.grid(
    X1 = seq(0.05, 0.95, length.out = 4L),
    X2 = seq(0.05, 0.95, length.out = 4L),
    X3 = seq(0.05, 0.95, length.out = 4L)
  ))
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * X1) + cos(2 * pi * X2) + X1 * X2 * X3
  )
  weights <- seq(0.5, 1.5, length.out = NROW(coordinates))
  weights[c(6L, 19L, 35L, 52L)] <- 0
  keep <- weights > 0
  upper.bounds <- rep(list(c(0.25, 0.5, 0.75, 1) + 1e-15), 3L)
  coordinate.ranges <- rep(list(c(0, 1)), 3L)

  weighted <- GeDS:::stageAOneStepND(
    coordinates = coordinates,
    response = response,
    upper.bounds = upper.bounds,
    coordinate.ranges = coordinate.ranges,
    placement.ranges = coordinate.ranges,
    weights = weights
  )
  dropped <- GeDS:::stageAOneStepND(
    coordinates = coordinates[keep, , drop = FALSE],
    response = response[keep],
    upper.bounds = upper.bounds,
    coordinate.ranges = coordinate.ranges,
    placement.ranges = coordinate.ranges,
    weights = weights[keep]
  )

  expect_equal(weighted$effective.nobs, sum(keep))
  expect_equal(weighted$coefficients, dropped$coefficients, tolerance = 1e-11)
  expect_equal(weighted$predicted[keep], dropped$predicted, tolerance = 1e-11)
  expect_equal(weighted$candidates, dropped$candidates, tolerance = 1e-11)
  expect_equal(weighted$selected, dropped$selected, tolerance = 1e-11)

  no_signal <- GeDS:::placeDimKnot(
    dim.index = 1L,
    intknots = NULL,
    coordinates = coordinates,
    residuals = rep(0, NROW(coordinates)),
    strip.id = rep(1L, NROW(coordinates)),
    beta = 0.5,
    dim.range = c(0, 1)
  )
  expect_true(no_signal$flagDim)
  expect_true(is.na(no_signal$Dim.newknot))
})

test_that("multivariate NGeDS rejects all-zero weights", {
  dat <- expand.grid(
    x0 = seq(0.05, 0.95, length.out = 3L),
    x1 = seq(0.05, 0.95, length.out = 3L),
    x2 = seq(0.05, 0.95, length.out = 3L)
  )
  dat$y <- with(dat, sin(2 * pi * x0) + x1 * x2)

  expect_error(
    NGeDS(
      y ~ f(x0, x1, x2), data = dat,
      weights = rep(0, NROW(dat)), max.intknots = 2L,
      stoptype = "RD", higher_order = FALSE
    ),
    "at least one positive value"
  )
})

test_that("a four-dimensional complete mesh uses separable least squares", {
  coordinates <- as.matrix(expand.grid(
    x0 = seq(0, 1, length.out = 4L),
    x1 = seq(0, 1, length.out = 4L),
    x2 = seq(0, 1, length.out = 4L),
    x3 = seq(0, 1, length.out = 4L)
  ))
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * x0) + x1 * x2 - x3 + x0 * x1 * x2 * x3
  )
  fit <- GeDS:::fitTensorSplineND(
    coordinates = coordinates,
    response = response,
    spline.order = 2L
  )

  expect_identical(fit$solver, "tensor-mesh-qr")
  expect_true(fit$mesh$complete)
  expect_null(fit$design)
  expect_equal(fit$design.dim, c(256L, 16L))
  expect_equal(fit$rank, 16L)
  expect_false(fit$rank.deficient)
  expect_equal(fit$effective.nobs, 256L)
  expect_equal(length(fit$coefficients), 16L)
  expect_equal(length(fit$predicted), 256L)
})

test_that("the multivariate fitter agrees with a bivariate NGeDSgam learner", {
  dat <- stats::na.omit(datasets::airquality)
  dat$Ozone <- dat$Ozone^(1 / 3)

  gam <- suppressMessages(suppressWarnings(
    NGeDSgam(
      Ozone ~ f(Wind, Temp),
      data = dat,
      phi = 0.8
    )
  ))
  multivariate <- GeDS:::MultivariateFitter(
    coordinates = as.matrix(dat[c("Wind", "Temp")]),
    response = dat$Ozone,
    phi = 0.8
  )

  gam_predictions <- unname(gam$predictions)
  multivariate_predictions <- lapply(
    2:4,
    function(spline_order) predict(multivariate, n = spline_order)
  )
  for (i in seq_along(gam_predictions)) {
    expect_equal(
      multivariate_predictions[[i]],
      gam_predictions[[i]],
      tolerance = 1e-10
    )
  }

  mse <- vapply(
    multivariate_predictions,
    function(prediction) mean((dat$Ozone - prediction)^2),
    numeric(1)
  )
  expect_equal(
    mse,
    c(0.268985116963, 0.204605651109, 0.198708755782),
    tolerance = 1e-10
  )
})

test_that("a joint fit captures a pure three-way interaction", {
  set.seed(123)
  n <- 800L
  # Each one-dimensional marginal is zero under independent Uniform(0, 1)
  # coordinates, so the signal is entirely a three-way interaction.
  truth_function <- function(x0, x1, x2) {
    sin(2 * pi * x0) * sin(2 * pi * x1) * sin(2 * pi * x2)
  }
  x0 <- runif(n)
  x1 <- runif(n)
  x2 <- runif(n)
  truth <- truth_function(x0, x1, x2)
  y <- rnorm(n, mean = truth, sd = 0.1)
  dat <- data.frame(y, x0, x1, x2)

  additive <- suppressMessages(suppressWarnings(
    NGeDSgam(y ~ f(x0) + f(x1) + f(x2), data = dat)
  ))
  joint <- NGeDS(y ~ f(x0, x1, x2), data = dat)
  expect_s3_class(joint, "GeDSfitND")
  expect_equal(joint$args$stoptype, "SR")

  set.seed(456)
  n_test <- 5000L
  newdata <- data.frame(
    x0 = runif(n_test, min(x0), max(x0)),
    x1 = runif(n_test, min(x1), max(x1)),
    x2 = runif(n_test, min(x2), max(x2))
  )
  test_truth <- with(newdata, truth_function(x0, x1, x2))
  additive_predictions <- lapply(
    2:4,
    function(spline_order) predict(additive, newdata, n = spline_order)
  )
  joint_predictions <- lapply(
    2:4,
    function(spline_order) predict(joint, newdata, n = spline_order)
  )
  truth_mse <- function(predictions) {
    vapply(
      predictions,
      function(prediction) mean((test_truth - prediction)^2),
      numeric(1)
    )
  }

  expect_equal(
    truth_mse(additive_predictions),
    c(0.126950594354684, 0.127284608146684, 0.128045883170027),
    tolerance = 1e-8
  )
  expect_equal(
    truth_mse(joint_predictions),
    c(0.0768980702441790, 0.0852003127970236, 0.00419966281976198),
    tolerance = 1e-8
  )
  expect_equal(
    lengths(additive$internal_knots$linear.int.knots),
    c(`f(x0)` = 0L, `f(x1)` = 0L, `f(x2)` = 0L)
  )
  expect_equal(joint$Nintknots, c(x0 = 2L, x1 = 1L, x2 = 1L))
})

test_that("public generalized and additive fitters support three-dimensional smoothers", {
  set.seed(204)
  dat <- data.frame(x1 = runif(70), x2 = runif(70), x3 = runif(70))
  eta <- with(dat, -0.4 + x1 - x2 + 0.8 * x1 * x2 * x3)
  dat$count <- rpois(nrow(dat), exp(eta))
  dat$y <- eta + rnorm(nrow(dat), sd = 0.08)

  generalized <- GGeDS(
    count ~ f(x1, x2, x3), data = dat, family = poisson(),
    max.intknots = 3L, higher_order = FALSE
  )
  expect_s3_class(generalized, "GeDSfitND")
  expect_true(all(is.finite(predict(generalized, dat[1:4, ], n = 2L))))

  additive <- NGeDSgam(
    y ~ f(x1, x2, x3), data = dat, min_iterations = 1L,
    max_iterations = 2L, internal_knots = 3L,
    higher_order = FALSE
  )
  expect_s3_class(additive, "GeDSgam")
  expect_true(all(is.finite(predict(additive, dat[1:4, ], n = 2L))))

  boosted <- NGeDSboost(
    y ~ f(x1, x2, x3), data = dat, min_iterations = 1L,
    max_iterations = 2L, int.knots_init = 2L, int.knots_boost = 1L,
    higher_order = FALSE
  )
  expect_s3_class(boosted, "GeDSboost")
  expect_true(all(is.finite(predict(boosted, dat[1:4, ], n = 2L))))
})

test_that("four-dimensional GGeDS supports weights, parametric terms, and offsets", {
  set.seed(205)
  dat <- as.data.frame(replicate(4, runif(90)))
  names(dat) <- paste0("x", 1:4)
  dat$z <- rnorm(nrow(dat))
  dat$exposure <- runif(nrow(dat), 0.7, 1.4)
  eta <- with(dat, -0.5 + x1 * x2 - x3 * x4 + 0.3 * z + log(exposure))
  dat$y <- rpois(nrow(dat), exp(eta))
  weights <- rep(1, nrow(dat)); weights[1:5] <- 0

  fit <- GGeDS(
    y ~ f(x1, x2, x3, x4) + z + offset(log(exposure)),
    data = dat, family = poisson(), weights = weights,
    max.intknots = 2L, higher_order = FALSE
  )
  expect_s3_class(fit, "GeDSfitND")
  expect_length(fit$dimensions, 4L)
  expect_true(all(is.finite(predict(fit, dat[1:6, ], n = 2L))))

  changed <- dat
  changed$y[weights == 0] <- changed$y[weights == 0] + 1000
  refit <- GGeDS(
    y ~ f(x1, x2, x3, x4) + z + offset(log(exposure)),
    data = changed, family = poisson(), weights = weights,
    max.intknots = 2L, higher_order = FALSE
  )
  expect_equal(fit$linear.intknots, refit$linear.intknots)
})

test_that("multivariate boosting memory retains supplied coordinate knots", {
  set.seed(206)
  coordinates <- matrix(runif(240), ncol = 3)
  colnames(coordinates) <- c("x1", "x2", "x3")
  initial <- list(x1 = 0.25, x2 = 0.5, x3 = NULL)
  fit <- GeDS:::MultivariateFitter(
    coordinates, rowSums(coordinates), intknots_init = initial,
    max.steps = 2L, stoptype = "RD", spline.orders = 2L, max.coef = 500L
  )
  expect_true(all(vapply(seq_along(initial), function(j) {
    all(initial[[j]] %in% fit$stageA$history[[1L]]$intknots[[j]])
  }, logical(1))))
})

test_that("multivariate prediction does not resolve f in the formula environment", {
  set.seed(207)
  dat <- data.frame(x1 = runif(70), x2 = runif(70), x3 = runif(70))
  dat$y <- with(dat, x1 * x2 * x3 + rnorm(nrow(dat), sd = 0.05))
  isolated <- new.env(parent = baseenv())
  formula <- stats::as.formula("y ~ f(x1, x2, x3)", env = isolated)
  fit <- NGeDS(
    formula, data = dat, max.intknots = 2L,
    higher_order = FALSE
  )
  expect_false(exists("f", envir = isolated, inherits = TRUE))
  expect_true(all(is.finite(predict(fit, dat[1:5, ], n = 2L))))
})

test_that("the Stage A loop stops before fitting a saturated tensor basis", {
  coordinates <- as.matrix(expand.grid(
    X1 = c(0.25, 0.75), X2 = c(0.25, 0.75), X3 = c(0.25, 0.75)
  ))
  loop <- GeDS:::stageALoopND(
    coordinates = coordinates,
    response = rowSums(coordinates),
    upper.bounds = rep(list(c(0.5, 1) + 1e-15), 3L),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    max.steps = 3L
  )

  expect_equal(loop$completed, 0L)
  expect_length(loop$history, 0L)
  expect_equal(loop$stop.reason, "linear tensor basis is saturated")

  coordinates_with_zero_weights <- as.matrix(expand.grid(
    X1 = c(0.1, 0.5, 0.9), X2 = c(0.25, 0.75), X3 = c(0.25, 0.75)
  ))
  weights <- c(rep(1, 7L), rep(0, 5L))
  effective_loop <- GeDS:::stageALoopND(
    coordinates = coordinates_with_zero_weights,
    response = rowSums(coordinates_with_zero_weights),
    upper.bounds = list(c(0.4, 0.7, 1) + 1e-15,
                        c(0.5, 1) + 1e-15,
                        c(0.5, 1) + 1e-15),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    weights = weights,
    max.steps = 3L
  )
  expect_equal(effective_loop$effective.nobs, 7L)
  expect_length(effective_loop$history, 0L)
  expect_equal(effective_loop$stop.reason, "linear tensor basis is saturated")

  expect_error(
    GeDS:::stageALoopND(
      coordinates = coordinates_with_zero_weights,
      response = rowSums(coordinates_with_zero_weights),
      upper.bounds = list(c(0.4, 0.7, 1) + 1e-15,
                          c(0.5, 1) + 1e-15,
                          c(0.5, 1) + 1e-15),
      weights = rep(0, NROW(coordinates_with_zero_weights)),
      max.steps = 3L
    ),
    "at least one positive value"
  )
})

test_that("tensor coefficient limits stop growth and omit oversized orders", {
  coordinates <- as.matrix(expand.grid(
    X1 = seq(0.1, 0.9, length.out = 4L),
    X2 = seq(0.1, 0.9, length.out = 4L),
    X3 = seq(0.1, 0.9, length.out = 4L)
  ))
  response <- with(
    as.data.frame(coordinates),
    sin(2 * pi * X1) + X2 * X3
  )
  loop <- GeDS:::stageALoopND(
    coordinates = coordinates,
    response = response,
    upper.bounds = rep(list(seq(0.25, 1, length.out = 4L) + 1e-15), 3L),
    coordinate.ranges = rep(list(c(0, 1)), 3L),
    max.steps = 5L,
    stop.rule = "none",
    max.coef = 8L
  )
  expect_equal(loop$stop.reason, "maximum coefficient limit reached")
  expect_equal(loop$basis.sizes, 8)
  expect_equal(loop$ncoef, 8L)
  expect_equal(loop$max.coef, 8)

  coordinates4 <- as.matrix(expand.grid(
    x0 = seq(0, 1, length.out = 4L),
    x1 = seq(0, 1, length.out = 4L),
    x2 = seq(0, 1, length.out = 4L),
    x3 = seq(0, 1, length.out = 4L)
  ))
  expect_warning(
    stage_b <- GeDS:::stageBND(
      coordinates = coordinates4,
      response = rowSums(coordinates4),
      intknots = setNames(rep(list(NULL), 4L), colnames(coordinates4)),
      spline.orders = 2:4,
      max.coef = 20L
    ),
    "Omitting Stage B quadratic and cubic fits"
  )
  expect_equal(stage_b$basis.sizes, c(linear = 16, quadratic = 81, cubic = 256))
  expect_type(stage_b$fits$linear, "list")
  expect_null(stage_b$fits$quadratic)
  expect_null(stage_b$fits$cubic)
  expect_true(is.na(stage_b$skipped[["linear"]]))
  expect_match(stage_b$skipped[["quadratic"]], "requires 81 coefficients")
  expect_match(stage_b$skipped[["cubic"]], "requires 256 coefficients")

  expect_error(
    GeDS:::MultivariateFitter(
      coordinates = coordinates4,
      response = rowSums(coordinates4),
      max.coef = 15L
    ),
    "initial linear tensor basis requires 16 coefficients"
  )
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
