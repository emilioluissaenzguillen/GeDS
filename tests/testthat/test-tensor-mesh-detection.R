test_that("scattered high-dimensional data do not construct overflowing grid indices", {
  set.seed(915)
  coordinates <- matrix(runif(500), ncol = 5)
  colnames(coordinates) <- paste0("x", 1:5)
  expect_warning(mesh <- GeDS:::detectTensorMeshND(coordinates), NA)
  expect_false(mesh$complete)
  expect_null(mesh$cell.id)
  expect_equal(unname(mesh$axis.lengths), rep(100L, 5))
  expect_identical(colnames(mesh$axis.index), colnames(coordinates))

  response <- rowSums(coordinates^2)
  bases <- lapply(1:5, function(j) {
    splines::splineDesign(rep(range(coordinates[, j]), 2),
                         coordinates[, j], ord = 2, outer.ok = TRUE)
  })
  expect_warning(fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates, response, bases, rep(1, nrow(coordinates))), NA)
  expect_identical(fit$solver, "dense")
  reference <- .lm.fit(GeDS:::tensorProdND(bases), response)
  expect_equal(fit$coefficients, unname(reference$coefficients), tolerance = 1e-11)
  expect_equal(fit$predicted,
               as.numeric(GeDS:::tensorProdND(bases) %*% reference$coefficients),
               tolerance = 1e-11)
})

test_that("shuffled five-dimensional grids retain the fast separable solver", {
  coordinates <- as.matrix(expand.grid(rep(list(c(0, .4, 1)), 5)))
  set.seed(916)
  coordinates <- coordinates[sample(nrow(coordinates)), , drop = FALSE]
  response <- rowSums(coordinates^2) + apply(coordinates, 1, prod)
  bases <- lapply(1:5, function(j) {
    splines::splineDesign(rep(c(0, 1), 2), coordinates[, j], ord = 2)
  })
  expect_warning(fit <- GeDS:::fitTensorLeastSquaresND(
    coordinates, response, bases, rep(1, nrow(coordinates))), NA)
  expect_true(fit$mesh$complete)
  expect_identical(sort(fit$mesh$cell.id), seq_len(nrow(coordinates)))
  expect_identical(fit$solver, "tensor-mesh-qr")
  expect_null(fit$design)
  reference <- .lm.fit(GeDS:::tensorProdND(bases), response)
  expect_equal(fit$coefficients, unname(reference$coefficients), tolerance = 1e-11)
  expect_equal(fit$predicted,
               as.numeric(GeDS:::tensorProdND(bases) %*% reference$coefficients),
               tolerance = 1e-11)
})

test_that("matching grid size alone does not accept duplicated or missing cells", {
  coordinates <- as.matrix(expand.grid(x = 0:1, y = 0:1))
  coordinates[4, ] <- coordinates[1, ]
  expect_warning(mesh <- GeDS:::detectTensorMeshND(coordinates), NA)
  expect_false(mesh$complete)
  expect_true(anyDuplicated(mesh$cell.id) > 0)
  expect_warning(incomplete <- GeDS:::detectTensorMeshND(coordinates[1:3, ]), NA)
  expect_false(incomplete$complete)
  expect_null(incomplete$cell.id)
})
