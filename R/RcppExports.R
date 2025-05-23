# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

whmx <- function(vector) {
    .Call('_GeDS_whmx', PACKAGE = 'GeDS', vector)
}

Knotnew <- function(weights, residuals, x, dcum, oldknots, tol) {
    .Call('_GeDS_Knotnew', PACKAGE = 'GeDS', weights, residuals, x, dcum, oldknots, tol)
}

makenewknots <- function(knots, degree) {
    .Call('_GeDS_makenewknots', PACKAGE = 'GeDS', knots, degree)
}

makeEpsilonsb <- function(data, Xs, Ys, degree) {
    .Call('_GeDS_makeEpsilonsb', PACKAGE = 'GeDS', data, Xs, Ys, degree)
}

makeRatSplines <- function(matrice, h) {
    .Call('_GeDS_makeRatSplines', PACKAGE = 'GeDS', matrice, h)
}

makeWeights <- function(x) {
    .Call('_GeDS_makeWeights', PACKAGE = 'GeDS', x)
}

tensorProd <- function(Xmat, Ymat) {
    .Call('_GeDS_tensorProd', PACKAGE = 'GeDS', Xmat, Ymat)
}

makeNewMatrCPP <- function(matrix, tab = NULL, by_row = FALSE) {
    .Call('_GeDS_makeNewMatrCPP', PACKAGE = 'GeDS', matrix, tab, by_row)
}

findNewDimKnot <- function(dcumFixedDim_Dim, Dim_weights, Dim_oldknots, matrFixedDim, Dim_index) {
    .Call('_GeDS_findNewDimKnot', PACKAGE = 'GeDS', dcumFixedDim_Dim, Dim_weights, Dim_oldknots, matrFixedDim, Dim_index)
}

