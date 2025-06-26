################################################################################
################################################################################
################################# GeDS-package #################################
################################################################################
################################################################################
#' @title Geometrically Designed Spline Regression
#' @name GeDS-package
#' @description
#' Geometrically designed spline (GeDS) regression is a non-parametric method
#' for fitting spline regression models with variable knots. The GeDS technique
#' is inspired by geometric principles and falls within the domain of
#' generalized non-linear models (GNM), which include generalized linear models
#' (GLM) as a special case. GeDS regression is fitted  based on a sample of
#' \eqn{N} observations of a response variable \eqn{y}, dependent on a set of
#' (currently up to two) covariates, assuming \eqn{y} has a distribution from
#' the exponential family. In addition, GeDS methodology is implemented both in
#' the context of generalized additive models (GAM) and functional gradient
#' boosting (FGB). On the one hand, GAM consist of an additive modeling
#' technique where the impact of the predictor variables is captured through
#' smooth (GeDS, in this case) functions. On the other hand, GeDS incorporates
#' gradient boosting machine learning technique by implementing functional
#' gradient descent algorithm to optimize general risk functions utilizing
#' component-wise GeDS estimates.
#'
#' @details
#' GeDS provides a novel solution to the spline regression problem and, in
#' particular, to the problem of estimating the number and position of the knots.
#' The GeDS estimation method is based on two stages: first, in stage A, a piecewise
#' linear fit (spline fit of order 2) capturing the underlying functional shape
#' determined by the data is constructed; second, in stage B, the latter fit is
#' approximated through shape preserving (variation diminishing) spline fits of
#' higher orders (\eqn{n = 3}, \eqn{n = 4},\eqn{\dots}, i.e., degrees 2, 3,\eqn{\dots}).
#' As a result, GeDS simultaneously produces a linear, a quadratic and
#' a cubic spline fit.
#'
#' The GeDS method was originally developed by Kaishev et al. (2016) for the univariate
#' case, assuming the response variable \eqn{y} to be normally distributed and a
#' corresponding \emph{Mathematica} code was provided.
#'
#' The GeDS method was extended by Dimitrova et al. (2023) to cover any
#' distribution from the exponential family. The \pkg{GeDS} \R package presented
#' here provides an enhanced implementation of the original normal GeDS
#' \emph{Mathematica} code, through the \code{\link{NGeDS}} function; it also
#' includes a generalization, \code{\link{GGeDS}}, which extends the method to
#' any distribution in the exponential family.
#'
#' The \pkg{GeDS} package allows also to fit two dimensional response surfaces
#' and to construct multivariate (predictor) models with a GeD spline component
#' and a parametric component (see the functions \code{\link{f}},
#' \code{\link[=formula.GeDS]{formula}}, \code{\link{NGeDS}} and
#' \code{\link{GGeDS}} for details).
#' 
#' Dimitrova et al. (2025) have recently made significant enhancements to the
#' \pkg{GeDS} methodology, by incorporating generalized additive models (GAM-GeDS)
#' and functional gradient boosting (FGB-GeDS). On the one hand, generalized additive
#' models are encompassed by implementing the \emph{local-scoring} algorithm
#' using normal GeD splines (i.e., \code{\link{NGeDS}}) as function smoothers
#' within the \emph{backfitting} iterations. This is implemented through the function
#' \code{\link{NGeDSgam}}. On the other hand, the \pkg{GeDS} package incorporates
#' functional gradient descent algorithm by utilizing normal GeD splines (i.e.,
#' \code{\link{NGeDS}}) as base learners within the boosting iterations. Unlike
#' typical boosting methods, the final FGB-GeDS model is expressed as a single
#' spline model rather than as a sum of base-learner fits. For this,
#' \code{\link{NGeDSboost}} leverages the piecewise polynomial representation of
#' B-splines, and, at each boosting iteration, performs a piecewise update of the
#' corresponding polynomial coefficients.
#'
#' The outputs of both \code{\link{NGeDS}} and \code{\link{GGeDS}} functions are
#' \code{"GeDS"} class objects, while the outputs of \code{\link{NGeDSgam}}
#' and \code{\link{NGeDSboost}} functions are\code{"GeDSgam"} class and
#' \code{"GeDSboost"} class objects, respectively. \code{"GeDS"} class,
#' \code{"GeDSgam"} class and \code{"GeDSboost"} class objects contain
#' second, third and fourth order spline fits. As described in 
#' Kaishev et al. (2016), Dimitrova et al. (2023) and  Dimitrova et al. (2025),
#' the "final" GeDS fit is the one minimizing the empirical deviance. Nevertheless,
#' the user can choose to use any of the available fits.
#'
#' The \pkg{GeDS} package also includes some datasets where GeDS regression
#' proves to be very efficient and some user friendly functions that are designed
#' to easily extract required information. Several methods are also provided to
#' handle GeDS, GAM-GeDS and FGB-GeDS output results (see \code{\link{NGeDS}}/\code{\link{GGeDS}},
#' \code{\link{NGeDSgam}} and \code{\link{NGeDSboost}}, respectively).
#'
#' Throughout this document, we use the terms GeDS predictor model, GeDS
#' regression and GeDS fit interchangeably.
#'
#' Please report any issue arising or bug in the code to
#' \email{Emilio.Saenz-Guillen@citystgeorges.ac.uk}.
#'
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S., & Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \doi{10.1007/s00180-015-0621-7}
#' 
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}
#' 
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. L. (2025).
#' \pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
#' Models and Functional Gradient Boosting, based on Geometrically Designed
#' (GeD) Splines. \emph{Manuscript submitted for publication.}
#'
#' @keywords package
#' @aliases GeDS
#'
#' @author  Dimitrina S. Dimitrova <D.Dimitrova@citystgeorges.ac.uk>,\cr 
#' Vladimir K. Kaishev <V.Kaishev@citystgeorges.ac.uk>,\cr
#' Andrea Lattuada <Andrea.Lattuada@hotmail.com>,\cr
#' Emilio L. Sáenz Guillén <Emilio.Saenz-Guillen@citystgeorges.ac.uk> and\cr
#' Richard J. Verrall <R.J.Verrall@citystgeorges.ac.uk>
#'
#' 
"_PACKAGE"
NULL

