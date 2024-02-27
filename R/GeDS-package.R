################################################################################
################################################################################
################################# GeDS-package #################################
################################################################################
################################################################################
#' @title GeDS
#' @name GeDS-package
#' @description
#' Geometrically Designed Splines (GeDS) regression is a non-parametric method
#' inspired by geometric principles, which is designed for fitting spline 
#' predictor models with variable knots. This method falls within the domain of
#' generalized non-linear models (GNM), which include generalized linear models
#' (GLM) as a special case. GeDS regression is fitted  based on a sample of
#' \eqn{N} observations of a response variable \eqn{y}, dependent on a set of
#' (currently up to two) covariates, assuming \eqn{y} has a distribution from
#' the exponential family. In addition, GeDS methodology is
#' implemented both in the context of Generalized Additive Models (GAM) and
#' Functional Gradient Boosting (FGB). On the one hand, GAM consist of an 
#' additive modeling technique where the impact of the predictor variables is
#' captured through smooth (GeDS, in this case) functions. On the other hand,
#' GeDS incorporates gradient boosting machine learning technique by
#' implementing functional gradient descent algorithm to optimize general risk
#' functions utilizing component-wise GeDS estimates.
#'
#' @details
#' The GeDS estimation method is based on: first constructing a piecewise linear
#' fit (spline fit of order 2) at stage A  which captures the shape of the data
#' and; second approximating this fit with shape preserving (variation
#' diminishing) spline fits of higher orders 3, 4,\eqn{\dots} (i.e., degrees 2,
#' 3,\eqn{\dots}) at stage B. As a result of this,  GeDS estimates the number
#' and location of the knots and the order of the spline fit in a fast and
#' efficient way.
#'
#' The GeDS method was originally developed by Kaishev et al. (2016) assuming
#' the response variable \eqn{y} to be normally distributed and a corresponding
#' \emph{Mathematica} code was provided.
#'
#' The GeDS method was extended by Dimitrova et al. (2023) to cover any
#' distribution from the exponential family. The \pkg{GeDS} \R package presented
#' here includes an enhanced \R implementation of the original Normal GeDS 
#' \emph{Mathematica} code due to Kaishev et al. (2016), implemented as the
#' \code{\link{NGeDS}} function and a generalization of it in the function
#' \code{\link{GGeDS}} which covers the case of any distribution from the
#' exponential family.
#'
#' The \pkg{GeDS} package allows also to fit two dimensional response surfaces
#' currently implemented only in the Normal case via the function
#' \code{\link{NGeDS}}. It also allows to construct multivariate (predictor)
#' models with a GeD spline component and a parametric component (see the
#' functions \code{\link{f}}, \code{\link[=formula.GeDS]{formula}},
#' \code{\link{NGeDS}} and \code{\link{GGeDS}} for details).
#' 
#' Dimitrova et al. (2024) have recently made significant enhancements to the
#' \pkg{GeDS} methodology, extending it towards the benchmark of generalized
#' additive models (GAM) and functional gradient boosting (FGB). On the one hand,
#' generalized additive models are encompassed by implementing the
#' \emph{local-scoring} algorithm using normal GeD splines (i.e., 
#' \code{\link{NGeDS}}) as function smoothers within the \emph{backfitting}
#' iterations. This is implemented via the function \code{\link{NGeDSgam}}. On
#' the other hand, \pkg{GeDS} package incorporates functional gradient descent
#' algorithm by utilizing normal GeD splines (i.e., \code{\link{NGeDS}}) as base
#' learners. This is implemented via the function \code{\link{NGeDSboost}}.
#'
#' The outputs of both \code{\link{NGeDS}} and \code{\link{GGeDS}} functions are
#' \code{\link{GeDS-class}} objects, while the outputs of \code{\link{NGeDSgam}}
#' and \code{\link{NGeDSboost}} are \code{\link{GeDSgam-class}} and
#' \code{\link{GeDSboost-class}} objects, respectively. As described in
#' Kaishev et al. (2016), Dimitrova et al. (2023) and  Dimitrova et al. (2024),
#' the final GeDS fit is the one whose order is chosen according to a strategy
#' described in stage B of the algorithm. However, \code{\link{GeDS-class}},
#' \code{\link{GeDSgam-class}} and \code{\link{GeDSboost-class}} objects contain
#' second, third and fourth order spline fits and the user has the possibility
#' to choose among them.
#'
#' This package also includes some datasets where GeDS regression proves to be
#' very efficient and some user friendly functions that are designed to easily
#' extract required information. Several methods are also provided to handle
#' GeDS, GAM-GeDS and FGB-GeDS output results (see \code{\link{GeDS-class}},
#' \code{\link{GeDSgam-class}} and \code{\link{GeDSboost-class}}, respectively).
#'
#' Throughout this document, we use the terms GeDS predictor model, GeDS
#' regression and GeDS fit interchangeably.
#'
#' Please report any issue arising or bug in the code to
#' \email{emilio.saenz-guillen@bayes.city.ac.uk}.
#'
#' \tabular{rl}{
#' Package: \tab GeDS\cr
#' Version: \tab 0.2.0 \cr
#' Date: \tab 2024-01-28\cr
#' License: \tab GPL-3 \cr
#' }
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
#' Dimitrova, D. S., Guillen, E. S. and Kaishev, V. K.  (2024).
#' \pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
#' Models and Functional Gradient Boosting, based on Geometrically Designed
#' (GeD) Splines.
#'
#' @aliases GeDS
#' @keywords package
#'
#' @author  Dimitrina S. Dimitrova <D.Dimitrova@city.ac.uk>,
#' Emilio S. Guillen <emilio.saenz-guillen@bayes.city.ac.uk>,
#' Vladimir K. Kaishev <V.Kaishev@city.ac.uk>,
#' Andrea Lattuada <Andrea.Lattuada@unicatt.it> and
#' Richard J. Verrall <R.J.Verrall@city.ac.uk>
#'
#' 
#' @docType package
NULL

