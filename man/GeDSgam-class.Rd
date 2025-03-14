% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeDSClass.R
\docType{class}
\name{GeDSgam-class}
\alias{GeDSgam-class}
\alias{GeDSgam-Class}
\title{GeDSgam Class}
\description{
A fitted GeDSgam object returned by the function \code{\link{NGeDSgam}},
inheriting the methods for class \code{"GeDSgam"}. Methods for functions
\code{coef}, \code{knots}, \code{plot}, \code{print} and \code{predict} are 
available.
}
\section{Slots}{

\describe{
\item{\code{extcall}}{call to the \code{\link{NGeDSgam}} function.}

\item{\code{formula}}{a formula object representing the model to be fitted.}

\item{\code{args}}{a list containing the arguments passed to the \code{\link{NGeDSgam}}
function. This list includes:
\itemize{
  \item \code{response}: \code{data.frame} containing the response variable
  observations.
  \item \code{predictors}: \code{data.frame} containing the corresponding 
  observations of the predictor variables included in the model.
  \item \code{base_learners}: description of the model's base learners
  ('smooth functions').
  \item \code{family}: the statistical family. The possible options are
  \itemize{
  \item \code{binomial(link = "logit", "probit", "cauchit", "log", "cloglog")}
  \item \code{gaussian(link = "identity", "log", "inverse")}
  \item \code{Gamma(link = "inverse", "identity", "log")}
  \item \code{inverse.gaussian(link = "1/mu^2", "inverse", "identity", "log")}
  \item \code{poisson(link = "log", "identity", "sqrt")}
  \item \code{quasi(link = "identity", variance = "constant")}
  \item \code{quasibinomial(link = "logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")}
  \item \code{quasipoisson(llink = "logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")}   
  }
  \item \code{normalize_data}: if \code{TRUE}, then response and predictors
  were standardized before running the local-scoring algorithm.
  \item \code{X_mean}: mean of the predictor variables (only if
  \code{normalize_data = TRUE}).
  \item \code{X_sd}: standard deviation of the predictors (only if
  \code{normalize_data = TRUE}, else is \code{NULL}).
  \item \code{Y_mean}: mean of the response variable (only if
  \code{normalize_data = TRUE}, else is \code{NULL}).
  \item \code{Y_sd}: standard deviation of the response variable (only if
  \code{normalize_data = TRUE}, else is \code{NULL}).
}}

\item{\code{final_model}}{A list detailing the final GeDSgam model selected after
running the local scoring algorithm. The chosen model minimizes deviance
across all models generated by each local-scoring iteration. This list
includes:
\itemize{
  \item \code{model_name}: local-scoring iteration that yielded the "best"
  model. Note that when \code{family = "gaussian"}, it will always correspond
  to \code{iter1}, as only one local-scoring iteration is conducted in this
  scenario. This occurs because, with \code{family = "gaussian"}, the
  algorithm is tantamount to directly implementing backfitting.
  \item \code{DEV}: the deviance for the fitted predictor model, defined as
  in Dimitrova et al. (2023), which for \code{family = "gaussian"} coincides
  with the Residual Sum of Squares.
  \item \code{Y_hat}: fitted values.
  \itemize{
     \item \code{eta}: additive predictor.
     \item \code{mu}: vector of means.
     \item \code{z}: adjusted dependent variable.      
  }
  \item \code{base_learners}: a list containing, for each base-learner, the
  corresponding linear fit piecewise polynomial coefficients. It includes the
  knots for each order fit, resulting from computing the averaging knot
  location. Although if the number of internal knots of the final linear fit
  is less than $n-1$, the averaging knot location is not computed.
  \item \code{Linear.Fit}: final model linear fit in B-spline form.
  See \code{\link{SplineReg}} for details.
  \item \code{Quadratic.Fit}: quadratic fit obtained via Schoenberg variation
  diminishing spline approximation. See \code{\link{SplineReg}} for details.
  \item \code{Cubic.Fit}: cubic fit obtained via Schoenberg variation
  diminishing spline approximation. See \code{\link{SplineReg}} for details.
}}

\item{\code{predictions}}{A list containing the predicted values obtained for each of
the fits (linear, quadratic, and cubic). Each of the predictions contains
both the additive predictor \code{eta} and the vector of means \code{mu}.}

\item{\code{internal_knots}}{A list detailing the internal knots obtained for the fits
of different order (linear, quadratic, and cubic).}
}}

\references{
Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
Geometrically designed variable knot splines in generalized (non-)linear
models.
\emph{Applied Mathematics and Computation}, \strong{436}. \cr
DOI: \doi{10.1016/j.amc.2022.127493}

Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. L. (2025).
\pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
Models and Functional Gradient Boosting, based on Geometrically Designed
(GeD) Splines. \emph{Manuscript submitted for publication.}
}
