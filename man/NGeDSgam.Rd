% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NGeDSgam.R
\name{NGeDSgam}
\alias{NGeDSgam}
\title{NGeDSgam: Local Scoring Algorithm with GeD Splines in Backfitting}
\usage{
NGeDSgam(
  formula,
  family = "gaussian",
  data,
  weights = NULL,
  normalize_data = FALSE,
  min_iterations,
  max_iterations,
  phi_gam_exit = 0.99,
  q_gam = 2L,
  beta = 0.5,
  phi = 0.99,
  internal_knots = 500L,
  q = 2L,
  higher_order = TRUE
)
}
\arguments{
\item{formula}{A description of the model structure to be fitted,
specifying both the dependent and independent variables. Unlike \code{\link{NGeDS}}
and \code{\link{GGeDS}}, this formula supports multiple additive (normal) GeD
spline regression components as well as linear components. For example, setting
\code{formula = Y ~ f(X1) + f(X2) + X3} implies using a normal linear GeD
spline as the smoother for \code{X1} and for \code{X2}, while for \code{X3} a
linear model would be used.}

\item{family}{A character string indicating the response variable distribution
and link function to be used. Default is \code{"gaussian"}. This should be a
character or a \code{\link[stats]{family}} object.}

\item{data}{A \code{data.frame} containing the variables referenced in the formula.}

\item{weights}{An optional vector of "prior weights" to be put on the
observations during the fitting process. It should be \code{NULL} or a numeric
vector of the same length as the response variable defined in the formula.}

\item{normalize_data}{A logical that defines whether the data should be
normalized (standardized) before fitting the baseline linear model, i.e.,
before running the local-scoring algorithm. Normalizing the data involves
scaling the predictor variables to have a mean of 0 and a standard deviation
of 1. This process alters the scale and interpretation of the knots and
coefficients estimated. Default is equal to \code{FALSE}.}

\item{min_iterations}{Optional parameter to manually set a minimum number of
local-scoring iterations to be run. If not specified, it defaults to 0L.}

\item{max_iterations}{Optional parameter to manually set the maximum number
of local-scoring iterations to be run. If not specified, it defaults to \code{100L}.
This setting serves as a fallback when the stopping rule, based on
consecutive deviances and tuned by \code{phi_gam_exit} and \code{q_gam},
does not trigger an earlier termination (see Dimitrova et al. (2025)).
Therefore, users can increase/decrease the number of local-scoring iterations,
by increasing/decreasing the value \code{phi_gam_exit} and/or \code{q_gam},
or directly specify \code{max_iterations}.}

\item{phi_gam_exit}{Convergence threshold for local-scoring and backfitting.
Both algorithms stop when the relative change in the deviance is below this
threshold. Default is \code{0.99}.}

\item{q_gam}{Numeric parameter which allows to fine-tune the stopping rule of
the local-scoring and backfitting iterations. By default equal to \code{2L}.}

\item{beta}{Numeric parameter in the interval \eqn{[0,1]}
tuning the knot placement in stage A of GeDS, for each of the GeD spline
components of the model. Default is equal to \code{0.5}.
See Details in \code{\link{NGeDS}}.}

\item{phi}{Numeric parameter in the interval \eqn{(0,1)} specifying the
threshold for the stopping rule  (model selector) in stage A of GeDS, for each
of the GeD spline components of the model. Default is equal to \code{0.99}.
See Details in \code{\link{NGeDS}}.}

\item{internal_knots}{The maximum number of internal knots that can be added
by the GeDS smoothers at each backfitting iteration, effectively setting the
value of \code{max.intknots} in \code{\link{NGeDS}} at each backfitting
iteration. Default is \code{500L}.}

\item{q}{Numeric parameter which allows to fine-tune the stopping rule of
stage A of GeDS, for each of the GeD spline components of the model. By
default equal to \code{2L}. See Details in \code{\link{NGeDS}}.}

\item{higher_order}{a logical that defines whether to compute the higher order
fits (quadratic and cubic) after the local-scoring algorithm is run. Default
is \code{TRUE}.}
}
\value{
An object of class \code{"GeDSgam"} (a named list) with components:
\describe{
  \item{extcall}{Call to the \code{\link{NGeDSgam}} function.}
  \item{formula}{A formula object representing the model to be fitted.}
  \item{args}{A list containing the arguments passed to the \code{\link{NGeDSgam}} function.
  This list includes:
    \describe{
      \item{\code{response}}{\code{data.frame} containing the response variable
      observations.}
      \item{\code{predictors}}{\code{data.frame} containing the corresponding
      observations of the predictor variables included in the model.}
      \item{\code{base_learners}}{Description of the model's base learners
      ("smooth functions").}
      \item{\code{family}}{The statistical family. The possible options are:
      \itemize{
      \item{\code{binomial(link = "logit", "probit", "cauchit", "log", "cloglog")},}
      \item{\code{gaussian(link = "identity", "log", "inverse")},}
      \item{\code{Gamma(link = "inverse", "identity", "log")},}
      \item{\code{inverse.gaussian(link = "1/mu^2", "inverse", "identity", "log")},}
      \item{\code{poisson(link = "log", "identity", "sqrt")},}
      \item{\code{quasi(link = "identity", variance = "constant")},}
      \item{\code{quasibinomial(link = "logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")} and}
      \item{\code{quasipoisson(link = "log", "identity", "sqrt")}.}
      }         
      }
      \item{\code{normalize_data}}{If \code{TRUE}, then response and predictors
      were standardized before running the local-scoring algorithm.}
      \item{\code{X_mean}}{Mean of the predictor variables (only if
      \code{normalize_data = TRUE}).}
      \item{\code{X_sd}}{Standard deviation of the predictors (only if
      \code{normalize_data = TRUE}, otherwise this is \code{NULL}).}
      \item{\code{Y_mean}}{Mean of the response variable (only if
      \code{normalize_data = TRUE}, otherwise this is \code{NULL}).}
      \item{\code{Y_sd}}{Standard deviation of the response variable (only if
      \code{normalize_data = TRUE}, otherwise this is \code{NULL}).}
    }
  }

  \item{final_model}{A list detailing the final \code{"GeDSgam"} model selected after
  running the local scoring algorithm. The chosen model minimizes deviance
  across all models generated by each local-scoring iteration. This list includes:
    \describe{
    \item{\code{model_name}}{Local-scoring iteration that yielded the "best"
    model. Note that when \code{family = "gaussian"}, it will always correspond
    to \code{iter1}, as only one local-scoring iteration is conducted in this
    scenario. This occurs because, when \code{family = "gaussian"}, the
    algorithm is tantamount to directly implementing backfitting.}
    \item{\code{dev}}{Deviance of the final model. For \code{family = "gaussian"}
    this coincides with the Residual Sum of Squares.}
    \item{\code{Y_hat}}{Fitted values, including:
      - \code{eta}: the additive predictor,
      - \code{mu}: the vector of means,
      - \code{z}: the adjusted dependent variable.}
    \item{\code{base_learners}}{A list containing, for each base-learner, the
    corresponding linear fit piecewise polynomial coefficients. It includes the
    knots for each order fit, resulting from computing the averaging knot
    location. Although if the number of internal knots of the final linear fit
    is less than \eqn{n-1}, the averaging knot location is not computed.}
    \item{\code{linear.fit}}{Final linear fit in B-spline form (see \code{\link{SplineReg}}).}
    \item{\code{quadratic.fit}}{Quadratic fit obtained via Schoenberg
    variation diminishing approximation (see \code{\link{SplineReg}}).}
    \item{\code{cubic.fit}}{Cubic fit obtained via via Schoenberg variation
    diminishing approximation (see \code{\link{SplineReg}}).}
    }
  }
  \item{predictions}{A list containing the predicted values obtained for each of
  the fits (linear, quadratic, and cubic). Each of the predictions contains
  both the additive predictor \code{eta} and the vector of means \code{mu}.}

  \item{internal_knots}{A list detailing the internal knots obtained for the
  fits of different order (linear, quadratic, and cubic).}
}
}
\description{
Implements the Local Scoring Algorithm (Hastie and Tibshirani
(1986)), applying normal linear GeD splines (i.e., \code{\link{NGeDS}}
function) to fit the targets within each backfitting iteration. Higher order 
fits are computed by pursuing stage B of GeDS after the local-scoring algorithm
is run.
}
\details{
The  \code{NGeDSgam} function employs the local scoring algorithm to fit a
generalized additive model (GAM). This algorithm iteratively fits weighted
additive models by backfitting. Normal linear GeD splines, as well as linear
learners, are supported as function smoothers within the backfitting
algorithm. The local-scoring algorithm ultimately produces a linear fit.
Higher order fits (quadratic and cubic) are then computed by calculating the
Schoenberg’s variation diminishing spline (VDS) approximation of the linear
fit.

On the one hand, \code{NGeDSgam} includes all the parameters of
\code{\link{NGeDS}}, which in this case tune the function smoother fit at each
backfitting iteration. On the other hand, \code{NGeDSgam} includes some
additional parameters proper to the local-scoring procedure. We describe
the main ones as follows. 

The \code{family} chosen determines the link function, adjusted dependent
variable and weights to be used in the local-scoring algorithm. The number of
local-scoring and backfitting iterations is controlled by a
\emph{Ratio of Deviances} stopping rule similar to the one presented for
\code{\link{NGeDS}}/\code{\link{GGeDS}}. In the same way \code{phi} and \code{q}
tune the stopping rule of \code{\link{NGeDS}}/\code{\link{GGeDS}},
\code{phi_gam_exit} and \code{q_gam} tune the stopping rule of \code{NGeDSgam}.
The user can also manually control the number of local-scoring iterations
through \code{min_iterations} and \code{max_iterations}.

A model term wrapped in \code{offset()} is treated as a known (fixed) component
and added directly to the linear predictor when fitting the model. In case
more than one covariate is fixed, the user should sum the corresponding
coordinates of the fixed covariates to produce one common \eqn{N}-vector of
coordinates. See \code{\link[=formula.GeDS]{formula}}.
}
\examples{

# Load package
library(GeDS) 

data(airquality) 
data = na.omit(airquality)
data$Ozone <- data$Ozone^(1/3)

formula = Ozone ~ f(Solar.R) + f(Wind, Temp)
Gmodgam <- NGeDSgam(formula = formula, data = data,
phi = 0.8)
MSE_Gmodgam_linear <- mean((data$Ozone - Gmodgam$predictions$pred_linear)^2)
MSE_Gmodgam_quadratic <- mean((data$Ozone - Gmodgam$predictions$pred_quadratic)^2)
MSE_Gmodgam_cubic <- mean((data$Ozone - Gmodgam$predictions$pred_cubic)^2)

cat("\n", "MEAN SQUARED ERROR", "\n",
"Linear NGeDSgam:", MSE_Gmodgam_linear, "\n",
"Quadratic NGeDSgam:", MSE_Gmodgam_quadratic, "\n",
"Cubic NGeDSgam:", MSE_Gmodgam_cubic, "\n")

## S3 methods for class 'GeDSgam'
# Print 
print(Gmodgam); summary(Gmodgam)
# Knots
knots(Gmodgam, n = 2)
knots(Gmodgam, n = 3)
knots(Gmodgam, n = 4)
# Coefficients
coef(Gmodgam, n = 2)
coef(Gmodgam, n = 3)
coef(Gmodgam, n = 4)
# Wald-type confidence intervals
confint(Gmodgam, n = 2)
confint(Gmodgam, n = 3)
confint(Gmodgam, n = 4)
# Deviances
deviance(Gmodgam, n = 2)
deviance(Gmodgam, n = 3)
deviance(Gmodgam, n = 4)

}
\references{
Hastie, T. and Tibshirani, R. (1986). Generalized Additive Models.
\emph{Statistical Science} \strong{1 (3)} 297 - 310. \cr
DOI: \doi{10.1214/ss/1177013604}

Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
Geometrically designed, variable knot regression splines.
\emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
DOI: \doi{10.1007/s00180-015-0621-7}

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
\seealso{
\code{\link{NGeDS}}; \code{\link{GGeDS}}; S3 methods such as
\code{\link[=coef.GeDSgam]{coef}}, \code{\link[=confint.GeDSgam]{confint}},
\code{\link{deviance.GeDSgam}}, \code{\link{family}}, \code{\link{formula}},
\code{\link[=knots.GeDSgam]{knots}}, \code{\link{logLik}},
\code{\link[=predict.GeDSgam]{predict}}, \code{\link[=print.GeDSgam]{print}},
\code{\link[=summary.GeDSgam]{summary}}.
}
