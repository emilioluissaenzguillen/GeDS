################################################################################
################################################################################
################################## GeDS Class ##################################
################################################################################
################################################################################
#' @title GeDS Class
#' @name GeDS-class
#' @description
#' A fitted GeDS object returned by the function \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}, inheriting the methods for class \code{"GeDS"}.
#' Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict},
#' \code{plot}, and \code{lines} are available.
#'
#' @slot Type Character string indicating the type of regression performed.
#' This can be \code{"LM - Univ"}/\code{"LM - Biv"} or
#' \code{"GLM - Univ"}/\code{"GLM - Biv"}, respectively corresponding to Normal
#' univariate/bivariate GeDS (implemented by \code{\link{NGeDS}}), and to
#' generalized (GNM-GLM) univariate/bivariate GeDS (implemented by
#' \code{\link{GGeDS}}).
#' @slot Linear.Knots vector containing the locations of the knots of the
#' second order GeD spline fit produced at stage A.
#' @slot Quadratic.Knots vector containing the locations of the knots of the
#' third order GeD spline fit produced in stage B.
#' @slot Cubic.knots  Vector containing the locations of the knots of the fourth
#' order GeD spline fit produced in stage B.
#' @slot Dev.Linear deviance of the second order GeD spline fit, produced in
#' stage A.
#' @slot Dev.Quadratic deviance of the third order GeD spline fit, produced in
#' stage B.
#' @slot Dev.Cubic deviance of the fourth order GeD spline fit, produced in
#' stage B.
#' @slot RSS vector containing the deviances of the second order spline
#' fits computed at each stage A's GeDS iteration. 
#' @slot Linear list containing the results from running \code{\link{SplineReg}}
#' function to fit the second order spline fit of stage A.
#' @slot Quadratic list containing the results from running \code{\link{SplineReg}}
#' function used to fit the third order spline fit in stage B.
#' @slot Cubic list containing the results from a \code{\link{SplineReg}}
#' function used to fit the fourth order spline fit in stage B.
#' @slot Stored Matrix containing the knot locations estimated at each iteration
#' of stage A.
#' @slot Args list containing the input arguments passed on the
#' \code{\link{Fitters}} functions.
#' @slot Call \code{call} to the \code{\link{Fitters}} functions.
#' @slot Nintknots the final number of internal knots of the second order GeD
#' spline fit produced in stage A.
#' @slot iters number of iterations performed during stage A of the GeDS fitting
#' procedure.
#' @slot Guesses initial values for the coefficients used at each iteration of
#' stage A in order to estimate the spline coefficients. Since the initial
#' values are used only in the IRLS procedure, this slot is empty if the object
#' is not created by \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}}
#' functions.
#' @slot Coefficients matrix containing the fitted coefficients of the GeD
#' spline regression component and the parametric component at each iteration
#' of stage A.
#' @slot deviance vector containing the deviances of the second order spline
#' fits computed at each IRLS iteration in stage A. Since the IRLS procedure is
#' used only in \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}}, this
#' slot is empty if the object is not created by one of these functions.
#' @slot iterIrls vector containing the numbers of IRLS iterations for all
#' iterations of stage A cumulatively. Since the IRLS procedure is used only in
#' \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}}, this slot is empty
#' if the object is not created by one of these functions.
#' @slot stopinfo list of values providing information related to the stopping
#' rule of stage A of GeDS. The sub-slots of \code{stopinfo} are \code{phis},
#' \code{phis_star}, \code{oldintc} and \code{oldslp}. The sub-slot \code{phis}
#' is a vector containing the values of the ratios of deviances (or the
#' difference of deviances if the \code{LR} stopping rule was chosen). The
#' sub-slots \code{phis_star}, \code{oldintc} and \code{oldslp} are non-empty
#' slots if the \code{SR} stopping rule was chosen. These respectively contain
#' the values at each iteration of stage A of \eqn{\hat{\phi}_{\kappa}},
#' \eqn{\hat{\gamma}_0} and \eqn{\hat{\gamma}_1}. See Dimitrova et al. (2023)
#' for further details on these parameters.
#' @slot Formula the model \code{\link[=formula.GeDS]{formula}}.
#' @slot extcall \code{call} to the \code{\link{NGeDS}} or \code{\link{GGeDS}}
#' functions.
#' @slot terms \code{terms} object containing information on the model frame.
#'
#' @aliases GeDS-Class GeDS-class
#' @rdname GeDS-class
#' 
#' @references
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}

setClass(
  "GeDS",
  representation(
    Type = "character", Linear.Knots = "numeric", Quadratic.Knots = "numeric",
    Cubic.Knots = "numeric", RMS.Linear  = "numeric", RMS.Quadratic = "numeric",
    RMS.Cubic  = "numeric", Knots  = "numeric", RSS = "numeric",
    Linear = "list", Quadratic = "list", Cubic = "list", Stored = "matrix",
    Args = "list", Call = "call", Nintknots = "numeric", iters = "numeric",
    Guesses = "matrix", Coefficients = "matrix", deviance = "numeric",
    iter = "numeric", stopinfo = "list", Formula = "formula", extcall = "call"
    )
  )


################################################################################
################################################################################
################################ GeDSboost Class ###############################
################################################################################
################################################################################
#' @title GeDSboost Class
#' @name GeDSboost-class
#' @description
#' A fitted GeDSboost object returned by the function \code{\link{NGeDSboost}},
#' inheriting the methods for class \code{"GeDSboost"}. Methods for functions
#' \code{coef}, \code{knots}, \code{plot}, \code{print}, \code{predict},
#' \code{visualize_boosting}, and \code{bl_imp} are available.
#' 
#' @slot extcall call to the \code{\link{NGeDSboost}} function.
#' @slot formula a formula object representing the model to be fitted.
#' @slot args 
#' a list containing the arguments passed to the \code{\link{NGeDSboost}}
#' function. This list includes:
#' \itemize{
#'   \item \code{response}: \code{data.frame} containing the response variable
#'   observations.
#'   \item \code{predictors}: \code{data.frame} containing the observations
#'   corresponding to the predictor variables included in the model.
#'   \item \code{base_learners}: description of model's base learners.
#'   \item \code{family}: the statistical family; the possible options are
#'   \itemize{
#'   \item \code{mboost::Binomial(type = c("adaboost", "glm"),
#'   link = c("logit", "probit", "cloglog", "cauchit", "log"), ...)}
#'   \item \code{mboost::Gaussian()}
#'   \item \code{mboost::Poisson()}
#'   \item \code{mboost::GammaReg(nuirange = c(0, 100))}
#'   }
#'   Other \code{mboost} families may be suitable; however, these have not yet
#'   been thoroughly tested and are therefore not recommended for use.
#'   \item \code{initial_learner}: if \code{TRUE} a \code{\link{NGeDS}} or
#'   \code{\link{GGeDS}} fit was used as the initial learner; otherwise, the
#'   empirical risk minimizer corresponding to the selected \code{family}
#'   was employed.
#'   \item \code{int.knots_init}: if \code{initial_learner = TRUE}, this
#'   corresponds to the maximum number of internal knots set in the
#'   \code{\link{NGeDS}}/\code{\link{GGeDS}} function before the initial learner
#'   fit.
#'   \item \code{shrinkage}: shrinkage/step-length/learning rate utilized
#'   throughout the boosting iterations.
#'   \item \code{normalize_data}: if \code{TRUE}, then response and predictors
#'   were standardized before running the FGB algorithm.
#'   \item \code{X_mean}: mean of the predictor variables (only if
#'   \code{normalize_data = TRUE}, otherwise this is \code{NULL}).
#'   \item \code{X_sd}: standard deviation of the predictors (only if
#'   \code{normalize_data = TRUE}, otherwise this is \code{NULL}).
#'   \item \code{Y_mean}: mean of the response variable (only if
#'   \code{normalize_data = TRUE}, otherwise this is \code{NULL}).
#'   \item \code{Y_sd}: standard deviation of the response variable (only if
#'   \code{normalize_data = TRUE}, otherwise this is \code{NULL}).
#'}
#' @slot models A list containing the 'model' generated at each boosting
#' iteration. Each of these \code{models} includes:
#' \itemize{
#'  \item \code{best_bl}: fit of the base learner that minimized the residual
#'  sum of squares (RSS) in fitting the gradient at the \emph{i}-th boosting
#'  iteration.
#'   \item \code{Y_hat}: model fitted values at the \emph{i}-th boosting
#'   iteration.
#'   \item \code{base_learners}: knots and polynomial coefficients for each of the
#'   base-learners at the \emph{i}-th boosting iteration.  
#' }
#' @slot final_model A list detailing the final GeDSboost model after the
#' gradient descent algorithm is run:
#' \itemize{
#' \item \code{model_name}: the boosting iteration corresponding to the final
#'  model.
#'  \item \code{DEV}: deviance of the final model.
#'  \item \code{Y_hat}: fitted values.
#'  \item \code{base_learners}: a list containing, for each base-learner, the
#'  intervals defined by the piecewise linear fit and its corresponding
#'  polynomial coefficients. It also includes the knots corresponding to each
#'  order fit, which result from computing the corresponding averaging knot
#'  location. See Kaishev et al. (2016) for details. If the number of internal
#'  knots of the final linear fit is less than $n-1$, the averaging knot location
#'  is not computed.
#'  \item \code{Linear.Fit}/\code{Quadratic.Fit}/\code{Cubic.Fit}: final linear,
#'  quadratic and cubic fits in B-spline form. These include the same elements
#'  as \code{Linear}, \code{Quadratic} and \code{Cubic} in a \code{\link{GeDS-class}}
#'  object (see \code{\link{SplineReg}} for details).
#'  }
#'   
#' @slot predictions a list containing the predicted values obtained for each of
#' the fits (linear, quadratic and cubic).
#'
#' @slot internal_knots a list detailing the internal knots obtained for each of
#' the different order fits (linear, quadratic, and cubic).
#'
#' @aliases GeDSboost-Class GeDSboost-class
#' @rdname GeDSboost-class
#' 
#' @references
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

setClass(
  "GeDSboost",
  representation(
    formula = "formula",
    args = "list",
    models = "list",
    final_model = "list",
    predictions = "list",
    internal_knots = "list"
  )
)


################################################################################
################################################################################
################################# GeDSgam Class ################################
################################################################################
################################################################################
#' @title GeDSgam Class
#' @name GeDSgam-class
#' @description
#' A fitted GeDSgam object returned by the function \code{\link{NGeDSgam}},
#' inheriting the methods for class \code{"GeDSgam"}. Methods for functions
#' \code{coef}, \code{knots}, \code{plot}, \code{print} and \code{predict} are 
#' available.
#' 
#' @slot extcall call to the \code{\link{NGeDSgam}} function.
#' @slot formula a formula object representing the model to be fitted.
#' @slot args 
#' a list containing the arguments passed to the \code{\link{NGeDSgam}}
#' function. This list includes:
#' \itemize{
#'   \item \code{response}: \code{data.frame} containing the response variable
#'   observations.
#'   \item \code{predictors}: \code{data.frame} containing the corresponding 
#'   observations of the predictor variables included in the model.
#'   \item \code{base_learners}: description of the model's base learners
#'   ('smooth functions').
#'   \item \code{family}: the statistical family. The possible options are
#'   \itemize{
#'   \item \code{binomial(link = "logit", "probit", "cauchit", "log", "cloglog")}
#'   \item \code{gaussian(link = "identity", "log", "inverse")}
#'   \item \code{Gamma(link = "inverse", "identity", "log")}
#'   \item \code{inverse.gaussian(link = "1/mu^2", "inverse", "identity", "log")}
#'   \item \code{poisson(link = "log", "identity", "sqrt")}
#'   \item \code{quasi(link = "identity", variance = "constant")}
#'   \item \code{quasibinomial(link = "logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")}
#'   \item \code{quasipoisson(llink = "logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")}   
#'   }
#'   \item \code{normalize_data}: if \code{TRUE}, then response and predictors
#'   were standardized before running the local-scoring algorithm.
#'   \item \code{X_mean}: mean of the predictor variables (only if
#'   \code{normalize_data = TRUE}).
#'   \item \code{X_sd}: standard deviation of the predictors (only if
#'   \code{normalize_data = TRUE}, else is \code{NULL}).
#'   \item \code{Y_mean}: mean of the response variable (only if
#'   \code{normalize_data = TRUE}, else is \code{NULL}).
#'   \item \code{Y_sd}: standard deviation of the response variable (only if
#'   \code{normalize_data = TRUE}, else is \code{NULL}).
#'}
#' @slot final_model A list detailing the final GeDSgam model selected after
#' running the local scoring algorithm. The chosen model minimizes deviance
#' across all models generated by each local-scoring iteration. This list
#' includes:
#' \itemize{
#'   \item \code{model_name}: local-scoring iteration that yielded the "best"
#'   model. Note that when \code{family = "gaussian"}, it will always correspond
#'   to \code{iter1}, as only one local-scoring iteration is conducted in this
#'   scenario. This occurs because, with \code{family = "gaussian"}, the
#'   algorithm is tantamount to directly implementing backfitting.
#'   \item \code{DEV}: the deviance for the fitted predictor model, defined as
#'   in Dimitrova et al. (2023), which for \code{family = "gaussian"} coincides
#'   with the Residual Sum of Squares.
#'   \item \code{Y_hat}: fitted values.
#'   \itemize{
#'      \item \code{eta}: additive predictor.
#'      \item \code{mu}: vector of means.
#'      \item \code{z}: adjusted dependent variable.      
#'   }
#'   \item \code{base_learners}: a list containing, for each base-learner, the
#'   corresponding linear fit piecewise polynomial coefficients. It includes the
#'   knots for each order fit, resulting from computing the averaging knot
#'   location. Although if the number of internal knots of the final linear fit
#'   is less than $n-1$, the averaging knot location is not computed.
#'   \item \code{Linear.Fit}: final model linear fit in B-spline form.
#'   See \code{\link{SplineReg}} for details.
#'   \item \code{Quadratic.Fit}: quadratic fit obtained via Schoenberg variation
#'   diminishing spline approximation. See \code{\link{SplineReg}} for details.
#'   \item \code{Cubic.Fit}: cubic fit obtained via Schoenberg variation
#'   diminishing spline approximation. See \code{\link{SplineReg}} for details.
#' } 
#' @slot predictions A list containing the predicted values obtained for each of
#' the fits (linear, quadratic, and cubic). Each of the predictions contains
#' both the additive predictor \code{eta} and the vector of means \code{mu}.
#'
#' @slot internal_knots A list detailing the internal knots obtained for the fits
#' of different order (linear, quadratic, and cubic).
#'
#' @aliases GeDSgam-Class GeDSgam-class
#' @rdname GeDSgam-class
#'
#' @references 
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

setClass(
  "GeDSgam",
  representation(
    formula = "formula",
    args = "list",
    final_model = "list",
    predictions = "list"
  )
)

