#' GeDS Class
#' @name GeDS-class
#' @rdname GeDS-class
#' @aliases GeDS-Class GeDS-class
#'
#' @description A fitted GeDS object returned by functions \code{\link{NGeDS}} or \code{\link{GGeDS}}
#' inheriting the methods from class \code{"GeDS"}.
#'  Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict}, \code{plot},
#' and  \code{lines} are available.
#'
#' @slot Type Character string indicating the type of the regression performed.
#' One of \code{"LM - Univ"}, \code{"LM - Biv"} or \code{"GLM - Univ"} corresponding to
#' the Normal univariate GeDS, the Normal bivariate GeDS performed by \code{\link{NGeDS}} and the
#' generalized (GNM-GLM) univariate GeDS performed by \code{\link{GGeDS}}.
#' @slot Linear.Knots  Vector containing the
#'  locations of the knots of the second order GeDS spline fit generated at stage A.
#' @slot Quadratic.Knots  Vector containing the locations
#' of the knots of the third order GeDS spline fitted in stage B.
#' @slot Cubic.knots  Vector containing the locations of
#' the knots of the fourth order GeDS spline fitted in stage B.
#' @slot Dev.Linear Deviance of the second order GeD spline fit of stage A.
#' @slot Dev.Quadratic Deviance of the third order GeD spline fit of stage B.
#' @slot Dev.Cubic Deviance of the fourth order GeD spline fit of stage B.
#' @slot Linear List containing the results from running a \code{\link{SplineReg}}
#' function used to fit the second order spline of stage A.
#' @slot Quadratic List containing the results from running \code{\link{SplineReg}}
#' function used to fit the third order spline in stage B.
#' @slot Cubic List containing the results from a \code{\link{SplineReg}}
#' function used to fit the fourth order spline in stage B.
#' @slot Stored Matrix containing the knot locations estimated at each step of stage A.
#' @slot Args List containing the input arguments passed on the \code{\link{Fitters}} functions.
#' @slot Call \code{call} to the \code{\link{Fitters}} functions.
#' @slot Nintknots The final number of internal knots of the second order GeD spline fit of stage A.
#' @slot iters Number of iterations performed in stage A  of the GeDS fitting procedure.
#' @slot Guesses Initial values for the coefficients used at each
#' iteration of stage A in order to estimate the spline coefficients.
#' Since the initial values are used only in the IRLS procedure,
#' this slot is empty if the object is not created by \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}}
#'  functions.
#' @slot Coefficients Matrix containing the fitted coefficients of the GeD spline regression  component and the
#' parametric component at each iteration of stage A.
#' @slot deviance Vector containing the deviances of the second order spline fits computed at each IRLS
#' iteration in stage A.  Since the IRLS procedure is used only in \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}},
#' this slot is empty if the object is not created by one of these functions.
#' @slot iterIrls Vector containing the numbers of IRLS iterations for all iterations of stage A cumulatively.
#' Since the IRLS procedure is used only in \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}},
#' this slot is empty if the object is not created by one of these functions.
#' @slot stopinfo List of values providing information related to
#' the stopping rule of stage A of GeDS. The sub-slots of \code{stopinfo} are \code{phis},  \code{phis_star},
#' \code{oldintc} and \code{oldslp}. The sub-slot \code{phis} is a vector containing the values
#' of the ratios of deviances (or the difference of deviances if the \code{LR} stopping
#' rule was chosen). The sub-slots \code{phis_star}, \code{oldintc} and \code{oldslp} are non-empty slots
#' if the \code{SR} stopping rule was chosen. They contain respectively \eqn{\hat{\phi}_{\kappa}}, \eqn{\hat{\gamma}_0} and
#' \eqn{\hat{\gamma}_1} computed at each iteration of stage A, see Dimitrova et al. (2017).
#' @slot Formula The model \code{\link[=formula.GeDS]{formula}}.
#' @slot extcall \code{call} to the \code{\link{NGeDS}} or \code{\link{GGeDS}} functions.
#' @slot terms \code{terms} object containing information on the model frame.
#'
#' @references
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{https://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#'
#'
setClass(
  "GeDS",
  representation(
    Type = "character", Linear.Knots = "numeric", Quadratic.Knots = "numeric",
    Cubic.Knots = "numeric", RMS.Linear  = "numeric", RMS.Quadratic = "numeric",
    RMS.Cubic  = "numeric", Knots  = "numeric", RSS = "numeric",
    Linear = "list", Quadratic = "list", Cubic = "list", Stored = "matrix",
    Args = "list", Call = "call", Nintknots = "numeric", iters = "numeric", Guesses = "matrix",
    Coefficients = "matrix", deviance = "numeric", iter = "numeric", stopinfo = "list", Formula = "formula", extcall = "call"
    )
  )


#' GeDSboost Class
#'
#' @name GeDSboost-class
#' @rdname GeDSboost-class
#' @aliases GeDSboost-Class GeDSboost-class
#'
#' @description 
#' A class representing a fitted GeDSboost object. This object encapsulates details about the
#' GeDS gradient boosting model, including its formula, arguments, models, the final model,
#' predictions, and internal knots. Objects of this class are typically returned by the \code{\link{NGeDSboost}} function.
#'
#' @slot formula A formula object representing the relationship between the variables used in the model.
#'
#' @slot args A list containing the arguments passed to the \code{\link{NGeDSboost}} function. 
#' It includes details about base learners, family, initial learner, knots for initialization, normalization options, outcome, 
#' predictors, shrinkage parameter, mean and standard deviation of predictors (X_mean, X_sd) and outcome (Y_mean, Y_sd).
#'
#' @slot models A list containing the models generated during the boosting process.
#'
#' @slot final_model A list detailing the final GeDS model after all boosting iterations. This includes 
#' base learners, the best base learner (best_bl), cubic fit, model name, quadratic fit, and predicted values (Y_hat).
#'
#' @slot predictions A list containing predictions obtained from the model for different fit types 
#' (linear, quadratic, and cubic).
#'
#' @slot internal_knots A list detailing the internal knots used for the different fit types 
#' (linear, quadratic, and cubic) during the boosting process.
#'
#' @references 
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#' @seealso \code{\link{NGeDSboost}} which generates objects of this class.
#'
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


#' GeDSgam Class
#'
#' @name GeDSgam-class
#' @rdname GeDSgam-class
#' @aliases GeDSgam-Class GeDSgam-class
#'
#' @description 
#' A class representing a fitted GeDSgam object returned by function \code{\link{NGeDSgam}}.
#'
#' @slot formula A formula object representing the relationship between the variables used in the model.
#'
#' @slot args A list containing:
#' \describe{
#'   \item{base_learners}{List of base learners used in the model.}
#'   \item{family}{Family of the model.}
#'   \item{normalize_data}{Indicator whether the data should be normalized.}
#'   \item{outcome}{Outcome variable of the model.}
#'   \item{predictors}{Predictors used in the model.}
#' }
#'
#' @slot final_model The final GeDSgam model object containing:
#' \describe{
#'   \item{base_learners}{List of base learners used in the model. Each base learner contains:
#'     \describe{
#'       \item{coefficients}{Coefficients of the base learner.}
#'       \item{cubic.int.knots}{Knots used in the cubic fit.}
#'       \item{knots}{General knots used in the model.}
#'       \item{quadratic.int.knots}{Knots used in the quadratic fit.}
#'     }
#'   }
#'   \item{Y_hat}{Predicted values.}
#'   \item{model_name}{Name of the model.}
#'   \item{Quadratic.Fit}{Results of the quadratic fit.}
#'   \item{Cubic.Fit}{Results of the cubic fit.}
#' }
#'
#' @slot predictions A list containing:
#' \describe{
#'   \item{pred_linear}{Predictions from the linear fit.}
#'   \item{pred_quadratic}{Predictions from the quadratic fit.}
#'   \item{pred_cubic}{Predictions from the cubic fit.}
#' }
#'
#' @references 
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#'@seealso \code{\link{NGeDSgam}} which generates objects of this class.
#'
setClass(
  "GeDSgam",
  representation(
    formula = "formula",
    args = "list",
    final_model = "list",
    predictions = "list"
  )
)
