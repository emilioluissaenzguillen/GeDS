################################################################################
##################################### COEF #####################################
################################################################################
#' @title Coef method for GeDS objects
#' @name coef.GeDS
#' @description
#' Methods for the functions \code{\link[stats]{coef}} that allow to extract the
#' estimated coefficients of a fitted GeDS regression from a \code{\link{GeDS-Class}}
#' object.
#' @param object the  \code{\link{GeDS-class}} object from which the
#' coefficients of the selected GeDS regression should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit whose coefficients should be
#' extracted. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param onlySpline logical variable specifying whether only the coefficients
#' for the GeDS  component of a fitted multivariate regression model should be
#' extracted or if, alternatively, also the coefficients of the parametric component
#' should also be extracted.
#' @param ... potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#' 
#' @return A named vector containing the required coefficients of the fitted
#' univariate or multivariate predictor model. The coefficients corresponding to
#' the variables that enter the parametric component of the fitted multivariate
#' predictor model are named as the variables themselves. The  coefficients of
#' the GeDS component are coded as "\code{N}" followed by the index of the
#' corresponding B-spline.
#'
#' @details
#' These are simple methods for the functions \code{\link[stats]{coef}} and
#' \code{\link[stats]{coefficients}}.
#'
#' As \code{\link{GeDS-class}} objects contain three different fits (linear,
#' quadratic and cubic), it is possible to specify the order of the fit for
#' which GeDS regression coefficients are required via the input argument
#' \code{n}.
#'
#' As mentioned in the details of \code{\link[=formula.GeDS]{formula}}, the
#' predictor model may be multivariate and it may include a (univariate or
#' bivariate) GeD spline component whereas the remaining variables may be part
#' of a parametric component. If the \code{onlySpline} argument is set to
#' \code{TRUE} (the default value), only the coefficients corresponding to the
#' GeD spline component of order \code{n} of the multivariate predictor model 
#' are extracted.
#' 
#' @examples
#' # Generate a data sample for the response variable
#' # and the covariates
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' Z <- runif(N)
#' # Specify a model for the mean of the response Y to be a superposition of
#' # a non-linear component f_1(X), a linear component 2*Z and a
#' # free term 1, i.e.
#' means <- f_1(X) + 2*Z + 1
#' # Add normal noise to the mean of y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit to this sample a predictor model of the form f(X) + Z, where
#' # f(X) is the GeDS component and Z is the linear (additive) component
#' # see ?formula.GeDS for details
#' (Gmod <- NGeDS(Y ~ f(X) + Z, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Extract the GeD spline regression coefficients
#' coef(Gmod, n = 3)
#'
#' # Extract all the coefficients, including the one for the linear component
#' coef(Gmod, onlySpline = FALSE, n = 3)
#' 
#' @seealso \code{\link[stats]{coef}} for the standard definition;
#' \code{\link{NGeDS}} for examples.
#'  
#' @export
#' 
#' @aliases coef.GeDS
#' @rdname coef

coef.GeDS <- function(object, n = 3L, onlySpline = TRUE, ...)
  {
  
  # CDF GeDS
  if (object$Type == "CDF") {
    coef <- object$coef
    return(coef)
  }
  
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")
  
  # Check if order was wrongly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # 1. Linear
  if(n == 2L) {
    theta <- object$Linear.Fit$Theta
    if (object$Type == "LM - Univ" || object$Type == "GLM - Univ") {
      nth <- length(object$Linear.Fit$Polygon$Kn)
      } else if (object$Type == "LM - Biv" || object$Type == "GLM - Biv") {
        nth <- NCOL(object$Linear.Fit$XBasis) * NCOL(object$Linear.Fit$YBasis)
        }
  # 2. Quadratic
    } else if (n == 3L) {
      theta <- object$Quadratic.Fit$Theta
      if (object$Type == "LM - Univ" || object$Type == "GLM - Univ") {
        nth <- length(object$Quadratic.Fit$Polygon$Kn)
        } else if (object$Type == "LM - Biv" || object$Type == "GLM - Biv") {
          nth <- NCOL(object$Quadratic.Fit$XBasis) * NCOL(object$Quadratic.Fit$YBasis)
          }
  # 3. Cubic
      } else if (n == 4L) {
        theta <- object$Cubic.Fit$Theta
        if (object$Type == "LM - Univ" || object$Type == "GLM - Univ") {
          nth <- length(object$Cubic.Fit$Polygon$Kn)
          } else if (object$Type == "LM - Biv" || object$Type == "GLM - Biv") {
            nth <- NCOL(object$Cubic.Fit$XBasis) * NCOL(object$Cubic.Fit$YBasis)
          }
      }
  
  if(!is.null(object$Args$Z) && !onlySpline){
    znames <- attr(object$terms,"term.labels")[-1]
    names(theta) <- c(paste0("N",1:nth),znames)
  } else {
    theta <- theta[1:nth]
    names(theta) <- paste0("N",1:nth)
  }
  return(theta)
}

################################################################################
################################### CONFINT ####################################
################################################################################
#' @title Confidence Intervals for GeDS Model Coefficients
#' @name confint.GeDS
#' @description
#' Method for \code{\link[stats]{confint.default}} to compute confidence intervals for
#' the coefficients of a fitted GeDS model stored in a \code{\link{GeDS-class}} object.
#'
#' @param object the \code{\link{GeDS-class}} object from which the confidence intervals 
#' for the selected GeDS regression should be extracted.
#' @param parm a specification of which parameters are to be given confidence intervals, 
#' either a vector of numbers or names; defaults to all parameters.
#' @param level the confidence level required (default is 0.95).
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit for which to compute confidence
#' intervals. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param ... additional arguments passed to \code{\link[stats]{confint}}.
#'
#' @return A matrix with columns giving lower and upper confidence limits for
#' each spline coefficient of the selected GeDS model (by default 2.5% and 97.5%).
#'
#' @seealso \code{\link[stats]{confint.default}}, \code{\link{GeDS-class}}
#'
#' @aliases confint.GeDS
#' @rdname confint.GeDS
#' @export

confint.GeDS <- function(object, parm, level = 0.95, n = 3L, ...) {
  # Validate order
  n <- as.integer(n)
  if (!(n %in% 2L:4L)) {
    warning("'n' must be 2, 3, or 4. Defaulting to n = 3.")
    n <- 3L
  }
  
  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "Linear.Fit",
                     "3" = "Quadratic.Fit",
                     "4" = "Cubic.Fit"
  )
  
  # Check model availability
  if (inherits(object, "GeDS")) {
    fit_obj <- object[[fit_name]]
  } else if (inherits(object, "GeDSboost") || inherits(object, "GeDSgam")){
    fit_obj <- object$final[[fit_name]]
  }
  if (is.null(fit_obj) || is.null(fit_obj$temporary)) {
    stop("The requested model (", fit_name, ") is not available in the GeDS object.")
  }
  
  # Call stats::confint on the lm/glm object
  if (is.numeric(fit_obj$Theta)) {
    stats::confint.default(fit_obj$temporary, parm = parm, level = level, ...)
    } else if (fit_obj$Theta == "When using bivariate base-learners, the 'single spline representation' (in pp form or B-spline form) of the boosted fit is not available.") {
      cat("Note:\n")
      cat(fit_obj$Theta, "\n")
      if(!inherits(fit_obj$temporary, "glm")) {
        cat("As a result, the intervals printed below are pointwise confidence intervals for the fitted values (not for the coefficients).\n\n")
        
        ci_mat <- as.matrix(cbind(fit_obj$NCI$Low, fit_obj$NCI$Upp))
        colnames(ci_mat) <- c("2.5 %", "97.5 %")
        print(ci_mat)
      }
      
    }
}

################################################################################
################################### DEVIANCE ###################################
################################################################################
#' @title Deviance method for GeDS, GeDSboost, GeDSgam
#' @name deviance.GeDS
#' @description
#' Method for the function \code{\link[stats]{deviance}} that allows the user to
#' extract  the value of the deviance corresponding to a selected GeDS, GeDSboost
#' or GeDSgam fit from a \code{\link{GeDS-Class}},
#' \code{\link{GeDSboost-Class}} or \code{\link{GeDSgam-Class}} object.
#' @param object the \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-class}} object from which the deviance should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit whose deviance should be
#' extracted. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#' 
#' @return A numeric value corresponding to the  deviance of the selected
#' GeDS/GeDSboost/GeDSgam fit.
#' 
#' @details
#' This is a method for the function \code{\link[stats]{deviance}}. As
#' \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} and
#' \code{\link{GeDSgam-class}} objects contain three different fits (linear,
#' quadratic and cubic), it is possible to specify the order of the GeDS fit
#' for which the deviance is required via the input argument \code{n}.
#' 
#' @seealso \code{\link[stats]{deviance}} for the standard definition;
#' \code{\link{GGeDS}} for examples.
#' 
#' @export
#' 
#' @aliases deviance.GeDS deviance.GeDSboost deviance.GeDSgam

deviance.GeDS <- function(object, n = 3L, ...)
  {
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object' and 'n' arguments will be considered")
  
  # Validate order
  n <- as.integer(n)
  if (!(n %in% 2L:4L)) {
    warning("'n' must be 2, 3, or 4. Defaulting to n = 3.")
    n <- 3L
  }
  
  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "Linear.Fit",
                     "3" = "Quadratic.Fit",
                     "4" = "Cubic.Fit"
  )
  
  # Check model availability
  fit_obj <- object[[fit_name]]
  if (is.null(fit_obj)) {
    stop("The requested model (", fit_name, ") is not available in the GeDS object.")
  }
  
  dev <- as.numeric(fit_obj$RSS)
  return(dev)
}

################################################################################
#################################### FAMILY ####################################
################################################################################
#' @title Extract Family from a GeDS Object
#' @name family.GeDS
#' @description
#' Method for \code{\link[stats]{family}} that returns the error distribution 
#' family used in the fitted GeDS model.
#'
#' @param object A \code{\link{GeDS-class}} object.
#' @param ... further arguments (ignored).
#'
#' @return An object of class \code{\link[stats]{family}} describing the 
#' distribution and link function used in the GeDS fit.
#'
#' @seealso \code{\link[stats]{family}}, \code{\link{GeDS-class}}
#'
#' @aliases family.GeDS
#' @rdname family.GeDS
#' @export

family.GeDS  <- function(object, ...) {
  if (object$Type == "LM - Univ" || object$Type == "LM - Biv") {
    stats::gaussian()
    
  } else if (object$Type == "GLM - Univ" || object$Type == "GLM - Biv") {
    object$Args$family
  }
}

################################################################################
################################### FORMULA ####################################
################################################################################
#' @title Formula for the predictor model
#' @name formula.GeDS
#' @description
#' A description of the structure of the predictor model fitted using
#' \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param x fitted \code{\link{GeDS-class}} object, produced by 
#' \code{\link{NGeDS}} or \code{\link{GGeDS}}, from which the predictor model
#' \code{\link[stats]{formula}} should be extracted.
#' @param ... unused in this case.
#'
#' @details
#' In GeDS GNM (GLM) regression, implemented with \code{\link{NGeDS}} and
#' \code{\link{GGeDS}}, the mean of the response variable, correspondingly
#' transformed through an appropriate link function, is modeled using a
#' potentially multivariate predictor model. The latter comprises two components:
#' a GeD variable-knot spline regression involving up to two of the independent
#' variables, and a parametric component for the remaining independent variables.
#' The formula defines the structure of this potentially multivariate predictor.
#'
#' The formulae that are input in \code{\link{NGeDS}} and \code{\link{GGeDS}}
#' are similar to those input in \code{\link[stats]{lm}} or
#' \code{\link[stats]{glm}} except that the function \code{\link{f}} should be
#' specified in order to identify which of the covariates enter the GeD spline
#' regression part of the predictor model. For example, if the predictor model
#' is univariate and it links the transformed mean of \code{y} to \code{x1},
#' the predictor has only a GeD spline component and the
#' \code{\link[stats]{formula}} should be in the form \code{y ~ f(x1)}.
#' 
#' As noted, there may be additional independent variables, \code{x2},
#' \code{x3}, ... which may enter linearly into the parametric component of the
#' predictor model and not be part of the GeD spline regression component. For
#' example one may use the formula \code{y ~ f(x1) + x2 + x3} which assumes a
#' spline regression only between the transformed mean of \code{y} and \code{x1},
#' while \code{x2} and \code{x3} enter the predictor model linearly.
#'
#' Both function \code{\link{NGeDS}} and function \code{\link{GGeDS}}, generate
#' bivariate GeDS regression models. Therefore, if the functional dependence of
#' the mean of the response variable \code{y} on \code{x1} and \code{x2} needs
#' to be jointly modeled and there are no other covariates, the formula for the
#' corresponding two dimensional predictor model should be specified as
#' \code{y ~ f(x1,x2)}.
#' 
#' Within the argument \code{formula}, similarly as in other \R functions, it is
#' possible to specify one or more offset variables, i.e. known terms with fixed
#' regression coefficients equal to 1. These terms should be identified via the
#' function \code{\link[stats]{offset}}.
#' 
#' For \code{\link{NGeDSboost}} and \code{\link{NGeDSgam}}, a formula of the form 
#' \code{y ~ f(x1) + f(x2, x3) + x4} should be provided, where \code{f()} denotes 
#' GeD spline-based (univariate or bivariate) regression base-learners/smoothing 
#' functions, and \code{x4} is included as a linear term in the predictor model.

#'
#' @export
#' 
#' @aliases formula.GeDS

formula.GeDS <- function(x,...)
{
  formula <- stats::formula(x$Formula)
  if(is.null(formula)) stop("Unable to extract the formula. \n
                            Please re-fit using 'NGeDS', 'GGeDS', 'NGeDSboost or 'NGeDSgam")
  return(formula)
}

################################################################################
##################################### KNOTS ####################################
################################################################################
#' @title Knots method for GeDS, GeDSboost, GeDSgam
#' @name knots.GeDS
#' @description
#' Method for the generic function \code{\link[stats]{knots}} that allows the
#' user to extract the vector of knots of a GeDS, GeDSboost or GeDSgam fit of a
#' specified order contained in a \code{\link{GeDS-class}},
#' \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object,
#' respectively.
#' @param Fn the \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-class}} object from which the vector of knots for the
#' specified GeDS, FGB-GeDS or GAM-GeDS fit should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit whose knots should be
#' extracted. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param options a character string specifying whether "\code{all}" knots,
#' including the left-most and the right-most limits of the interval embedding
#' the observations (the default) or only the "\code{internal}" knots should be
#' extracted.
#' @param ... potentially further arguments (required for compatibility with the
#' definition of the generic function). Currently ignored, but with a warning.
#'
#' @return A vector in which each element represents a knot of the GeDS/FGB-GeDS/GAM-GeDS fit of the required order.
#'
#' @details
#' This is a method for the function \code{\link[stats]{knots}} in the
#' \pkg{stats} package.
#'
#' As \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} and
#' \code{\link{GeDSgam-class}} objects contain three different fits (linear,
#' quadratic and cubic), it is possible to specify the order of the GeDS fit
#' whose knots are required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{knots}} for the definition of the generic function; \code{\link{NGeDS}}, \code{\link{GGeDS}},
#' \code{\link{NGeDSboost}} and \code{\link{NGeDSgam}} for examples.
#'
#' @export
#'
#' @aliases knots.GeDS knots.GeDSboost, knots.GeDSgam
#' @rdname knots

knots.GeDS <- function(Fn, n = 3L, options = c("all","internal"), ...)
  {
  
  # CDF GeDS
  if (Fn$Type == "CDF") {
    kn <- Fn$knots
    return(kn)
  }
  
  # Handle additional arguments
  if(!missing(...)) warning("Arguments other than 'Fn', 'n' and 'options' currenly igored. \n Please check if the input parameters have been correctly specified.")
  
  options <- match.arg(options)
  n <- as.integer(n)
  
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # 1. Linear
  if(n == 2L) {
    kn <- Fn$Linear.IntKnots
  # 2. Quadratic
    } else if (n == 3L) {
      kn <- Fn$Quadratic.IntKnots
  # 3. Cubic
      } else if (n == 4L) {
        kn <- Fn$Cubic.IntKnots
      }
  
  if (options == "all") {
    if(Fn$Type == "LM - Univ" || Fn$Type == "GLM - Univ"){
      kn <- sort(c(rep(Fn$Args$extr,n), kn))
    } else if (Fn$Type =="LM - Biv" || Fn$Type =="GLM - Biv") {
      kn$Xk <- sort(c(rep(Fn$Args$Xextr,n), kn$Xk))
      kn$Yk <- sort(c(rep(Fn$Args$Yextr,n), kn$Yk))
    }
  }
  return(kn)
}

################################################################################
#################################### LOGLIK ####################################
################################################################################
#' @title Extract Log-Likelihood from a GeDS Object
#' @name logLik.GeDS
#' @description
#' Method for \code{\link[stats]{logLik}} that returns the log-likelihood of 
#' the selected GeDS regression model.
#'
#' @param object a \code{\link{GeDS-class}} object.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit whose loglikelihood should be
#' extracted. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param ... additional arguments passed to \code{\link[stats]{logLik}}.
#'
#' @return An object of class \code{\link[stats]{logLik}}.
#'
#' @seealso \code{\link[stats]{logLik}}, \code{\link{GeDS-class}}
#'
#' @aliases logLik.GeDS
#' @rdname logLik.GeDS
#' @export

logLik.GeDS <- function(object, n = 3L, ...) {
  
  # Validate order
  n <- as.integer(n)
  if (!(n %in% 2L:4L)) {
    warning("'n' must be 2, 3, or 4. Defaulting to n = 3.")
    n <- 3L
  }
  
  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "Linear.Fit",
                     "3" = "Quadratic.Fit",
                     "4" = "Cubic.Fit"
  )
  
  # Check model availability
  if (inherits(object, "GeDS")) {
    fit_obj <- object[[fit_name]]
  } else if (inherits(object, "GeDSboost") || inherits(object, "GeDSgam")){
    fit_obj <- object$final[[fit_name]]
  }
  
  if (is.null(fit_obj) || is.null(fit_obj$temporary)) {
    stop("The requested model (", fit_name, ") is not available in the GeDS object.")
  }
  
  # Call stats::logLik on the lm/glm object
  stats::logLik(fit_obj$temporary, ...)
  
}

################################################################################
#################################### PREDICT ###################################
################################################################################
#' @title Predict method for GeDS objects
#' @name predict.GeDS
#' @description
#' This is a user friendly method to compute predictions from GeDS objects.
#' @param object the \code{\link{GeDS-class}} object for which the
#' computation of the predicted values is required.
#' @param newdata an optional \code{data.frame}, \code{list} or
#' \code{environment} containing values of the independent variables for which
#' predicted values of the predictor model (including the GeDS and the
#' parametric components) should be computed. If left empty the values are
#' extracted from the object \code{x} itself.
#' @param type character string indicating the type of prediction required. By
#' default it is equal to \code{"response"}, i.e. the result is on the scale of
#' the response variable. See details for the other options.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS/GeDSboost/GeDSgam fit whose predicted values should be
#' computed. By default equal to \code{3L}; non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the
#' generic function). They are ignored, but with a warning.
#'
#' @return A numeric vector corresponding to the predicted values (if
#' \code{type = "link"} or \code{type = "response"}). If \code{type = "terms"} a
#' numeric matrix with a column per term.
#'
#' @details
#' This is a method for the function \code{\link[stats]{predict}} that allows
#' the user to handle \code{\link{GeDS-Class}} objects.
#'
#' In analogy with the function \code{\link[stats]{predict.glm}} in the
#' \pkg{stats} package, the user can specify the scale on which the predictions
#' should be computed through the argument \code{type}. If the predictions are
#' required to be on the scale of the response variable, the user should set
#' \code{type = "response"}, which is the default. Alternatively if one wants
#' the predictions to be on the predictor scale, it is necessary to set
#' \code{type = "link"}.

#' By specifying \code{type = "terms"}, it is possible to inspect the predicted
#' values separately for each single independent variable which enter either the
#' GeD spline component or the parametric component of the predictor model. In
#' this case the returned result is a matrix whose columns correspond to the
#' terms supplied via \code{newdata} or extracted from the \code{object}.
#'
#' As GeDS objects contain three different fits (linear, quadratic and cubic),
#' it is possible to specify the order for which GeDS predictions are required
#' via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{predict}} for the standard definition;
#' \code{\link{GGeDS}} for examples.
#'
#' @export
#'
#' @aliases predict.GeDS

predict.GeDS <- function(object, newdata,
                         type = c("response", "link", "terms"), n = 3L, ...)
  {
  # Handle additional arguments
  if (!missing(...)) warning("Only 'object', 'newdata, 'type' and 'n' arguments will be considered")
  
  # Check if object is of class "GeDS"
  if (!inherits(object, "GeDS"))
    warning("calling predict.GeDS(<fake-GeDS-object>) ...")
  
  # Extract order and model terms
  n <- as.integer(n)
  mt <- object$terms
  
  # 1. Univariate
  if (object$Type == "LM - Univ" || object$Type == "GLM - Univ") {
    
    # If newdata was not provided
    if (missing(newdata) || is.null(newdata)) {
      X <- object$Args$X
      Z <- object$Args$Z
      offset <- object$Args$offset
      if (is.null(offset)) offset <- rep(0, NROW(X))
    # If newdata was provided
      } else {
        mt <- delete.response(mt)
        newdata <- as.list(newdata)
        newdata$f <- f
        mm <- model.matrix(mt,newdata)
        mf <- model.frame(mt,newdata)
        spec <- attr(mt,"specials")$f
        X <- mf[,spec]
        if(ncol(mm) > ncol(X)) {
          Z <- mf[, -c(spec, attr(mt,"response")), drop = T]
          } else {
            Z <- NULL
          }
      offset <- rep(0, NROW(X))
      if (!is.null(off.num <- attr(mt, "offset")))
        for (i in off.num) offset <- offset + eval(attr(mt, "variables")[[i + 1]], newdata)
      if (!is.null(object$call$offset))
        offset <- offset + eval(object$call$offset, newdata)
      }
    
    # Knots
    kn <- knots(object, n = n, options = "all")
    if (min(X) < min(kn) || max(X) > max(kn)) warning("Input values out of the boundary knots")
    # Design matrix
    basisMatrix <- splineDesign(knots = kn, derivs = rep(0,length(X)), x = X, ord = n, outer.ok = T)
    
    type <- match.arg(type)
    
    # (i) Response or Link
    if (type != "terms") {
      coefs <- coef(object,n=n, onlySpline = FALSE)
      basisMatrix2 <- cbind(basisMatrix,Z)
      predicted <- basisMatrix2 %*% coefs + offset
      
      if(type=="response" & !is.null(object$Args$family)) {
        predicted <- object$Args$family$linkinv(predicted)
      }
      
    # (ii) Terms
      } else {
        coefs <- coef(object, n = n, onlySpline = TRUE)
        coefs1 <- coef(object, n = n, onlySpline = FALSE)
        predicted <- basisMatrix %*% coefs
        colnames(predicted) <- "Spline"
        predicted1 <- if(!is.null(Z)) {
          Z*matrix(coefs1[-c(1:length(coefs))],
                   ncol = length(coefs1[-c(1:length(coefs))]), nrow=NROW(Z))
          } else NULL
        
        if(!is.null(predicted1)) {
          predicted1 <- as.matrix(predicted1)
          colnames(predicted1) <- object$znames
        }
        
    predicted <- cbind(predicted, predicted1)
      }
    
  # 2. Bivariate
  } else if (object$Type == "LM - Biv" || object$Type == "GLM - Biv") {
    
    # If newdata was not provided
    if (missing(newdata) || is.null(newdata)) {
      X <- object$Args$X
      Y <- object$Args$Y
      W <- object$Args$W
      
    # If newdata was provided
    } else {
      mt <- delete.response(mt)
      newdata <- as.list(newdata)
      newdata$f <- f
      mm <- model.matrix(mt, newdata)
      mf <- model.frame(mt, newdata)
      spec <- attr(mt,"specials")$f
      X <- mf[,spec][, 1]
      Y <- mf[,spec][, 2]
      if(ncol(mm) > ncol(mf[,spec])) {
        W <- mf[, -c(spec, attr(mt, "response")), drop = T]
      } else {
        W <- NULL
      }
    }
    
    # Knots
    kn <- knots(object, n = n, options = "all")
    if(min(X) < min(kn$Xk) || max(X) > max(kn$Xk) || min(Y) < min(kn$Yk) | max(Y) > max(kn$Yk))
      warning("Input values out of the boundary knots")
    # Design matrix
    basisMatrixX <- splineDesign(knots = kn$Xk, derivs = rep(0, length(X)), x = X, ord = n, outer.ok = T)
    basisMatrixY <- splineDesign(knots = kn$Yk, derivs = rep(0, length(Y)), x = Y, ord = n, outer.ok = T)

    basisMatrixbiv <- tensorProd(basisMatrixX, basisMatrixY)

    type <- match.arg(type)
    
    # (i) Response or Link
    if (type != "terms") {
      coefs <- coef(object, n = n, onlySpline = FALSE)
      coefs[is.na(coefs)] <- 0
      basisMatrixbiv2 <- cbind(basisMatrixbiv,W)
      predicted <- basisMatrixbiv2 %*% coefs
      
      if(type == "response" & !is.null(object$Args$family)) {
        predicted <- object$Args$family$linkinv(predicted)
      }
      
    # (ii) Terms
    } else {
      coefs <- coef(object, n = n, onlySpline = TRUE)
      coefs[is.na(coefs)] <- 0
      coefs1 <- coef(object, n = n, onlySpline = FALSE)
      coefs1[is.na(coefs1)] <- 0
      
      # Spline prediction
      predicted <- basisMatrixbiv %*% coefs
      colnames(predicted) <- "Spline"
      # Linear prediction
      predicted1 <- if(!is.null(W)) {
        W*matrix(coefs1[-c(1:length(coefs))],
                 ncol = length(coefs1[-c(1:length(coefs))]), nrow=NROW(W))
        } else NULL
      
      if(!is.null(predicted1)) {
        predicted1 <- as.matrix(predicted1)
        colnames(predicted1) <- object$znames
      }
      
      predicted <- cbind(predicted, predicted1)
    }
    
  }
  return(as.numeric(predicted))
}

################################################################################
##################################### PRINT ####################################
################################################################################
#' @title Print method for GeDS, GeDSboost, GeDSgam
#' @name print.GeDS
#' @description
#' Method for the generic function \code{\link[base]{print}} that allows to
#' print on screen the main information related to a fitted \code{\link{GeDS-class}},
#' \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} model.
#' @param x the \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-class}} object for which the main information should be
#' printed on screen.
#' @param digits number of digits to be printed.
#' @param ... potentially further arguments (required by the definition of the
#' generic function).
#'
#' @return This function returns (invisibly) the same input object, but adding
#' the slot \code{Print} that contains the three sub-slots:
#' \item{Nknots}{ the number of internal knots of the linear
#' GeDS/FGB-GeDS/GAM-GeDS fit}
#' \item{Deviances}{ the deviances of the three (linear, quadratic and cubic)
#' GeDS/FGB-GeDS/GAM-GeDS fits}
#' \item{Call}{ the \code{call} to the function that produced the \code{x}
#' object}
#' 
#' @details
#' This method allows to print on screen basic information related to the fitted
#' predictor model such as the function \code{call}, the number of internal
#' knots for the linear GeDS/FGB-GeDS/GAM-GeDS fit and the deviances for the
#' three (linear, quadratic and cubic) fitted predictor models embedded in the
#' \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object.
#'
#' @seealso \code{\link[base]{print}} for the standard definition.
#' 
#' @export
#'
#' @aliases print.GeDS print.GeDSboost print.GeDSgam

print.GeDS <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  kn <- knots(x, n=2, options="int")
  names(kn) <- NULL
  DEVs <- numeric(3)
  names(DEVs) <- c("Order 2","Order 3","Order 4")
  
  # 1) Univariate
  if(x$Type == "LM - Univ" || x$Type == "GLM - Univ") {
    
    if (length(kn)) {
      
      cat(paste0("Number of internal knots of the second order (linear) spline: ", length(kn)))
      cat("\n")
      
      DEVs[1] <- x$Dev.Linear
      DEVs[2] <- if(!is.null(x$Dev.Quadratic)) x$Dev.Quadratic else NA
      DEVs[3] <- if(!is.null(x$Dev.Cubic)) x$Dev.Cubic else NA
      cat("Deviances:\n")
      
      print.default(format(DEVs, digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
    else {
      cat("No internal knots found\n")
    }
    cat("\n")
    print <- list("Nknots" = length(kn), "Deviances" = DEVs, "Call" = x$extcall)
    
    
  # 2) Bivariate
  } else if (x$Type == "LM - Biv" || x$Type == "GLM - Biv") {
    
    if(length(kn[[1]])||length(kn[[2]])) {
      
      cat(paste0("Number of internal knots of the second order (linear) spline in the X direction: ", length(kn[[1]])))
      cat("\n")
      cat(paste0("Number of internal knots of the second order (linear) spline in the Y direction: ", length(kn[[2]])))
      cat("\n")
      cat("\n")
      
      DEVs[1] <- x$Dev.Linear
      DEVs[2] <- if(!is.null(x$Dev.Quadratic)) x$Dev.Quadratic else NA
      DEVs[3] <- if(!is.null(x$Dev.Cubic)) x$Dev.Cubic else NA
      cat("\nDeviances:\n")
      
      print.default(format(DEVs, digits = digits), print.gap = 2L,
                    quote = FALSE)
      
    } else {
      cat("No internal knots found\n")
      }
    cat("\n")
    print <- list("Nknots" = c(length(kn[[1]]),length(kn[[2]])), "Deviances" = DEVs, "Call" = x$extcall)
  }
  
  x$Print <- print
  invisible(x)
}

################################################################################
################################### SUMMARY ####################################
################################################################################
#' @title Summary method for GeDS, GeDSboost, GeDSgam
#' @name summary.GeDS
#' @description
#' Method for the generic function \code{\link[base]{summary}} that allows you to
#' print on screen the main information related to a fitted \code{\link{GeDS-class}},
#' \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} model.
#' Similar to \code{\link{print.GeDS}} but with some extra detail.
#' @param object the \code{\link{GeDS-class}}, \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-class}} object for which the main information should be
#' printed on screen.
#' @param ... potentially further arguments (required by the definition of the
#' generic function).
#' @seealso \code{\link{print.GeDS}}
#' @export
#' @aliases summary.GeDS summary.GeDSboost summary.GeDSgam

summary.GeDS <- function(object, ...)
{
  
  # 1) Header
  cat("\nGeometrically Desined (GeD) Spine Regression Model\n\n")
    
  # 2) Core print
  print.GeDS(object, ...)
  
  # 3) Additional details
  if (object$Type == "LM - Univ" || object$Type == "LM - Biv") {
    print(gaussian())
    
  } else if (object$Type == "GLM - Univ" || object$Type == "GLM - Biv") {
    print(object$Args$family)
  }
  
  cat("phi = ", object$Args$phi, "and q = ", object$Args$q, "(stopping rule parameters);\n")
  cat("beta = ", object$Args$beta)
  
  if (object$Args$beta == 0.5) {
    cat(", meaning that the within-cluster mean residual and the cluster range were considered equally important when placing the knots.\n")
  } else if (object$Args$beta > 0.5) {
    cat(", meaning that more weight was given to the within-cluster mean residual than to the cluster range when placing the knots.\n")
  } else if (object$Args$beta < 0.5) {
    cat(", meaning that more weight was given to the cluster range than to the within-cluster mean residual when placing the knots.\n")
  }
  
  
  invisible(object)
}


