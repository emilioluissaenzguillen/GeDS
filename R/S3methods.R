################################################################################
##################################### COEF #####################################
################################################################################
#' @title Coef Method for GeDS and GeDSfitND Objects
#' @name coef.GeDS
#' @description
#' Method for the function \code{\link[stats]{coef}} that allows to extract the
#' estimated coefficients of a fitted GeDS regression model from a \code{"GeDS"} class
#' object.
#'
#' @param object The \code{"GeDS"} or \code{"GeDSfitND"} class object from which the
#' coefficients of the selected GeDS regression model should be extracted.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"}
#' fit whose coefficients should be extracted. By default equal to \code{3L};
#' non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param onlySpline Logical variable specifying whether only the coefficients
#' for the GeDS  component of a fitted multivariate regression model should be
#' extracted or whether the coefficients of both the GeDS and the parametric
#' components should be returned.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#'
#' @details
#' Simple method for the function \code{\link[stats]{coef}}.
#'
#' As \code{"GeDS"} class objects contain three different fits (linear,
#' quadratic and cubic), the argument \code{n} can be used to specify the order
#' of the GeDS fit for which regression coefficients are required.
#'
#' As mentioned in the Details of \code{\link[=formula.GeDS]{formula}}, the
#' predictor model may be multivariate and it may include a (univariate or
#' bivariate) GeD spline component, plus a parametric component involving the
#' remaining variables. If the \code{onlySpline} argument is set to
#' \code{TRUE} (the default value), only the coefficients corresponding to the
#' GeD spline component of order \code{n} of the multivariate predictor model
#' are extracted.
#'
#' @return A named vector containing the required coefficients of the fitted
#' univariate or multivariate predictor model. The coefficients corresponding to
#' the variables that enter the parametric component of the fitted multivariate
#' predictor model are named as the variables themselves. The  coefficients of
#' the GeDS component are coded as "\code{N}" followed by the index of the
#' corresponding B-spline. For \code{"GeDSfitND"} objects, coefficient names
#' identify the basis index in every tensor dimension and the same indices are
#' available as a data frame in the \code{"basis.index"} attribute.
#'
#' @seealso \code{\link[stats]{coef}} for the standard definition;
#' \code{\link{NGeDS}} for more examples.
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
#' # Add normal noise to the mean of Y
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
#' @rdname coef
#' @aliases coef.GeDS coef.GeDSfitND
#' @export

coef.GeDS <- function(object, n = 3L, onlySpline = TRUE, ...)
  {

  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")

  # Validate order
  n <- validate_GeDS_order(n)

  # 1. Linear
  if(n == 2L) {
    theta <- object$linear.fit$theta
    if (object$type == "LM - Univ" || object$type == "GLM - Univ") {
      nth <- length(object$linear.fit$polygon$kn)
      } else if (object$type == "LM - Biv" || object$type == "GLM - Biv") {
        nth <- NCOL(object$linear.fit$Xbasis) * NCOL(object$linear.fit$Ybasis)
        }
  # 2. Quadratic
    } else if (n == 3L) {
      theta <- object$quadratic.fit$theta
      if (object$type == "LM - Univ" || object$type == "GLM - Univ") {
        nth <- length(object$quadratic.fit$polygon$kn)
        } else if (object$type == "LM - Biv" || object$type == "GLM - Biv") {
          nth <- NCOL(object$quadratic.fit$Xbasis) * NCOL(object$quadratic.fit$Ybasis)
          }
  # 3. Cubic
      } else if (n == 4L) {
        theta <- object$cubic.fit$theta
        if (object$type == "LM - Univ" || object$type == "GLM - Univ") {
          nth <- length(object$cubic.fit$polygon$kn)
          } else if (object$type == "LM - Biv" || object$type == "GLM - Biv") {
            nth <- NCOL(object$cubic.fit$Xbasis) * NCOL(object$cubic.fit$Ybasis)
          }
      }

  if(!is.null(object$args$Z) && !onlySpline){
    znames <- attr(object$terms,"term.labels")[-1]
    names(theta) <- c(paste0("N",1:nth),znames)
  } else {
    theta <- theta[1:nth]
    names(theta) <- paste0("N",1:nth)
  }
  return(theta)
}

#' @rdname coef
#' @export
coef.GeDSfitND <- function(object, n = 3L, onlySpline = TRUE, ...)
{
  if (!missing(...)) {
    warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")
  }
  if (length(onlySpline) != 1L || is.na(onlySpline) || !is.logical(onlySpline)) {
    stop("'onlySpline' must be TRUE or FALSE.", call. = FALSE)
  }

  n <- validate_GeDS_order(n)
  fit.name <- switch(
    as.character(n),
    "2" = "linear.fit",
    "3" = "quadratic.fit",
    "4" = "cubic.fit"
  )
  fit <- object[[fit.name]]
  if (is.null(fit)) {
    stop("The requested model (", fit.name,
         ") is not available in the GeDS object.")
  }

  basis.counts <- vapply(fit$basis.matrices, NCOL, integer(1))
  dimensions <- names(basis.counts)
  reversed.grid <- do.call(
    expand.grid,
    c(
      setNames(lapply(rev(basis.counts), seq_len), rev(dimensions)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  basis.index <- reversed.grid[, rev(seq_along(dimensions)), drop = FALSE]
  names(basis.index) <- dimensions

  coefficients <- fit$coefficients
  names(coefficients) <- apply(
    basis.index,
    1L,
    function(index) {
      paste0("B[", paste0(dimensions, "=", index, collapse = ","), "]")
    }
  )
  attr(coefficients, "basis.index") <- basis.index
  coefficients
}

################################################################################
################################### CONFINT ####################################
################################################################################
#' @title Confidence Intervals for GeDS Models Coefficients
#' @name confint.GeDS
#' @description
#' Method for \code{\link[stats]{confint.default}} to compute confidence intervals for
#' the coefficients of a fitted GeDS model stored in a \code{"GeDS"}, \code{"GeDSgam"}
#' or \code{"GeDSboost"} class object.
#'
#' @param object The \code{"GeDS"}/\code{"GeDSgam"}/\code{"GeDSboost"} class object
#' from which the confidence intervals for the selected order \code{n} should be extracted.
#' @param parm A specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or names; defaults to all parameters.
#' @param level The confidence level required (default is 0.95).
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} fit
#' for which to compute confidence intervals. By default equal to \code{3L};
#' non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... Additional arguments passed to \code{\link[stats]{confint.default}}.
#'
#' @return A matrix with columns giving lower and upper confidence limits for
#' each spline coefficient of the selected GeDS model (by default 2.5\% and 97.5\%).
#'
#' @seealso \code{\link[stats]{confint.default}}, \code{\link{NGeDS}}, \code{\link{GGeDS}},
#' \code{\link{NGeDSgam}}, \code{\link{NGeDSboost}}
#'
#' @rdname confint.GeDS
#' @aliases confint.GeDS
#' @importFrom stats confint.default
#' @export

confint.GeDS <- function(object, parm, level = 0.95, n = 3L, ...) {
  # Validate order
  n <- validate_GeDS_order(n)

  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "linear.fit",
                     "3" = "quadratic.fit",
                     "4" = "cubic.fit"
  )

  # Check model availability
  if (inherits(object, "GeDS")) {
    fit_obj <- object[[fit_name]]
    # names
    if (object$type == "LM - Univ" || object$type == "GLM - Univ") {
      nth <- NCOL(fit_obj$basis)
    } else if (object$type == "LM - Biv" || object$type == "GLM - Biv") {
      nth <- NCOL(fit_obj$Xbasis) * NCOL(fit_obj$Ybasis)
    }
    if (!is.null(object$args$Z)) {
      znames <- attr(object$terms, "term.labels")[-1]
      names <- c(paste0("N", 1:nth), znames)
    } else {
      names <- paste0("N", 1:nth)
    }

  } else if (inherits(object, "GeDSboost") || inherits(object, "GeDSgam")){
    fit_obj <- object$final_model[[fit_name]]
    names <- names(fit_obj$theta)
  }
  if (is.null(fit_obj) || is.null(fit_obj$temporary)) {
    stop("The requested model (", fit_name, ") is not available in the GeDS object.")
  }

  # Call stats::confint on the lm/glm object
  if (is.numeric(fit_obj$theta)) {
    ci <- confint.default(fit_obj$temporary, parm = parm, level = level, ...)
    rownames(ci) <- names

    ci

    } else if (fit_obj$theta == "When using bivariate base-learners, the 'single spline representation' (in pp form or B-spline form) of the boosted fit is not available.") {
      cat("Note:\n")
      cat(fit_obj$theta, "\n")
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
#' @title Deviance Method for GeDS, GeDSfitND, GeDSgam, GeDSboost
#' @name deviance.GeDS
#' @description
#' Method for the function \code{\link[stats]{deviance}} that allows the user to
#' extract the value of the deviance corresponding to a selected GeDS, GeDSboost
#' or GeDSgam fit typically returned by \code{\link{NGeDS}}/\code{\link{GGeDS}},
#' \code{\link{NGeDSgam}} or \code{\link{NGeDSboost}}.
#'
#' @param object The \code{"GeDS"}, \code{"GeDSfitND"}, \code{"GeDSgam"} or
#' \code{"GeDSboost"} class object from which the deviance should be extracted.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} fit
#' whose deviance should be extracted. By default equal to \code{3L}; non-integer
#' values will be passed to the function \code{\link{as.integer}}.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#'
#' @return A numeric value corresponding to the  deviance of the selected
#' \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} fit.
#'
#' @details
#' This is a method for the function \code{\link[stats]{deviance}} in the
#' \pkg{stats} package. As \code{"GeDS"}, \code{"GeDSgam"} and \code{"GeDSboost"}
#' class objects contain three different fits (linear, quadratic and cubic), it
#' is possible to specify the order of the GeDS fit for which the deviance is
#' required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{deviance}} for the standard definition;
#' \code{\link{NGeDS}}, \code{\link{GGeDS}}, \code{\link{NGeDSgam}},
#' \code{\link{NGeDSboost}} for examples.
#'
#' @aliases deviance.GeDS deviance.GeDSfitND deviance.GeDSboost deviance.GeDSgam
#' @importFrom stats deviance
#' @export

deviance.GeDS <- function(object, n = 3L, ...)
  {
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object' and 'n' arguments will be considered")

  # Validate order
  n <- validate_GeDS_order(n)

  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "linear.fit",
                     "3" = "quadratic.fit",
                     "4" = "cubic.fit"
  )

  # Check model availability
  fit_obj <- object[[fit_name]]
  if (is.null(fit_obj)) {
    stop("The requested model (", fit_name, ") is not available in the GeDS object.")
  }

  dev <- as.numeric(fit_obj$rss)
  return(dev)
}

################################################################################
#################################### FAMILY ####################################
################################################################################
#' @title Extract Family from a GeDS, GeDSgam, GeDSboost Object
#' @name family.GeDS
#' @description
#' Method for \code{\link[stats]{family}} that returns the error distribution
#' family used in the fitted GeDS model.
#'
#' @param object A \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} class object.
#' @param ... Further arguments (ignored).
#'
#' @return An object of class \code{\link[stats]{family}} describing the
#' distribution and link function used in the GeDS fit.
#'
#' @seealso \code{\link[stats]{family}}
#'
#' @rdname family.GeDS
#' @aliases family.GeDS
#' @importFrom stats family gaussian
#' @export

family.GeDS  <- function(object, ...) {
  if (object$type == "LM - Univ" || object$type == "LM - Biv") {
    gaussian()

  } else if (object$type == "GLM - Univ" || object$type == "GLM - Biv") {
    object$args$family
  }
}

################################################################################
################################### formula ####################################
################################################################################
#' @title Formula for the Predictor Model
#' @name formula.GeDS
#' @description
#' A description of the structure of the predictor model fitted using
#' \code{\link{NGeDS}}, \code{\link{GGeDS}}, \code{\link{NGeDSgam}} or
#' \code{\link{NGeDSboost}}.
#' @param x Fitted \code{"GeDS"}, \code{"GeDSfitND"}, \code{"GeDSgam"} or
#' \code{"GeDSboost"} class
#' object, produced by \code{\link{NGeDS}}, \code{\link{GGeDS}},
#' \code{\link{NGeDSgam}} or \code{\link{NGeDSboost}} from which the predictor model
#' \code{\link[stats]{formula}} should be extracted.
#' @param ... Unused in this case.
#'
#' @details
#' In GeDS GNM (GLM) regression (implemented through \code{\link{NGeDS}} and
#' \code{\link{GGeDS}}) the mean of the response variable, correspondingly
#' transformed through an appropriate link function, is modeled using a
#' potentially multivariate predictor model. Established \code{NGeDS} and
#' \code{GGeDS} fits allow one or two covariates in a GeD spline component.
#' Experimental Normal \code{NGeDS} fits allow more than two covariates in one
#' joint tensor-product spline. \code{GGeDS} remains limited to at most two.
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
#' As noted, there may be additional independent variables \code{x2},
#' \code{x3}, ... which may enter linearly into the parametric component of the
#' predictor model and not be part of the GeD spline regression component. For
#' example one may use the formula \code{y ~ f(x1) + x2 + x3} which assumes a
#' spline regression only between the transformed mean of \code{y} and \code{x1},
#' while \code{x2} and \code{x3} enter the predictor model linearly.
#'
#' Both \code{\link{NGeDS}} and \code{\link{GGeDS}} functions, generate
#' bivariate GeDS regression models. Therefore, if the functional dependence of
#' the mean of the response variable \code{y} on \code{x1} and \code{x2} needs
#' to be jointly modeled and there are no other covariates, the formula for the
#' corresponding two dimensional predictor model should be specified as
#' \code{y ~ f(x1,x2)}.
#'
#' For an experimental joint Normal fit involving three or more covariates,
#' \code{NGeDS} accepts, for example, \code{y ~ f(x1,x2,x3)}. This route does
#' not yet support parametric terms, offsets, initial knots or fits, or the LR
#' stopping rule.
#'
#' Within the argument \code{formula}, similarly as in other \R functions, it is
#' possible to specify one or more offset variables, i.e., known terms with fixed
#' regression coefficients equal to 1. These terms should be identified via the
#' function \code{\link[stats]{offset}}.
#'
#' For \code{\link{NGeDSgam}} and \code{\link{NGeDSboost}}, more than one GeD spline
#' component can be included in the formula, e.g., \code{y ~ f(x1) + f(x2,x3) + x4},
#' where \code{f()} denotes GeD spline-based (univariate or bivariate) regression smoothing
#' functions/base-learners, and \code{x4} is included as a linear term in the
#' predictor model. Offset terms are not supported by \code{\link{NGeDSboost}} and
#' will be ignored if included in the formula. Known additive components can
#' instead be manually incorporated into the response variable prior to fitting the model.
#'
#'
#' @aliases formula.GeDS
#' @importFrom stats formula
#' @export

formula.GeDS <- function(x, ...)
{
  formula <- formula(x$formula)
  if(is.null(formula)) stop("Unable to extract the formula. \n
                            Please re-fit using 'NGeDS', 'GGeDS', 'NGeDSboost or 'NGeDSgam")
  return(formula)
}

#' @rdname deviance.GeDS
#' @export
deviance.GeDSfitND <- deviance.GeDS

#' @rdname formula.GeDS
#' @export
formula.GeDSfitND <- formula.GeDS

################################################################################
##################################### KNOTS ####################################
################################################################################
#' @title Knots Method for GeDS, GeDSfitND, GeDSgam, GeDSboost
#' @name knots.GeDS
#' @description
#' Method for the generic function \code{\link[stats]{knots}} that allows the
#' user to extract the vector of knots of a GeDS, GAM-GeDS or FGB-GeDS fit of a
#' specified order contained in a \code{"GeDS"}, \code{"GeDSgam"} or
#' \code{"GeDSboost"} class, respectively.
#' @param Fn The \code{"GeDS"}, \code{"GeDSfitND"}, \code{"GeDSgam"} or
#' \code{"GeDSboost"} class object from which the vector of knots for the
#' specified GeDS, GAM-GeDS or FGB-GeDS fit should be extracted.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} fit
#' whose knots should be extracted. By default equal to \code{3L}; non-integer
#' values will be passed to the function \code{\link{as.integer}}.
#' @param options A character string specifying whether "\code{all}" knots,
#' including the left-most and the right-most limits of the interval embedding
#' the observations (the default) or only the "\code{internal}" knots should be
#' extracted.
#' @param ... Potentially further arguments (required for compatibility with the
#' definition of the generic function). Currently ignored, but with a warning.
#'
#' @return For univariate fits, a vector in which each element represents a
#' knot. For bivariate and multivariate fits, a named list containing one knot
#' vector per spline coordinate.
#'
#' @details
#' This is a method for the function \code{\link[stats]{knots}} in the
#' \pkg{stats} package.
#'
#' As \code{"GeDS"} class, \code{\link{NGeDSgam}} and
#' \code{\link{NGeDSboost}} objects contain three different fits (linear,
#' quadratic and cubic), it is possible to specify the order of the GeDS fit
#' whose knots are required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{knots}} for the definition of the generic function; \code{\link{NGeDS}}, \code{\link{GGeDS}},
#' \code{\link{NGeDSboost}} and \code{\link{NGeDSgam}} for examples.
#'
#' @rdname knots
#' @aliases knots.GeDS knots.GeDSfitND knots.GeDSgam knots.GeDSboost
#' @importFrom stats knots
#' @export

knots.GeDS <- function(Fn, n = 3L, options = c("all","internal"), ...)
  {


  # Handle additional arguments
  if(!missing(...)) warning("Arguments other than 'Fn', 'n' and 'options' currenly igored. \n Please check if the input parameters have been correctly specified.")

  options <- match.arg(options)

  # Validate order
  n <- validate_GeDS_order(n)

  # 1. Linear
  if(n == 2L) {
    kn <- Fn$linear.intknots
  # 2. Quadratic
    } else if (n == 3L) {
      kn <- Fn$quadratic.intknots
  # 3. Cubic
      } else if (n == 4L) {
        kn <- Fn$cubic.intknots
      }

  if (options == "all") {
    if(Fn$type == "LM - Univ" || Fn$type == "GLM - Univ"){
      kn <- sort(c(rep(Fn$args$extr,n), kn))
    } else if (Fn$type =="LM - Biv" || Fn$type =="GLM - Biv") {
      kn$Xk <- sort(c(rep(Fn$args$Xextr,n), kn$Xk))
      kn$Yk <- sort(c(rep(Fn$args$Yextr,n), kn$Yk))
    }
  }
  return(kn)
}

#' @rdname knots
#' @export
knots.GeDSfitND <- function(Fn, n = 3L,
                            options = c("all", "internal"), ...)
{
  if (!missing(...)) {
    warning("Only 'Fn', 'n' and 'options' arguments will be considered")
  }

  n <- validate_GeDS_order(n)
  options <- match.arg(options)
  fit.name <- switch(
    as.character(n),
    "2" = "linear.fit",
    "3" = "quadratic.fit",
    "4" = "cubic.fit"
  )
  fit <- Fn[[fit.name]]
  if (is.null(fit)) {
    stop("The requested model (", fit.name,
         ") is not available in the GeDS object.")
  }

  if (options == "internal") fit$intknots else fit$full.knots
}

################################################################################
#################################### LOGLIK ####################################
################################################################################
#' @title Extract Log-Likelihood from a GeDS Object
#' @name logLik.GeDS
#' @description
#' Method for \code{\link[stats]{logLik}} that returns the log-likelihood of
#' the selected GeDS, GeDS-GAM or FGB-GeDS model.
#'
#' @param object A \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} class object.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} fit
#' whose loglikelihood should be extracted. By default equal to \code{3L};
#' non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... Additional arguments passed to \code{\link[stats]{logLik}}.
#'
#' @return An object of class \code{\link[stats]{logLik}}.
#'
#' @seealso \code{\link[stats]{logLik}}
#'
#' @rdname logLik.GeDS
#' @aliases logLik.GeDS
#' @importFrom stats logLik
#' @export

logLik.GeDS <- function(object, n = 3L, ...) {

  # Validate order
  n <- validate_GeDS_order(n)

  # Select model based on order
  fit_name <- switch(as.character(n),
                     "2" = "linear.fit",
                     "3" = "quadratic.fit",
                     "4" = "cubic.fit"
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
  logLik(fit_obj$temporary, ...)

}

################################################################################
#################################### PREDICT ###################################
################################################################################
#' @title Predict Method for GeDS Objects
#' @name predict.GeDS
#' @description
#' This is a user friendly method to compute predictions from GeDS objects.
#'
#' @param object The \code{"GeDS"} class object for which the
#' computation of the predicted values is required.
#' @param newdata An optional \code{data.frame}, \code{list} or
#' \code{environment} containing values of the independent variables for which
#' predicted values of the predictor model (including the GeDS and the
#' parametric components) should be computed. If left empty the values are
#' extracted from the object \code{x} itself.
#' @param type Character string specifying the type of prediction to return. The
#' default is \code{"response"}, which gives predictions on the scale of the
#' response variable. See Details for other available options.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"}
#' fit whose predicted values should be computed. By default equal to \code{3L};
#' non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). They are ignored, but with a warning.
#'
#' @details
#' This is a method for the function \code{\link[stats]{predict}} in the
#' \pkg{stats} package, that allows the user to handle \code{"GeDS"} class objects.
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
#' @return A numeric vector corresponding to the predicted values (if
#' \code{type = "link"} or \code{type = "response"}). If \code{type = "terms"}, a
#' numeric matrix with a column per term.
#'
#' @seealso \code{\link[stats]{predict}} for the standard definition;
#' \code{\link{GGeDS}} for examples.
#'
#' @aliases predict.GeDS
#' @importFrom splines splineDesign
#' @importFrom stats delete.response model.frame model.matrix
#' @export

predict.GeDS <- function(object, newdata,
                         type = c("response", "link", "terms"), n = 3L, ...)
  {
  # Handle additional arguments
  if (!missing(...)) warning("Only 'object', 'newdata, 'type' and 'n' arguments will be considered")

  # Check if object is of class "GeDS"
  if (!inherits(object, "GeDS"))
    warning("calling predict.GeDS(<fake-GeDS-object>) ...")

  # Validate order
  n <- validate_GeDS_order(n)

  # Extract order and model terms
  n <- as.integer(n)
  mt <- object$terms

  # 1. Univariate
  if (object$type == "LM - Univ" || object$type == "GLM - Univ") {

    # If newdata was not provided
    if (missing(newdata) || is.null(newdata)) {
      X <- object$args$X
      Z <- object$args$Z
      offset <- object$args$offset
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

      if(type=="response" & !is.null(object$args$family)) {
        predicted <- object$args$family$linkinv(predicted)
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
  } else if (object$type == "LM - Biv" || object$type == "GLM - Biv") {

    # If newdata was not provided
    if (missing(newdata) || is.null(newdata)) {
      X <- object$args$X
      Y <- object$args$Y
      W <- object$args$W

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

      if(type == "response" & !is.null(object$args$family)) {
        predicted <- object$args$family$linkinv(predicted)
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
##################################### print ####################################
################################################################################
#' @title Print Method for GeDS, GeDSgam, GeDSboost
#' @name print.GeDS
#' @description
#' Method for the generic function \code{\link[base]{print}} that allows to
#' print on screen the main information related to a fitted \code{"GeDS"},
#' \code{"GeDSgam"}, \code{"GeDSboost"} class model.
#'
#' @param x The \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} class
#' object for which the main information should be printed on screen.
#' @param digits Number of digits to be printed.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function).
#'
#' @details
#' This method allows to print on screen basic information related to the fitted
#' predictor model such as the function \code{call}, the number of internal
#' knots for the linear GeDS/FGB-GeDS/GAM-GeDS fit and the deviances for the
#' three (linear, quadratic and cubic) fitted predictor models embedded in the
#' \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} object.
#'
#' @return This function returns (invisibly) the same input object, but adding
#' the slot \code{print} that contains the three sub-slots:
#' \describe{
#' \item{Nintknots}{ the number of internal knots of the linear
#' GeDS/FGB-GeDS/GAM-GeDS fit}
#' \item{deviances}{ the deviances of the three (linear, quadratic and cubic)
#' GeDS/FGB-GeDS/GAM-GeDS fits}
#' \item{call}{ the \code{call} to the function that produced the \code{x}
#' object}
#' }
#'
#' @seealso \code{\link[base]{print}} for the standard definition.
#'
#' @aliases print.GeDS print.GeDSboost print.GeDSgam
#'
#' @export

print.GeDS <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  int.knt <- knots(x, n = 2, options = "int")
  names(int.knt) <- NULL
  devs <- numeric(3)
  names(devs) <- c("Order 2","Order 3","Order 4")

  # 1) Univariate
  if(x$type == "LM - Univ" || x$type == "GLM - Univ") {

    if (length(int.knt)) {

      cat(paste0("Number of internal knots of the second order (linear) spline: ", length(int.knt)))
      cat("\n")

      devs[1] <- x$dev.linear
      devs[2] <- if(!is.null(x$dev.quadratic)) x$dev.quadratic else NA
      devs[3] <- if(!is.null(x$dev.cubic)) x$dev.cubic else NA
      cat("Deviances:\n")

      print.default(format(devs, digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
    else {
      cat("No internal knots found\n")
    }
    cat("\n")
    print <- list("Nintknots" = length(int.knt), "Deviances" = devs, "Call" = x$extcall)


  # 2) Bivariate
  } else if (x$type == "LM - Biv" || x$type == "GLM - Biv") {

    if(length(int.knt[[1]])||length(int.knt[[2]])) {

      cat(paste0("Number of internal knots of the second order (linear) spline in the X direction: ", length(int.knt[[1]])))
      cat("\n")
      cat(paste0("Number of internal knots of the second order (linear) spline in the Y direction: ", length(int.knt[[2]])))
      cat("\n")
      cat("\n")

      devs[1] <- x$dev.linear
      devs[2] <- if(!is.null(x$dev.quadratic)) x$dev.quadratic else NA
      devs[3] <- if(!is.null(x$dev.cubic)) x$dev.cubic else NA
      cat("\nDeviances:\n")

      print.default(format(devs, digits = digits), print.gap = 2L,
                    quote = FALSE)

    } else {
      cat("No internal knots found\n")
      }
    cat("\n")
    print <- list("Nintknots" = c(length(int.knt[[1]]),length(int.knt[[2]])), "deviances" = devs, "call" = x$extcall)
  }

  x$print <- print
  invisible(x)
}

#' Print an Experimental Multivariate GeDS Fit
#'
#' Prints the call, dimensions, selected knot counts, residual sums of squares,
#' and least-squares solver used for each available spline order.
#'
#' @param x A \code{"GeDSfitND"} object returned by \code{NGeDS()}.
#' @param digits Number of significant digits used for the residual sums of
#'   squares.
#' @param ... Further arguments, currently ignored.
#'
#' @return The input object, invisibly.
#' @seealso \code{\link{NGeDS}}, \code{\link{summary.GeDSfitND}}
#' @export
print.GeDSfitND <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!missing(...)) warning("Arguments in '...' are currently ignored.")

  fit.call <- if (!is.null(x$extcall)) x$extcall else x$Call
  cat("\nExperimental multivariate Normal GeDS fit\n")
  if (!is.null(fit.call)) {
    cat("\nCall:\n", paste(deparse(fit.call), collapse = "\n"), "\n", sep = "")
  }
  cat("\nDimensions (", length(x$dimensions), "): ",
      paste(x$dimensions, collapse = ", "), "\n", sep = "")

  knot.counts <- x$Nintknots
  if (is.null(knot.counts)) knot.counts <- setNames(integer(0), character(0))
  cat("Linear internal knots: ",
      if (length(knot.counts)) {
        paste0(names(knot.counts), "=", knot.counts, collapse = ", ")
      } else {
        "none"
      },
      "\n", sep = "")

  rss <- c(
    `Order 2 (linear)` = x$dev.linear,
    `Order 3 (quadratic)` = x$dev.quadratic,
    `Order 4 (cubic)` = x$dev.cubic
  )
  if (length(rss)) {
    cat("\nResidual sum of squares:\n")
    print.default(format(rss, digits = digits), quote = FALSE)
  }
  solvers <- c(
    `Order 2` = if (!is.null(x$linear.fit)) x$linear.fit$solver,
    `Order 3` = if (!is.null(x$quadratic.fit)) x$quadratic.fit$solver,
    `Order 4` = if (!is.null(x$cubic.fit)) x$cubic.fit$solver
  )
  if (length(solvers)) {
    cat("Least-squares solver: ",
        paste0(names(solvers), "=", solvers, collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

################################################################################
################################### SUMMARY ####################################
################################################################################
#' @title Summary Method for GeDS, GeDSgam, GeDSboost
#' @name summary.GeDS
#' @description
#' Method for the generic function \code{\link[base]{summary}} that allows you to
#' print on screen the main information related to a fitted \code{"GeDS"},
#' \code{"GeDSgam"} or \code{"GeDSboost"} model.
#' Similar to \code{\link{print.GeDS}} but with some extra detail.
#'
#' @param object The \code{"GeDS"}, \code{"GeDSgam"} or \code{"GeDSboost"} object
#' for which the main information should be printed on screen.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function).
#'
#' @seealso \code{\link{print.GeDS}}
#'
#' @aliases summary.GeDS summary.GeDSboost summary.GeDSgam
#' @export

summary.GeDS <- function(object, ...)
{

  # 1) Header
  cat("\nGeometrically Desined (GeD) Spine Regression Model\n\n")

  # 2) Core print
  print.GeDS(object, ...)

  # 3) Additional details
  if (object$type == "LM - Univ" || object$type == "LM - Biv") {
    print(gaussian())

  } else if (object$type == "GLM - Univ" || object$type == "GLM - Biv") {
    print(object$args$family)
  }

  cat("phi = ", object$args$phi, "and q = ", object$args$q, "(stopping rule parameters);\n")
  cat("beta = ", object$args$beta)

  if (object$args$beta == 0.5) {
    cat(", meaning that the within-cluster mean residual and the cluster range were considered equally important when placing the knots.\n")
  } else if (object$args$beta > 0.5) {
    cat(", meaning that more weight was given to the within-cluster mean residual than to the cluster range when placing the knots.\n")
  } else if (object$args$beta < 0.5) {
    cat(", meaning that more weight was given to the cluster range than to the within-cluster mean residual when placing the knots.\n")
  }


  invisible(object)
}

#' Summarize an Experimental Multivariate GeDS Fit
#'
#' Prints the model dimensions, selected linear knot counts, residual sums of
#' squares, details of the Stage A stopping rule and tensor-basis size, and any
#' Stage B fits omitted by the coefficient limit.
#'
#' @param object A \code{"GeDSfitND"} object returned by \code{NGeDS()}.
#' @param digits Number of significant digits used for numerical output.
#' @param ... Further arguments passed to \code{print.GeDSfitND()}.
#'
#' @return The input object, invisibly.
#' @seealso \code{\link{print.GeDSfitND}}, \code{\link{NGeDS}}
#' @export
summary.GeDSfitND <- function(
    object, digits = max(3L, getOption("digits") - 3L), ...)
{
  print.GeDSfitND(object, digits = digits, ...)
  cat("\nStage A:\n")
  cat("  Stopping rule: ", object$args$stoptype,
      " (phi = ", format(object$args$phi, digits = digits),
      ", q = ", object$args$q, ")\n", sep = "")
  cat("  Knot-placement beta: ",
      format(object$args$beta, digits = digits), "\n", sep = "")
  cat("  Fitted iterations: ", object$iters, "\n", sep = "")
  cat("  Selected iteration: ", object$selected.iteration, "\n", sep = "")
  cat("  Maximum coefficients: ",
      format(object$args$max.coef, scientific = FALSE), "\n", sep = "")
  if (length(object$stageA$basis.sizes)) {
    cat("  Largest fitted tensor basis: ",
        format(max(object$stageA$basis.sizes), scientific = FALSE),
        " coefficients\n", sep = "")
  }
  if (!is.null(object$stageA$stop.reason)) {
    cat("  Termination: ", object$stageA$stop.reason, "\n", sep = "")
  }
  skipped <- object$stageB$skipped
  if (length(skipped) && any(!is.na(skipped))) {
    cat("\nStage B fits omitted by max.coef:\n")
    for (fit.name in names(skipped)[!is.na(skipped)]) {
      cat("  ", fit.name, ": ", skipped[[fit.name]], "\n", sep = "")
    }
  }
  invisible(object)
}

################################################################################
############################### SHAPECONSTRAIN #################################
################################################################################
.shape_constrained_refit <- function(X, Y, Z = NULL, weights = NULL, offset = NULL,
                                     full_knots, n, shape_constraint, eps, ridge,
                                     center_basis = FALSE) {

  if (is.null(weights)) weights <- rep(1, length(Y))
  if (is.null(offset)) offset <- rep(0, length(Y))
  if (!is.null(Z)) Z <- as.matrix(Z)

  basisMatrix <- splineDesign(
    knots = full_knots,
    x = X,
    ord = n,
    derivs = rep(0, length(X)),
    outer.ok = TRUE
  )

  basisFit <- if (isTRUE(center_basis)) {
    sweep(basisMatrix, 2L, colMeans(basisMatrix), check.margin = FALSE)
  } else {
    basisMatrix
  }

  basisMatrix2 <- cbind(basisFit, Z)
  p_spline <- NCOL(basisMatrix)

  Y0 <- Y - offset

  shape_cons <- make_shape_constraint(
    knots = full_knots,
    n = n,
    p = p_spline,
    shape_constraint = shape_constraint,
    eps = eps
  )

  C <- shape_cons$C
  b <- shape_cons$b

  n_extra <- NCOL(basisMatrix2) - p_spline

  if (n_extra > 0) {
    C <- cbind(C, matrix(0, nrow = nrow(C), ncol = n_extra))
  }

  theta <- constrained_wls(
    X = basisMatrix2,
    y = Y0,
    weights = weights,
    C = C,
    b = b,
    ridge = ridge
  )

  predicted <- as.numeric(basisMatrix2 %*% theta + offset)
  residuals <- Y - predicted
  rss <- .weighted_rss(residuals, weights)

  list(
    theta = theta,
    predicted = predicted,
    residuals = residuals,
    rss = rss,
    basis = basisMatrix,
    p_spline = p_spline
  )
}

#' @title Apply Shape Constraints to GeDS Fits
#' @name shapeConstrain
#' @description
#' Generic function for imposing shape constraints on fitted model objects.
#'
#' @param object A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' This generic function is currently implemented for selected \code{"GeDS"},
#' \code{"GeDSboost"} and \code{"GeDSgam"} class objects.
#'
#' @return An object of the same general class as \code{object}, with the selected
#' fit updated according to the requested shape constraint.
#'
#' @seealso \code{\link{NGeDS}}, \code{\link{coef.GeDS}},
#' \code{\link{knots.GeDS}}, \code{\link{predict.GeDS}}
#'
#' @aliases shapeConstrain
#' @importFrom splines splineDesign
#' @export

shapeConstrain <- function(object, ...) {
  UseMethod("shapeConstrain")
}

#' @title Apply Shape Constraints to a Fitted GeDS Model
#' @name shapeConstrain.GeDS
#' @description
#' Method for \code{\link{shapeConstrain}} that imposes monotonicity or
#' convexity/concavity constraints on a fitted univariate \code{"GeDS"} model.
#'
#' @param object A fitted \code{"GeDS"} class object produced by
#' \code{\link{NGeDS}}.
#' @param n Integer value (2, 3 or 4) specifying the order
#' (\eqn{=} degree \eqn{+ 1}) of the GeDS fit to be shape-constrained.
#' By default equal to \code{3L}; non-integer values will be passed to
#' \code{\link{as.integer}}.
#' @param shape_constraint Character vector specifying the shape constraint to
#' impose. Possible values are \code{"increasing"}, \code{"decreasing"},
#' \code{"convex"} and \code{"concave"}. Monotonicity and curvature constraints
#' may be combined, for example \code{c("increasing", "convex")}.
#' @param eps Numeric value specifying the lower bound used in the linear
#' inequality constraints. By default equal to \code{0}.
#' @param ridge Small positive numeric value added to the diagonal of the
#' quadratic programming matrix for numerical stability. By default equal to
#' \code{1e-8}.
#' @param ... Further arguments. Currently ignored.
#'
#' @details
#' The method leaves the knot vector selected by the GeDS algorithm unchanged
#' and re-estimates the B-spline coefficients conditional on those knots by
#' solving a linearly constrained weighted least-squares problem.
#'
#' Monotonicity constraints are imposed through linear inequality constraints on
#' first divided differences of the B-spline coefficients. Convexity and
#' concavity constraints are imposed analogously through second divided
#' differences.
#'
#' At present, this method is implemented only for univariate GeDS fits with
#' Gaussian response, i.e. objects with type \code{"LM - Univ"}. Standard
#' \code{lm}-based confidence intervals are not computed for the constrained fit,
#' since the coefficient estimates are obtained from a constrained least-squares
#' problem.
#'
#' @return A \code{"GeDS"} class object with the selected fit of order \code{n}
#' updated after re-estimating the coefficients under the requested shape
#' constraint. The returned object also contains a \code{shape_constraint}
#' component recording the imposed constraint.
#'
#' @seealso \code{\link{shapeConstrain}}, \code{\link{NGeDS}},
#' \code{\link{coef.GeDS}}, \code{\link{knots.GeDS}},
#' \code{\link{predict.GeDS}}
#'
#' @examples
#' if (requireNamespace("quadprog", quietly = TRUE)) {
#'   set.seed(123)
#'   N <- 300
#'   X <- sort(runif(N, 0, 1))
#'   f_1 <- function(x) exp(2 * x)
#'   Y <- f_1(X) + rnorm(N, sd = 0.2)
#'
#'   Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995)
#'   Gmod_sc <- shapeConstrain(Gmod, n = 4,
#'                             shape_constraint = c("increasing", "convex"))
#'
#'   plot(X, Y, pch = 16, col = "grey70")
#'   lines(X, f_1(X), lwd = 2)
#'   lines(X, predict(Gmod, n = 4), lwd = 2, col = "steelblue4")
#'   lines(X, predict(Gmod_sc, n = 4), lwd = 2, col = "darkred")
#' }
#'
#' @rdname shapeConstrain.GeDS
#' @aliases shapeConstrain.GeDS
#' @method shapeConstrain GeDS
#' @export

shapeConstrain.GeDS <- function(object,
                                n = 3L,
                                shape_constraint = "increasing",
                                eps = 0,
                                ridge = 1e-8, ...) {

  if (!missing(...)) {
    warning("Arguments in '...' are currently ignored.")
  }

  n <- validate_GeDS_order(n)

  if (!identical(object$type, "LM - Univ")) {
    stop("Shape constraints are currently implemented only for univariate LM GeDS fits.",
         call. = FALSE)
  }

  shape_constraint <- match.arg(
    shape_constraint,
    choices = c("increasing", "decreasing", "convex", "concave"),
    several.ok = TRUE
  )

  if (all(c("increasing", "decreasing") %in% shape_constraint)) {
    stop("Cannot impose both increasing and decreasing constraints.",
         call. = FALSE)
  }

  if (all(c("convex", "concave") %in% shape_constraint)) {
    stop("Cannot impose both convex and concave constraints.",
         call. = FALSE)
  }

  fit_name <- switch(as.character(n),
                     "2" = "linear.fit",
                     "3" = "quadratic.fit",
                     "4" = "cubic.fit")

  dev_name <- switch(as.character(n),
                     "2" = "dev.linear",
                     "3" = "dev.quadratic",
                     "4" = "dev.cubic")

  refit <- .shape_constrained_refit(
    X = object$args$X,
    Y = object$args$Y,
    Z = object$args$Z,
    weights = object$args$weights,
    offset = object$args$offset,
    full_knots = knots(object, n = n, options = "all"),
    n = n,
    shape_constraint = shape_constraint,
    eps = eps,
    ridge = ridge
  )

  object[[fit_name]]$theta <- refit$theta
  object[[fit_name]]$predicted <- refit$predicted
  object[[fit_name]]$residuals <- refit$residuals
  object[[fit_name]]$rss <- refit$rss
  object[[fit_name]]$basis <- refit$basis

  # Usual lm-based confidence intervals are not valid after constrained LS
  object[[fit_name]]$temporary <- NULL
  object[[fit_name]]$nci <- NULL
  object[[fit_name]]$aci <- NULL

  if (!is.null(object[[fit_name]]$polygon)) {
    object[[fit_name]]$polygon$thetas <- refit$theta[seq_len(refit$p_spline)]
  }

  object[[dev_name]] <- refit$rss

  object$shape_constraint <- list(
    n = n,
    constraint = shape_constraint,
    eps = eps,
    note = "Coefficients refitted under shape constraints conditional on the selected GeDS knots."
  )

  class(object) <- unique(c("shapeConstrainedGeDS", class(object)))

  object
}

# ###############################################################################
# ############################### SEPARATE LINEAR ###############################
# ###############################################################################
# # Greville abscissae for an order-n basis on `knots`: B %*% rep(1,p) == 1 and
# # B %*% .greville(knots,n) == x, i.e. the coefficient-space images of {1, x}.
# .greville <- function(knots, n) {
#   p <- length(knots) - n
#   if (n == 1L) return(knots[seq_len(p)])
#   vapply(seq_len(p), function(i) mean(knots[(i + 1):(i + n - 1)]), numeric(1))
# }
# 
# # Project a fitted spline curve `eta` (predictor/link scale) onto {1, x} in the
# # (weighted) observation inner product. Returns identifiable intercept/slope and
# # the orthogonal non-linear remainder f.
# .separate_linear_curve <- function(eta, x, w = NULL) {
#   if (is.null(w)) w <- rep(1, length(eta))
#   Tm <- cbind(1, x)
#   cf <- stats::lm.wfit(Tm, eta, w)$coefficients
#   if (anyNA(cf))
#     stop("Could not separate the linear component (is the covariate constant?).",
#          call. = FALSE)
#   lin <- as.numeric(Tm %*% cf)
#   list(beta0 = unname(cf[1]), beta1 = unname(cf[2]),
#        linear = lin, f = as.numeric(eta - lin))
# }
# 
# #' @title Separate the Linear Component from a GeDS Spline Fit
# #' @name separateLinear
# #' @description
# #' Re-expresses the GeD spline component(s) of a fitted model as the sum of an
# #' identifiable linear term and a purely non-linear ("curvature") term, i.e.
# #' \eqn{\eta \approx \beta_0 + \beta_1 X + f(X)} with \eqn{f} orthogonal to
# #' \eqn{\{1, X\}}. The fitted curve is unchanged; only its parameterization is.
# #' @param object A fitted \code{"GeDS"}, \code{"GeDSboost"} or \code{"GeDSgam"}
# #' object.
# #' @param ... Further arguments passed to methods.
# #' @details
# #' A GeD spline basis reproduces constant and linear functions exactly, so a
# #' linear term in the spline variable is not separately identifiable from the
# #' spline. This generic resolves that by orthogonalizing the fitted spline curve
# #' against \eqn{\{1, X\}}. For \code{"GeDS"} fits the result applies to the
# #' single spline component; for \code{"GeDSboost"}/\code{"GeDSgam"} it is applied
# #' per univariate GeD spline base-learner, yielding a linear + additive
# #' non-linear decomposition. For non-Gaussian families the slope is on the
# #' \emph{link} scale (reported in the \code{scale} field).
# #' @return
# #' An object of class \code{"separateLinear"}. For a \code{"GeDS"} fit, a list with:
# #' \describe{
# #'   \item{\code{beta0}, \code{beta1}}{Intercept and slope of the identifiable
# #'   linear component (\code{beta1} is the weighted least-squares slope on the
# #'   scale given by \code{scale}).}
# #'   \item{\code{f}}{Non-linear component at the observed covariate values,
# #'   orthogonal to \eqn{\{1, X\}}.}
# #'   \item{\code{theta_f}}{Coefficients of \code{f} in the original B-spline basis
# #'   (same knots), so \code{f} can be evaluated at new covariate values.}
# #'   \item{\code{knots}, \code{n}}{Full knot vector and order of the decomposed fit.}
# #'   \item{\code{scale}}{\code{"response"} or \code{"link"} -- the scale of
# #'   \code{beta1} and \code{f}.}
# #'   \item{\code{fitted_spline}}{The spline contribution
# #'   \eqn{\beta_0 + \beta_1 X + f(X)}; identical to the original fit.}
# #' }
# #' For \code{"GeDSboost"}/\code{"GeDSgam"} fits, a list with \code{learners} (one
# #' such decomposition -- \code{variable}, \code{beta0}, \code{beta1}, \code{f},
# #' \code{theta_f}, \code{knots}, \code{n} -- per univariate GeD spline
# #' base-learner), plus \code{n} and \code{scale}.
# #'
# #' @examples
# #' set.seed(123)
# #' N <- 300
# #' X <- sort(runif(N, -2, 2))
# #' Y <- 2 * X + cos(2 * X) + rnorm(N, sd = 0.2)
# #'
# #' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995)
# #' dec  <- separateLinear(Gmod, n = 3)
# #'
# #' dec$beta1
# #' all.equal(as.numeric(predict(Gmod, n = 3)), dec$fitted_spline) # fit unchanged
# #'
# #' linear_part <- dec$beta0 + dec$beta1 * X
# #'
# #' oldpar <- par(no.readonly = TRUE)
# #' par(mfrow = c(1, 2))
# #'
# #' plot(X, Y, pch = 16, col = "grey75",
# #'      main = "Linear component",
# #'      xlab = "X", ylab = "Y")
# #' lines(X, predict(Gmod, n = 3), lwd = 2)
# #' lines(X, linear_part, lwd = 2, col = "red")
# #' legend("topleft",
# #'        legend = c("GeDS fit", "Linear component"),
# #'        col = c("black", "red"),
# #'        lwd = 2, bty = "n")
# #'
# #' plot(X, dec$f, type = "l", lwd = 2,
# #'      main = "Non-linear component",
# #'      xlab = "X", ylab = "f(X)")
# #' abline(h = 0, lty = 3)
# #'
# #' par(oldpar)
# #' @seealso \code{\link{shapeConstrain}}, \code{\link{NGeDS}},
# #' \code{\link{NGeDSboost}}, \code{\link{NGeDSgam}}
# #' @importFrom splines splineDesign
# #' @export
# separateLinear <- function(object, ...) UseMethod("separateLinear")
# 
# #' @rdname separateLinear
# #' @param n Integer (2, 3 or 4) order of the fit to decompose. Default \code{3L}.
# #' @method separateLinear GeDS
# #' @importFrom splines splineDesign
# #' @importFrom stats lm.wfit
# #' @export
# separateLinear.GeDS <- function(object, n = 3L, ...) {
#   if (!missing(...)) warning("Arguments in '...' are currently ignored.")
#   n <- validate_GeDS_order(n)
#   if (!object$type %in% c("LM - Univ", "GLM - Univ"))
#     stop("separateLinear() is implemented only for univariate GeDS fits ",
#          "(type 'LM - Univ' or 'GLM - Univ').", call. = FALSE)
# 
#   X <- as.numeric(object$args$X)
#   w <- object$args$weights; if (is.null(w)) w <- rep(1, length(X))
# 
#   kn    <- knots(object, n = n, options = "all")
#   theta <- as.numeric(coef(object, n = n, onlySpline = TRUE))
#   B     <- splineDesign(kn, x = X, ord = n,
#                         derivs = rep(0, length(X)), outer.ok = TRUE)
#   eta_spline <- as.numeric(B %*% theta)          # spline part, predictor scale
# 
#   dec <- .separate_linear_curve(eta_spline, X, w)
#   xi  <- .greville(kn, n)
#   scale <- if (object$type == "GLM - Univ") "link" else "response"
# 
#   structure(
#     list(beta0 = dec$beta0, beta1 = dec$beta1, f = dec$f,
#          theta_f = theta - dec$beta0 - dec$beta1 * xi,   # f in original B-basis
#          knots = kn, n = n, scale = scale, fitted_spline = eta_spline),
#     class = "separateLinear")
# }
# 
# #' @method print separateLinear
# #' @export
# print.separateLinear <- function(x, digits = 4L, ...) {
#   cat("Linear / non-linear decomposition of a GeDS fit\n")
#   cat("  order:", x$n, " scale:", x$scale, "\n\n")
#   if (!is.null(x$learners)) {
#     for (nm in names(x$learners)) {
#       l <- x$learners[[nm]]
#       cat(sprintf("  %-16s slope = % .*f    ||f|| = %.*f\n",
#                   nm, digits, l$beta1, digits, sqrt(sum(l$f^2))))
#     }
#   } else {
#     cat(sprintf("  intercept = %.*f   slope = %.*f   ||f|| = %.*f\n",
#                 digits, x$beta0, digits, x$beta1, digits, sqrt(sum(x$f^2))))
#   }
#   invisible(x)
# }

