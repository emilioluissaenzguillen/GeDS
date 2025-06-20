################################################################################
################################################################################
############################## lines.GeDS method ###############################
################################################################################
################################################################################
#' @title Lines Method for GeDS Objects
#' @name lines.GeDS
#' @description
#' Lines method for GeDS objects. Adds a GeDS curve to an existing plot. 
#' @param x A GeDS object, as returned by \code{NGeDS()} or \code{GGeDS()}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS fit that should be plotted. By default equal to
#' \code{3L}. Non-integer values will be passed to the function
#' \code{\link{as.integer}}.
#' @param transform A function that can be used to transform the scale of the Y
#' axis. Typically it can be the inverse of the link function if the plot is on
#' the scale of the response variable.
#' @param onlySpline Logical variable specifying whether only the spline
#' component of the fitted GeDS predictor model  should be plotted or
#' alternatively also the parametric component (see
#' \code{\link[=formula.GeDS]{formula}}) should be plotted.
#' @param data An optional \code{data.frame}, \code{list}, or \code{environment}.
#' It should contain values of the independent variables for which the GeDS predicted
#' values should be plotted. If left empty, the values are extracted from the object \code{x}.
#' 
#' @param ... Further arguments to be passed to the default
#' \code{\link[graphics]{lines}} function.
#'
#' @details
#' This method can be used to add a curve corresponding to a particular GeDS fit
#' to an active plot.
#'
#' As GeDS objects contain three different fits (linear, quadratic and cubic),
#' it is possible to specify the order of the GeDS regression to be plotted via
#' the input argument \code{n}.
#'
#' @examples
#'
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a GeDS regression model using NGeDS
#' (Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Plot the GeDS third order fit (the quadratic one)
#' # without its corresponding Polygon
#' plot(Gmod, type = "none")
#'
#' # Add a curve corresponding to the second order fit (the linear one)
#' lines(Gmod, n = 2, col = "green", lwd = 2, lty = 3)
#' 
#' @seealso \code{\link[graphics]{lines}} for the definition of the generic
#' function; \code{\link{NGeDS}} and \code{\link{GGeDS}} for examples.
#' @importFrom graphics lines lines.default
#' @importFrom splines splineDesign
#' @export
#' @method lines GeDS

lines.GeDS <- function(x , n=3L, transform = function(x) x,
                       onlySpline = TRUE, data = data.frame(), ...)
{
  object <- x
  if(object$Type == "LM - Biv" || object$Type == "GLM - Biv") stop("Works only with univariate spline objects")
  
  extr <- object$Args$extr
  
  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # Extract fit
  if(n == 2L) {
    temp <- object$Linear.Fit
  } else if (n == 3L) {
    temp <- object$Quadratic.Fit
  } else if(n == 4L) {
    temp <- object$Cubic.Fit
  }
  
  kn <- knots.GeDS(Fn = object, n = n, options= "internal")
  fitters <- F
  if(is.null(object$terms)) fitters <- T
  if (fitters) {
    Predicted <- if (x$Type == "LM - Univ") temp$Predicted else if (x$Type == "GLM - Univ") x$Args$family$linkfun(temp$Predicted)
    Xvalues <- object$Args$X
    
  } else {
    
    if(!missing(data)) {
      dati2 <- read.formula(object$Formula,data)
      Xvalues <- dati2$X
      
      mm <- splineDesign(knots = sort(c(kn,rep(extr,n))), derivs = rep(0,length(Xvalues)),
                         x = Xvalues, ord = n, outer.ok = T)
      if (!onlySpline && !is.null(object$Args$Z)) mm <- cbind(mm,dati2$Z)
      
      offset <- if(!onlySpline && !is.null(object$Args$offset)) {
        dati2$offset
      } else {
        rep(0,length(Xvalues))
      }
      
    } else {
      Xvalues <- object$Args$X
      if(onlySpline && length(unique(Xvalues)) < 1000) {
        step <- (range(extr)[2] - range(extr)[1])/(1000)
        step <- rep(step,(1000))
        Xvalues <- min(extr)+c(0,cumsum(step))
      }
      
      mm <- splineDesign(knots = sort(c(kn,rep(extr,n))), derivs = rep(0,length(Xvalues)),
                         x = Xvalues, ord = n, outer.ok = T)
      if(!onlySpline && !is.null(object$Args$Z)) mm <- cbind(mm, object$Args$Z)
      
      offset <- if(!onlySpline && !is.null(object$Args$offset)) {
        object$Args$offset
      } else {
        rep(0,length(Xvalues))
      }
    }
    
    th <- coef.GeDS(object, n=n, onlySpline = onlySpline)
    Predicted <- mm%*%th + offset
  }
  
  Predicted <- transform(Predicted)
  lines.default(Xvalues, Predicted,...)
}


