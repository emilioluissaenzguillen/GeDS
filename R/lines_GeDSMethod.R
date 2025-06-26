################################################################################
################################################################################
############################## lines.GeDS method ###############################
################################################################################
################################################################################
#' @title Lines Method for GeDS Objects
#' @name lines.GeDS
#' @description
#' Lines method for GeDS objects. Adds a GeDS curve to an existing plot. 
#' @param x A \code{"GeDS"} class object, as returned by \code{NGeDS()} or \code{GGeDS()}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS fit that should be plotted. By default equal to
#' \code{3L}. Non-integer values will be passed to the function
#' \code{\link{as.integer}}.
#' @param transform A function that can be used to transform the scale of the
#' \eqn{y}-axis. Typically it can be the inverse of the link function if the plot
#' is on the scale of the response variable.
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

lines.GeDS <- function(x , n = 3L, transform = function(x) x,
                       onlySpline = TRUE, data = data.frame(), ...)
{
  object <- x
  if(object$type == "LM - Biv" || object$type == "GLM - Biv") stop("Works only with univariate spline objects")
  
  extr <- object$args$extr
  
  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # Extract fit
  if(n == 2L) {
    temp <- object$linear.fit
  } else if (n == 3L) {
    temp <- object$quadratic.fit
  } else if(n == 4L) {
    temp <- object$cubic.fit
  }
  
  kn <- knots.GeDS(Fn = object, n = n, options= "internal")
  fitters <- F
  if(is.null(object$terms)) fitters <- T
  if (fitters) {
    predicted <- if (x$type == "LM - Univ") temp$predicted else if (x$type == "GLM - Univ") x$args$family$linkfun(temp$predicted)
    Xvalues <- object$args$X
    
  } else {
    
    if(!missing(data)) {
      dati2 <- read.formula(object$Formula,data)
      Xvalues <- dati2$X
      
      mm <- splineDesign(knots = sort(c(kn,rep(extr,n))), derivs = rep(0,length(Xvalues)),
                         x = Xvalues, ord = n, outer.ok = T)
      if (!onlySpline && !is.null(object$args$Z)) mm <- cbind(mm,dati2$Z)
      
      offset <- if(!onlySpline && !is.null(object$args$offset)) {
        dati2$offset
      } else {
        rep(0,length(Xvalues))
      }
      
    } else {
      Xvalues <- object$args$X
      if(onlySpline && length(unique(Xvalues)) < 1000) {
        step <- (range(extr)[2] - range(extr)[1])/(1000)
        step <- rep(step,(1000))
        Xvalues <- min(extr)+c(0,cumsum(step))
      }
      
      mm <- splineDesign(knots = sort(c(kn,rep(extr,n))), derivs = rep(0,length(Xvalues)),
                         x = Xvalues, ord = n, outer.ok = T)
      if(!onlySpline && !is.null(object$args$Z)) mm <- cbind(mm, object$args$Z)
      
      offset <- if(!onlySpline && !is.null(object$args$offset)) {
        object$args$offset
      } else {
        rep(0,length(Xvalues))
      }
    }
    
    th <- coef.GeDS(object, n=n, onlySpline = onlySpline)
    predicted <- mm%*%th + offset
  }
  
  predicted <- transform(predicted)
  lines.default(Xvalues, predicted,...)
}


