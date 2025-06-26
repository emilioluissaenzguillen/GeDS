################################################################################
################################################################################
################################### Integrate ##################################
################################################################################
################################################################################
#' @title Defined Integral of GeDS Objects
#' @name Integrate 
#' @description
#' This function computes defined integrals of a fitted GeDS regression model.
#' 
#' @param object An object of class \code{"GeDS"} containing the  GeDS fit
#' which should be integrated. It should be the result of fitting a univariate
#' GeDS regression via \code{\link{NGeDS}} or \code{\link{GGeDS}}. If this is
#' provided, the \code{knots} and \code{coef} parameters will be automatically
#' extracted from the \code{GeDS} object. If \code{object} is \code{NULL}, the
#' user must provide the \code{knots} and \code{coef} vectors explicitly.
#' @param knots A numeric vector of knots. This is required if \code{object} is 
#'   \code{NULL}. If a \code{GeDS} object is provided, this parameter is ignored.
#' @param coef A numeric vector of coefficients. This is required if \code{object} is 
#'   \code{NULL}. If a \code{GeDS} object is provided, this parameter is ignored
#' @param from Optional numeric vector containing the lower limit(s) of
#' integration. It should be either of size one or of the same size as the
#' argument \code{to}. If left unspecified, by default it is set to the left-most
#' limit of the interval embedding the observations of the independent variable.
#' @param to Numeric vector containing the upper limit(s) of integration.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{ + 1}) of the GeDS fit to be integrated. By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#'
#' @details
#' The function relies on the well known property (c.f. De Boor, 2001, Chapter
#' X, formula (33)) that the integral of a linear combination of appropriately
#' normalized B-splines (i.e., the standard representation of a
#' GeDS regression model) is equal to the sum of its corresponding coefficients.
#'
#' Since the function is based on this property, it is designed to work only on
#' the predictor scale in the GNM (GLM) framework.
#'
#' If the argument \code{from} is a single value, then it is taken as the lower
#' limit of integration for all the defined integrals required, whereas the upper
#' limits of integration are the values contained in the argument \code{to}. If
#' the arguments \code{from} and \code{to} are of same size, the integrals
#' (as many as the size) are computed by sequentially taking the pairs of values
#' in the \code{from} and \code{to} vectors as limits of integration.
#'
#' @examples
#' 
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' # see Dimitrova et al. (2023), section 4.1
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#' # Fit GeDS regression using NGeDS
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = .995, Xextr = c(-2,2))
#' # Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
#' # and $\int_{1}^{1} f(x)dx$
#' # $f$ being the quadratic fit
#' Integrate(Gmod, from = 1, to = c(-1,1), n = 3)
#' # Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
#' # and $\int_{-1}^{1} f(x)dx$
#' # $f$ being the quadratic fit
#' Integrate(Gmod, from = c(1,-1), to = c(-1,1), n = 3)
#' 
#' # Compute $\int_{-\infty}^{x} f(s)ds$
#' Integrate(Gmod, from = rep(-Inf, N), to = X, n = 3)
#' 
#' @export
#'
#' @references
#' De Boor, C. (2001). \emph{A Practical Guide to Splines (Revised Edition)}.
#' Springer, New York.

Integrate <- function(object = NULL, knots = NULL, coef = NULL, from, to, n = 3L)
{
  # If object is not NULL, check if it's a GeDS object
  if (!is.null(object)) {
    if (!inherits(object, "GeDS")) stop("Incorrect object class")
    if(!is.null(knots) || !is.null(coef)) {
      warning("object and knots/coefficients were providely simultaneously, integral will be computed from object")
    }
    
    # Extract knots and coefficients from the GeDS object
    kn <- knots(object, n=n, options="all")
    lastkn <- length(kn)
    theta <- coef(object, n = n)
  } else {
    # If object is NULL, ensure knots and coefficients are provided
    if (is.null(knots) || is.null(coef)) {
      stop("Either provide a GeDS object or specify knots and coefficients")
    }
    if ( length(coef) != length(knots)-n ) stop("length(coef) should be equal to length(knots)-n!")
    kn <- knots; lastkn <- length(kn); theta <- coef
  }
  
  n <- as.integer(n)

  to <- as.numeric(to)
  if (missing(from)) {
    from = min(kn)
  } else {
    from <- as.numeric(from)
  }
  
  if (!length(from) %in% c(1,length(to))) stop("length of argument 'from' must be either 1 or length(to)")
  
  # coefficents
  p <- length(theta)
  newtheta <- numeric(p)
  # knots: (t_{k+n} - t_{k})/n
  knnew <- (kn[-(1:n)]-kn[-((lastkn-n+1):lastkn)])/n
  
  # \sum_{j=1}^{s-1}\sum_{k=1}^j\theta_k\frac{\bar{\tau}_{k+n}-\bar{\tau}_{k}}{n} * N_{j,n+1}(x)
  # newtheta[1] <- theta[1]*knnew[1]
  # for(i in 2:p) {
  #   newtheta[i] <- newtheta[i-1] + theta[i]*knnew[i]
  # }
  newtheta <- cumsum(theta * knnew)
  
  resFrom <- sapply(from, gedsint, knts = kn, coefs = newtheta, n = n)
  resTo <- sapply(to, gedsint, knts = kn, coefs = newtheta, n = n)
  
  res <- rowSums(cbind(resTo,-resFrom))
  return(res)
}

#' @importFrom splines splineDesign
gedsint <- function(val, knts, coefs, n){
  if (val <= min(knts)) return(0)
  if (val >= max(knts))  return(1)
  pos <- min(which(knts>=val))
  basisMatrix <- splineDesign(knots = c(knts,max(knts)),
                              derivs = rep(0,length(val)),
                              x = val, ord = n+1, outer.ok = T)
  
  ris <- as.numeric(basisMatrix[,1:(pos-1)]%*%coefs[1:(pos-1)])
  return(ris)
}

################################################################################
################################################################################
#################################### Derive ####################################
################################################################################
################################################################################
#' @title Derivative of GeDS Objects
#' @name Derive
#' @description
#' This function computes derivatives of a fitted GeDS regression model.
#' @param object An object of class \code{"GeDS"} containing the GeDS fit
#' which should be differentiated. It should be the result of fitting a 
#' univariate GeDS regression via \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param order Integer value indicating the order of differentiation required
#' (e.g. first, second or higher derivatives). Note that \code{order} should be
#' lower than \code{n} and that non-integer values will be passed to the
#' function \code{\link{as.integer}}.
#' @param x Numeric vector containing values of the independent variable at
#' which the derivatives of order \code{order} should be computed.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{ + 1}) of the GeDS fit to be differentiated. By default equal to
#' \code{3L}.
#' @details
#' This function relies on the \code{\link[splines]{splineDesign}} function to compute 
#' the exact derivatives of the GeDS fit. Specifically, it leverages the well-known 
#' property that the \eqn{m}-th derivative of a spline (for \eqn{m = 1, 2, \ldots}) 
#' can be obtained by differentiating its B-spline basis functions. This property is 
#' detailed, e.g., in De Boor (2001, Chapter X, formula (15)).
#'
#' Note that the GeDS fit is a B-spline representation of the predictor. Consequently, 
#' the derivative is computed with respect to the predictor scale and not the response 
#' scale. This implies that, in the GNM(GLM) framework, the function 
#' does not return derivatives of the conditional mean on the response scale, but rather 
#' of the underlying linear predictor scale.
#' 
#' @examples
#' 
#' # Generate a data sample for the response variable
#' # Y and the covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit GeDS regression with only a spline component in the predictor model
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2))
#'
#' # Compute the second derivative of the cubic GeDS fit
#' # at the points 0, -1 and 1
#' Derive(Gmod, x = c(0, -1, 1), order = 2, n = 4)
#' @importFrom splines splineDesign
#' @export
#'
#' @references De Boor, C. (2001). \emph{A Practical Guide to Splines (Revised Edition)}.
#' Springer, New York.

Derive <- function(object, order = 1L, x, n = 3L)
  {
  if (!inherits(object, "GeDS")) stop("incorrect object class")
  if (!(object$type %in% c("LM - Univ","GLM - Univ"))) stop("Implemented only for the univariate case")
  x <- as.numeric(x)
  l <- length(x)
  if (order >= n) stop("'order' must be less than 'n'")
  if (length(n)!=1) stop("'n' must have length 1")
  n <- as.integer(n)
  if (length(order)!=1) stop("'order' must have length 1")
  order <- as.integer(order)
  kn <- knots(object, options = "all", n = n)
  thetas <- coef(object, n = n)
  basis <- splineDesign(knots = kn, x = x, ord = n, derivs = rep(order,l), outer.ok = TRUE)
  der <- as.numeric(basis%*%thetas)
  der
}

