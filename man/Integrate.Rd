% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntegrateDerive.R
\name{Integrate}
\alias{Integrate}
\title{Defined Integral of GeDS Objects}
\usage{
Integrate(object = NULL, knots = NULL, coef = NULL, from, to, n = 3L)
}
\arguments{
\item{object}{An object of class \code{"GeDS"} containing the  GeDS fit
which should be integrated. It should be the result of fitting a univariate
GeDS regression via \code{\link{NGeDS}} or \code{\link{GGeDS}}. If this is
provided, the \code{knots} and \code{coef} parameters will be automatically
extracted from the \code{GeDS} object. If \code{object} is \code{NULL}, the
user must provide the \code{knots} and \code{coef} vectors explicitly.}

\item{knots}{A numeric vector of knots. This is required if \code{object} is 
\code{NULL}. If a \code{GeDS} object is provided, this parameter is ignored.}

\item{coef}{A numeric vector of coefficients. This is required if \code{object} is 
\code{NULL}. If a \code{GeDS} object is provided, this parameter is ignored}

\item{from}{Optional numeric vector containing the lower limit(s) of
integration. It should be either of size one or of the same size as the
argument \code{to}. If left unspecified, by default it is set to the left-most
limit of the interval embedding the observations of the independent variable.}

\item{to}{Numeric vector containing the upper limit(s) of integration.}

\item{n}{Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
\eqn{ + 1}) of the GeDS fit to be integrated. By default equal to \code{3L}.
Non-integer values will be passed to the function \code{\link{as.integer}}.}
}
\description{
This function computes defined integrals of a fitted GeDS regression model.
}
\details{
The function relies on the well known property (c.f. De Boor, 2001, Chapter
X, formula (33)) that the integral of a linear combination of appropriately
normalized B-splines (i.e., the standard representation of a
GeDS regression model) is equal to the sum of its corresponding coefficients.

Since the function is based on this property, it is designed to work only on
the predictor scale in the GNM (GLM) framework.

If the argument \code{from} is a single value, then it is taken as the lower
limit of integration for all the defined integrals required, whereas the upper
limits of integration are the values contained in the argument \code{to}. If
the arguments \code{from} and \code{to} are of same size, the integrals
(as many as the size) are computed by sequentially taking the pairs of values
in the \code{from} and \code{to} vectors as limits of integration.
}
\examples{

# Generate a data sample for the response variable
# Y and the single covariate X
# see Dimitrova et al. (2023), section 4.1
set.seed(123)
N <- 500
f_1 <- function(x) (10*x/(1+100*x^2))*4+4
X <- sort(runif(N, min = -2, max = 2))
# Specify a model for the mean of Y to include only
# a component non-linear in X, defined by the function f_1
means <- f_1(X)
# Add (Normal) noise to the mean of Y
Y <- rnorm(N, means, sd = 0.1)
# Fit GeDS regression using NGeDS
Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = .995, Xextr = c(-2,2))
# Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
# and $\int_{1}^{1} f(x)dx$
# $f$ being the quadratic fit
Integrate(Gmod, from = 1, to = c(-1,1), n = 3)
# Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
# and $\int_{-1}^{1} f(x)dx$
# $f$ being the quadratic fit
Integrate(Gmod, from = c(1,-1), to = c(-1,1), n = 3)

# Compute $\int_{-\infty}^{x} f(s)ds$
Integrate(Gmod, from = rep(-Inf, N), to = X, n = 3)

}
\references{
De Boor, C. (2001). \emph{A Practical Guide to Splines (Revised Edition)}.
Springer, New York.
}
