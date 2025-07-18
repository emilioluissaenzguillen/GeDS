% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.formula.R
\name{f}
\alias{f}
\title{Defining the Covariates for the Spline Component in a GeDS Formula}
\usage{
f(x, xx = NULL, ...)
}
\arguments{
\item{x}{Numeric vector containing \eqn{N} sample values of the covariate
chosen to enter the spline
regression component of the predictor model.}

\item{xx}{Numeric vector containing \eqn{N} sample values for the second
covariate (in case \code{\link{NGeDS}}/\code{\link{GGeDS}} is run for two
dimensions). It has to be either \code{NULL} (the default) or a vector of size
\eqn{N}, same as \code{x}.}

\item{...}{Further arguments. As GeDS currently allows for up to two
covariates, specification of further arguments will return an error.}
}
\description{
In general the GeDS predictor model may include a GeD spline regression
component with respect to one or two independent variables and a parametric
component in which the remaining covariates may enter as additive terms.
GAM-GeDS and FGB-GeDS models may include more than one GeD spline regression
component.

The function \code{f} is to be used in the
\code{\link[=formula.GeDS]{formula}} argument of \code{\link{NGeDS}},
\code{\link{GGeDS}}, \code{\link{NGeDSgam}} or \code{\link{NGeDSboost}} in 
order to specify which independent variables (covariates) should be included
in the GeD spline regression component of the predictor model.
}
\note{
This function is intended to be used only as part of the
\code{\link[=formula.GeDS]{formula}} in a GeDS model via
\code{\link{NGeDS}}, \code{\link{GGeDS}}, \code{\link{NGeDSgam}} or
\code{\link{NGeDSboost}} and not to be called in other cases by the user.
}
\examples{
# Generate a data sample for the response variable Y and
# the covariates X, reg1, reg2 and off
set.seed(123)
N <- 500
f_1 <- function(x) (10*x/(1+100*x^2))*4+4
X <- sort(runif(N ,min = -2, max = 2))
reg1 <- runif(500, min = -0.1, max = 0.1)
reg2 <- runif(500, min = -0.2, max = 0.2)
off <- runif(500, min = -1, max = 1)
# Specify a model for the mean of Y to include a component non linear
# in X defined by the function f_1 and a linear one in the other covariates
means <- f_1(X) + 2*reg1 + 0.5*reg2 + off
# Add Normal noise to the mean of Y
Y <- rnorm(N, means, sd = 0.1)

# Specify a formula that will be used to model Y as a
# function of X, reg1, reg2 and off.
# The covariate X is for the spline component modeled as GeDS,
# reg1 and reg2 enter linearly, off is an offset, i.e. no coefficient
# will be estimated for it
formula <- Y ~ f(X) + reg1 + reg2 + offset(off)

# Fit a GeDS model specified in formula using NGeDS
(Gmod <- NGeDS(formula, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))


}
\seealso{
\code{\link[=formula.GeDS]{formula}}; \link{NGeDS}; \link{GGeDS};
\link{NGeDSgam}; \link{NGeDSboost}
}
