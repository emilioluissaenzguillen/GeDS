% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NGeDS.R
\name{NGeDS}
\alias{NGeDS}
\title{Geometrically Designed Spline Regression Estimation}
\usage{
NGeDS(
  formula,
  data,
  weights,
  beta = 0.5,
  phi = 0.99,
  min.intknots,
  max.intknots,
  q = 2L,
  Xextr = NULL,
  Yextr = NULL,
  show.iters = FALSE,
  stoptype = "RD",
  higher_order = TRUE,
  intknots_init = NULL,
  fit_init = NULL,
  only_pred = FALSE
)
}
\arguments{
\item{formula}{A description of the structure of the model to be fitted,
including the dependent and independent variables. See
\code{\link[=formula.GeDS]{formula}} for details.}

\item{data}{An optional \code{data.frame}, \code{list} or \code{environment}
containing the variables of the model. If not found in \code{data}, the variables
are taken from \code{environment(formula)}, typically the environment from which
\code{NGeDS} is called.}

\item{weights}{An optional vector of "prior weights" to be put on the
observations in the fitting process in case the user requires weighted GeDS
fitting. It should be \code{NULL} or a numeric vector of the same length as
the response variable in the argument \code{\link[=formula.GeDS]{formula}}.}

\item{beta}{Numeric parameter in the interval \eqn{[0,1]} tuning the knot
placement in stage A of GeDS. See Details.}

\item{phi}{Numeric parameter in the interval \eqn{(0,1)} specifying the
threshold for the stopping rule  (model selector) in stage A of GeDS. See
also \code{stoptype} and Details below.}

\item{min.intknots}{Optional parameter allowing the user to set a minimum
number of internal knots required. By default equal to zero.}

\item{max.intknots}{Optional parameter allowing the user to set a maximum
number of internal knots to be added by the GeDS estimation algorithm. By
default equal to the number of knots for the saturated GeDS model
(i.e., \eqn{\kappa = N - 2}).}

\item{q}{Numeric parameter which allows to fine-tune the stopping rule of
stage A of GeDS, by default equal to \code{2L}. See Details.}

\item{Xextr}{Numeric vector of 2 elements representing the left-most and
right-most limits of the interval embedding the observations of the first
independent variable. See Details.}

\item{Yextr}{Numeric vector of 2 elements representing the left-most and
right-most limits of the interval embedding the observations of the second
independent variable (if bivariate GeDS is run). See Details.}

\item{show.iters}{Logical variable indicating whether or not to print 
information about the fitting at each step.}

\item{stoptype}{A character string indicating the type of GeDS stopping rule
to be used. It should be either one of \code{"SR"}, \code{"RD"} or 
\code{"LR"}, partial match allowed. See Details.}

\item{higher_order}{Logical. If \code{TRUE}, after completing stage A, the
function proceeds to stage B and fits higher-order models (quadratic and cubic). 
Default is \code{TRUE}.}

\item{intknots_init}{A vector of initial internal knots for stage A. The default
is \code{NULL}, in which case stage A starts with a straight-line fit (i.e.,
no internal knots).}

\item{fit_init}{A list containing fitted values, \code{pred}, along with
corresponding internal knots, \code{intknots}, and coefficients, \code{coef},
representing the initial fit from which to begin stage A GeDS iteration (i.e.
departing from step 2). See Details.}

\item{only_pred}{Logical. If \code{TRUE} only predictions are computed.}
}
\value{
An object of class \code{"GeDS"} (a named list) with components:
\describe{
\item{type}{Character string indicating the type of regression performed.
This can be \code{"LM - Univ"}/\code{"LM - Biv"} respectively corresponding
to normal univariate/bivariate GeDS (implemented by \code{\link{NGeDS}}).}
\item{linear.intknots}{Vector containing the locations of the internal knots
of the second order GeD spline fit produced at stage A.}
\item{quadratic.intknots}{Vector containing the locations of the internal knots
of the third order GeD spline fit produced in stage B.}
\item{cubic.intknots}{Vector containing the locations of the internal knots
of the fourth order GeD spline fit produced in stage B.}
\item{dev.linear}{Deviance of the second order GeD spline fit, produced in
stage A.}
\item{dev.quadratic}{Deviance of the third order GeD spline fit, produced in
stage B.}
\item{dev.cubic}{Deviance of the fourth order GeD spline fit, produced in
stage B.}
\item{rss}{Vector containing the deviances of the second-order spline fits
computed at each iteration of stage A of GeDS.}
\item{linear.fit}{List containing the results from running \code{\link{SplineReg}}
function to fit the second order spline fit of stage A.}
\item{quadratic.fit}{List containing the results from running \code{\link{SplineReg}}
function to fit the third order spline fit of stage B.}
\item{cubic.fit}{List containing the results from running \code{\link{SplineReg}}
function to fit the fourth order spline fit of stage B.}
\item{stored}{Matrix containing the knot locations estimated at each iteration
of stage A.}
\item{args}{List containing the input arguments passed on the
\code{\link{Fitters}} functions.}
\item{call}{\code{call} to the \code{\link{Fitters}} functions.}
\item{Nintknots}{The final number of internal knots of the second order GeD
spline fit produced in stage A.}
\item{iters}{Number of iterations performed during stage A of the GeDS fitting
procedure.}   
\item{coefficients}{Matrix containing the fitted coefficients of the GeD
spline regression component and the parametric component at each iteration
of stage A.}   
\item{stopinfo}{List of values providing information related to the stopping
rule of stage A of GeDS. The sub-slots of \code{stopinfo} are \code{phis},
\code{phis_star}, \code{oldintc} and \code{oldslp}. The sub-slot \code{phis}
is a vector containing the values of the ratios of deviances (or the
difference of deviances if the \code{LR} stopping rule was chosen). The
sub-slots \code{phis_star}, \code{oldintc} and \code{oldslp} are non-empty
slots if the \code{SR} stopping rule was chosen. These respectively contain
the values at each iteration of stage A of \eqn{\hat{\phi}_{\kappa}},
\eqn{\hat{\gamma}_0} and \eqn{\hat{\gamma}_1}. See Dimitrova et al. (2023)
for further details on these parameters.}
\item{formula}{The model \code{\link[=formula.GeDS]{formula}}.}
\item{extcall}{\code{call} to the \code{\link{NGeDS}} functions.}
\item{terms}{\code{terms} object containing information on the model frame.}
}
}
\description{
\code{NGeDS} constructs a geometrically designed variable knot spline
regression model  referred to as a GeDS model, for a response having a normal
distribution.
}
\details{
The  \code{NGeDS} function implements the GeDS methodology, developed by
Kaishev et al. (2016) and extended in the \code{\link{GGeDS}} function for
the more general GNM (GLM) context, allowing for the response
to have any distribution from the exponential family. Under the GeDS approach
the (non-)linear predictor is viewed as a spline with variable knots which
are estimated along with the regression coefficients and the order of the
spline, using a two stage algorithm. In stage A, a linear variable-knot
spline is fitted to the data applying iteratively least squares  regression
(see \code{\link[stats]{lm}} function). In stage B, a Schoenberg variation
diminishing spline approximation to the fit from stage A is constructed, thus
simultaneously producing spline fits of order 2, 3 and 4, all of which are
included in the output (an object of class \code{"GeDS"}).

As noted in \code{\link[=formula.GeDS]{formula}}, the argument \code{formula}
allows the user to specify models with two components, a spline regression
(non-parametric) component involving part of the independent variables
identified through the function \code{f} and an optional  parametric
component involving the remaining independent variables. For \code{NGeDS} one
or two independent variables are allowed for the spline component and
arbitrary many independent variables for the parametric component. Failure to
specify the independent variable for the  spline regression component through
the function \code{f} will return an error. See
\code{\link[=formula.GeDS]{formula}}.

Within the argument \code{formula}, similarly as in other \R functions, it is
possible to specify one or more offset variables, i.e. known terms with fixed
regression coefficients equal to 1. These terms should be identified via the
function \code{\link[stats]{offset}}.

The parameter \code{beta} tunes the placement of a new knot in stage A of the
algorithm. Once a current second-order  spline is fitted to the data the
regression residuals are computed and grouped by their sign. A new knot is
placed at the residual cluster for which a certain measure
attains its maximum. The latter measure is defined as a weighted linear
combination of the range of each cluster and the mean of the absolute
residuals within it. The parameter \code{beta} determines the weights in this
measure correspondingly as \code{beta} and \code{1 - beta}. The  higher it
is, the more weight is put to the mean of the residuals and the less to the
range of their corresponding x-values. The default value of \code{beta} is
\code{0.5}.

The argument \code{stoptype} allows to choose between three alternative
stopping rules for the knot selection in stage A of GeDS: \code{"RD"},
that stands for \emph{Ratio of Deviances}, \code{"SR"}, that stands for
\emph{Smoothed Ratio} of deviances and \code{"LR"}, that stands for
\emph{Likelihood Ratio}. The latter is based on the difference of deviances
rather than on their ratio as in the case of \code{"RD"} and \code{"SR"}.
Therefore \code{"LR"} can be viewed as a log likelihood ratio test performed
at each iteration of the knot placement. In each of these cases the
corresponding stopping criterion is compared with a threshold value
\code{phi}.

The argument \code{phi} provides a threshold value required for the stopping
rule to exit the knot placement in stage A of GeDS. The higher the value of
\code{phi}, the more knots are added under the \code{"RD"} and \code{"SR"}
stopping rules contrary to the case of the stopping rule \code{"LR"} where
the lower \code{phi} is, more knots are included in the spline regression.
Further details for each of the three alternative stopping rules can be found
in Dimitrova et al. (2023).

The argument \code{q} is an input parameter that allows to fine-tune the
stopping rule in stage A. It identifies the number of consecutive iterations
over which the deviance should exhibit stable convergence so that the knot
placement in stage A is terminated. More precisely, under any of the rules
\code{"RD"}, \code{"SR"}, or \code{"LR"}, the deviance at the current
iteration is compared to the deviance computed \code{q} iterations before,
i.e., before selecting the last \code{q} knots. Setting a higher \code{q}
will lead to more knots being added before exiting stage A of GeDS.
}
\examples{
###################################################
# Generate a data sample for the response variable
# Y and the single covariate X
set.seed(123)
N <- 500
f_1 <- function(x) (10*x/(1+100*x^2))*4+4
X <- sort(runif(N, min = -2, max = 2))
# Specify a model for the mean of Y to include only a component
# non-linear in X, defined by the function f_1
means <- f_1(X)
# Add (Normal) noise to the mean of Y
Y <- rnorm(N, means, sd = 0.1)

# Fit a Normal GeDS regression using NGeDS
(Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2)))

# Apply some of the available methods, e.g.
# coefficients, knots and deviance extractions for the
# quadratic GeDS fit
# Note that the first call to the function knots returns
# also the left and right limits of the interval containing
# the data
coef(Gmod, n = 3)
confint(Gmod, n = 3)
knots(Gmod, n = 3)
knots(Gmod, n = 3, options = "internal")
deviance(Gmod, n = 3)

# Add a covariate, Z, that enters linearly
Z <- runif(N)
Y2 <- Y + 2*Z + 1
# Re-fit the data using NGeDS
(Gmod2 <- NGeDS(Y2 ~ f(X) + Z, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
coef(Gmod2, n = 3)
coef(Gmod2, onlySpline = FALSE, n = 3)

\dontrun{
##########################################
# Real data example
# See Kaishev et al. (2016), section 4.2
data('BaFe2As2')
(Gmod2 <- NGeDS(intensity ~ f(angle), data = BaFe2As2, beta = 0.6, phi = 0.99, q = 3))
plot(Gmod2)
}

#########################################
# bivariate example
# See Dimitrova et al. (2023), section 5

# Generate a data sample for the response variable
# Z and the covariates X and Y assuming Normal noise
set.seed(123)
doublesin <- function(x){
 sin(2*x[,1])*sin(2*x[,2])
}

X <- (round(runif(400, min = 0, max = 3),2))
Y <- (round(runif(400, min = 0, max = 3),2))
Z <- doublesin(cbind(X,Y))
Z <- Z+rnorm(400, 0, sd = 0.1)
# Fit a two dimensional GeDS model using NGeDS
(BivGeDS <- NGeDS(Z ~ f(X, Y), phi = 0.9))

# Extract quadratic coefficients/knots/deviance
coef(BivGeDS, n = 3)
confint(BivGeDS, n = 3)
knots(BivGeDS, n = 3)
deviance(BivGeDS, n = 3)

# Surface plot of the generating function (doublesin)
plot(BivGeDS, f = doublesin)
# Surface plot of the fitted model
plot(BivGeDS)

}
\references{
Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
Geometrically designed, variable knot regression splines.
\emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
DOI: \doi{10.1007/s00180-015-0621-7}

Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
Geometrically designed variable knot splines in generalized (non-)linear
models.
\emph{Applied Mathematics and Computation}, \strong{436}. \cr
DOI: \doi{10.1016/j.amc.2022.127493}
}
\seealso{
\link{GGeDS}; S3 methods such as \code{\link{coef.GeDS}},
\code{\link{confint.GeDS}}, \code{\link{deviance.GeDS}}, \code{\link{family}},
\code{\link{formula}}, \code{\link{knots.GeDS}}, \code{\link{lines.GeDS}},
\code{\link{logLik}}, \code{\link{plot.GeDS}}, \code{\link{predict.GeDS}},
\code{\link{print.GeDS}}, \code{\link{summary.GeDS}}; \link{Integrate} and
\link{Derive}; \link{PPolyRep}.
}
