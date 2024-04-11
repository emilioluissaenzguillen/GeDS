################################################################################
################################################################################
################################### SplineReg ##################################
################################################################################
################################################################################
#' @title Estimation of the coefficients of a predictor model with spline and
#' possibly parametric components.
#' @name SplineReg
#' @description
#' Functions that estimate the coefficients of a predictor model involving a
#' spline component and possibly a parametric component applying (Iteratively
#' Re-weighted) Least Squares (IR)LS iteration.
#' @param X a numeric vector containing \eqn{N} sample values of the covariate
#' chosen to enter the spline regression component of the predictor model.
#' @param Y a vector of size \eqn{N} containing the observed values of the
#' response variable \eqn{y}.
#' @param Z a design matrix with \eqn{N} rows containing other covariates
#' selected to enter the parametric component of the predictor model (see
#' \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected,
#' it is set to \code{NULL} by default.
#' @param offset a vector of size \eqn{N} that can be used to specify a fixed
#' covariate to be included in the predictor model avoiding the estimation of
#' its corresponding regression coefficient. In case  more than one covariate is
#' fixed, the user should sum the corresponding coordinates of the fixed
#' covariates to produce one common \eqn{N}-vector of coordinates. The argument
#' \code{offset} is particularly useful in \code{Splinereg_GLM} if the link
#' function used is not the identity.
#' @param weights an optional vector of `prior weights' to be put on the
#' observations in the fitting process in case the user requires weighted
#' fitting. It is a vector of 1s by default.
#' @param InterKnots a numeric vector containing the locations of the internal
#' knots necessary to compute the B-splines. In GeDS these are the internal
#' knots in a current iteration of stage A.
#' @param n integer value specifying  the order of the spline to be evaluated.
#' It should be 2 (linear spline), 3 (quadratic spline) or 4 (cubic spline).
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param extr optional numeric vector of 2 elements representing the left-most
#' and right-most limits of the interval embedding the sample values of
#' \code{X}. By default equal correspondingly to the smallest and largest values
#' of \code{X}.
#' @param family a description of the error distribution and link function to be
#' used in the model. This can be a character string naming a family function, a
#' family function or the result of a call to a family function. See
#' \code{\link[stats]{family}} for details of family functions.
#' @param inits a numeric vector of length
#' \code{length(InterKnots) + n + NCOL(Z)} providing initial values for the
#' coefficients, to be used in the IRLS estimation (alternative to providing the
#' \code{mustart} vector).
#' @param mustart initial values for the vector of means in the IRLS estimation.
#' Must be a vector of length \eqn{N}.
#' @param etastart initial values for the predictor in the IRLS estimation
#' (alternative to providing either \code{inits} or \code{mustart}). Must be a
#' vector of length \eqn{N}.
#' @param prob the confidence level to be used for the confidence bands in the
#' \code{SplineReg_LM} fit. See details below.
#' @param coefficients optional vector of spline coefficients. If provided,
#' \code{SplineReg} computes only the corresponding predicted values.
#' 
#' @return A \code{list} containing:
#' \item{Theta}{ a vector containing the fitted coefficients of the spline
#' regression component and the parametric component of the predictor model.}
#' \item{Predicted}{ a vector of \eqn{N} predicted mean values of the response
#' variable computed at the sample values of the covariate(s).}
#' \item{Residuals}{ a vector containing the normal regression residuals if
#' \code{SplineReg_LM} is called or the residuals described in Details if
#' \code{SplineReg_GLM} is called.}
#' \item{RSS}{ the deviance for the fitted predictor model, defined as in
#' Dimitrova et al. (2023), which for \code{SplineReg_LM} coincides with the
#' Residual Sum of Squares.}
#' \item{NCI}{ a list containing the lower (\code{Low}) and upper (\code{Upp})
#' limits of the approximate confidence intervals computed at the sample values
#' of the covariate(s). See  details above.}
#' \item{Basis}{ the matrix of B-spline regression functions and the covariates
#' of the parametric part evaluated at the sample values of the covariate(s).}
#' \item{Polygon}{ a list containing x-y coordinates ("\code{Kn}" and
#' "\code{Thetas}") of the vertices of the Control Polygon, see
#' Dimitrova et al. (2023).}
#' \item{deviance}{ a vector containing deviances computed at each IRLS step
#' (computed only with the \code{SplineReg_GLM}).}
#' \item{temporary}{ the result of the function \code{\link[stats]{lm}} if
#' \code{SplineReg_LM} is used or the output of the function
#' \code{\link{IRLSfit}} (which is similar to the output from
#' \code{\link[stats]{glm.fit}}), if \code{SplineReg_GLM} is used.}
#'
#' @details
#' The functions estimate the coefficients of a predictor model with a spline
#' component (and possibly a parametric component) for a given, fixed order and
#' vector of knots of the spline and a specified distribution of the response
#' variable (from the Exponential Family). The functions \code{SplineReg_LM} and
#' \code{SplineReg_GLM} are based correspondingly on LS and IRLS and used
#' correspondingly in \code{\link{NGeDS}} and \code{\link{GGeDS}}, to estimate
#' the coefficients of the final GeDS fits of stage B, after their knots have
#' been positioned to coincide with the Greville abscissas of the knots of the
#' linear fit from stage A (see Dimitrova et al. 2023). Additional inference
#' related quantities are also computed (see Value below). The function
#' \code{SplineReg_GLM} is also used to estimate the coefficients of the linear
#' GeDS fit of stage A within \code{\link{GGeDS}}, whereas in
#' \code{\link{NGeDS}} this estimation is performed internally leading to faster
#' R code.
#'
#' In addition \code{SplineReg_LM} computes some useful quantities among which
#' confidence intervals and the Control Polygon (see Section 2 of
#' Kaishev et al. 2016).
#'
#' The confidence intervals contained in the output slot \code{NCI} are
#' approximate local bands obtained considering the knots as fixed constants.
#' Hence the columns of the design matrix are seen as covariates and standard
#' methodology relying on the \code{se.fit} option of \code{predict.lm} or
#' \code{predict.glm} is used. In the \code{ACI} slot, asymptotic confidence
#' intervals are provided, following Kaishev et al (2006).
#'
#' As mentioned, \code{SplineReg_GLM} is intensively used in Stage A of the GeDS
#' algorithm implemented in \code{\link{GGeDS}} and in order to make it as fast
#' as possible input data validation is mild. Hence it is expected that the user
#' checks carefully the input parameters before using \code{SplineReg_GLM}. The
#' "\code{Residuals}" in the output of this function are similar to the so
#' called ``working residuals" in the \code{\link[stats]{glm}} function. 
#' "\code{Residuals}"  are the residuals \eqn{r_i} used in the knot placement
#' procedure, i.e. \deqn{r_i= (y_i - \hat{\mu}_i){d \mu_i \over d \eta_i },} but
#' in contrast to \code{\link[stats]{glm}} ``working residuals", they are
#' computed using the final IRLS fitted \eqn{\hat{\mu}_i}. "\code{Residuals}"
#' are then used in locating the knots of the linear spline fit of Stage A.
#' 
#' In \code{SplineReg_GLM} confidence intervals are not computed.
#'
#' @seealso \code{\link{NGeDS}}, \code{\link{GGeDS}}, \code{\link{Fitters}},
#' \code{\link{IRLSfit}}, \code{\link[stats]{lm}} and
#' \code{\link[stats]{glm.fit}}.
#' 
#' @export 
#' 
#' @references
#' Kaishev, V. K., Dimitrova, D. S., Haberman, S. & Verrall, R. J. (2006).
#' Geometrically designed, variable know regression splines: asymptotics and
#' inference \emph{(Statistical Research Paper No. 28)}.
#' London, UK: Faculty of Actuarial Science & Insurance, City University London. \cr
#' URL: \href{https://openaccess.city.ac.uk/id/eprint/2372}{openaccess.city.ac.uk}
#'
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S., & Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \doi{10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models. \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}
#'
#' @aliases SplineReg
#' @rdname SplineReg

SplineReg_LM <- function(X, Y, Z = NULL, offset = rep(0,NROW(X)), weights = rep(1,length(Y)),
                         InterKnots, n, extr = range(X), prob = 0.95,
                         coefficients = NULL)
  {
  # Convert spline order to integer
  n <- as.integer(n)
  # Create spline basis matrix using specified knots, evaluation points and order
  matrice <- splineDesign(knots = sort(c(InterKnots,rep(extr,n))), x = X, ord = n,
                          derivs = rep(0,length(X)), outer.ok = T)
  # Combine spline basis with parametric design matrix (if provided)
  matrice2 <- cbind(matrice,Z)
  
  # 1) If coefficients are NOT provided estimate the corresponding regression model
  if (is.null(coefficients)) {
    # Substract offset (if any) from Y
    Y0 <- Y - offset
    # Fit linear model without intercept, using weights
    tmp <- lm(Y0 ~ -1 + matrice2, weights = as.numeric(weights))
    # Extract fitted coefficients
    theta <- coef(tmp)
    # Compute predicted values
    predicted <- matrice2 %*% theta + offset
    # Knots for control polygon
    nodes <- sort(c(InterKnots,rep(extr,n)))[-c(1,NCOL(matrice)+1)]
    polyknots <- makenewknots(nodes, degree = n)
    # Calculate residuals
    resid <- Y - predicted
    # Residual standard error
    df <- tmp$df
    sigma_hat <- sqrt(sum(resid^2)/df)
    # Recalculate the probability for a two-tailed test
    prob <- 1-.5*(1-prob)
    #  CI_j =\hat{y_j} ± t_{α/2,df}*\hat{σ}*\sqrt{H_{jj}}; H = X(X'X)^{−1}X'
    band <- qt(prob,df) * sigma_hat * influence(tmp)$hat^.5
    
    # Huang's method for confidence band width (alternative approach)
    n <- length(Y)
    N <- NCOL(matrice)
    matcb <- matrix(0,N,N)
    # to be written in cpp
    for(i in 1:n){
      matcb <- matcb + matrice[i,] %*% t(matrice[i,])
    }
    matcb <- matcb/n
    matcbinv <- solve(matcb)
    band_width_huang <- qnorm(prob) * n^(-.5) * diag(sigma_hat*matrice%*%matcbinv%*%t(matrice))
  
  # 2) If coefficients are provided, use them to compute predicted values directly
  } else {
    theta <- coefficients
    predicted <- matrice2 %*% theta + offset
    # Set these variables to NA/NULL as they are not computed in this branch
    polyknots <- resid <- band <- band_width_huang <- NA; tmp <- NULL
  }
  
  out <- list("Theta" = theta, "Predicted" = predicted, "Residuals" = resid, "RSS" = t(resid)%*%resid,
              "NCI" = list("Upp" = predicted + band, "Low" = predicted - band),
              "Basis" = matrice, "Polygon" = list("Kn" = polyknots,
                                                  "Thetas" = theta[1:NCOL(matrice)]),
              "temporary" = tmp, "ACI"=list("Upp" = predicted + band_width_huang,
                                            "Low" = predicted - band_width_huang))
  return(out)
}

 
#' @export
#' @rdname SplineReg
 
SplineReg_GLM <- function(X, Y, Z, offset = rep(0,nobs), weights = rep(1,length(Y)),
                          InterKnots, n, extr = range(X), family, mustart,
                          inits = NULL, etastart = NULL)
  {
  # Ensure X, Y, Z, InterKnots are numeric, and n and extr are integers
  X <- as.numeric(X); Y <- as.numeric(Y); Z <- as.numeric(Z)
  InterKnots <- as.numeric(InterKnots); n <- as.integer(n); extr <- as.numeric(extr)
  
  # Check that 'n' (spline order) has length 1
  if(length(n) != 1) stop("'n' must have length 1")
  ord <- n # to avoid problem in use of family$initialize e.g. binomial()
  # Check X and Y have the same length
  if(length(X) != length(Y)) stop("length of 'X' must match length of 'Y'")
  # Ensure weights are set
  if(is.null(weights)) weights <- rep(1,length(Y))
  # Set required environment variables for family$initialize and IRLSfit
  y <- Y; nobs <- NROW(Y)
  
  # Create spline basis matrix using specified knots, order, and evaluation points
  matrice <- splineDesign(knots = sort(c(InterKnots,rep(extr,n))), x = X, ord = n,
                          derivs = rep(0,length(X)), outer.ok = T)
  # Combine spline basis with parametric design matrix (if provided)
  matrice2 <- cbind(matrice,Z)
  
  # Initialize mustart based on input or defaults
  if (missing(mustart) || is.null(mustart)) {
    
    if (is.null(inits)) {
      # Set environment to parent frame
      env <- parent.frame()
      eval(family$initialize)
      mustart <- env$mustart
    } else {
      # Validate length of 'inits'
      if (length(inits)!= NCOL(matrice2)) stop("'inits' must be of length length(InterKnots) + n + NCOL(Z)")
      # Calculate initial mustart based on 'inits' (initial value for spline coefficients)
      mustart <- family$linkinv(matrice2 %*% inits)
    }
  }
  
  tmp <- IRLSfit(matrice2, Y, offset = offset,
                 family = family, mustart = mustart, weights = weights)
  # Extract fitted coefficients
  theta <- coef(tmp)
  # Compute predicted mean values of the response variable
  predicted <- family$linkinv(matrice2%*%theta + offset)
  # Knots for control polygon
  nodes <- sort(c(InterKnots,rep(extr,ord)))[-c(1,NCOL(matrice)+1)]
  polyknots <- makenewknots(nodes, degree = ord)
  # Extract residuals
  resid <- tmp$res2
  
  out <- list("Theta" = theta, "Predicted" = predicted, "Residuals" = resid,
              "RSS" = tmp$lastdeviance, "Basis" = matrice, 
              "Polygon" = list("Kn" = polyknots,
                               "Thetas" = theta[1:NCOL(matrice)]),
              "temporary" = tmp, "deviance" = tmp$deviance)
  return(out)
}


# this is an unexported function.
# it may be useful in case knots are seleced according to some other algorithms (e.g. package earth?)
# and then one want to use GeDS to get higher order spline fits
SplineReg_GLM2 <- function(X,Y,Z,offset=rep(0,nobs),
                          weights=rep(1,length(X)),InterKnots,n,extr=range(X),family,
                          inits = NULL,mustart){
  X <- as.numeric(X)
  Y <- as.numeric(Y)
  #Z <- as.numeric(Z)
  nobs <- NROW(Y)
  y <- Y
  if(length(X) != length(Y)) stop("length of 'X' must match length of 'Y'")
  InterKnots <- as.numeric(InterKnots)
  if(length(n) != 1) stop("'n' must have length 1")
  n <- as.integer(n)
  extr <- as.numeric(extr)
  matrice <- splineDesign(knots=sort(c(InterKnots,rep(extr,n))),
                           derivs=rep(0,length(X)),x=X,ord=n,outer.ok = T)
  matrice2 <- cbind(matrice,Z)
  ord <- n #may cause problem the use of family$initialize e.g. binomial()
  if(missing(mustart)||is.null(mustart)){

    env <- parent.frame()
    if(is.null(weights)) weights <- rep(1,length(y))
    if (is.null(inits)) {
      eval(family$initialize)

      mustart <- env$mustart
    } else {
      if(length(inits)!= NCOL(matrice2)) stop("'inits' must be of length length(InterKnots) + n + NCOL(Z)")
      mustart <- family$linkinv(matrice2%*%inits)
    }
  }

  tmp <- glm.fit(matrice2, Y,
                 family=family)
  theta <- coef(tmp)
  nodes<-sort(c(InterKnots,rep(extr,ord)))[-c(1,NCOL(matrice)+1)]
  polyknots <- makenewknots(nodes,degree=ord)
  predicted <- matrice2%*%theta+offset
  resid <- tmp$res2
  out <- list("Theta"=theta,"Predicted"=predicted,
              "Residuals"=resid,"RSS"=tmp$lastdeviance,
              "Basis"= matrice,"Polygon"=list("Kn"=polyknots,"Thetas"=theta[1:NCOL(matrice)]),
              "temporary"=tmp, "deviance"= tmp$deviance)
  return(out)
}

