################################################################################
################################################################################
############################## Univariate Fitters ##############################
################################################################################
################################################################################
#' @title Functions used to fit GeDS objects w/univariate spline regression
#' component
#' @name UnivariateFitters
#' @aliases Fitters UnivariateFitter GenUnivariateFitter
#' @description
#' These are computing engines called by \code{\link{NGeDS}} and
#' \code{\link{GGeDS}}, needed for the underlying fitting procedures.
#' @param X a numeric vector containing \eqn{N} sample values of the covariate
#' chosen to enter the spline regression component of the predictor model.
#' @param Y a vector of size \eqn{N} containing the observed values of the
#' response variable \eqn{y}.
#' @param Z a design matrix with \eqn{N} rows containing other covariates
#' selected to enter the parametric component of the predictor model (see
#' \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected, it
#' is set to \code{NULL} by default.
#' @param family a description of the error distribution and link function to be
#' used in the model. This can be a character string naming a family function
#' (e.g. \code{"gaussian"}), the family function itself (e.g.
#' \code{\link[stats]{gaussian}}) or the result of a call to a family function
#' (e.g. \code{gaussian()}). See \link[stats]{family} for details on family
#' functions.
#' @param weights an optional vector of size \eqn{N} of `prior weights' to be
#' put on the observations in the fitting process in case the user requires
#' weighted GeDS fitting. It is \code{NULL} by default.
#' @param beta numeric parameter in the interval \eqn{[0,1]} tuning the knot
#' placement in stage A of GeDS. See the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the
#' threshold for the stopping rule  (model selector) in stage A of GeDS. See
#' also \code{stoptype} and details in the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param min.intknots optional parameter allowing the user to set a minimum
#' number of internal knots required. By default equal to zero.
#' @param max.intknots optional parameter allowing the user to set a maximum
#' number of internal knots to be added by the GeDS estimation algorithm. By
#' default equal to the number of internal knots \eqn{\kappa} for the saturated
#' GeDS model (i.e. \eqn{\kappa=N-2}).
#' @param q numeric parameter which allows to fine-tune the stopping rule of
#' stage A of GeDS, by default equal to 2. See details in the description of
#' \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param extr numeric vector of 2 elements representing the left-most and
#' right-most limits of the interval embedding the sample values of \code{X}. By
#' default equal correspondingly to the smallest and largest values of \code{X}.
#' @param show.iters logical variable indicating whether or not to print 
#' information at each step. By default equal to \code{FALSE}.
#' @param stoptype a character string indicating the type of GeDS stopping rule
#' to be used. It should be either \code{"SR"}, \code{"RD"} or \code{"LR"},
#' partial match allowed. See details of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param offset a vector of size \eqn{N} that can be used to specify a fixed
#' covariate to be included in the predictor model  avoiding the estimation of
#' its corresponding regression coefficient. In case  more than one covariate is
#' fixed, the user should sum the corresponding coordinates of the fixed
#' covariates to produce one common \eqn{N}-vector of coordinates. The
#' \code{offset} argument is particularly useful when using 
#' \code{GenUnivariateFitter} if the link function used is not the identity.
#' @param tol numeric value indicating the tolerance to be used in the knot
#' placement steps in stage A. By default equal to 1E-12. See details below.
#' 
#' @return A \code{\link{GeDS-Class}} object, but without the \code{Formula},
#' \code{extcall}, \code{terms} and \code{znames} slots.
#' 
#' @details
#' The functions \code{UnivariateFitter} and \code{GenUnivariateFitter} are in
#' general not intended to be used directly, they should be called through
#' \code{\link{NGeDS}} and \code{\link{GGeDS}}. However, in case there is a need
#' for multiple GeDS fitting (as may be the case e.g. in Monte Carlo simulations)
#' it may be efficient to use the fitters outside the main functions.
#' 
#' The argument \code{tol} is used in the knot placement procedure of stage A of
#' the GeDS algorithm in order to check whether the current knot \eqn{\delta^*} 
#' is set at an acceptable location or not. If there exists a knot \eqn{\delta_i}
#' such that \eqn{|\delta^* - \delta_i| < }\code{tol}, \eqn{\delta^*}, then the
#' new knot is considered to be coalescent with an existing one, it is discarded
#' and the algorithm seeks alternative knot locations. By default it is equal to
#' 1e-12.
#'
#' See \code{\link{NGeDS}} and \code{\link{GGeDS}}, Kaishev et al. (2016) and
#' Dimitrova et al. (2023) for further details.
#'
#' @examples
#' # Examples similar to the ones
#' # presented in NGeDS and in GGeDS
#'
#' # Generate a data sample for the response variable
#' # Y and the covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression model using the fitter function
#' (Gmod <- UnivariateFitter(X, Y, beta = 0.6, phi = 0.995,
#'            extr = c(-2,2)))
#'
#' ##############################################################
#' # second: very similar example, but based on Poisson data
#' set.seed(123)
#' X <- sort(runif(N , min = -2, max = 2))
#' means <- exp(f_1(X))
#' Y <- rpois(N,means)
#' (Gmod2 <- GenUnivariateFitter(X, Y, beta = 0.2,
#'             phi = 0.995, family = poisson(), extr = c(-2,2)))
#'
#' # a plot showing quadratic and cubic fits,
#' # in the predictor scale
#' plot(X,log(Y), xlab = "x", ylab = expression(f[1](x)))
#' lines(Gmod2, n = 3, col = "red")
#' lines(Gmod2, n = 4, col = "blue", lty = 2)
#' legend("topleft", c("Quadratic","Cubic"),
#'      col = c("red","blue"), lty = c(1,2))
#'
#' @seealso  \code{\link{NGeDS}} and \code{\link{GGeDS}}.
#' 
#' @export
#' @rdname UnivariateFitters
#' 
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S., & Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \doi{10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}

UnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0,NROW(Y)),
                             weights = rep(1,length(X)), beta=.5, phi = 0.5,
                             min.intknots = 0, max.intknots = 300, q = 2,
                             extr = range(X), show.iters=FALSE,
                             tol = as.double(1e-12), stoptype = c("SR","RD","LR"))
  {
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "weights" = weights,
               "beta" = beta, "phi" = phi, "min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "extr" = extr, "tol" = tol)
  
  # Initialize RSS and phis
  RSSnew <- numeric()
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL;  oldslp <- NULL
  
  # Stop type, Indicator, distinctX
  stoptype <- match.arg(stoptype)
  Indicator <- table(X)
  distinctX <- unique(X)
 
  # Initialize knots and coefficients matrices and internal knots
  previous <- matrix(nrow = max.intknots + 1,     # maximum number of GeDS iterations (if max.intknots + 1 =< length(Y) - 2, max j = max.intknots + 1
                                                  #                                    o.w. max j = length(Y) - 2 < max.intknots + 1 )
                     ncol = max.intknots + 4)     # the n (= 2) left and right most knots are assumed to be coalescent (l + 4 knots)
  nz <- if(!is.null(Z)) NCOL(Z) else 0            # number of linear covariates
  oldcoef <- matrix(nrow = max.intknots + 1,
                    ncol = max.intknots + 2 + nz) # number of B-splines is p = l + 2
  
  intknots <- NULL
  
  
  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  for(j in 1:min(max.intknots + 1, length(Y) - 2)) {
    
    #############################################################
    ## STEP 1/STEP 8: Find the least squares linear spline fit ##
    #############################################################
    first.deg <- SplineReg_fast_weighted_zed(X = X, Y = Y, Z = Z, weights = weights, offset = offset,
                                           extr = extr, InterKnots = intknots, n = 2) #first regression
    # Store knots and coefficients
    previous[j, 1:(j+3)] <- sort(c(intknots, rep(extr, 2)))
    oldcoef[j, 1:(j+1+nz)] <- first.deg$Theta
    
    #####################################
    ## STEP 9: Store residuals and RSS ##
    ####################################
    res.tmp <- first.deg$Residuals
    RSSnew <- c(RSSnew, first.deg$RSS)
    
    ###########################
    ## STEP 10: Stopping Rule ##
    ############################
    # (I) Smoothed Ratio of deviances
    if(stoptype == "SR") {
      if(j > q  && length(phis) >= 3) {
        phis <- c(phis, RSSnew[j]/RSSnew[j-q])
        if(j - q > min.intknots){
          phismod <- log(1-phis)
          ccc <- .lm.fit(cbind(1, (q+1):j), phismod)$coef
          phis_star <- c(phis_star, 1-exp(ccc[1])*exp(ccc[2]*j))
          oldintc <- c(oldintc, ccc[1])
          oldslp <- c(oldslp, ccc[2])
          prnt <- paste0(", phi_hat = ",round(1-exp(ccc[1])*exp(ccc[2]*j),3))

          if(1-exp(ccc[1])*exp(ccc[2]*j) >= phi) {
            break
          }
        }
      }
    }
    # (II) Ratio of Deviances
    if(stoptype == "RD" | (stoptype == "SR" & length(phis) < 3)) {
      phis <- c(phis, RSSnew[j]/RSSnew[j-q])
      if(j > q && (j-q > min.intknots)) {
        prnt <- paste0(", phi = ", round(RSSnew[j]/RSSnew[j-q],3))
        if(RSSnew[j] / RSSnew[j-q] >= phi) {
          break
        }
      }
    }
    # (III) Likelihood Ratio
    if(stoptype == "LR") {
      phis <- c(phis,RSSnew[j-q]-RSSnew[j])
      if(j > q && (j-q > min.intknots)) {
        prnt <- paste0(", p = ",round(pchisq(-(RSSnew[j]-RSSnew[j-q]), df = q),3))
        if(-(RSSnew[j]-RSSnew[j-q]) < qchisq(phi,df=q)) {     
          break
        }
      }
    }
    
    ############
    ## STEP 2 ##
    ############
    d <- numeric()
    if(any(Indicator > 1)){
      # Average res.tmp*weights for repeated values of X
      res.weighted <- makeNewRes(resold = res.tmp*weights, recurr = as.numeric(Indicator)) # already fast - useless to do in C++
    } else {
      res.weighted <- res.tmp*weights
    }
    # Group the consecutive residuals into clusters by their sign
    signs <- sign(res.weighted)
    for(i in 1:length(distinctX)) {
      # If all residual values have the same sign, count the entire set as one cluster
      if (all(signs == signs[1])) {
        d[i] <- length(signs)
        break
      } else {
        # If signs change, identify the first change to split the cluster
        d[i] <- min(which(signs!=signs[1])) - 1 # number of consecutive residuals with same sign
        # Update signs to exclude the identified cluster
        signs <- signs[-(1:d[i])]               # extract cluster from signs
      }
    }
    
    ####################################################################
    ## STEP 3: within residual cluster means +  within-cluster ranges ##
    ####################################################################
    u <-length(d)
    dcum <- cumsum(d)
    # initialize means and wc.range
    means <- wc.range <- numeric(u)
    means[1] <- abs(mean(res.weighted[1:dcum[1]]))
    wc.range[1] <- distinctX[dcum[1]] - X[1]
    for (i in 2:u) { # embed in C++ useless
      means[i] <- abs(mean((res.weighted[(dcum[i-1] + 1):dcum[i]])))
      wc.range[i] <- distinctX[dcum[i]] - distinctX[dcum[i-1] + 1]
    }
    
    ######################################################################
    ## STEP 4: calculate the normalized within-cluster means and ranges ##
    ######################################################################
    means <- means/max(means)
    wc.range <- wc.range/max(wc.range)
    
    ###########################################
    ## STEP 5: calculate the cluster weights ##
    ###########################################
    w <- beta*means + (1 - beta)*wc.range
    
    ###############
    ## STEPS 6/7 ##
    ###############
    newknot <- Knotnew(wht = w, restmp = res.weighted, x = distinctX, dcm = dcum,
                       oldknots = c(rep(extr, 2 + 1), intknots), tol = tol)[1]
    intknots <- c(intknots, newknot)
    
    # Print iteration
    if(show.iters) {
      indent <- rep(" ", nchar(options()$prompt))
      indent <- paste(indent, collapse="")
      if (j > q) {
        toprint <- paste0(indent, "Iteration ",j,": New Knot = ", round(newknot, 3),
                          ", RSS = " , round(RSSnew[j],3), prnt,"\n")
      } else {
          toprint <- paste0(indent,"Iteration ",j,": New Knot = ", round(newknot, 3),
                            ", RSS = " , round(RSSnew[j],3), "\n")
      }
      cat(toprint)
    }
    
  }
  
  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  if (j == max.intknots + 1) warning("Maximum number of iterations exceeded")
  if (j <= max.intknots) {
    # Eliminate NAs in previous
    previous <- previous[-((j+1):(max.intknots+1)), ] # eliminate all the NA rows
    previous <- previous[ ,-((j+4):(max.intknots+4))] #  k = j - 1 internal knots (i.e. k + 4 = j + 3 knots) in the last iteration
    # Eliminate NAs in oldcoef
    oldcoef  <- oldcoef[-((j+1):(max.intknots+1)), ]  # eliminate all the NA rows
    oldcoef  <- oldcoef[ ,1:(j+1+nz)]                 # p = n + k = j + 1  B-splines
  }
  
  # 1. LINEAR
  if(j - q < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = ll, n = 2)
    } else {
      ik <- previous[j-q, -c(1, 2, (j + 2 - q):(j + 3))] # keep the internal knots in the "j-q"th row
      # Stage B.1 (averaging knot location)
      ll <- makenewknots(ik, 2)
      # Stage B.2
      lin <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = ll, n = 2)
    }
  # 2. QUADRATIC
  if(j - q < 3) {
    warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qq <- NULL
    squ <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = qq, n = 3)
    } else {
      # Stage B.1 (averaging knot location)
      qq <- makenewknots(ik, 3)
      # Stage B.2
      squ <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = qq, n = 3)
    }
  # 3. CUBIC
  if(j-q < 4) {
    warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    cc <- NULL
    cub <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = cc, n = 4)
    } else {
      # Stage B.1 (averaging knot location)
      cc <- makenewknots(ik, 4)
      # Stage B.2
      cub <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr, InterKnots = cc, n = 4)
    }
  
  out <- list("Type" = "LM - Univ", "Linear.Knots" = ll, "Quadratic.Knots" = qq, "Cubic.Knots" = cc,
              "Dev.Linear" = lin$RSS, "Dev.Quadratic" = squ$RSS, "Dev.Cubic" = cub$RSS,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previous,
              "Args"= args, "Call"= save, "Nintknots" = j - q - 1, "iters" = j, "Guesses" = NULL,
              "Coefficients" = oldcoef, stopinfo = list("phis" = phis, "phis_star" = phis_star, "oldintc" = oldintc, "oldslp" = oldslp))

  class(out) <- "GeDS"
  return(out)
}


################################################################################
############################# GenUnivariateFitter ##############################
################################################################################
#' @rdname UnivariateFitters
#' @export
GenUnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0, NROW(Y)),
                                weights = rep(1, length(X)), family = gaussian(),
                                beta = 0.5, phi = 0.5, min.intknots = 0,
                                max.intknots = 300, q = 2, extr = range(X),
                                show.iters = F, tol = as.double(1e-12),
                                stoptype = c("SR","RD","LR"))
{
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "weights"=weights,
               "beta" = beta, "phi" = phi, "family"=family, "min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "extr" = extr)
  
  # Initialize RSS and phis
  RSSnew <- numeric()
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  
  # Initialize knots and coefficients matrices and internal knots
  previous <- matrix(nrow = max.intknots + 1,     # maximum number of GeDS iterations (if max.intknots + 1 =< length(Y) - 2, max j = max.intknots + 1
                                                  # o.w. max j = length(Y) - 2 < max.intknots + 1 )
                     ncol = max.intknots + 4)     # the n (= 2) left and right most knots are assumed to be coalescent (l + 4 knots)
  nz <- if(!is.null(Z)) NCOL(Z) else 0            # number of linear covariates
  oldcoef <- matrix(nrow = max.intknots + 1,
                    ncol = max.intknots + 2 + nz) # number of B-splines is p = l + 2
  intknots <- NULL
  # Initial values for the coefficients used at each iteration of stage A in order to estimate the spline coefficients
  oldguess <- matrix(nrow = max.intknots + 1,
                     ncol = max.intknots + 2)    # number of B-splines is p = l + 2 (max number of coef)
  
  # Accumulated number of IRLS iterations/ Vector of IRLS deviances
  irlsAccumIterCount <- devianceTracking <- NULL
  # Control basis matrix singularity
  flag <- FALSE
  
  # Stop type, Indicator, distinctX
  stoptype <- match.arg(stoptype)
  Indicator <- table(X)
  distinctX <- unique(X)
  
  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  for(j in 1:min(max.intknots + 1, length(Y) - 2)) {
    
    if(flag) j <- j - 2
    
    # Sort internal knots and update the oldguess matrix; create new guess vector
    if(j > 1) {
      intknots <- sort(intknots)
      oldguess[j, 1:(j+1)] <- guess
      guess <- c(guess, guess_z)
    } else if (j == 1) {
      guess <- NULL
    }
    
    #############################################
    ## STEP 1: Find the IRLS linear spline fit ##
    #############################################
    # Linear spline regression using specified parameters
    first.deg <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                               InterKnots = intknots, n = 2, extr = extr,
                               family = family, inits = guess)
    
    # 1. Check for NA values in the Theta vector to handle potential singularities
    if(anyNA(first.deg$Theta)) {
      # Calculate the rank of the basis matrix and check for singularities
      rank.basis <- Matrix::rankMatrix(first.deg$Basis)
      cols <- NCOL(first.deg$Basis)
      
      # (i) Handle the case when the basis matrix is singular
      if(rank.basis < cols) {
        # (a) If Basis was singular for second consecutive time, break loop
        if(flag) {
          warning("Matrix singular for the second time. Breaking the loop.")
          break
        }
        # (b) Identify and remove problematic knots and guesses based on NA positions
        check <- which(is.na(first.deg$Theta))
        intknots <- intknots[-(check+2)]
        guess <- guess[1:length(first.deg$Theta)][-check]
        toprint <- paste0("Basis Matrix singular, deleting one knot")
        print(toprint)
        flag  <- T 
        # (c) Check if the number of knots equals the number of unique X values and issue a warning if true
        if(cols == length(distinctX)) {
          warning("Number of knots equal to number of unique Xs. Breaking the loop.")
          break
          } else {
            # Continue to the next iteration otherwise
            next  
          }
        
      # (ii) NA values in the Theta vector, but basis matrix is not singular (i.e. other issues)
        } else {
          stop("NA(s) in the coefficients")
        }
      
    # 2. If no NAs, update guess (coefficients initial value in the next iteration)
    } else {
      guess <- first.deg$Theta[1:(j+1)]
    }
    
    # Vector accumulating the IRLS iterations deviances obtained at each GeDS iteration
    devianceTracking <- c(devianceTracking, first.deg$deviance)
    # Accumulated number of IRLS iterations at each GeDS iteration
    irlsAccumIterCount <- c(irlsAccumIterCount, length(devianceTracking))
    
    # Store knots and coefficients
    previous[j,1:(j+3)] <- sort(c(intknots,rep(extr,2)))
    oldcoef[j,1:(j+1+nz)] <- first.deg$Theta
    guess_z <- if(nz > 0) first.deg$Theta[(j+2):(j+1+nz)] else NULL
    
    # Store residuals and deviance 
    res.tmp <- first.deg$Residuals
    RSS.tmp <- first.deg$temporary$lastdev
    RSSnew <- c(RSSnew, RSS.tmp)
    # Working weights (weights in the final iteration of the IRLS fit)
    working.weights <- first.deg$temporary$weights  
    
    ###########################
    ## STEP 2: Stopping Rule ##
    ###########################
    # (I) Smoothed Ratio of deviances
    if(stoptype =="SR") {
      if(j > q  && length(phis) >= 3) {
        phis <- c(phis, RSSnew[j]/RSSnew[j-q])
        if (j - q > min.intknots) {
          phismod <- log(1-phis)
          ccc <- .lm.fit(cbind(1,(q+1):j), phismod)$coef
          phis_star <- c(phis_star, 1-exp(ccc[1])*exp(ccc[2]*j))
          oldintc <- c(oldintc, ccc[1])
          oldslp <- c(oldslp, ccc[2])
          prnt <- paste0(", phi_hat = ", round(1-exp(ccc[1])*exp(ccc[2]*j),3))
          if(1-exp(ccc[1])*exp(ccc[2]*j) >= phi) {
            break
          }
        }
      }
    }
    # (II) Ratio of Deviances
    if(stoptype == "RD" | (stoptype == "SR" & length(phis) < 3)) {
      phis <- c(phis, RSSnew[j]/RSSnew[j-q])
      if(j > q && (j - q > min.intknots)) {
        prnt <- paste0(", phi = ", round(RSSnew[j]/RSSnew[j-q],3))
        if(RSSnew[j]/RSSnew[j-q] >= phi) {
          break
        }
      }
    }
    # (III) Likelihood Ratio
    if(stoptype=="LR"){
      if(j > q && (j-q > min.intknots)) {
        phis <- c(phis, RSSnew[j-q]-RSSnew[j])
        prnt <- paste0(", p = ", round(pchisq(-(RSSnew[j]-RSSnew[j-q]), df=q),3))
        if(-(RSSnew[j]-RSSnew[j-q]) < qchisq(phi,df=q)) {
          break
        }
      }
    }
    
    ###############################################
    ## STEPS 3-8 (i.e. steps 2-7 of Normal GeDS) ##
    ###############################################
    d <- numeric()
    if(any(Indicator > 1)) {
      # Average res.tmp*weights for repeated values of X
      res.weighted <- makeNewRes2(resold = res.tmp, weights = weights*working.weights, recurr = as.numeric(Indicator)) #already fast - useless to do in c++
    } else {
      res.weighted <- res.tmp*weights*working.weights
    }
    # Group the consecutive residuals into clusters by their sign
    signs <- sign(res.weighted)
    for (i in 1:length(distinctX)) {
      # If all residual values have the same sign, count the entire set as one cluster
      if (all(signs==signs[1])) {
        d[i]<-length(signs)
        break
      } else {
        # If signs change, identify the first change to split the cluster
        d[i] <- min(which(signs!=signs[1]) - 1) # number of consecutive residuals with same sign
        # Update signs to exclude the identified cluster
        signs <- signs[-(1:d[i])]               # extract cluster from signs
      }
    }
    ####################################################################
    ## STEP 4: within residual cluster means +  within-cluster ranges ##
    ####################################################################
    u <- length(d)
    dcum <- cumsum(d)
    # initialize means and wc.range
    means <- wc.range <- numeric(u)
    means[1] <- abs(mean(res.weighted[1:dcum[1]]))
    wc.range[1] <- distinctX[dcum[1]]-X[1]
    for(i in 2:u){ # embed in C++ useless
      means[i] <- abs(mean((res.weighted[(dcum[i-1]+1):dcum[i]])))
      wc.range[i] <- distinctX[dcum[i]]-distinctX[dcum[i-1]+1]
    }
    ######################################################################
    ## STEP 5: calculate the normalized within-cluster means and ranges ##
    ######################################################################
    means    <- means/max(means)
    wc.range <- wc.range/max(wc.range)
    
    ###########################################
    ## STEP 6: calculate the cluster weights ##
    ###########################################
    w <- beta*means + (1-beta)*wc.range+1e-8
    
    ###############
    ## STEPS 7/8 ##
    ###############
    # i. Calculate new knot
    newknot <- Knotnew(wht = w, restmp = res.weighted, x = distinctX, dcm = dcum,
                       oldknots = c(rep(extr,2 + 1), intknots), tol = tol)[1]
    
    # ii. Calculate guess-coefficient for newknot
    guess <- newknot.guess(intknots, extr, guess, newknot)
    
    # iii. Update internal knots vector
    intknots <- c(intknots, newknot)
    
    # Print iteration
    if (show.iters) {
      indent <- rep(" ", nchar(options()$prompt))
      indent <- paste(indent, collapse="")
      if(j > q) {
        toprint <- paste0(indent, "Iteration ",j,": New Knot = ", round(newknot,3),
                          ", DEV = " , round(RSSnew[j],3), prnt,"\n")
      } else {
        toprint <- paste0(indent, "Iteration ", j, ": New Knot = ", round(newknot,3),
                          ", DEV = " , round(RSSnew[j],3), "\n")
      }
      cat(toprint)
    }
    
  }
  
  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  if (j == max.intknots + 1) warning("Maximum number of iterations exceeded")
  if (j <= max.intknots) {
    # Eliminate NAs in previous
    previous <- previous[-((j+1):(max.intknots+1)), ] # eliminate all the NA rows
    previous <- previous[ ,-((j+4):(max.intknots+4))] #  k = j - 1 internal knots (i.e. k + 4 = j + 3 knots) in the last iteration
    # Eliminate NAs in oldcoef
    oldcoef  <- oldcoef[-((j+1):(max.intknots+1)), ]  # eliminate all the NA rows
    oldcoef  <- oldcoef[ ,1:(j+1+nz)]                 # p = n + k = j + 1  B-splines
  }
  
  # 1. LINEAR
  if (j - q < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = ll, n = 2, family = family,
                         inits = c(oldguess[j-q, 1:(j-q+1)], guess_z))
  } else {
    ik <- previous[j-q,-c(1,2,(j+2-q):(j+3))]
    # Stage B.1 (averaging knot location)
    ll <- makenewknots(ik, 2)
    # Stage B.2
    lin <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = ll, n = 2, family = family,
                         inits = c(oldguess[j-q, 1:(j-q+1)], guess_z))
  }
  # 2. QUADRATIC
  if (j - q < 3) {
    warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qq <- NULL
    squ <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = qq, n = 3, family = family,
                         mustart = lin$Predicted)
  } else {
    # Stage B.1 (averaging knot location)
    qq <- makenewknots(ik, 3)
    # Stage B.2
    squ <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = qq, n = 3, family = family,
                         mustart = lin$Predicted)
  }
  # 3. CUBIC
  if (j - q < 4) {
    warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    cc <- NULL
    cub <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = cc, n = 4, family = family,
                         mustart = squ$Predicted)
  } else {
    # Stage B.1 (averaging knot location)
    cc <- makenewknots(ik, 4)
    # Stage B.2
    cub <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = cc, n = 4, family = family,
                         mustart = squ$Predicted)
  }
  
  out <- list("Type" = "GLM - Univ", "Linear.Knots" = ll, "Quadratic.Knots" = qq, "Cubic.Knots" = cc,
              "Dev.Linear" = lin$RSS, "Dev.Quadratic" = squ$RSS, "Dev.Cubic" = cub$RSS, "Knots" = intknots,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previous,
              "Args" = args, "Call" = save, "Nintknots" = j - q - 1, "iters" = j, "Guesses" = oldguess,
              "Coefficients" = oldcoef, "deviance" = devianceTracking, "iterIrls" = irlsAccumIterCount,
              stopinfo = list("phis" = phis,"phis_star" = phis_star, "oldintc" = oldintc, "oldslp" = oldslp))
  
  class(out) <- "GeDS"
  return(out)
}

