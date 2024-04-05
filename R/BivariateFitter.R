################################################################################
################################################################################
############################## Bivariate Fitters ###############################
################################################################################
################################################################################
#' @title Fitter function for GeD Spline Regression for bivariate data
#' @name BivariateFitters
#' @aliases BivariateFitters BivariateFitter
#' @description
#' These are computing engines called by \code{\link{NGeDS}}, needed for the
#' underlying fitting procedures.
#' @param X a numeric vector containing \eqn{N} sample values of the first
#' independent variable chosen to enter the spline regression component of the
#' predictor model.
#' @param Y  a numeric vector containing \eqn{N} sample values of the second
#' independent variable chosen to enter the spline regression component of the
#' predictor model.
#' @param Z a vector of size \eqn{N} containing the observed values of the
#' response variable.
#' @param W a design matrix with \eqn{N} rows containing other covariates
#' selected to enter the parametric component of the predictor model (see
#' \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected, it
#' is set to \code{NULL} by default.
#' @param weights an optional vector of size \eqn{N} of `prior weights' to be
#' put on the observations in the fitting process in case the user requires
#' weighted GeDS fitting. It is \code{NULL} by default.
#' @param family a description of the error distribution and link function to be
#' used in the model. This can be a character string naming a family function
#' (e.g. \code{"gaussian"}), the family function itself (e.g.
#' \code{\link[stats]{gaussian}}) or the result of a call to a family function
#' (e.g. \code{gaussian()}). See \link[stats]{family} for details on family
#' functions.
#' @param Indicator contingency table of \code{X} and \code{Y}.
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
#' @param Xextr boundary knots in the \code{X} direction. By default equal to
#' the range of \code{X}.
#' @param Yextr boundary knots in the \code{Y} direction. By default equal to
#' the range of \code{Y}.
#' @param show.iters logical variable indicating whether or not to print 
#' information at each step. By default equal to \code{FALSE}.
#' @param stoptype a character string indicating the type of GeDS stopping rule
#' to be used. It should be either \code{"SR"}, \code{"RD"} or \code{"LR"},
#' partial match allowed. See details of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param tol numeric value indicating the tolerance to be used in checking
#' whether two knots should be considered different during the knot placement
#' steps in stage A.
#' 
#' @return A \code{\link{GeDS-Class}} object, but without the \code{Formula},
#' \code{extcall}, \code{terms} and \code{znames} slots.
#' 
#' @seealso \code{\link{NGeDS}} and \code{\link{UnivariateFitters}}.
#' 
#' @export
#' @rdname BivariateFitters
#' 
#' @references
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}

BivariateFitter <- function(X, Y, Z, W, weights=rep(1,length(X)), Indicator,
                            beta = 0.5, phi = 0.99, min.intknots = 0,
                            max.intknots = 300, q = 2, Xextr = range(X),
                            Yextr = range(Y), show.iters = TRUE,
                            tol = as.double(1e-12))
  {
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol)
  
  # Initialize RSS
  RSSnew <- numeric()
  
  # Initialize knots matrix
  previousX <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  previousY <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  # Initialize coefficients matrix
  nw <- if(!is.null(W)) NCOL(W) else 0
  oldcoef   <- matrix(nrow = max.intknots + 1,
                      ncol = round((max.intknots/2 + 2)^2) + nw) # max number of coef; comes from maximizing f(x) = (x + 2)(max.intknots - x + 2)
  # Initialize internal knots
  Xintknots <- Yintknots <- NULL
  # Matrix for X, Y and residuals
  ordX <- order(X, Y); ordY <- order(Y, X)
  matr <- matrix(ncol = 3, nrow = length(Z))
  
  ##################################################################################
  ## STEP 1: Divide the sample space D into M_1/M_2 rectangular strips in X_1/X_2 ##
  ##################################################################################
  
  # Set the number of intervals for dividing the X_1 and X_2 dimensions (M_1 and M_2)
  nintX <- nintY <- 10
  
  # D_{1j} = [a_1 + (j - 1)(b_1 - a_1)/M_1, a_1 + j(b_1 - a_1)/M_1] \times [a_2, b_2], j = 1, ..., M_1
  # upperX = a_1 + j(b_1 - a_1)/M_1, i.e., the interval upper bound
  upperX <- seq(from = Xextr[1], to = Xextr[2], length = nintX + 1)[-1]
  # Divide the X domain in nintX  intervals (= M1)
  dX <- numeric(nintX)
  upperX <- upperX+1e-15
  # Count unique data points in the first X_1 interval
  dX[1] <- sum(unique(X) <= upperX[1])
  # Count unique data points in subsequent X_1 intervals, avoiding double counting
  for(i in 2:nintX) {
    dX[i] <- sum(unique(X) <= upperX[i]) - sum(dX)
  }
  # Identify intervals with no data points
  zeroesX <- dX == 0
  # Calculate the cumulative sum of counts for X_1 intervals
  dcumX <- cumsum(dX)
  
  # D_{2j} = [a_1, b_1] \times [a_2 + (j - 1)(b_2 - a_2)/M_2, a_2 + j(b_2 - a_2)/M_2], j = 1, ..., M_2
  # upperY = a_2 + j(b_2 - a_2)/M_2, i.e., the interval upper bound
  upperY <- seq(from = Yextr[1], to = Yextr[2], length = nintY + 1)[-1]
  # Divide the Y domain in nintY  intervals (= M2)
  dY <- numeric(nintY)
  upperY <- upperY + 1e-15
  # Count unique data points in the first X_2 interval
  dY[1] <- sum(unique(Y) <= upperY[1])
  # Count unique data points in subsequent X_1 intervals, avoiding double counting
  for(i in 2:nintY) {
    dY[i] <- sum(unique(Y) <= upperY[i]) - sum(dY)
  }
  # Identify intervals with no data points
  zeroesY <- dY == 0
  # Calculate the cumulative sum of counts for X_2 intervals
  dcumY <- cumsum(dY)
  
  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  
  Xctrl <- Yctrl <- FALSE # Initialize control flag indicating a new X/Y knot was added

  for (j in 1:(max.intknots + 1)) {
    
    if (j > 1) {
      # Sort internal knots vector if new X/Y intknot was added on previous iteration
      if(Xctrl)  Xintknots <- sort(Xintknots)
      if(Yctrl)  Yintknots <- sort(Yintknots)
    }
    
    ########################################################################
    ## STEP 2: Apply the IRLS procedure to find a bivariate ML spline fit ##
    ########################################################################
    first.deg <- SplineReg_fast_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                                    InterKnotsX = Xintknots, InterKnotsY = Yintknots,
                                    Xextr = Xextr, Yextr = Yextr, n = 2)
    
    # Store knots and coefficients
    previousX[j, 1:(length(Xintknots)+4)] <- sort(c(Xintknots, rep(Xextr,2)))
    previousY[j, 1:(length(Yintknots)+4)] <- sort(c(Yintknots, rep(Yextr,2)))
    oldcoef[j, 1:length(first.deg$Theta)] <- first.deg$Theta
    # Store weighted residuals
    matr <- cbind(X, Y, first.deg$Residuals*weights)
    
    ###########################
    ## STEP 3: Stopping rule ##
    ###########################
    RSS.tmp <- first.deg$RSS
    RSSnew <- c(RSSnew, RSS.tmp)
    if(j > q && (j - q > min.intknots)) {
      if(RSSnew[j]/RSSnew[j-q] >= phi) {
        break
      }
    }
    
    ###################################
    ## STEP 4. (i) X knot placement ##
    ##################################
    placeXKnot <- placeKnot(Dim = "X", Dim.intknots = Xintknots, matr = matr, Indicator = Indicator,
                            FixedDim = Y, ordFixedDim = ordY, nintFixedDim = nintY, zeroesFixedDim = zeroesY,
                            dcumFixedDim = dcumY, beta = beta) 
    
    
    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim

    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(Dim = "Y", Dim.intknots = Yintknots, matr = matr, Indicator = Indicator,
                            FixedDim = X, ordFixedDim = ordX, nintFixedDim = nintX, zeroesFixedDim = zeroesX,
                            dcumFixedDim = dcumX, beta = beta) 
    
    Ynewknot = placeYKnot$Dim.newknot; weightY = placeYKnot$weightDim; flagY = placeYKnot$flagDim
    
    
    # Check if both X and Y dimensions have flags indicating no valid knots could be found
    if(flagX && flagY) {
      print("Unable to find other knots satisfying required conditions")
      break # Exit the loop since no further knots can be added
    } else {
      # Adjust weights if only one dimension has no valid knot
      if (flagX) {
        # If no valid X knot, then set new knot on Y dimension
        weightX <- 0; weightY <- 1
      } else {
        if(flagY) {
          # If no valid Y knot, then set new knot on X dimension
          weightY <- 0; weightX <- 1
        }
      }
    }
    
    #############################################################################################
    ## STEP 4. (iii): if \omega_1^* => \omega_2^* a new knot \delta_1^* is added and viceversa ##
    #############################################################################################
    # A. If weight for X is greater, then add a new X knot
    if(weightX > weightY) {
      Ynewknot <- NULL
      Xctrl <- TRUE # Control flag indicating an X knot is to be added
      
      # Print iteration details if show.iters = TRUE
      if (show.iters) {
        if (j > q) {
          toprint <- paste0("Iteration ", j,": New X Knot = ", round(Xnewknot, 3), ", RSS = " ,
                            round(RSSnew[j],3), ", phi = ", round(RSSnew[j]/RSSnew[j-q], 3))
          } else {
            toprint <- paste0("Iteration ", j,": New X Knot = ", round(Xnewknot, 3),", RSS = " ,
                              round(RSSnew[j], 3))
          }
        print(toprint)
      }
      
    # B. If weight for Y is greater or equal, then add a new Y knot
      } else {
        Xnewknot <- NULL
        Yctrl <- TRUE # Control flag indicating a Y knot is to be added
        
        # Print iteration details if show.iters = TRUE
        if(show.iters) {
          if (j > q) {
            toprint <- paste0("Iteration ", j,": New Y Knot = ", round(Ynewknot,3), ", RSS = ",
                              round(RSSnew[j],3), ", phi = ", round(RSSnew[j]/RSSnew[j-q], 3))
            } else {
              toprint <- paste0("Iteration ", j,": New Y Knot = ", round(Ynewknot, 3),", RSS = ",
                                round(RSSnew[j], 3))
              }
          print(toprint)
        }
      }
    
    # Update knots vectors
    Yintknots <- c(Yintknots, Ynewknot)
    Xintknots <- c(Xintknots, Xnewknot)
    
    # Check if the total number of knots exceeds a threshold based on the length of the response
    if((length(Yintknots)+3)*(length(Xintknots)+3)>=length(Z)) {
      warning("Exiting stage A: Too many knots found")
      break # Exit the loop to avoid adding too many knots (prevent overfitting)
      }
  }
  
  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  
  # Keep the non-NA columns from the "j"th row
  # Keep the non-NA columns from the "j"th row
  toBeSaved <- sum(!is.na(previousX[j,]))
  previousX <- previousX[ ,-((toBeSaved + 1):(max.intknots + 4))]
  toBeSaved <- sum(!is.na(previousY[j,]))
  previousY <- previousY[ ,-((toBeSaved + 1):(max.intknots + 4))]
  
  # Keep the corresponding (intknotsX + 2) * (intknotsY + 2) coefficients
  oldcoef <- oldcoef[, 1:(NCOL(previousX) - 4 + 2) * (NCOL(previousY) - 4 + 2)]
  
  if (j == max.intknots + 1) {
    warning("Maximum number of iterations exceeded")
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    iter <- j
    } else {
      # Delete from the "j+1th" row until the "max.intknots+1th" row (i.e. keep the j first rows)
      previousX <- previousX[-((j+1):(max.intknots+1)), ] 
      previousY <- previousY[-((j+1):(max.intknots+1)), ]
      oldcoef   <- oldcoef[-((j+1):(max.intknots+1)),]
      
      lastXknots <- sum(!is.na(previousX[j-q, ]))
      lastYknots <- sum(!is.na(previousY[j-q, ]))
      iter <- j - q
    }
  
  # 1. LINEAR
  if (iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    llX <- llY <- NULL
    lin <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr, n = 2)
    } else {
      ikX <- previousX[iter, 3:(lastXknots - 2)]
      ikY <- previousY[iter, 3:(lastYknots - 2)]
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr, n = 2)
      }
  # 2. QUADRATIC
  if (iter < 3) {
    warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qqX <- qqY <- NULL
    squ <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr, n = 3)
    } else {
      # Stage B.1 (averaging knot location)
      qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
      qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
      # Stage B.2
      squ <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr, n = 3)
      }
  # 3. CUBIC
  if(iter < 4) {
    warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ccX <- ccY <- NULL
    cub <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr, n = 4)
    } else {
      # Stage B.1 (averaging knot location)
      ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
      ccY <- if (length(ikY) < 3) NULL else makenewknots(ikX, 4)
      # Stage B.2
      cub <- SplineReg_biv(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr, n = 4)
      }
  
  out <- list("Type" = "LM - Biv", "Linear.IntKnots" = list("Xk" = llX, "Yk" = llY), "Quadratic.IntKnots" = list("Xk" = qqX, "Yk" = qqY),
              "Cubic.IntKnots" = list("Xk" = ccX,"Yk" = ccY),"Dev.Linear" = lin$RSS, "Dev.Quadratic" = squ$RSS, "Dev.Cubic" = cub$RSS,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = list("previousX" = previousX, "previousY" = previousY),
              "Args"= args, "Call"= save, "Nintknots"= list("X"= length(llX), "Y"= length(llY)), "iters" = j, "Guesses" = NULL,
              "Coefficients" = oldcoef)
  class(out) <- "GeDS"
  return(out)
}


################################################################################
############################## GenBivariateFitter ##############################
################################################################################
#' @rdname BivariateFitters
#' @export

GenBivariateFitter <- function(X, Y, Z, W, family = family, weights = rep(1,length(X)),
                               Indicator, beta = 0.5, phi = 0.5, min.intknots = 0, 
                               max.intknots = 300, q = 2, Xextr=range(X), Yextr=range(Y),
                               show.iters=TRUE, tol = as.double(1e-12), stoptype = "RD")
{
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol, family = family)
  
  # Initialize RSS and phis
  RSSnew <- numeric()
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  
  # Initialize knots matrix
  previousX <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  previousY <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  # Initialize coefficients matrix
  nw <- if(!is.null(W)) NCOL(W) else 0
  oldcoef <- matrix(nrow = max.intknots + 1,
                    ncol = round((max.intknots/2 + 2)^2) + nw) # max number of coef; comes from maximizing f(x) = (x + 2)(max.intknots - x + 2)
  # Initialize internal knots
  Xintknots <- Yintknots <- NULL
  # Initial values for the coefficients used at each iteration of stage A in order to estimate the spline coefficients
  oldguess <- matrix(nrow = max.intknots + 1, ncol = NROW(Z))
  # oldguess <- matrix(nrow = max.intknots + 1,
  #                    ncol = round((max.intknots/2 + 2)^2)) # max number of coef; comes from maximizing f(x) = (x + 2)(max.intknots - x + 2)
  # Matrix for X, Y and residuals
  ordX <- order(X, Y); ordY <- order(Y, X)
  matr <- matrix(ncol = 3, nrow = length(Z))
  
  ##################################################################################
  ## STEP 1: Divide the sample space D into M_1/M_2 rectangular strips in X_1/X_2 ##
  ##################################################################################
  
  # Set the number of intervals for dividing the X_1 and X_2 dimensions (M_1 and M_2)
  nintX <- nintY <- 10
  
  # D_{1j} = [a_1 + (j - 1)(b_1 - a_1)/M_1, a_1 + j(b_1 - a_1)/M_1] \times [a_2, b_2], j = 1, ..., M_1
  # upperX = a_1 + j(b_1 - a_1)/M_1, i.e., the interval upper bound
  upperX <- seq(from = Xextr[1], to = Xextr[2], length = nintX + 1)[-1]
  # Divide the X domain in nintX  intervals (= M1)
  dX <- numeric(nintX)
  upperX <- upperX+1e-15
  # Count unique data points in the first X_1 interval
  dX[1] <- sum(unique(X) <= upperX[1])
  # Count unique data points in subsequent X_1 intervals, avoiding double counting
  for(i in 2:nintX) {
    dX[i] <- sum(unique(X) <= upperX[i]) - sum(dX)
  }
  # Identify intervals with no data points
  zeroesX <- dX == 0
  # Calculate the cumulative sum of counts for X_1 intervals
  dcumX <- cumsum(dX)
  
  # D_{2j} = [a_1, b_1] \times [a_2 + (j - 1)(b_2 - a_2)/M_2, a_2 + j(b_2 - a_2)/M_2], j = 1, ..., M_2
  # upperY = a_2 + j(b_2 - a_2)/M_2, i.e., the interval upper bound
  upperY <- seq(from = Yextr[1], to = Yextr[2], length = nintY + 1)[-1]
  # Divide the Y domain in nintY  intervals (= M2)
  dY <- numeric(nintY)
  upperY <- upperY + 1e-15
  # Count unique data points in the first X_2 interval
  dY[1] <- sum(unique(Y) <= upperY[1])
  # Count unique data points in subsequent X_1 intervals, avoiding double counting
  for(i in 2:nintY) {
    dY[i] <- sum(unique(Y) <= upperY[i]) - sum(dY)
  }
  # Identify intervals with no data points
  zeroesY <- dY == 0
  # Calculate the cumulative sum of counts for X_2 intervals
  dcumY <- cumsum(dY)
  
  # Stop type, min.X/Yintknots
  stoptype <- match.arg(stoptype)
  min.Xintknots <- min.intknots
  min.Yintknots <- min.intknots
  
  guess <- devi <- iter <- ncoef <- NULL
  
  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  
  Xctrl <- Yctrl <- FALSE # Initialize control flag indicating a new X/Y knot was added
  
  for (j in 1:(max.intknots + 1)) {
    
    if (j > 1)  {
      # Sort internal knots if new intknot was added on previous iteration and update the oldguess matrix
      if(Xctrl)  Xintknots <- sort(Xintknots)
      if(Yctrl)  Yintknots <- sort(Yintknots)
      oldguess[j,] <- guess
      # oldguess[j, 1:(lth-nw)] <- guess
      # guess <- c(guess, guess_w)
    }
    
    
    ########################################################################
    ## STEP 2: Apply the IRLS procedure to find a bivariate ML spline fit ##
    ########################################################################
    first.deg <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, W = W, weights = weights,
                                   InterKnotsX = Xintknots, InterKnotsY = Yintknots,
                                   n = 2, Xextr = Xextr, Yextr = Yextr, 
                                   family = family, mustart = guess)
    
    # 1. Check for NA values in the Theta vector to handle potential singularities
    lth <- length(first.deg$Theta)
    if (anyNA(first.deg$Theta)) {
      rank.basis <- Matrix::rankMatrix(first.deg$Basisbiv)
      cols <- NCOL(first.deg$Basisbiv)
      # (i) Handle the case when the basis matrix is singular
      if(rank.basis < cols) {
        warning("Matrix singular for the second time. Breaking the loop.")
        break
      # (ii) NA values in the Theta vector, but basis matrix is not singular (i.e. other issues)
        } else {
          stop("NA(s) in the coefficients")
        }
      
    # 2. If no NAs, update guess (=mustart in the next iteration)
    } else {
      guess <- first.deg$Predicted
    }
    
    # Vector accumulating the IRLS iterations deviances obtained at each GeDS iteration
    devi <- c(devi, first.deg$deviance)
    # Accumulated number of IRLS iterations at each GeDS iteration
    iter <- c(iter,length(devi))
    
    
    # Store knots and coefficients
    previousX[j, 1:(length(Xintknots)+4)] <- sort(c(Xintknots, rep(Xextr,2)))
    previousY[j, 1:(length(Yintknots)+4)] <- sort(c(Yintknots, rep(Yextr,2)))
    oldcoef[j, 1:lth] <- first.deg$Theta
    # guess_w <- if(nw > 0) first.deg$Theta[-(1:(lth-nw))] else NULL
    
    # Store residuals and deviance 
    res.tmp <- first.deg$Residuals
    RSS.tmp <- first.deg$temporary$lastdev
    RSSnew <- c(RSSnew, RSS.tmp)
    # Working weights (weights in the final iteration of the IRLS fit)
    working.weights <- first.deg$temporary$weights  
    # Store weighted residuals
    matr <- cbind(X, Y, first.deg$Residuals*working.weights*weights)
    
    ###########################
    ## STEP 2: Stopping Rule ##
    ###########################
    xl <- length(Xintknots)
    yl <- length(Yintknots)
    ncoef <- c(ncoef, lth)
    
    # (I) Smoothed Ratio of deviances
    if (stoptype == "SR") {
      if (j>q && (xl-q > min.Xintknots) && (yl-q > min.Yintknots) && length(phis) >= 3) {
        phnew     <- (RSSnew[j]/RSSnew[j-q])^(1/(ncoef[j]-ncoef[j-q]))
        phis      <- c(phis, phnew)
        phismod   <- log(1-phis)
        ccc       <- .lm.fit(cbind(1,(q+1):j),phismod)$coef
        phis_star <- c(phis_star,1-exp(ccc[1])*exp(ccc[2]*j))
        oldintc   <- c(oldintc,ccc[1])
        oldslp    <- c(oldslp,ccc[2])
        prnt      <- paste0(", phi_hat = ", round(1-exp(ccc[1])*exp(ccc[2]*j),3), ", ",
                            ncoef[j], " coefficients")
        
        if(1-exp(ccc[1])*exp(ccc[2]*j) >= phi) {
          break
        }
      }
    }
    # (II) Ratio of Deviances
    if (stoptype=="RD" | (stoptype=="SR" & length(phis) < 3)) {
      phnew <- (RSSnew[j]/RSSnew[j-q])^(1/(ncoef[j]-ncoef[j-q]))
      phis <- c(phis,phnew)
      if(j>q) prnt <- paste0(", phi = ",round(phnew,3), ", ",
                             ncoef[j]," coefficients")
      if(j>q && (xl-q > min.Xintknots)&& (yl-q > min.Yintknots)  ){
        if(RSSnew[j]/RSSnew[j-q] >= phi^(ncoef[j]-ncoef[j-q])) {      # stop rule
          break
        }
      }
    }
    # (III) Likelihood Ratio
    if (stoptype=="LR") {
      phis <- c(phis,RSSnew[j-q]-RSSnew[j])
      if (j > q) prnt <- paste0(", p = ", round(pchisq(-(RSSnew[j]-RSSnew[j-q]), df = (ncoef[j]-ncoef[j-q])),3), ", ",
                                 ncoef[j]," coefficients")
      if (j > q && (xl-q > min.Xintknots)&& (yl-q > min.Yintknots)) {
        if(-(RSSnew[j]-RSSnew[j-q]) < qchisq(phi,df=(ncoef[j]-ncoef[j-q]))) {
          break
        }
      }
    }
    
    ###################################
    ## STEP 4. (i) X knot placement ##
    ##################################
    placeXKnot <- placeKnot(Dim = "X", Dim.intknots = Xintknots, matr = matr, Indicator = Indicator,
                            FixedDim = Y, ordFixedDim = ordY, nintFixedDim = nintY, zeroesFixedDim = zeroesY,
                            dcumFixedDim = dcumY, beta = beta) 
    
    
    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim
    
    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(Dim = "Y", Dim.intknots = Yintknots, matr = matr, Indicator = Indicator,
                            FixedDim = X, ordFixedDim = ordX, nintFixedDim = nintX, zeroesFixedDim = zeroesX,
                            dcumFixedDim = dcumX, beta = beta) 
    
    Ynewknot = placeYKnot$Dim.newknot; weightY = placeYKnot$weightDim; flagY = placeYKnot$flagDim
    
    
    # Check if both X and Y dimensions have flags indicating no valid knots could be found
    if(flagX && flagY) {
      print("Unable to find other knots satisfying required conditions")
      break # Exit the loop since no further knots can be added
    } else {
      # Adjust weights if only one dimension has no valid knot
      if (flagX) {
        # If no valid X knot, then set new knot on Y dimension
        weightX <- 0; weightY <- 1
      } else {
        if(flagY) {
          # If no valid Y knot, then set new knot on X dimension
          weightY <- 0; weightX <- 1
        }
      }
    }
    
    #############################################################################################
    ## STEP 4. (iii): if \omega_1^* => \omega_2^* a new knot \delta_1^* is added and viceversa ##
    #############################################################################################
    # A. If weight for X is greater, then add a new X knot
    if (weightX > weightY) {
      Ynewknot <- NULL
      previousX <- rbind(previousX, c(Xnewknot,j))
      Xctrl <- TRUE # Control flag indicating an X knot is to be added
      
      # Print iteration details if show.iters = TRUE
      if(show.iters) {
        if (j > q) {
          toprint <- paste0("Iteration ",j,": New X Knot = ", round(Xnewknot,3),
                            ", RSS = " , round(RSSnew[j],3), prnt)
        } else {
            toprint <- paste0("Iteration ", j, ": New X Knot = ", round(Xnewknot,3), prnt)
        }
        print(toprint)
      }
      
      # B. If weight for Y is greater or equal, then add a new Y knot
      } else {
        Xnewknot <- NULL
        previousY <- rbind(previousY,c(Ynewknot,j))
        Yctrl <- TRUE
        
        # Print iteration details if show.iters = TRUE
        if (show.iters) {
          if (j > q) {
            toprint <- paste0("Iteration ",j,": New Y Knot = ", round(Ynewknot,3),
                              ", RSS = " , round(RSSnew[j],3), prnt)
          } else {
              toprint <- paste0("Iteration ",j,": New Y Knot = ", round(Ynewknot,3) , prnt)
              }
        print(toprint)
        }
      }
    
    # Calculate guess-coefficients for newknot
    # guess <- newknot.guess_biv(X, Y, Xintknots, Yintknots, Xextr, Yextr, guess, Xnewknot)
    
    # Update knots vectors
    Yintknots <- c(Yintknots,Ynewknot)
    Xintknots <- c(Xintknots,Xnewknot)
    
    # Check if the total number of knots exceeds a threshold based on the length of the response
    if((length(Yintknots)+3)*(length(Xintknots)+3)>=length(Z)) {
      warning("Exiting stage A: Too many knots found")
      break # Exit the loop to avoid adding too many knots (prevent overfitting)
    }
  }
  
  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  
  # Keep the non-NA columns from the "j"th row
  toBeSaved <- sum(!is.na(previousX[j,]))
  previousX <- previousX[ ,-((toBeSaved + 1):(max.intknots + 4))]
  toBeSaved <- sum(!is.na(previousY[j,]))
  previousY <- previousY[ ,-((toBeSaved + 1):(max.intknots + 4))]
  
  # Keep the corresponding (intknotsX + 2) * (intknotsY + 2) coefficients
  oldcoef <- oldcoef[, 1:(NCOL(previousX) - 4 + 2) * (NCOL(previousY) - 4 + 2)]
  
  if (j == max.intknots + 1) {
    warning("Maximum number of iterations exceeded")
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    iter <- j
  } else {
    # Delete from the "j+1th" row until the "max.intknots+1th" row (i.e. keep the j first rows)
    previousX <- previousX[-((j+1):(max.intknots+1)), ] 
    previousY <- previousY[-((j+1):(max.intknots+1)), ]
    oldcoef   <- oldcoef[-((j+1):(max.intknots+1)),]
    
    lastXknots <- sum(!is.na(previousX[j-q, ]))
    lastYknots <- sum(!is.na(previousY[j-q, ]))
    iter <- j - q
  }
  
  # 1. LINEAR
  if(iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    llX <- llY <- NULL
    lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
                             n = 2, family = family, mustart = oldguess[iter,])
    } else {
      ikX <- previousX[iter, 3:(lastXknots-2)]
      ikY <- previousY[iter, 3:(lastYknots-2)]
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
                               n = 2, family = family, mustart = oldguess[iter,])
    }
  # 2. QUADRATIC
  if (iter < 3) {
    warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qqX <- qqY <- NULL
    guess_lin <- lin$Predicted
    squ <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                             n = 3, family = family, mustart = guess_lin)
    } else {
      # Stage B.1 (averaging knot location)
      qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
      qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
      # Stage B.2
      guess_lin <- lin$Predicted
      squ <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                               n = 3, family = family, mustart = guess_lin)
    }
  # 3. CUBIC
  if (iter < 4) {
    warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ccX <- ccY <- NULL
    guess_sq <- squ$Predicted
    cub <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
                             n = 4, family = family, mustart = guess_sq)
    } else {
      # Stage B.1 (averaging knot location)
      ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
      ccY <- if (length(ikY) < 3) NULL else makenewknots(ikX, 4)
      # Stage B.2
      guess_sq <- squ$Predicted
      cub <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
                               n = 4, family = family, mustart = guess_sq)
      }
    
  out <- list("Type" = "GLM - Biv","Linear.IntKnots"=list("Xk" = llX,"Yk" = llY),"Quadratic.IntKnots"=list("Xk" = qqX,"Yk" = qqY),
              "Cubic.IntKnots"=list("Xk" = ccX,"Yk" = ccY),"Dev.Linear" = lin$RSS, "Dev.Quadratic" = squ$RSS,"Dev.Cubic" = cub$RSS,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previousX, "Args"= args,
              "Call"= save, "Nintknots"= list("X"= length(llX), "Y"= length(llY)),"iters" = j, "Guesses" = NULL,
              "Coefficients" = oldcoef)
  class(out) <- "GeDS"
  return(out)
}


#########################################
## STEP 4. (i)/(ii) X/Y knot placement ##
#########################################
placeKnot <- function(Dim, Dim.intknots, matr, Indicator, FixedDim, ordFixedDim, nintFixedDim, zeroesFixedDim, dcumFixedDim, beta)
{
  if (Dim == "X") {
    Dim.index <- 1
    by.row <- TRUE
  } else if (Dim == "Y") {
    Dim.index <- 2
    by.row <- FALSE
  }
  # Order matrFixedDim as (FixedDim, Dim)
  matrFixedDim <- matr[ordFixedDim,]
  # Clean matrFixedDim in case there are repeated observations
  matrFixedDim <- makeNewMatr(matrFixedDim, Indicator, by.row = by.row)
  
  # Initialize empty vectors for storing mean Dim values, Dim interval widths, counts, and distances for cluster formation
  Dim.mean <- Dim.width <- Dim.num <- dFixedDim.Dim <- numeric()
  # Initialize counter for clusters
  kk <- 1
  
  # Loop through each FixedDim interval to form clusters and determine Dim knot placement
  for (i in 1:nintFixedDim) {
    # Check if the current FixedDim interval contains data points
    if (!zeroesFixedDim[i]) {
      # i. Handle the first FixedDim interval separately
      if (i == 1) {
        # Extract the Dim values corresponding to the first FixedDim interval
        tmpDim <- matrFixedDim[1:dcumFixedDim[i], Dim.index]
        # If the interval contains more than one data point, sort the matrix by Dim values
        if (length(tmpDim) > 1) {
          matrFixedDim[1:dcumFixedDim[i],] <- matrFixedDim[1:dcumFixedDim[i],][order(tmpDim),]
        }
        # Extract the residuals within the first FixedDim interval
        tmpR <- matrFixedDim[1:dcumFixedDim[i], 3]
        
        # ii. Rest of FixedDim intervals
      } else {
        # Extract the Dim values corresponding to the current FixedDim interval
        tmpDim <- matrFixedDim[(dcumFixedDim[i - 1] + 1):dcumFixedDim[i], Dim.index]
        # If the interval contains more than one data point, sort the matrix by Dim values
        if (length(tmpDim) > 1) {
          matrFixedDim[(dcumFixedDim[i - 1] + 1):dcumFixedDim[i],] <- matrFixedDim[(dcumFixedDim[i - 1] + 1):dcumFixedDim[i],][order(tmpDim),]
        }
        # Extract the residuals within the current FixedDim interval
        tmpR <- matrFixedDim[(dcumFixedDim[i - 1] + 1):dcumFixedDim[i], 3]
      }
      # Group the consecutive residuals into clusters by their sign
      signs <- sign(tmpR)
      for (jj in 1:length(tmpR)) {
        # If all residual values have the same sign, count the entire set as one cluster
        if (all(signs == signs[1])) {
          dFixedDim.Dim[kk] <- length(signs) 
          kk = kk + 1
          break
        } else {
          # If signs change, identify the first change to split the cluster
          dFixedDim.Dim[kk] <- min(which(signs != signs[1]) - 1) # number of consecutive residuals with same sign
          # Update signs to exclude the identified cluster
          signs <- signs[-(1:dFixedDim.Dim[kk])]                 # extract cluster from signs
          kk = kk + 1
        }
      }
    }
  }
  
  # (Step 3 - UnivariateFitter) Within residual cluster means +  within-cluster ranges 
  dcumFixedDim.Dim <- cumsum(dFixedDim.Dim)
  Dim.mean  <- Dim.width <- numeric(length(dFixedDim.Dim))
  # Calculate the mean absolute residual and Dim-width for the first cluster
  Dim.mean[1]  <- abs(mean(matrFixedDim[1:dcumFixedDim.Dim[1], 3]))
  Dim.width[1] <- diff(range(matrFixedDim[1:dcumFixedDim.Dim[1], Dim.index]))
  # Loop to calculate the mean absolute residual and FixedDim-width for subsequent clusters
  for (i in 2:length(dFixedDim.Dim)) {
    Dim.mean[i]  <- abs(mean(matrFixedDim[(dcumFixedDim.Dim[i - 1] + 1):dcumFixedDim.Dim[i], 3]))
    Dim.width[i] <- diff(range(matrFixedDim[(dcumFixedDim.Dim[i - 1] + 1):dcumFixedDim.Dim[i], Dim.index]))
  }
  # (Step 4 - UnivariateFitter) Calculate the normalized within-cluster means and ranges 
  Dim.mean  <- Dim.mean/max(Dim.mean)
  Dim.width <- Dim.width/max(Dim.width)
  # Calculate the cluster weights (Step 5 - UnivariateFitter)
  Dim.weights <- beta*Dim.mean + (1 - beta)*Dim.width
  
  # Loop through each cluster to find the optimal placement for a new Dim knot
  u     <- length(dFixedDim.Dim) # total number of clusters
  flagDim <- F                  # flag to handle cases where all calculated weights are non-positive
  for (i in 1:u) {
    
    if (all(Dim.weights < 0)) {
      flagDim <- TRUE # Set the flagDim = TRUE if all weights are non-positive, indicating no valid knot can be found
      break
    }
    
    # Find the index of the cluster with the highest weight
    indice <- which.max(Dim.weights)
    # Determine the (index) boundaries of the cluster with the highest weight
    if (indice == 1) {dcumInf = 1} else {dcumInf = dcumFixedDim.Dim[indice - 1] + 1}
    dcumSup <- dcumFixedDim.Dim[indice]
    # Calculate the superior and inferior Dim-bounds
    sup <- matrFixedDim[dcumSup, Dim.index]
    inf <- matrFixedDim[dcumInf, Dim.index]
    
    # (Step 7 - UnivariateFitter) Compute the new Dim knot as a weighted average of Dim values
    # within the selected cluster, weighted by their residuals 
    Dim.newknot <- matrFixedDim[dcumSup:dcumInf, 3]%*%matrFixedDim[dcumSup:dcumInf, Dim.index]/sum(matrFixedDim[dcumSup:dcumInf, 3])
    
    # Check conditions to ensure the new knot is valid and does not conflict with existing knots
    # This involves ensuring there are no existing knots within the bounds of the selected cluster
    if (
      (((dcumSup - dcumInf) != 0) && (!any((Dim.intknots >= inf) * (Dim.intknots <= sup))))  ||
      ((dcumSup - dcumInf) == 0 && (dcumInf == 1 || dcumSup == length(FixedDim))) # for the case in which the entire set is within one cluster
    ) {
      break # If conditions are met, exit the loop as a valid knot has been found
    } else {
      Dim.weights[indice] <- -Inf # Invalidate the current cluster by setting its weight to negative infinity and continue the search
    }
  }
  
  weightDim <- Dim.weights[indice] # Store the weight of the selected cluster for further use
  
  # Return the new Dim knot and its weight, along with the flag indicating if a valid knot was found
  return(list(Dim.newknot = as.numeric(Dim.newknot), weightDim = as.numeric(weightDim), flagDim = flagDim))
}



