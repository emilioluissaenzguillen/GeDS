################################################################################
################################################################################
############################## Bivariate Fitters ###############################
################################################################################
################################################################################
#' @title Fitter Function for GeD Spline Regression for Bivariate Data
#' @name BivariateFitters
#' @aliases BivariateFitters BivariateFitter
#' @description
#' These are computing engines called by \code{\link{NGeDS}} and
#' \code{\link{GGeDS}}, needed for the underlying fitting procedures.
#'
#' @param X A numeric vector containing \eqn{N} sample values of the first
#' independent variable chosen to enter the spline regression component of the
#' predictor model.
#' @param Y A numeric vector containing \eqn{N} sample values of the second
#' independent variable chosen to enter the spline regression component of the
#' predictor model.
#' @param Z A vector of size \eqn{N} containing the observed values of the
#' response variable.
#' @param W A design matrix with \eqn{N} rows containing other covariates
#' selected to enter the parametric component of the predictor model (see
#' \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected, it
#' is set to \code{NULL} by default.
#' @param weights An optional vector of size \eqn{N} of `prior weights' to be
#' put on the observations in the fitting process in case the user requires
#' weighted GeDS fitting. It is \code{NULL} by default.
#' @param family A description of the error distribution and link function to be
#' used in the model. This can be a character string naming a family function
#' (e.g. \code{"gaussian"}), the family function itself (e.g.
#' \code{\link[stats]{gaussian}}) or the result of a call to a family function
#' (e.g. \code{gaussian()}). See \link[stats]{family} for details on family
#' functions. Note that this argument applies only to \code{GenBivariateFitter}.
#' @param indicator A contingency table (i.e., frequency of observations) for the
#' independent variables \code{X} and \code{Y}.
#' @param beta Numeric parameter in the interval \eqn{[0,1]} tuning the knot
#' placement in stage A of GeDS. See the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param phi Numeric parameter in the interval \eqn{(0,1)} specifying the
#' threshold for the stopping rule  (model selector) in stage A of GeDS. See
#' also \code{stoptype} and Details in the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param min.intknots Optional parameter specifying the minimum number of
#' internal knots required in Stage A's fit. Default is \code{0L}.
#' @param max.intknots Optional parameter allowing the user to set a maximum
#' number of internal knots to be added in Stage A by the GeDS estimation
#' algorithm. Default equals \code{300L}.
#' @param q Numeric parameter which allows to fine-tune the stopping rule of
#' stage A of GeDS, by default equal to 2. See Details in the description of
#' \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param Xextr Boundary knots in the \code{X} direction. By default equal to
#' the range of \code{X}.
#' @param Yextr Boundary knots in the \code{Y} direction. By default equal to
#' the range of \code{Y}.
#' @param show.iters Logical variable indicating whether or not to print fitting
#' information at each step. Default is \code{FALSE}.
#' @param stoptype A character string indicating the type of GeDS stopping rule
#' to be used. It should be either \code{"SR"}, \code{"RD"} or \code{"LR"},
#' partial match allowed. See details of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param tol Numeric value indicating the tolerance to be used in checking
#' whether two knots should be considered different during the knot placement
#' steps in stage A.
#' @param higher_order A logical defining whether to compute the higher
#' order fits (quadratic and cubic) after stage A is run. Default is
#' \code{TRUE}.
#' @param Xintknots A vector of starting internal knots in the \code{X} direction.
#' Allows the user to begin Stage A's GeDS algorithm with a linear (least-squares)
#' spline fit using a predefined vector of internal \code{X} knots, instead of
#' starting with a straight line fit (i.e., with zero internal knots). Note that
#' this is not available for \code{GenBivariateFitter}. Default is \code{NULL}.
#' @param Yintknots A vector of starting internal knots in the \code{Y} direction.
#' Allows the user to begin Stage A's GeDS algorithm with a linear (least-squares)
#' spline fit using a predefined vector of internal \code{Y} knots, instead of
#' starting with a straight line fit (i.e., with zero internal knots). Note that
#' this is not available for \code{GenBivariateFitter}. Default is \code{NULL}.
#'
#' @return A \code{"GeDS"} class object, but without the \code{formula},
#' \code{extcall}, \code{terms} and \code{znames} slots.
#'
#' @references
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}
#'
#' @seealso \code{\link{NGeDS}}, \code{\link{GGeDS}} and \code{\link{UnivariateFitters}}.
#'
#' @rdname BivariateFitters
#' @importFrom stats .lm.fit qchisq pchisq
#' @export

BivariateFitter <- function(X, Y, Z, W, weights = rep(1,length(X)), indicator,
                            beta = 0.5, phi = 0.99, min.intknots = 0L,
                            max.intknots = 300L, q = 2L, Xextr = range(X),
                            Yextr = range(Y), show.iters = TRUE,
                            tol = as.double(1e-12), stoptype = c("SR","RD","LR"),
                            higher_order = TRUE, Xintknots = NULL, Yintknots = NULL)
  {
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol)

  # Initialize rss and phis
  n_starting_intknots <- length(Xintknots) + length(Yintknots)
  rssnew <- numeric(n_starting_intknots)
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  # Stop type
  stoptype <- match.arg(stoptype)

  # Initialize knots matrix
  previousX <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  previousY <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  # Initialize coefficients matrix
  nw <- if(!is.null(W)) NCOL(W) else 0
  oldcoef   <- matrix(nrow = max.intknots + 1,
                      ncol = round((max.intknots/2 + 2)^2) + nw) # max number of coef; comes from maximizing f(x) = (x + 2)(max.intknots - x + 2)

  # Matrix for X, Y and residuals
  ordX <- order(X, Y); ordY <- order(Y, X)
  matr <- matrix(ncol = 3, nrow = length(Z))

  ##################################################################################
  ## STEP 1: Divide the sample space D into M_1/M_2 rectangular strips in X_1/X_2 ##
  ##################################################################################

  # Set the number of intervals for dividing the X_1 and X_2 dimensions (M_1 and M_2)
  nintX <- nintY <- as.integer(sqrt(length(Z)))

  # D_{1j} = [a_1 + (j - 1)(b_1 - a_1)/M_1, a_1 + j(b_1 - a_1)/M_1] \times [a_2, b_2], j = 1, ..., M_1
  # upperX = a_1 + j(b_1 - a_1)/M_1, i.e., the interval upper bound
  upperX <- seq(from = Xextr[1], to = Xextr[2], length = nintX + 1)[-1]
  upperX <- upperX + 1e-15

  # D_{2j} = [a_1, b_1] \times [a_2 + (j - 1)(b_2 - a_2)/M_2, a_2 + j(b_2 - a_2)/M_2], j = 1, ..., M_2
  # upperY = a_2 + j(b_2 - a_2)/M_2, i.e., the interval upper bound
  upperY <- seq(from = Yextr[1], to = Yextr[2], length = nintY + 1)[-1]
  upperY <- upperY + 1e-15

  # Initialized iter and ncoef
  iter <- ncoef <- NULL

  # GeDS iterations start by j = n_starting_intknots + 1
  init.iter <- if (is.null(Xintknots) && is.null(Yintknots)) 1 else  n_starting_intknots + 1

  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################

  Xctrl <- Yctrl <- FALSE # Initialize control flag indicating a new X/Y knot was added

  for (j in init.iter:(max.intknots + 1)) {

    if (j > 1) {
      # Sort internal knots vector if new X/Y intknot was added on previous iteration
      if(Xctrl)  Xintknots <- sort(Xintknots)
      if(Yctrl)  Yintknots <- sort(Yintknots)
    }

    ########################################################################
    ## STEP 2: Apply the IRLS procedure to find a bivariate ML spline fit ##
    ########################################################################
    first.deg <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                               InterKnotsX = Xintknots, InterKnotsY = Yintknots,
                               Xextr = Xextr, Yextr = Yextr, n = 2, fast = TRUE)

    # Store knots and coefficients
    previousX[j, 1:(length(Xintknots)+4)] <- sort(c(Xintknots, rep(Xextr,2)))
    previousY[j, 1:(length(Yintknots)+4)] <- sort(c(Yintknots, rep(Yextr,2)))
    lth <- length(first.deg$theta); ncoef <- c(ncoef, lth)
    oldcoef[j, 1:lth] <- first.deg$theta

    # Store weighted residuals
    matr <- cbind(X, Y, first.deg$residuals*weights)
    # Store rss
    rss.tmp <- first.deg$rss
    rssnew <- c(rssnew, rss.tmp)

    ###########################
    ## STEP 3: Stopping rule ##
    ###########################
    if (j > q + n_starting_intknots) {

      if (rssnew[j]/rssnew[j-q] > 1) break

      # Adding the current ratio of deviances to the 'phis' vector
      if (stoptype == "LR") {
        phis <- c(phis, rssnew[j-q]-rssnew[j])
        } else {
          phnew <- (rssnew[j]/rssnew[j-q])^(1/(ncoef[j]-ncoef[j-q]))
          phis <- c(phis, phnew)
          }

      if (j - q > min.intknots) {

        # (I) Smoothed Ratio of deviances
        if (stoptype == "SR") {
          # \hat{φ}_κ = 1 - exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # 1-\hat{φ}_κ = exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # ln(1-\hat{φ}_κ) = \hat{γ}_0 + \hat{γ}_1*κ

          # Fit a linear model ln(1-φ) ~ \hat{γ}_0 + \hat{γ}_1*κ to the sample {φ_h, h}^κ_{h=q}
          phismod <- log(1-phis); kappa <- length(Xintknots) + length(Yintknots)
          gamma <- .lm.fit(cbind(1,(q+1):j),phismod)$coef
          # Calculate \hat{φ}_κ based on the estimated coefficients
          phi_kappa <- 1 - exp(gamma[1])*exp(gamma[2]*kappa)
          # Store \hat{φ}_κ and the estimated coefficients \hat{γ}_0 and \hat{γ}_1
          phis_star <- c(phis_star, phi_kappa)
          oldintc   <- c(oldintc, gamma[1]); oldslp <- c(oldslp, gamma[2])
          # Creating a print statement that shows the current adjusted phi value
          prnt      <- paste0(", phi_hat = ", round(phi_kappa, 3), ", ",
                              ncoef[j], " coefficients")
          # Check if \hat{φ}_κ ≥ φ_{exit}
          if(phi_kappa >= phi)  break

          # (II) Ratio of Deviances
          } else if (stoptype == "RD") {
            prnt <- paste0(", phi = ",round(phnew,3), ", ",
                           ncoef[j]," coefficients")
            # if (rssnew[j]/rssnew[j-q] >= phi^(ncoef[j]-ncoef[j-q])) break
            if (rssnew[j]/rssnew[j-q] >= phi) break

            # (III) Likelihood Ratio
            } else if (stoptype == "LR") {
              prnt <- paste0(", p = ",
                             round(pchisq(-(rssnew[j]-rssnew[j-q]), df = (ncoef[j]-ncoef[j-q])),3),
                             ", ", ncoef[j]," coefficients")
              if(-(rssnew[j]-rssnew[j-q]) < qchisq(phi,df=(ncoef[j]-ncoef[j-q]))) break
            }
      }
    }

    ###################################
    ## STEP 4. (i) X knot placement ##
    ##################################
    placeXKnot <- placeKnot(
      Dim = "X", Dim.intknots = Xintknots, matr = matr,
      indicator = indicator, FixedDim = Y, ordFixedDim = ordY,
      nintFixedDim = nintY, upperFixedDim = upperY, beta = beta
    )

    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim

    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(
      Dim = "Y", Dim.intknots = Yintknots, matr = matr,
      indicator = indicator, FixedDim = X, ordFixedDim = ordX,
      nintFixedDim = nintX, upperFixedDim = upperX, beta = beta
    )

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
    if (weightX > weightY || (weightX == weightY && length(Xintknots) > length(Yintknots))) {
      Xctrl <- TRUE # Control flag indicating an X knot is to be added
      Ynewknot <- NULL
      knottype <- "X"
      knotValue <- Xnewknot
    # B. If weight for Y is greater or equal, then add a new Y knot
    } else {
      Yctrl <- TRUE
      Xnewknot <- NULL
      knottype <- "Y"
      knotValue <- Ynewknot
    }

    # Print iteration details if show.iters is TRUE
    if (show.iters) {
      if (j > q) {
        toprint <- paste0("Iteration ", j, ": New ", knottype, " Knot = ", round(knotValue, 3),
                          ", rss = ", round(rssnew[j], 3), prnt)
      } else {
        toprint <- paste0("Iteration ", j, ": New ", knottype, " Knot = ", round(knotValue, 3), prnt)
      }
      print(toprint)
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
  toBeSaved <- sum(!is.na(previousX[j,]))
  previousX <- previousX[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]
  toBeSaved <- sum(!is.na(previousY[j,]))
  previousY <- previousY[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]

  # Keep the corresponding (intknotsX + 2) * (intknotsY + 2) coefficients
  oldcoef <- oldcoef[, 1:((NCOL(previousX) - 4 + 2) * (NCOL(previousY) - 4 + 2)), drop = FALSE]

  if (j == max.intknots + 1) {
    warning("Maximum number of iterations exceeded")
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    iter <- j
    } else {
      # Delete from the "j+1th" row until the "max.intknots+1th" row (i.e. keep the j first rows)
      previousX <- previousX[-((j+1):(max.intknots+1)), , drop = FALSE]
      previousY <- previousY[-((j+1):(max.intknots+1)), , drop = FALSE]
      oldcoef   <- oldcoef[-((j+1):(max.intknots+1)), , drop = FALSE]

      lastXknots <- sum(!is.na(previousX[j-q, ]))
      lastYknots <- sum(!is.na(previousY[j-q, ]))
      iter <- j - q
    }

  # 1. LINEAR
  if (iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    llX <- llY <- NULL
    lin <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                         InterKnotsX = llX, InterKnotsY = llY,
                         Xextr = Xextr, Yextr = Yextr, n = 2)
    } else {
      ikX <- if (lastXknots > 4) previousX[iter, 3:(lastXknots - 2)] else NULL
      ikY <- if (lastYknots > 4) previousY[iter, 3:(lastYknots - 2)] else NULL
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                           InterKnotsX = llX, InterKnotsY = llY,
                           Xextr = Xextr, Yextr = Yextr, n = 2)
    }
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
  # 2. QUADRATIC
  if (iter < 3) {
    warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qqX <- qqY <- NULL
    squ <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                         InterKnotsX = qqX, InterKnotsY = qqY,
                         Xextr = Xextr, Yextr = Yextr, n = 3)
    } else {
      # Stage B.1 (averaging knot location)
      qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
      qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
      # Stage B.2
      squ <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                           InterKnotsX = qqX, InterKnotsY = qqY,
                           Xextr = Xextr, Yextr = Yextr, n = 3)
      }
  # 3. CUBIC
  if(iter < 4) {
    warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ccX <- ccY <- NULL
    cub <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                         InterKnotsX = ccX, InterKnotsY = ccY,
                         Xextr = Xextr, Yextr = Yextr, n = 4)
    } else {
      # Stage B.1 (averaging knot location)
      ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
      ccY <- if (length(ikY) < 3) NULL else makenewknots(ikY, 4)
      # Stage B.2
      cub <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights,
                           InterKnotsX = ccX, InterKnotsY = ccY,
                           Xextr = Xextr, Yextr = Yextr, n = 4)
    }
    } else {
      qqX <- qqY <- squ <- ccX <- ccY <- cub <- NULL
    }

  out <- list("type" = "LM - Biv", "linear.intknots" = list("Xk" = llX, "Yk" = llY), "quadratic.intknots" = list("Xk" = qqX, "Yk" = qqY),
              "cubic.intknots" = list("Xk" = ccX,"Yk" = ccY),"dev.linear" = lin$rss, "dev.quadratic" = squ$rss, "dev.cubic" = cub$rss,
              "rss" = rssnew, "linear.fit" = lin, "quadratic.fit" = squ, "cubic.fit" = cub, "stored" = list("previousX" = previousX, "previousY" = previousY),
              "args" = args, "Call" = save, "Nintknots" = list("X" = length(llX), "Y" = length(llY)), "iters" = j, "guesses" = NULL,
              "coefficients" = oldcoef)
  class(out) <- "GeDS"
  return(out)
}


################################################################################
############################## GenBivariateFitter ##############################
################################################################################
#' @rdname BivariateFitters
#' @importFrom Matrix rankMatrix
#' @importFrom stats .lm.fit pchisq qchisq
#' @export

GenBivariateFitter <- function(X, Y, Z, W, family = family, weights = rep(1,length(X)),
                               indicator, beta = 0.5, phi = 0.5,
                               min.intknots = 0L, max.intknots = 300L, q = 2L,
                               Xextr = range(X), Yextr=range(Y),
                               show.iters=TRUE, tol = as.double(1e-12),
                               stoptype = c("SR","RD","LR"), higher_order = TRUE)
{
  # Capture the function call
  save <- match.call()
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol, family = family)

  # Initialize rss and phis
  rssnew <- numeric()
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  # Stop type
  stoptype <- match.arg(stoptype)

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
  nintX <- nintY <- as.integer(sqrt(length(Z)))

  # D_{1j} = [a_1 + (j - 1)(b_1 - a_1)/M_1, a_1 + j(b_1 - a_1)/M_1] \times [a_2, b_2], j = 1, ..., M_1
  # upperX = a_1 + j(b_1 - a_1)/M_1, i.e., the interval upper bound
  upperX <- seq(from = Xextr[1], to = Xextr[2], length = nintX + 1)[-1]
  upperX <- upperX + 1e-15

  # D_{2j} = [a_1, b_1] \times [a_2 + (j - 1)(b_2 - a_2)/M_2, a_2 + j(b_2 - a_2)/M_2], j = 1, ..., M_2
  # upperY = a_2 + j(b_2 - a_2)/M_2, i.e., the interval upper bound
  upperY <- seq(from = Yextr[1], to = Yextr[2], length = nintY + 1)[-1]
  upperY <- upperY + 1e-15

  # Stop type, min.X/Yintknots
  stoptype <- match.arg(stoptype)
  min.Xintknots <- min.intknots
  min.Yintknots <- min.intknots

  guess <- irlsAccumIterCount <- ncoef <- NULL

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
                                   family = family, mustart = guess, fast = TRUE)

    # 1. Check for NA values in the theta vector to handle potential singularities
    if (anyNA(first.deg$theta)) {
      basis.biv <- tensorProd(first.deg$Xbasis, first.deg$Ybasis)
      basis.full <- cbind(basis.biv, W)

      rank.basis <- rankMatrix(basis.full)
      cols <- NCOL(basis.full)
      # (i) Handle the case when the basis matrix is singular
      if(rank.basis < cols) {
        warning("Matrix singular for the second time. Breaking the loop.")
        break
      # (ii) NA values in the theta vector, but basis matrix is not singular (i.e. other issues)
        } else {
          stop("NA(s) in the coefficients")
        }

    # 2. If no NAs, update guess (=mustart in the next iteration)
    } else {
      guess <- first.deg$predicted
    }

    # Accumulated number of IRLS iterations at each GeDS iteration
    irlsAccumIterCount <- c(irlsAccumIterCount, first.deg$temporary$iter)


    # Store knots and coefficients
    previousX[j, 1:(length(Xintknots)+4)] <- sort(c(Xintknots, rep(Xextr,2)))
    previousY[j, 1:(length(Yintknots)+4)] <- sort(c(Yintknots, rep(Yextr,2)))
    lth <- length(first.deg$theta); ncoef <- c(ncoef, lth)
    oldcoef[j, 1:lth] <- first.deg$theta
    # guess_w <- if(nw > 0) first.deg$theta[-(1:(lth-nw))] else NULL

    # Store residuals and deviance
    res.tmp <- first.deg$residuals
    rss.tmp <- first.deg$rss
    rssnew <- c(rssnew, rss.tmp)
    # Working weights (weights in the final iteration of the IRLS fit)
    working.weights <- first.deg$temporary$weights
    # Store weighted residuals
    matr <- cbind(X, Y, first.deg$residuals*working.weights*weights)

    ###########################
    ## STEP 2: Stopping Rule ##
    ###########################
    if (j > q) {

      if (rssnew[j]/rssnew[j-q] > 1) break

      # Adding the current ratio of deviances to the 'phis' vector
      if (stoptype == "LR") {
        phis <- c(phis, rssnew[j-q]-rssnew[j])
      } else {
        phnew <- (rssnew[j]/rssnew[j-q])^(1/(ncoef[j]-ncoef[j-q]))
        phis <- c(phis, phnew)
      }

      if (j - q > min.intknots) {
        # (I) Smoothed Ratio of deviances
        if (stoptype == "SR") {
          # \hat{φ}_κ = 1 - exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # 1-\hat{φ}_κ = exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # ln(1-\hat{φ}_κ) = \hat{γ}_0 + \hat{γ}_1*κ

          # Fit a linear model ln(1-φ) ~ \hat{γ}_0 + \hat{γ}_1*κ to the sample {φ_h, h}^κ_{h=q}
          phismod <- log(1-phis); kappa <- kappa <- length(Xintknots) + length(Yintknots)
          gamma <- .lm.fit(cbind(1,(q+1):j),phismod)$coef
          # Calculate \hat{φ}_κ based on the estimated coefficients
          phi_kappa <- 1 - exp(gamma[1])*exp(gamma[2]*kappa)
          # Store \hat{φ}_κ and the estimated coefficients \hat{γ}_0 and \hat{γ}_1
          phis_star <- c(phis_star, phi_kappa)
          oldintc   <- c(oldintc, gamma[1]); oldslp <- c(oldslp, gamma[2])
          # Creating a print statement that shows the current adjusted phi value
          prnt      <- paste0(", phi_hat = ", round(phi_kappa, 3), ", ",
                              ncoef[j], " coefficients")
          # Check if \hat{φ}_κ ≥ φ_{exit}
          if(phi_kappa >= phi)  break
          # (II) Ratio of Deviances
        } else if (stoptype == "RD") {
          prnt <- paste0(", phi = ",round(phnew,3), ", ",
                         ncoef[j]," coefficients")
          if (rssnew[j]/rssnew[j-q] >= phi^(ncoef[j]-ncoef[j-q])) break
          # (III) Likelihood Ratio
        } else if (stoptype == "LR") {
          prnt <- paste0(", p = ",
                         round(pchisq(-(rssnew[j]-rssnew[j-q]), df = (ncoef[j]-ncoef[j-q])),3),
                         ", ", ncoef[j]," coefficients")
          if(-(rssnew[j]-rssnew[j-q]) < qchisq(phi,df=(ncoef[j]-ncoef[j-q]))) break
        }
      }
    }

    ###################################
    ## STEP 4. (i) X knot placement ##
    ##################################
    placeXKnot <- placeKnot(
      Dim = "X", Dim.intknots = Xintknots, matr = matr,
      indicator = indicator, FixedDim = Y, ordFixedDim = ordY,
      nintFixedDim = nintY, upperFixedDim = upperY, beta = beta
    )

    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim

    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(
      Dim = "Y", Dim.intknots = Yintknots, matr = matr,
      indicator = indicator, FixedDim = X, ordFixedDim = ordX,
      nintFixedDim = nintX, upperFixedDim = upperX, beta = beta
    )

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
    if (weightX > weightY || (weightX == weightY && length(Xintknots) > length(Yintknots))) {
      Xctrl <- TRUE # Control flag indicating an X knot is to be added
      Ynewknot <- NULL
      knottype <- "X"
      knotValue <- Xnewknot
    # B. If weight for Y is greater or equal, then add a new Y knot
    } else {
      Yctrl <- TRUE
      Xnewknot <- NULL
      knottype <- "Y"
      knotValue <- Ynewknot
    }

    # Print iteration details if show.iters is TRUE
    if (show.iters) {
      if (j > q) {
        toprint <- paste0("Iteration ", j, ": New ", knottype, " Knot = ", round(knotValue, 3),
                          ", rss = ", round(rssnew[j], 3), prnt)
      } else {
        toprint <- paste0("Iteration ", j, ": New ", knottype, " Knot = ", round(knotValue, 3), prnt)
      }
      print(toprint)
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
  previousX <- previousX[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]
  toBeSaved <- sum(!is.na(previousY[j,]))
  previousY <- previousY[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]

  # Keep the corresponding (intknotsX + 2) * (intknotsY + 2) coefficients
  oldcoef <- oldcoef[, 1:((NCOL(previousX) - 4 + 2) * (NCOL(previousY) - 4 + 2)), drop = FALSE]

  if (j == max.intknots + 1) {
    warning("Maximum number of iterations exceeded")
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    iter <- j
  } else {
    # Delete from the "j+1th" row until the "max.intknots+1th" row (i.e. keep the j first rows)
    previousX <- previousX[-((j+1):(max.intknots+1)), , drop = FALSE]
    previousY <- previousY[-((j+1):(max.intknots+1)), , drop = FALSE]
    oldcoef   <- oldcoef[-((j+1):(max.intknots+1)), , drop = FALSE]

    lastXknots <- sum(!is.na(previousX[j-q, ]))
    lastYknots <- sum(!is.na(previousY[j-q, ]))
    iter <- j - q
  }

  # If model selected is from first iteration
  if (iter == 1) {
    mustart <- NULL
    } else {
      mustart <- oldguess[iter,]
    }

  # 1. LINEAR
  if(iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    llX <- llY <- NULL
    lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
                             n = 2, family = family, mustart = mustart)
    } else {
      ikX <- if (lastXknots > 4) previousX[iter, 3:(lastXknots-2)] else NULL
      ikY <- if (lastYknots > 4) previousY[iter, 3:(lastYknots-2)] else NULL
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
                               n = 2, family = family, mustart = mustart)
    }
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
    # 2. QUADRATIC
    if (iter < 3) {
      warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      qqX <- qqY <- NULL
      guess_lin <- lin$predicted
      squ <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                               n = 3, family = family, mustart = guess_lin)
      } else {
        # Stage B.1 (averaging knot location)
        qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
        qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
        # Stage B.2
        guess_lin <- lin$predicted
        squ <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                                 n = 3, family = family, mustart = guess_lin)
        }
    # 3. CUBIC
    if (iter < 4) {
      warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      ccX <- ccY <- NULL
      guess_sq <- squ$predicted
      cub <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
                               n = 4, family = family, mustart = guess_sq)
      } else {
        # Stage B.1 (averaging knot location)
        ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
        ccY <- if (length(ikY) < 3) NULL else makenewknots(ikY, 4)
        # Stage B.2
        guess_sq <- squ$predicted
        cub <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
                                 n = 4, family = family, mustart = guess_sq)
      }
    } else {
      qqX <- qqY <- squ <- ccX <- ccY <- cub <- NULL
      }

  out <- list("type" = "GLM - Biv", "linear.intknots" = list("Xk" = llX, "Yk" = llY), "quadratic.intknots" = list("Xk" = qqX, "Yk" = qqY),
              "cubic.intknots" = list("Xk" = ccX, "Yk" = ccY), "dev.linear" = lin$rss, "dev.quadratic" = squ$rss, "dev.cubic" = cub$rss,
              "rss" = rssnew, "linear.fit" = lin, "quadratic.fit" = squ, "cubic.fit" = cub, "stored" = list("previousX" = previousX, "previousY" = previousY),
              "args"= args, "Call" = save, "Nintknots" = list("X"= length(llX), "Y"= length(llY)), "iters" = j, "guesses" = NULL,
              "coefficients" = oldcoef, "iterIrls" = irlsAccumIterCount)
  class(out) <- "GeDS"
  return(out)
}


##########################################################
## Dimension-indexed core for Stage A knot placement    ##
##########################################################
placeDimKnot <- function(dim.index, intknots, coordinates, residuals,
                         strip.id, beta, dim.range = range(coordinates[, dim.index]))
{
  coordinates <- as.matrix(coordinates)
  dim.index <- as.integer(dim.index)

  if (length(dim.index) != 1L || is.na(dim.index) ||
      dim.index < 1L || dim.index > NCOL(coordinates)) {
    stop("'dim.index' must identify one column of 'coordinates'.", call. = FALSE)
  }
  if (NROW(coordinates) != length(residuals) ||
      NROW(coordinates) != length(strip.id)) {
    stop("'coordinates', 'residuals' and 'strip.id' must have compatible lengths.",
         call. = FALSE)
  }
  if (length(beta) != 1L || is.na(beta) || beta < 0 || beta > 1) {
    stop("'beta' must be a single value in [0, 1].", call. = FALSE)
  }

  # NA strip IDs mark observations belonging to strips that should be dropped.
  keep <- !is.na(strip.id)
  Dim.values <- coordinates[keep, dim.index]
  Dim.residuals <- residuals[keep]
  kept.strips <- strip.id[keep]

  if (!length(Dim.values)) {
    return(list(Dim.newknot = NA_real_, weightDim = NA_real_, flagDim = TRUE))
  }

  # Order observations within every fixed-coordinate strip by the candidate
  # coordinate. This is the only ordering required by the knot-selection core.
  ord <- order(kept.strips, Dim.values, method = "radix")
  Dim.values <- Dim.values[ord]
  Dim.residuals <- Dim.residuals[ord]
  kept.strips <- kept.strips[ord]

  # Form consecutive residual-sign clusters, restarting at strip boundaries.
  signs <- sign(Dim.residuals)
  breaks <- c(TRUE,
              kept.strips[-1L] != kept.strips[-length(kept.strips)] |
                signs[-1L] != signs[-length(signs)])
  cluster.starts <- which(breaks)
  cluster.ends <- c(cluster.starts[-1L] - 1L, length(signs))

  Dim.mean <- vapply(
    seq_along(cluster.starts),
    function(i) abs(mean(Dim.residuals[cluster.starts[i]:cluster.ends[i]])),
    numeric(1)
  )
  Dim.width <- vapply(
    seq_along(cluster.starts),
    function(i) diff(range(Dim.values[cluster.starts[i]:cluster.ends[i]])),
    numeric(1)
  )

  # Normalize cluster height and width as in the bivariate implementation.
  height.scale <- max(Dim.mean)
  if (!is.finite(height.scale) || height.scale == 0) {
    return(list(Dim.newknot = NA_real_, weightDim = NA_real_, flagDim = TRUE))
  }
  Dim.mean <- Dim.mean / height.scale
  if (max(Dim.width) != 0) Dim.width <- Dim.width / max(Dim.width)
  Dim.weights <- beta * Dim.mean + (1 - beta) * Dim.width

  if (is.null(intknots)) intknots <- NA_real_
  oldknots <- sort(c(intknots, dim.range))

  # R debugging implementation:
  # x <- findNewDimKnot_R(
  #   cluster.ends, Dim.weights, oldknots, Dim.values, Dim.residuals
  # )
  xx <- findNewDimKnot(
    cluster.ends, Dim.weights, oldknots, Dim.values, Dim.residuals
  )
  # if (isTRUE(all.equal(x, xx, tolerance = 1e-12))) print("Both equal!")

  list(
    Dim.newknot = as.numeric(xx$Dim.newknot),
    weightDim = as.numeric(xx$weightDim),
    flagDim = xx$flagDim
  )
}


##########################################################
## Fixed-coordinate cells for dimension-wise placement ##
##########################################################
fixedDimCellId <- function(fixed.coordinates, upper.bounds)
{
  fixed.coordinates <- as.matrix(fixed.coordinates)
  if (!NCOL(fixed.coordinates)) {
    stop("'fixed.coordinates' must contain at least one column.",
         call. = FALSE)
  }
  if (!is.list(upper.bounds)) upper.bounds <- list(upper.bounds)
  if (length(upper.bounds) != NCOL(fixed.coordinates)) {
    stop("'upper.bounds' must contain one vector per fixed-coordinate column.",
         call. = FALSE)
  }
  valid.bounds <- vapply(
    upper.bounds,
    function(bounds) {
      length(bounds) > 0L && !is.unsorted(bounds, strictly = FALSE)
    },
    logical(1)
  )
  if (!all(valid.bounds)) {
    stop("Each element of 'upper.bounds' must be non-empty and non-decreasing.",
         call. = FALSE)
  }

  # The bounds already include the tolerance used during strip construction.
  # Values on a boundary remain in the lower interval of that coordinate.
  coordinate.strips <- vapply(
    seq_len(NCOL(fixed.coordinates)),
    function(j) {
      id <- findInterval(
        fixed.coordinates[, j], upper.bounds[[j]], left.open = TRUE
      ) + 1L
      pmin.int(id, length(upper.bounds[[j]]))
    },
    integer(NROW(fixed.coordinates))
  )
  coordinate.strips <- matrix(
    coordinate.strips,
    nrow = NROW(fixed.coordinates),
    ncol = NCOL(fixed.coordinates)
  )

  # Mixed-radix encoding gives one ID to each Cartesian cell. The first fixed
  # coordinate changes fastest; with one fixed coordinate the IDs are unchanged.
  cells.per.coordinate <- lengths(upper.bounds)
  total.cells <- prod(cells.per.coordinate)
  if (total.cells > .Machine$integer.max) {
    stop("The fixed-coordinate grid contains too many cells.", call. = FALSE)
  }
  multipliers <- c(1, cumprod(cells.per.coordinate)[-length(cells.per.coordinate)])
  as.integer(
    1 + rowSums(sweep(coordinate.strips - 1L, 2L, multipliers, `*`))
  )
}


fixedDimStripId <- function(fixed.values, upper.bounds)
{
  fixedDimCellId(
    fixed.coordinates = matrix(fixed.values, ncol = 1L),
    upper.bounds = list(upper.bounds)
  )
}


##########################################################
## Dimension-independent knot-placement adapter        ##
##########################################################
placeKnotND <- function(target.index, intknots, coordinates, residuals,
                        fixed.bounds, beta, dim.range = NULL)
{
  coordinates <- as.matrix(coordinates)
  target.index <- as.integer(target.index)

  if (NCOL(coordinates) < 2L) {
    stop("'coordinates' must contain at least two predictor columns.",
         call. = FALSE)
  }
  if (length(target.index) != 1L || is.na(target.index) ||
      target.index < 1L || target.index > NCOL(coordinates)) {
    stop("'target.index' must identify one column of 'coordinates'.",
         call. = FALSE)
  }

  fixed.index <- setdiff(seq_len(NCOL(coordinates)), target.index)
  if (!is.list(fixed.bounds) && length(fixed.index) == 1L) {
    fixed.bounds <- list(fixed.bounds)
  }
  if (!is.list(fixed.bounds) || length(fixed.bounds) != length(fixed.index)) {
    stop("'fixed.bounds' must contain one vector per non-target coordinate.",
         call. = FALSE)
  }
  if (is.null(dim.range)) dim.range <- range(coordinates[, target.index])

  cell.id <- fixedDimCellId(
    fixed.coordinates = coordinates[, fixed.index, drop = FALSE],
    upper.bounds = fixed.bounds
  )

  placeDimKnot(
    dim.index = target.index,
    intknots = intknots,
    coordinates = coordinates,
    residuals = residuals,
    strip.id = cell.id,
    beta = beta,
    dim.range = dim.range
  )
}


##########################################################
## One dimension-independent Normal Stage A iteration  ##
##########################################################
tensorProdND <- function(basis.matrices)
{
  if (!is.list(basis.matrices) || !length(basis.matrices)) {
    stop("'basis.matrices' must be a non-empty list.", call. = FALSE)
  }
  basis.matrices <- lapply(basis.matrices, as.matrix)
  nrows <- vapply(basis.matrices, NROW, integer(1))
  if (length(unique(nrows)) != 1L) {
    stop("All tensor-product basis matrices must have the same number of rows.",
         call. = FALSE)
  }
  if (length(basis.matrices) == 1L) return(basis.matrices[[1L]])

  Reduce(tensorProd, basis.matrices)
}


detectTensorMeshND <- function(coordinates)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  axes <- lapply(seq_len(ndim), function(j) sort(unique(coordinates[, j])))
  axis.index <- vapply(
    seq_len(ndim),
    function(j) match(coordinates[, j], axes[[j]]),
    integer(NROW(coordinates))
  )
  if (ndim == 1L) axis.index <- matrix(axis.index, ncol = 1L)
  axis.lengths <- lengths(axes)
  mesh.size <- prod(axis.lengths)
  # Scattered observations can imply an enormous hypothetical Cartesian grid.
  # Only construct cell identifiers when a complete mesh is possible; otherwise
  # unused identifiers can overflow R's integer range during grid detection.
  cell.id <- NULL
  complete <- FALSE
  if (mesh.size == NROW(coordinates)) {
    multipliers <- cumprod(c(1, utils::head(axis.lengths, -1L)))
    cell.id <- as.integer(1 + (axis.index - 1L) %*% multipliers)
    complete <- !anyDuplicated(cell.id) &&
      identical(sort(cell.id), seq_len(NROW(coordinates)))
  }

  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(axes) <- names(axis.lengths) <- dimension.names
  colnames(axis.index) <- dimension.names

  list(
    complete = complete,
    axes = axes,
    axis.index = axis.index,
    axis.lengths = axis.lengths,
    cell.id = cell.id
  )
}


modeMultiplyND <- function(array, matrix, mode)
{
  dimensions <- dim(array)
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  permuted <- aperm(array, permutation)
  permuted.dimensions <- dim(permuted)
  unfolded <- matrix(permuted, nrow = permuted.dimensions[1L])
  multiplied <- matrix %*% unfolded
  result <- array(
    multiplied,
    dim = c(NROW(matrix), permuted.dimensions[-1L])
  )
  aperm(result, order(permutation))
}


modeLeastSquaresND <- function(array, design, mode)
{
  dimensions <- dim(array)
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  permuted <- aperm(array, permutation)
  permuted.dimensions <- dim(permuted)
  unfolded <- matrix(permuted, nrow = permuted.dimensions[1L])
  solved <- qr.solve(design, unfolded)
  solved <- matrix(solved, nrow = NCOL(design))
  result <- array(solved, dim = c(NCOL(design), permuted.dimensions[-1L]))
  aperm(result, order(permutation))
}


fitTensorLeastSquaresND <- function(coordinates, response, basis.matrices,
                                    weights)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  if (length(response) != nobs || length(weights) != nobs) {
    stop("'response' and 'weights' must have one value per observation.",
         call. = FALSE)
  }
  if (!is.list(basis.matrices) || length(basis.matrices) != ndim ||
      any(vapply(basis.matrices, NROW, integer(1)) != nobs)) {
    stop("'basis.matrices' must contain one observation-level matrix per coordinate.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (!any(weights > 0)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  effective.nobs <- sum(weights > 0)
  mesh <- detectTensorMeshND(coordinates)
  equal.positive.weights <- weights[1L] > 0 && all(weights == weights[1L])

  if (mesh$complete && equal.positive.weights) {
    axis.bases <- lapply(
      seq_len(ndim),
      function(j) {
        first.rows <- match(seq_along(mesh$axes[[j]]), mesh$axis.index[, j])
        basis.matrices[[j]][first.rows, , drop = FALSE]
      }
    )
    full.rank <- all(vapply(
      axis.bases,
      function(basis) NROW(basis) >= NCOL(basis) &&
        qr(basis)$rank == NCOL(basis),
      logical(1)
    ))

    if (full.rank) {
      mesh.response <- response[match(seq_len(nobs), mesh$cell.id)]
      coefficient.array <- array(mesh.response, dim = mesh$axis.lengths)
      for (j in seq_len(ndim)) {
        coefficient.array <- modeLeastSquaresND(
          coefficient.array, axis.bases[[j]], j
        )
      }
      fitted.array <- coefficient.array
      for (j in seq_len(ndim)) {
        fitted.array <- modeMultiplyND(fitted.array, axis.bases[[j]], j)
      }
      fitted.mesh <- as.vector(fitted.array)
      coefficients <- as.vector(
        aperm(coefficient.array, rev(seq_len(ndim)))
      )
      predicted <- fitted.mesh[mesh$cell.id]
      return(list(
        coefficients = coefficients,
        predicted = predicted,
        design = NULL,
        design.dim = c(
          as.integer(nobs),
          as.integer(prod(vapply(basis.matrices, NCOL, integer(1))))
        ),
        rank = as.integer(length(coefficients)),
        rank.deficient = FALSE,
        effective.nobs = as.integer(effective.nobs),
        solver = "tensor-mesh-qr",
        mesh = mesh
      ))
    }
  }

  design <- tensorProdND(basis.matrices)
  if (all(weights == 1)) {
    fit <- .lm.fit(design, response)
    if (fit$rank < NCOL(design)) fit <- lm.fit(design, response)
  } else {
    fit <- lm.wfit.light(design, response, weights)
    if (fit$rank < NCOL(design)) fit <- lm.wfit(design, response, weights)
  }
  coefficients <- as.numeric(coef(fit))
  coefficients[is.na(coefficients)] <- 0
  list(
    coefficients = coefficients,
    predicted = as.numeric(design %*% coefficients),
    design = design,
    design.dim = dim(design),
    rank = as.integer(fit$rank),
    rank.deficient = fit$rank < NCOL(design),
    effective.nobs = as.integer(effective.nobs),
    solver = "dense",
    mesh = mesh
  )
}


stageAOneStepND <- function(coordinates, response, upper.bounds,
                            intknots = NULL, coordinate.ranges = NULL,
                            placement.ranges = NULL,
                            weights = rep(1, NROW(coordinates)), beta = 0.5,
                            max.coef = Inf)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)

  if (ndim < 2L) {
    stop("'coordinates' must contain at least two predictor columns.",
         call. = FALSE)
  }
  if (length(response) != nobs || length(weights) != nobs) {
    stop("'response' and 'weights' must have one value per observation.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (!any(weights > 0)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  if (is.null(intknots)) intknots <- rep(list(NULL), ndim)
  if (!is.list(intknots) || length(intknots) != ndim) {
    stop("'intknots' must contain one element per coordinate.", call. = FALSE)
  }
  if (!is.list(upper.bounds) || length(upper.bounds) != ndim) {
    stop("'upper.bounds' must contain one vector per coordinate.",
         call. = FALSE)
  }
  if (is.null(coordinate.ranges)) {
    coordinate.ranges <- lapply(
      seq_len(ndim), function(j) range(coordinates[, j])
    )
  }
  valid.ranges <- is.list(coordinate.ranges) &&
    length(coordinate.ranges) == ndim &&
    all(vapply(coordinate.ranges, length, integer(1)) == 2L)
  if (!valid.ranges) {
    stop("'coordinate.ranges' must contain one two-value range per coordinate.",
         call. = FALSE)
  }
  if (is.null(placement.ranges)) {
    placement.ranges <- lapply(
      seq_len(ndim), function(j) range(coordinates[, j])
    )
  }
  valid.placement.ranges <- is.list(placement.ranges) &&
    length(placement.ranges) == ndim &&
    all(vapply(placement.ranges, length, integer(1)) == 2L)
  if (!valid.placement.ranges) {
    stop("'placement.ranges' must contain one two-value range per coordinate.",
         call. = FALSE)
  }

  basis.matrices <- lapply(
    seq_len(ndim),
    function(j) {
      splines::splineDesign(
        knots = sort(c(intknots[[j]], rep(coordinate.ranges[[j]], 2L))),
        x = coordinates[, j],
        ord = 2L,
        derivs = rep(0L, nobs),
        outer.ok = TRUE
      )
    }
  )
  basis.size <- prod(as.double(vapply(basis.matrices, NCOL, integer(1))))
  if (basis.size > max.coef) {
    stop("Tensor basis requires ", format(basis.size, scientific = FALSE),
         " coefficients, exceeding 'max.coef' = ",
         format(max.coef, scientific = FALSE), ".", call. = FALSE)
  }
  least.squares <- fitTensorLeastSquaresND(
    coordinates = coordinates,
    response = response,
    basis.matrices = basis.matrices,
    weights = weights
  )
  design <- least.squares$design
  coefficients <- least.squares$coefficients
  predicted <- least.squares$predicted
  residuals <- response - predicted
  placement.residuals <- residuals * weights
  positive.weight <- weights > 0

  candidates <- lapply(
    seq_len(ndim),
    function(target.index) {
      placeKnotND(
        target.index = target.index,
        intknots = intknots[[target.index]],
        coordinates = coordinates[positive.weight, , drop = FALSE],
        residuals = placement.residuals[positive.weight],
        fixed.bounds = upper.bounds[-target.index],
        beta = beta,
        dim.range = placement.ranges[[target.index]]
      )
    }
  )
  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(candidates) <- dimension.names

  candidate.weights <- vapply(candidates, `[[`, numeric(1), "weightDim")
  candidate.flags <- vapply(candidates, `[[`, logical(1), "flagDim")
  admissible <- which(!candidate.flags & !is.na(candidate.weights))

  selected.index <- NA_integer_
  if (length(admissible)) {
    maximum.weight <- max(candidate.weights[admissible])
    tie.tolerance <- 1e-12 * max(1, abs(maximum.weight))
    tied <- admissible[
      abs(candidate.weights[admissible] - maximum.weight) <= tie.tolerance
    ]
    knot.counts <- lengths(intknots)
    tied <- tied[knot.counts[tied] == max(knot.counts[tied])]
    # This reproduces the bivariate tie rule: with equal weights and equal
    # knot counts, the later coordinate is selected.
    selected.index <- max(tied)
  }

  selected <- if (is.na(selected.index)) {
    list(index = NA_integer_, dimension = NA_character_, newknot = NA_real_,
         weight = NA_real_, flag = TRUE)
  } else {
    list(
      index = selected.index,
      dimension = dimension.names[selected.index],
      newknot = candidates[[selected.index]]$Dim.newknot,
      weight = candidates[[selected.index]]$weightDim,
      flag = FALSE
    )
  }

  list(
    basis.matrices = basis.matrices,
    design = design,
    design.dim = least.squares$design.dim,
    rank = least.squares$rank,
    rank.deficient = least.squares$rank.deficient,
    effective.nobs = least.squares$effective.nobs,
    solver = least.squares$solver,
    mesh = least.squares$mesh,
    basis.size = basis.size,
    coefficients = coefficients,
    predicted = predicted,
    residuals = residuals,
    rss = .weighted_rss(residuals, weights),
    candidates = candidates,
    selected = selected
  )
}


##########################################################
## Bounded dimension-independent Normal Stage A loop   ##
##########################################################
stageALoopND <- function(coordinates, response, upper.bounds,
                         intknots = NULL, coordinate.ranges = NULL,
                         placement.ranges = NULL,
                         weights = rep(1, NROW(coordinates)), beta = 0.5,
                         max.steps = 10L, stop.rule = c("none", "SR", "RD"),
                         phi = 0.99, q = 2L, min.intknots = 0L,
                         max.coef = 100000L)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  max.steps <- as.integer(max.steps)
  stop.rule <- match.arg(stop.rule)
  q <- as.integer(q)
  min.intknots <- as.integer(min.intknots)
  max.coef <- as.numeric(max.coef)

  if (ndim < 2L) {
    stop("'coordinates' must contain at least two predictor columns.",
         call. = FALSE)
  }
  if (length(response) != nobs || length(weights) != nobs) {
    stop("'response' and 'weights' must have one value per observation.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (!any(weights > 0)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  effective.nobs <- sum(weights > 0)
  if (length(max.steps) != 1L || is.na(max.steps) || max.steps < 1L) {
    stop("'max.steps' must be a positive integer.", call. = FALSE)
  }
  if (length(max.coef) != 1L || is.na(max.coef) || max.coef < 1) {
    stop("'max.coef' must be a positive number.", call. = FALSE)
  }
  if (is.null(intknots)) intknots <- rep(list(NULL), ndim)
  if (!is.list(intknots) || length(intknots) != ndim) {
    stop("'intknots' must contain one element per coordinate.", call. = FALSE)
  }
  if (length(q) != 1L || is.na(q) || q < 1L) {
    stop("'q' must be a positive integer.", call. = FALSE)
  }
  if (length(min.intknots) != 1L || is.na(min.intknots) ||
      min.intknots < 0L) {
    stop("'min.intknots' must be a non-negative integer.", call. = FALSE)
  }
  if (length(phi) != 1L || is.na(phi) || phi <= 0 || phi >= 1) {
    stop("'phi' must be a single value in (0, 1).", call. = FALSE)
  }
  if (stop.rule == "SR" && sum(lengths(intknots)) != 0L) {
    stop("The initial SR implementation currently requires NULL internal knots.",
         call. = FALSE)
  }

  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(intknots) <- dimension.names

  history <- vector("list", max.steps)
  insertions <- vector("list", max.steps)
  completed <- 0L
  stop.reason <- "maximum steps reached"
  phis <- phis.star <- numeric()
  gamma0 <- gamma1 <- numeric()
  selected.iteration <- NA_integer_

  for (iteration in seq_len(max.steps)) {
    linear.basis.size <- prod(as.double(lengths(intknots) + 2L))
    if (linear.basis.size > max.coef) {
      stop.reason <- "maximum coefficient limit reached"
      break
    }
    if (linear.basis.size >= effective.nobs) {
      stop.reason <- "linear tensor basis is saturated"
      break
    }

    step <- stageAOneStepND(
      coordinates = coordinates,
      response = response,
      upper.bounds = upper.bounds,
      intknots = intknots,
      coordinate.ranges = coordinate.ranges,
      placement.ranges = placement.ranges,
      weights = weights,
      beta = beta,
      max.coef = max.coef
    )
    step$iteration <- iteration
    step$intknots <- lapply(intknots, identity)
    history[[iteration]] <- step

    if (stop.rule != "none" && iteration > q) {
      earlier <- history[[iteration - q]]
      rss.ratio <- step$rss / earlier$rss
      if (rss.ratio > 1) {
        stop.reason <- "RSS increased"
        selected.iteration <- iteration - q
        break
      }

      coefficient.increment <- length(step$coefficients) -
        length(earlier$coefficients)
      phnew <- rss.ratio^(1 / coefficient.increment)
      phis <- c(phis, phnew)

      if (iteration - q > min.intknots) {
        if (stop.rule == "RD" && isTRUE(rss.ratio >= phi)) {
          stop.reason <- "RD threshold reached"
          selected.iteration <- iteration - q
          break
        }
        if (stop.rule == "SR") {
          phismod <- log(1 - phis)
          gamma <- .lm.fit(
            cbind(1, (q + 1L):iteration),
            phismod
          )$coefficients
          kappa <- sum(lengths(intknots))
          phi.kappa <- 1 - exp(gamma[1L]) * exp(gamma[2L] * kappa)
          phis.star <- c(phis.star, phi.kappa)
          gamma0 <- c(gamma0, gamma[1L])
          gamma1 <- c(gamma1, gamma[2L])

          if (isTRUE(phi.kappa >= phi)) {
            stop.reason <- "SR threshold reached"
            selected.iteration <- iteration - q
            break
          }
        }
      }
    }

    if (isTRUE(step$selected$flag) || is.na(step$selected$index)) {
      stop.reason <- "no admissible knot"
      break
    }

    target.index <- step$selected$index
    newknot <- step$selected$newknot
    if (!is.finite(newknot)) {
      stop.reason <- "selected knot is not finite"
      break
    }
    existing <- intknots[[target.index]]
    if (length(existing) && any(abs(existing - newknot) < 1e-12)) {
      stop.reason <- "selected knot duplicates an existing knot"
      break
    }

    intknots[[target.index]] <- sort(c(existing, newknot))
    completed <- completed + 1L
    insertions[[completed]] <- list(
      iteration = iteration,
      index = target.index,
      dimension = dimension.names[target.index],
      knot = newknot,
      weight = step$selected$weight
    )
  }

  fitted.steps <- sum(!vapply(history, is.null, logical(1)))
  history <- history[seq_len(fitted.steps)]
  insertions <- insertions[seq_len(completed)]
  if (is.na(selected.iteration) && fitted.steps > 0L) {
    selected.iteration <- fitted.steps
  }
  selected.intknots <- if (is.na(selected.iteration)) {
    setNames(rep(list(NULL), ndim), dimension.names)
  } else {
    history[[selected.iteration]]$intknots
  }

  list(
    intknots = intknots,
    selected.intknots = selected.intknots,
    selected.iteration = selected.iteration,
    history = history,
    insertions = insertions,
    rss = vapply(history, `[[`, numeric(1), "rss"),
    ncoef = vapply(history, function(step) length(step$coefficients), integer(1)),
    basis.sizes = vapply(history, `[[`, numeric(1), "basis.size"),
    effective.nobs = as.integer(effective.nobs),
    max.coef = max.coef,
    phis = phis,
    phis.star = phis.star,
    gamma = list(intercept = gamma0, slope = gamma1),
    completed = completed,
    stop.reason = stop.reason
  )
}


##########################################################
## Dimension-independent Normal Stage B                ##
##########################################################
fitTensorSplineND <- function(coordinates, response, intknots = NULL,
                              spline.order = 2L,
                              coordinate.ranges = NULL,
                              weights = rep(1, NROW(coordinates)),
                              max.coef = 100000L)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  spline.order <- as.integer(spline.order)

  if (ndim < 2L) {
    stop("'coordinates' must contain at least two predictor columns.",
         call. = FALSE)
  }
  if (length(response) != nobs || length(weights) != nobs) {
    stop("'response' and 'weights' must have one value per observation.",
         call. = FALSE)
  }
  if (length(spline.order) != 1L || is.na(spline.order) ||
      spline.order < 2L) {
    stop("'spline.order' must be an integer greater than or equal to two.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (!any(weights > 0)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  max.coef <- as.numeric(max.coef)
  if (length(max.coef) != 1L || is.na(max.coef) || max.coef < 1) {
    stop("'max.coef' must be a positive number.", call. = FALSE)
  }
  if (is.null(intknots)) intknots <- rep(list(NULL), ndim)
  if (!is.list(intknots) || length(intknots) != ndim) {
    stop("'intknots' must contain one element per coordinate.", call. = FALSE)
  }
  if (is.null(coordinate.ranges)) {
    coordinate.ranges <- lapply(
      seq_len(ndim), function(j) range(coordinates[, j])
    )
  }
  valid.ranges <- is.list(coordinate.ranges) &&
    length(coordinate.ranges) == ndim &&
    all(vapply(coordinate.ranges, length, integer(1)) == 2L)
  if (!valid.ranges) {
    stop("'coordinate.ranges' must contain one two-value range per coordinate.",
         call. = FALSE)
  }

  full.knots <- lapply(
    seq_len(ndim),
    function(j) {
      sort(c(intknots[[j]], rep(coordinate.ranges[[j]], spline.order)))
    }
  )
  basis.matrices <- lapply(
    seq_len(ndim),
    function(j) {
      splines::splineDesign(
        knots = full.knots[[j]],
        x = coordinates[, j],
        ord = spline.order,
        derivs = rep(0L, nobs),
        outer.ok = TRUE
      )
    }
  )
  basis.size <- prod(as.double(vapply(basis.matrices, NCOL, integer(1))))
  if (basis.size > max.coef) {
    stop("Tensor basis requires ", format(basis.size, scientific = FALSE),
         " coefficients, exceeding 'max.coef' = ",
         format(max.coef, scientific = FALSE), ".", call. = FALSE)
  }
  least.squares <- fitTensorLeastSquaresND(
    coordinates = coordinates,
    response = response,
    basis.matrices = basis.matrices,
    weights = weights
  )
  design <- least.squares$design
  coefficients <- least.squares$coefficients
  predicted <- least.squares$predicted
  residuals <- response - predicted

  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(intknots) <- names(full.knots) <- names(basis.matrices) <- dimension.names

  list(
    spline.order = spline.order,
    intknots = intknots,
    full.knots = full.knots,
    basis.matrices = basis.matrices,
    design = design,
    design.dim = least.squares$design.dim,
    rank = least.squares$rank,
    rank.deficient = least.squares$rank.deficient,
    effective.nobs = least.squares$effective.nobs,
    solver = least.squares$solver,
    mesh = least.squares$mesh,
    basis.size = basis.size,
    coefficients = coefficients,
    predicted = predicted,
    residuals = residuals,
    rss = .weighted_rss(residuals, weights)
  )
}


##########################################################
## Dimension-independent generalized tensor fitting    ##
##########################################################
fitTensorGLMND <- function(coordinates, response, intknots = NULL,
                           spline.order = 2L,
                           coordinate.ranges = NULL,
                           weights = rep(1, NROW(coordinates)),
                           family = stats::gaussian(),
                           offset = rep(0, NROW(coordinates)),
                           parametric = NULL, mustart = NULL,
                           max.coef = 100000L)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  spline.order <- as.integer(spline.order)
  if (ndim < 2L) {
    stop("'coordinates' must contain at least two predictor columns.",
         call. = FALSE)
  }
  if (length(response) != nobs || length(weights) != nobs ||
      length(offset) != nobs) {
    stop("'response', 'weights', and 'offset' must have one value per observation.",
         call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      !any(weights > 0)) {
    stop("'weights' must be finite, non-negative, and contain a positive value.",
         call. = FALSE)
  }
  if (length(spline.order) != 1L || is.na(spline.order) || spline.order < 2L) {
    stop("'spline.order' must be an integer greater than or equal to two.",
         call. = FALSE)
  }
  if (is.null(intknots)) intknots <- rep(list(NULL), ndim)
  if (!is.list(intknots) || length(intknots) != ndim) {
    stop("'intknots' must contain one element per coordinate.", call. = FALSE)
  }
  if (is.null(coordinate.ranges)) {
    coordinate.ranges <- lapply(seq_len(ndim), function(j) range(coordinates[, j]))
  }
  if (!is.list(coordinate.ranges) || length(coordinate.ranges) != ndim ||
      any(lengths(coordinate.ranges) != 2L)) {
    stop("'coordinate.ranges' must contain one two-value range per coordinate.",
         call. = FALSE)
  }
  max.coef <- as.numeric(max.coef)
  full.knots <- lapply(seq_len(ndim), function(j) {
    sort(c(intknots[[j]], rep(coordinate.ranges[[j]], spline.order)))
  })
  basis.matrices <- lapply(seq_len(ndim), function(j) {
    splines::splineDesign(
      knots = full.knots[[j]], x = coordinates[, j], ord = spline.order,
      derivs = rep(0L, nobs), outer.ok = TRUE
    )
  })
  basis.size <- prod(as.double(vapply(basis.matrices, NCOL, integer(1))))
  if (basis.size > max.coef) {
    stop("Tensor basis requires ", format(basis.size, scientific = FALSE),
         " coefficients, exceeding 'max.coef' = ",
         format(max.coef, scientific = FALSE), ".", call. = FALSE)
  }
  design <- tensorProdND(basis.matrices)
  spline.basis.size <- NCOL(design)
  if (!is.null(parametric)) {
    parametric <- as.matrix(parametric)
    if (NROW(parametric) != nobs || !is.numeric(parametric)) {
      stop("'parametric' must be a numeric matrix with one row per observation.", call. = FALSE)
    }
    design <- cbind(design, parametric)
  }
  basis.size <- NCOL(design)
  if (basis.size > max.coef) {
    stop("Combined design requires ", basis.size,
         " coefficients, exceeding 'max.coef' = ", max.coef, ".",
         call. = FALSE)
  }
  fit <- stats::glm.fit(
    x = design, y = response, weights = weights, offset = offset,
    family = family, intercept = FALSE, mustart = mustart
  )
  if (anyNA(fit$coefficients)) {
    stop("The generalized tensor basis is rank deficient.", call. = FALSE)
  }
  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(intknots) <- names(full.knots) <- names(basis.matrices) <- dimension.names
  list(
    spline.order = spline.order,
    intknots = intknots,
    full.knots = full.knots,
    basis.matrices = basis.matrices,
    design = design,
    design.dim = dim(design),
    rank = fit$rank,
    rank.deficient = fit$rank < NCOL(design),
    effective.nobs = sum(weights > 0),
    solver = "IRLS",
    mesh = detectTensorMeshND(coordinates),
    basis.size = basis.size,
    spline.basis.size = spline.basis.size,
    parametric.names = colnames(parametric),
    coefficients = fit$coefficients,
    theta = fit$coefficients,
    predicted = as.numeric(fit$fitted.values),
    linear.predictors = as.numeric(fit$linear.predictors),
    residuals = as.numeric(fit$residuals),
    working.weights = as.numeric(fit$weights),
    rss = fit$deviance,
    deviance = fit$deviance,
    family = family,
    offset = offset,
    temporary = fit
  )
}


stageAGLMLoopND <- function(coordinates, response, upper.bounds,
                            family, offset = rep(0, NROW(coordinates)),
                            parametric = NULL,
                            intknots = NULL, coordinate.ranges = NULL,
                            placement.ranges = NULL,
                            weights = rep(1, NROW(coordinates)), beta = 0.5,
                            max.steps = 10L,
                            stop.rule = c("none", "SR", "RD", "LR"),
                            phi = 0.99, q = 2L, min.intknots = 0L,
                            max.coef = 100000L)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  stop.rule <- match.arg(stop.rule)
  max.steps <- as.integer(max.steps)
  q <- as.integer(q)
  min.intknots <- as.integer(min.intknots)
  if (is.null(intknots)) intknots <- rep(list(NULL), ndim)
  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(intknots) <- dimension.names
  effective.nobs <- sum(weights > 0)
  history <- vector("list", max.steps)
  insertions <- vector("list", max.steps)
  phis <- phis.star <- gamma0 <- gamma1 <- numeric()
  completed <- 0L
  selected.iteration <- NA_integer_
  stop.reason <- "maximum steps reached"
  mustart <- NULL

  for (iteration in seq_len(max.steps)) {
    basis.size <- prod(as.double(lengths(intknots) + 2L)) + NCOL(parametric)
    if (basis.size > max.coef) {
      stop.reason <- "maximum coefficient limit reached"
      break
    }
    if (basis.size >= effective.nobs) {
      stop.reason <- "linear tensor basis is saturated"
      break
    }
    fit <- fitTensorGLMND(
      coordinates = coordinates, response = response, intknots = intknots,
      spline.order = 2L, coordinate.ranges = coordinate.ranges,
      weights = weights, family = family, offset = offset,
      parametric = parametric, mustart = mustart, max.coef = max.coef
    )
    mustart <- fit$predicted
    placement.residuals <- fit$residuals * fit$working.weights * weights
    positive <- weights > 0 & is.finite(placement.residuals)
    candidates <- lapply(seq_len(ndim), function(target.index) {
      placeKnotND(
        target.index = target.index,
        intknots = intknots[[target.index]],
        coordinates = coordinates[positive, , drop = FALSE],
        residuals = placement.residuals[positive],
        fixed.bounds = upper.bounds[-target.index], beta = beta,
        dim.range = placement.ranges[[target.index]]
      )
    })
    candidate.weights <- vapply(candidates, `[[`, numeric(1), "weightDim")
    candidate.flags <- vapply(candidates, `[[`, logical(1), "flagDim")
    admissible <- which(!candidate.flags & is.finite(candidate.weights))
    selected.index <- NA_integer_
    if (length(admissible)) {
      maximum.weight <- max(candidate.weights[admissible])
      tied <- admissible[abs(candidate.weights[admissible] - maximum.weight) <=
                           1e-12 * max(1, abs(maximum.weight))]
      tied <- tied[lengths(intknots)[tied] == max(lengths(intknots)[tied])]
      selected.index <- max(tied)
    }
    selected <- list(
      index = selected.index,
      dimension = if (is.na(selected.index)) NA_character_ else dimension.names[selected.index],
      newknot = if (is.na(selected.index)) NA_real_ else candidates[[selected.index]]$Dim.newknot,
      weight = if (is.na(selected.index)) NA_real_ else candidate.weights[selected.index],
      flag = is.na(selected.index)
    )
    history[[iteration]] <- c(fit, list(
      iteration = iteration, intknots = lapply(intknots, identity),
      candidates = candidates, selected = selected
    ))

    if (stop.rule != "none" && iteration > q) {
      earlier <- history[[iteration - q]]
      if (fit$rss > earlier$rss) {
        stop.reason <- "deviance increased"
        selected.iteration <- iteration - q
        break
      }
      coefficient.increment <- fit$basis.size - earlier$basis.size
      if (stop.rule == "LR") {
        phnew <- earlier$rss - fit$rss
      } else {
        phnew <- (fit$rss / earlier$rss)^(1 / coefficient.increment)
      }
      phis <- c(phis, phnew)
      if (iteration - q > min.intknots) {
        should.stop <- FALSE
        if (stop.rule == "RD") should.stop <- fit$rss / earlier$rss >= phi^coefficient.increment
        if (stop.rule == "LR") should.stop <- earlier$rss - fit$rss < stats::qchisq(phi, coefficient.increment)
        if (stop.rule == "SR") {
          gamma <- stats::.lm.fit(
            cbind(1, (q + 1L):iteration), log(1 - phis)
          )$coefficients
          kappa <- sum(lengths(intknots))
          phi.kappa <- 1 - exp(gamma[1L]) * exp(gamma[2L] * kappa)
          phis.star <- c(phis.star, phi.kappa)
          gamma0 <- c(gamma0, gamma[1L]); gamma1 <- c(gamma1, gamma[2L])
          should.stop <- isTRUE(phi.kappa >= phi)
        }
        if (isTRUE(should.stop)) {
          stop.reason <- paste0(stop.rule, " threshold reached")
          selected.iteration <- iteration - q
          break
        }
      }
    }
    if (selected$flag || !is.finite(selected$newknot)) {
      stop.reason <- "no admissible knot"
      break
    }
    intknots[[selected.index]] <- sort(c(
      intknots[[selected.index]], selected$newknot
    ))
    completed <- completed + 1L
    insertions[[completed]] <- list(
      iteration = iteration, index = selected.index,
      dimension = selected$dimension, knot = selected$newknot,
      weight = selected$weight
    )
  }
  fitted.steps <- sum(!vapply(history, is.null, logical(1)))
  history <- history[seq_len(fitted.steps)]
  insertions <- insertions[seq_len(completed)]
  if (is.na(selected.iteration) && fitted.steps) selected.iteration <- fitted.steps
  selected.intknots <- if (is.na(selected.iteration)) {
    setNames(rep(list(NULL), ndim), dimension.names)
  } else history[[selected.iteration]]$intknots
  list(
    intknots = intknots, selected.intknots = selected.intknots,
    selected.iteration = selected.iteration, history = history,
    insertions = insertions,
    rss = vapply(history, `[[`, numeric(1), "rss"),
    ncoef = vapply(history, `[[`, numeric(1), "basis.size"),
    basis.sizes = vapply(history, `[[`, numeric(1), "basis.size"),
    effective.nobs = as.integer(effective.nobs), max.coef = max.coef,
    phis = phis, phis.star = phis.star,
    gamma = list(intercept = gamma0, slope = gamma1),
    completed = completed, stop.reason = stop.reason
  )
}


stageBGLMND <- function(coordinates, response, intknots, family,
                        offset = rep(0, NROW(coordinates)),
                        parametric = NULL,
                        coordinate.ranges = NULL,
                        weights = rep(1, NROW(coordinates)),
                        spline.orders = 2:4, max.coef = 100000L)
{
  spline.orders <- as.integer(spline.orders)
  order.names <- vapply(spline.orders, function(n) {
    if (n == 2L) "linear" else if (n == 3L) "quadratic" else if (n == 4L) "cubic" else paste0("order", n)
  }, character(1))
  transformed.knots <- lapply(spline.orders, function(n) stageBKnotsND(intknots, n))
  basis.sizes <- vapply(seq_along(spline.orders), function(i) {
    prod(as.double(lengths(transformed.knots[[i]]) + spline.orders[i])) + NCOL(parametric)
  }, numeric(1))
  fits <- vector("list", length(spline.orders))
  mustart <- NULL
  for (i in seq_along(spline.orders)) {
    if (basis.sizes[i] <= max.coef) {
      fits[[i]] <- fitTensorGLMND(
        coordinates, response, transformed.knots[[i]], spline.orders[i],
        coordinate.ranges, weights, family, offset, parametric, mustart, max.coef
      )
      mustart <- fits[[i]]$predicted
    }
  }
  names(transformed.knots) <- names(fits) <- names(basis.sizes) <- order.names
  skipped <- ifelse(basis.sizes > max.coef,
                    paste0("requires ", basis.sizes, " coefficients; max.coef = ", max.coef),
                    NA_character_)
  names(skipped) <- order.names
  omitted <- names(skipped)[!is.na(skipped)]
  if (length(omitted)) {
    warning("Omitting Stage B ", paste(omitted, collapse = " and "),
            " fit(s) because the tensor coefficient limit would be exceeded.",
            call. = FALSE)
  }
  list(stageA.intknots = intknots, intknots = transformed.knots,
       fits = fits, basis.sizes = basis.sizes, skipped = skipped,
       max.coef = max.coef)
}


stageBKnotsND <- function(intknots, spline.order)
{
  spline.order <- as.integer(spline.order)
  if (!is.list(intknots) || !length(intknots)) {
    stop("'intknots' must be a non-empty list.", call. = FALSE)
  }
  if (length(spline.order) != 1L || is.na(spline.order) ||
      spline.order < 2L) {
    stop("'spline.order' must be an integer greater than or equal to two.",
         call. = FALSE)
  }

  required.knots <- spline.order - 1L
  transformed <- lapply(
    intknots,
    function(knots) {
      if (length(knots) < required.knots) return(NULL)
      makenewknots(sort(knots), spline.order)
    }
  )
  names(transformed) <- names(intknots)
  transformed
}


stageBND <- function(coordinates, response, intknots,
                     coordinate.ranges = NULL,
                     weights = rep(1, NROW(coordinates)),
                     spline.orders = 2:4, max.coef = 100000L)
{
  spline.orders <- as.integer(spline.orders)
  max.coef <- as.numeric(max.coef)
  if (!length(spline.orders) || anyNA(spline.orders) ||
      any(spline.orders < 2L) || anyDuplicated(spline.orders)) {
    stop("'spline.orders' must contain distinct integers of at least two.",
         call. = FALSE)
  }
  if (length(max.coef) != 1L || is.na(max.coef) || max.coef < 1) {
    stop("'max.coef' must be a positive number.", call. = FALSE)
  }

  order.names <- vapply(
    spline.orders,
    function(spline.order) {
      if (spline.order == 2L) return("linear")
      if (spline.order == 3L) return("quadratic")
      if (spline.order == 4L) return("cubic")
      paste0("order", spline.order)
    },
    character(1)
  )
  transformed.knots <- lapply(
    spline.orders,
    function(spline.order) stageBKnotsND(intknots, spline.order)
  )
  basis.sizes <- vapply(
    seq_along(spline.orders),
    function(i) prod(as.double(lengths(transformed.knots[[i]]) +
                                  spline.orders[i])),
    numeric(1)
  )
  fits <- lapply(
    seq_along(spline.orders),
    function(i) {
      if (basis.sizes[i] > max.coef) return(NULL)
      fitTensorSplineND(
        coordinates = coordinates,
        response = response,
        intknots = transformed.knots[[i]],
        spline.order = spline.orders[i],
        coordinate.ranges = coordinate.ranges,
        weights = weights,
        max.coef = max.coef
      )
    }
  )
  names(transformed.knots) <- names(fits) <- order.names
  names(basis.sizes) <- order.names
  skipped <- ifelse(
    basis.sizes > max.coef,
    paste0("requires ", format(basis.sizes, scientific = FALSE, trim = TRUE),
           " coefficients; max.coef = ",
           format(max.coef, scientific = FALSE, trim = TRUE)),
    NA_character_
  )
  names(skipped) <- order.names
  omitted <- names(skipped)[!is.na(skipped)]
  if (length(omitted)) {
    warning(
      "Omitting Stage B ",
      paste(omitted, collapse = " and "),
      if (length(omitted) == 1L) " fit" else " fits",
      " because ",
      paste0(omitted, " ", skipped[omitted], collapse = "; "),
      ". Increase 'max.coef' only if sufficient memory is available.",
      call. = FALSE
    )
  }

  list(
    stageA.intknots = intknots,
    intknots = transformed.knots,
    fits = fits,
    basis.sizes = basis.sizes,
    skipped = skipped,
    max.coef = max.coef
  )
}


##########################################################
## Experimental Normal multivariate GeDS fitter        ##
##########################################################
makeGridBoundsND <- function(coordinates, coordinate.ranges = NULL,
                             nint = NULL, tolerance = 1e-15)
{
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)

  if (ndim < 2L || nobs < 1L) {
    stop("'coordinates' must contain observations on at least two dimensions.",
         call. = FALSE)
  }
  if (is.null(coordinate.ranges)) {
    coordinate.ranges <- lapply(
      seq_len(ndim), function(j) range(coordinates[, j])
    )
  }
  valid.ranges <- is.list(coordinate.ranges) &&
    length(coordinate.ranges) == ndim &&
    all(vapply(coordinate.ranges, length, integer(1)) == 2L)
  if (!valid.ranges) {
    stop("'coordinate.ranges' must contain one two-value range per coordinate.",
         call. = FALSE)
  }

  if (is.null(nint)) {
    root <- nobs^(1 / ndim)
    nearest.integer <- round(root)
    common.nint <- if (abs(root - nearest.integer) < 1e-10) {
      nearest.integer
    } else {
      floor(root)
    }
    nint <- rep(common.nint, ndim)
  } else if (length(nint) == 1L) {
    nint <- rep(nint, ndim)
  }
  nint <- as.integer(nint)
  if (length(nint) != ndim || anyNA(nint) || any(nint < 1L)) {
    stop("'nint' must provide a positive number of intervals per coordinate.",
         call. = FALSE)
  }
  if (length(tolerance) != 1L || is.na(tolerance) || tolerance < 0) {
    stop("'tolerance' must be a single non-negative value.", call. = FALSE)
  }

  upper.bounds <- lapply(
    seq_len(ndim),
    function(j) {
      seq(
        from = coordinate.ranges[[j]][1L],
        to = coordinate.ranges[[j]][2L],
        length.out = nint[j] + 1L
      )[-1L] + tolerance
    }
  )
  dimension.names <- colnames(coordinates)
  if (is.null(dimension.names)) dimension.names <- paste0("X", seq_len(ndim))
  names(coordinate.ranges) <- names(upper.bounds) <- names(nint) <- dimension.names

  list(
    nint = nint,
    coordinate.ranges = coordinate.ranges,
    upper.bounds = upper.bounds,
    tolerance = tolerance
  )
}


MultivariateFitter <- function(coordinates, response,
                               weights = rep(1, NROW(coordinates)),
                               beta = 0.5, phi = 0.99, q = 2L,
                               min.intknots = 0L, max.steps = 300L,
                               stoptype = c("SR", "RD"),
                               coordinate.ranges = NULL,
                               placement.ranges = NULL,
                               nint = NULL, intknots_init = NULL,
                               spline.orders = 2:4,
                               max.coef = 100000L)
{
  call <- match.call()
  stoptype <- match.arg(stoptype)
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  max.coef <- as.numeric(max.coef)
  if (length(max.coef) != 1L || is.na(max.coef) || max.coef < 1) {
    stop("'max.coef' must be a positive number.", call. = FALSE)
  }
  if (is.null(colnames(coordinates))) {
    colnames(coordinates) <- paste0("X", seq_len(ndim))
  }
  if (is.null(intknots_init)) intknots_init <- rep(list(NULL), ndim)
  if (!is.list(intknots_init) || length(intknots_init) != ndim) {
    stop("'intknots_init' must contain one component per coordinate.", call. = FALSE)
  }
  names(intknots_init) <- colnames(coordinates)
  initial.basis.size <- prod(as.double(lengths(intknots_init) + 2L))
  if (initial.basis.size > max.coef) {
    stop("The initial linear tensor basis requires ",
         format(initial.basis.size, scientific = FALSE),
         " coefficients, exceeding 'max.coef' = ",
         format(max.coef, scientific = FALSE), ".", call. = FALSE)
  }

  grid <- makeGridBoundsND(
    coordinates = coordinates,
    coordinate.ranges = coordinate.ranges,
    nint = nint
  )
  if (is.null(placement.ranges)) {
    placement.ranges <- lapply(
      seq_len(ndim), function(j) range(coordinates[, j])
    )
  }

  stageA <- stageALoopND(
    coordinates = coordinates,
    response = response,
    upper.bounds = grid$upper.bounds,
    intknots = intknots_init,
    coordinate.ranges = grid$coordinate.ranges,
    placement.ranges = placement.ranges,
    weights = weights,
    beta = beta,
    max.steps = max.steps,
    stop.rule = stoptype,
    phi = phi,
    q = q,
    min.intknots = min.intknots,
    max.coef = max.coef
  )
  stageB <- stageBND(
    coordinates = coordinates,
    response = response,
    intknots = stageA$selected.intknots,
    coordinate.ranges = grid$coordinate.ranges,
    weights = weights,
    spline.orders = spline.orders,
    max.coef = max.coef
  )

  fit.or.null <- function(name) {
    if (name %in% names(stageB$fits)) stageB$fits[[name]] else NULL
  }
  knots.or.null <- function(name) {
    if (name %in% names(stageB$intknots)) stageB$intknots[[name]] else NULL
  }
  linear.fit <- fit.or.null("linear")
  quadratic.fit <- fit.or.null("quadratic")
  cubic.fit <- fit.or.null("cubic")
  linear.intknots <- knots.or.null("linear")

  out <- list(
    type = "LM - Multivariate Tensor",
    dimensions = colnames(coordinates),
    linear.intknots = linear.intknots,
    quadratic.intknots = knots.or.null("quadratic"),
    cubic.intknots = knots.or.null("cubic"),
    dev.linear = if (is.null(linear.fit)) NULL else linear.fit$rss,
    dev.quadratic = if (is.null(quadratic.fit)) NULL else quadratic.fit$rss,
    dev.cubic = if (is.null(cubic.fit)) NULL else cubic.fit$rss,
    linear.fit = linear.fit,
    quadratic.fit = quadratic.fit,
    cubic.fit = cubic.fit,
    stageA = stageA,
    stageB = stageB,
    grid = grid,
    Nintknots = if (is.null(linear.intknots)) NULL else lengths(linear.intknots),
    iters = length(stageA$history),
    selected.iteration = stageA$selected.iteration,
    args = list(
      coordinates = coordinates,
      response = response,
      weights = weights,
      beta = beta,
      phi = phi,
      q = q,
      min.intknots = min.intknots,
      max.steps = max.steps,
      stoptype = stoptype,
      coordinate.ranges = grid$coordinate.ranges,
      placement.ranges = placement.ranges,
      nint = grid$nint,
      intknots_init = intknots_init,
      spline.orders = spline.orders,
      max.coef = max.coef
    ),
    Call = call
  )
  class(out) <- c("GeDSfitND", "list")
  out
}


GenMultivariateFitter <- function(coordinates, response, family,
                                  weights = rep(1, NROW(coordinates)),
                                  offset = rep(0, NROW(coordinates)),
                                  parametric = NULL,
                                  beta = 0.5, phi = 0.99, q = 2L,
                                  min.intknots = 0L, max.steps = 300L,
                                  stoptype = c("SR", "RD", "LR"),
                                  coordinate.ranges = NULL,
                                  placement.ranges = NULL,
                                  nint = NULL, spline.orders = 2:4,
                                  max.coef = 100000L)
{
  call <- match.call()
  stoptype <- match.arg(stoptype)
  coordinates <- as.matrix(coordinates)
  ndim <- NCOL(coordinates)
  if (is.null(colnames(coordinates))) colnames(coordinates) <- paste0("X", seq_len(ndim))
  max.coef <- as.numeric(max.coef)
  if (2^ndim + NCOL(parametric) > max.coef) {
    stop("The initial linear tensor basis requires ", 2^ndim,
         " coefficients, exceeding 'max.coef' = ", max.coef, ".",
         call. = FALSE)
  }
  grid <- makeGridBoundsND(coordinates, coordinate.ranges, nint)
  if (is.null(placement.ranges)) {
    placement.ranges <- lapply(seq_len(ndim), function(j) {
      range(coordinates[weights > 0, j])
    })
  }
  stageA <- stageAGLMLoopND(
    coordinates, response, grid$upper.bounds, family, offset, parametric,
    intknots = rep(list(NULL), ndim),
    coordinate.ranges = grid$coordinate.ranges,
    placement.ranges = placement.ranges, weights = weights, beta = beta,
    max.steps = max.steps, stop.rule = stoptype, phi = phi, q = q,
    min.intknots = min.intknots, max.coef = max.coef
  )
  stageB <- stageBGLMND(
    coordinates, response, stageA$selected.intknots, family, offset, parametric,
    grid$coordinate.ranges, weights, spline.orders, max.coef
  )
  fit.or.null <- function(name) if (name %in% names(stageB$fits)) stageB$fits[[name]] else NULL
  knots.or.null <- function(name) if (name %in% names(stageB$intknots)) stageB$intknots[[name]] else NULL
  linear.fit <- fit.or.null("linear")
  quadratic.fit <- fit.or.null("quadratic")
  cubic.fit <- fit.or.null("cubic")
  linear.intknots <- knots.or.null("linear")
  out <- list(
    type = "GLM - Multivariate Tensor", dimensions = colnames(coordinates),
    linear.intknots = linear.intknots,
    quadratic.intknots = knots.or.null("quadratic"),
    cubic.intknots = knots.or.null("cubic"),
    dev.linear = if (is.null(linear.fit)) NULL else linear.fit$rss,
    dev.quadratic = if (is.null(quadratic.fit)) NULL else quadratic.fit$rss,
    dev.cubic = if (is.null(cubic.fit)) NULL else cubic.fit$rss,
    linear.fit = linear.fit, quadratic.fit = quadratic.fit, cubic.fit = cubic.fit,
    stageA = stageA, stageB = stageB, grid = grid,
    Nintknots = if (is.null(linear.intknots)) NULL else lengths(linear.intknots),
    iters = length(stageA$history), selected.iteration = stageA$selected.iteration,
    args = list(
      coordinates = coordinates, response = response, weights = weights,
      offset = offset, parametric = parametric, family = family,
      beta = beta, phi = phi, q = q,
      min.intknots = min.intknots, max.steps = max.steps, stoptype = stoptype,
      coordinate.ranges = grid$coordinate.ranges,
      placement.ranges = placement.ranges, nint = grid$nint,
      spline.orders = spline.orders, max.coef = max.coef
    ),
    Call = call
  )
  class(out) <- c("GeDSfitND", "list")
  out
}


#' Predict from an experimental multivariate tensor GeDS fit
#'
#' @param object A fitted object returned by the internal
#'   \code{MultivariateFitter()} function.
#' @param newdata Optional matrix or data frame containing the predictor
#'   coordinates.
#' @param n Spline order: 2 (linear), 3 (quadratic), or 4 (cubic).
#' @param ... Currently unused.
#'
#' @return A numeric vector of fitted values.
#' @export
predict.GeDSfitND <- function(object, newdata = NULL, n = 3L, ...)
{
  if (!inherits(object, "GeDSfitND")) {
    warning("calling predict.GeDSfitND() with a non-GeDSfitND object")
  }
  if (length(list(...))) {
    warning("Additional arguments are currently ignored.")
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n)) {
    stop("'n' must identify one fitted spline order.", call. = FALSE)
  }

  fitted.orders <- vapply(
    object$stageB$fits, `[[`, integer(1), "spline.order"
  )
  fit.index <- match(n, fitted.orders)
  if (is.na(fit.index)) {
    stop("Spline order ", n, " was not fitted.", call. = FALSE)
  }
  stage.fit <- object$stageB$fits[[fit.index]]

  if (is.null(newdata)) {
    coordinates <- object$args$coordinates
    parametric <- object$args$parametric
    prediction.offset <- object$args$offset
    if (is.null(prediction.offset)) prediction.offset <- rep(0, NROW(coordinates))
  } else if (is.atomic(newdata) && is.null(dim(newdata))) {
    if (!is.null(object$args$parametric)) {
      stop("A named data frame is required when the model has parametric terms.", call. = FALSE)
    }
    if (length(newdata) != length(object$dimensions)) {
      stop("A newdata vector must provide one value per fitted dimension.",
           call. = FALSE)
    }
    coordinates <- matrix(newdata, nrow = 1L)
    colnames(coordinates) <- names(newdata)
    parametric <- NULL
    prediction.offset <- rep(0, NROW(coordinates))
  } else {
    coordinates <- as.matrix(newdata)
    parametric <- NULL
    prediction.terms <- NULL
    if (!is.null(object$terms) && inherits(object$terms, "terms")) {
      prediction.terms <- stats::delete.response(object$terms)
      environment(prediction.terms) <- asNamespace("GeDS")
    }
    if (!is.null(object$args$parametric)) {
      mm <- stats::model.matrix(prediction.terms, newdata)
      missing.parametric <- setdiff(object$znames, colnames(mm))
      if (length(missing.parametric)) {
        stop("'newdata' is missing parametric design columns: ",
             paste(missing.parametric, collapse = ", "), ".", call. = FALSE)
      }
      parametric <- mm[, object$znames, drop = FALSE]
    }
    if (!is.null(prediction.terms)) {
      mf <- stats::model.frame(prediction.terms, newdata)
      prediction.offset <- stats::model.offset(mf)
    } else {
      prediction.offset <- NULL
    }
    if (is.null(prediction.offset)) prediction.offset <- rep(0, NROW(coordinates))
  }

  if (!is.null(colnames(coordinates))) {
    missing.dimensions <- setdiff(object$dimensions, colnames(coordinates))
    if (length(missing.dimensions)) {
      stop("'newdata' is missing fitted dimensions: ",
           paste(missing.dimensions, collapse = ", "), ".", call. = FALSE)
    }
    coordinates <- coordinates[, object$dimensions, drop = FALSE]
  } else if (NCOL(coordinates) != length(object$dimensions)) {
    stop("'newdata' must have one column per fitted dimension.",
         call. = FALSE)
  }
  if (!is.numeric(coordinates)) {
    stop("All prediction coordinates must be numeric.", call. = FALSE)
  }

  outside <- vapply(
    seq_along(object$dimensions),
    function(j) {
      boundary <- range(stage.fit$full.knots[[j]])
      any(coordinates[, j] < boundary[1L] |
            coordinates[, j] > boundary[2L])
    },
    logical(1)
  )
  if (any(outside)) {
    warning("Some prediction coordinates are outside the boundary knots.")
  }

  basis.matrices <- lapply(
    seq_along(object$dimensions),
    function(j) {
      splines::splineDesign(
        knots = stage.fit$full.knots[[j]],
        x = coordinates[, j],
        ord = n,
        derivs = rep(0L, NROW(coordinates)),
        outer.ok = TRUE
      )
    }
  )
  design <- tensorProdND(basis.matrices)
  if (!is.null(parametric)) design <- cbind(design, parametric)
  eta <- as.numeric(design %*% stage.fit$coefficients) + prediction.offset
  if (!is.null(stage.fit$family)) stage.fit$family$linkinv(eta) else eta
}


##########################################################
## Bivariate adapter used by the existing fitters       ##
##########################################################
placeKnot <- function(Dim, Dim.intknots, matr, indicator, FixedDim,
                      ordFixedDim, nintFixedDim, upperFixedDim, beta)
{
  if (Dim == "X") {
    Dim.index <- 1L
    by.row <- TRUE
  } else if (Dim == "Y") {
    Dim.index <- 2L
    by.row <- FALSE
  } else {
    stop("'Dim' must be either 'X' or 'Y'.", call. = FALSE)
  }
  if (length(FixedDim) != NROW(matr)) {
    stop("'FixedDim' must have one value per row of 'matr'.", call. = FALSE)
  }
  if (length(upperFixedDim) != nintFixedDim) {
    stop("'upperFixedDim' must have one bound per fixed-dimension strip.",
         call. = FALSE)
  }

  # Preserve the current bivariate ordering and duplicate aggregation exactly.
  matrFixedDim <- matr[ordFixedDim, , drop = FALSE]
  matrFixedDim <- makeNewMatrCPP(matrFixedDim, indicator, by.row)

  # The bivariate fitter is the one-fixed-coordinate case of the general
  # dimension-wise knot-placement adapter.
  placeKnotND(
    target.index = Dim.index,
    intknots = Dim.intknots,
    coordinates = matrFixedDim[, 1:2, drop = FALSE],
    residuals = matrFixedDim[, 3],
    fixed.bounds = list(upperFixedDim),
    beta = beta,
    dim.range = range(matr[, Dim])
  )
}


##########################################################
## Unmodified safe copy of the pre-fix bivariate code   ##
##########################################################
placeKnot_biv_original <- function(Dim, Dim.intknots, matr, indicator,
                                   FixedDim, ordFixedDim, nintFixedDim,
                                   zeroesFixedDim, dcumFixedDim, beta)
{
  if (Dim == "X") {
    Dim.index <- 1
    by.row <- TRUE
  } else if (Dim == "Y") {
    Dim.index <- 2
    by.row <- FALSE
  }

  matrFixedDim <- matr[ordFixedDim, ]
  matrFixedDim <- makeNewMatrCPP(matrFixedDim, indicator, by.row)

  Dim.mean <- Dim.width <- dFixedDim.Dim <- numeric()
  strip <- rep.int(seq_len(nintFixedDim), diff(c(0L, dcumFixedDim)))
  keep_rows <- !zeroesFixedDim[strip]
  idx <- which(keep_rows)
  ord_idx <- idx[order(strip[idx], matrFixedDim[idx, Dim.index], method = "radix")]
  matrFixedDim[idx, ] <- matrFixedDim[ord_idx, , drop = FALSE]

  st <- strip[idx]
  s <- sign(matrFixedDim[idx, 3])
  breaks <- c(TRUE, (st[-1] != st[-length(st)]) | (s[-1] != s[-length(s)]))
  lens <- diff(c(which(breaks), length(s) + 1L))
  dFixedDim.Dim[seq_along(lens)] <- lens

  dcumFixedDim.Dim <- cumsum(dFixedDim.Dim)
  Dim.mean <- Dim.width <- numeric(length(dFixedDim.Dim))
  Dim.mean[1] <- abs(mean(matrFixedDim[1:dcumFixedDim.Dim[1], 3]))
  Dim.width[1] <- diff(range(matrFixedDim[1:dcumFixedDim.Dim[1], Dim.index]))
  for (i in 2:length(dFixedDim.Dim)) {
    Dim.mean[i] <- abs(mean(
      matrFixedDim[(dcumFixedDim.Dim[i - 1] + 1):dcumFixedDim.Dim[i], 3]
    ))
    Dim.width[i] <- diff(range(
      matrFixedDim[(dcumFixedDim.Dim[i - 1] + 1):dcumFixedDim.Dim[i], Dim.index]
    ))
  }
  Dim.mean <- Dim.mean / max(Dim.mean)
  if (max(Dim.width) != 0) Dim.width <- Dim.width / max(Dim.width)
  Dim.weights <- beta * Dim.mean + (1 - beta) * Dim.width

  if (is.null(Dim.intknots)) Dim.intknots <- NA_real_
  xx <- findNewDimKnot(
    dcumFixedDim.Dim,
    Dim.weights,
    sort(c(Dim.intknots, range(matr[, Dim]))),
    matrFixedDim[, Dim.index],
    matrFixedDim[, 3]
  )

  list(
    Dim.newknot = as.numeric(xx$Dim.newknot),
    weightDim = as.numeric(xx$weightDim),
    flagDim = xx$flagDim
  )
}


##########################################################
## Corrected reference bivariate implementation         ##
##########################################################
placeKnot_biv_legacy <- function(Dim, Dim.intknots, matr, indicator,
                                 FixedDim, ordFixedDim, nintFixedDim,
                                 upperFixedDim, beta)
{
  if (Dim == "X") {
    Dim.index <- 1L
    FixedDim.index <- 2L
    by.row <- TRUE
  } else if (Dim == "Y") {
    Dim.index <- 2L
    FixedDim.index <- 1L
    by.row <- FALSE
  } else {
    stop("'Dim' must be either 'X' or 'Y'.", call. = FALSE)
  }
  if (length(FixedDim) != NROW(matr)) {
    stop("'FixedDim' must have one value per row of 'matr'.", call. = FALSE)
  }
  if (length(upperFixedDim) != nintFixedDim) {
    stop("'upperFixedDim' must have one bound per fixed-dimension strip.",
         call. = FALSE)
  }

  matrFixedDim <- matr[ordFixedDim, , drop = FALSE]
  matrFixedDim <- makeNewMatrCPP(matrFixedDim, indicator, by.row)

  # Unlike placeKnot_biv_original(), this produces one strip ID per row.
  strip <- fixedDimCellId(
    matrFixedDim[, FixedDim.index, drop = FALSE],
    list(upperFixedDim)
  )

  Dim.mean <- Dim.width <- dFixedDim.Dim <- numeric()
  idx <- seq_along(strip)
  ord_idx <- idx[order(strip, matrFixedDim[, Dim.index], method = "radix")]
  matrFixedDim[idx, ] <- matrFixedDim[ord_idx, , drop = FALSE]

  st <- strip[ord_idx]
  s <- sign(matrFixedDim[, 3])
  breaks <- c(TRUE, (st[-1] != st[-length(st)]) | (s[-1] != s[-length(s)]))
  lens <- diff(c(which(breaks), length(s) + 1L))
  dFixedDim.Dim[seq_along(lens)] <- lens

  dcumFixedDim.Dim <- cumsum(dFixedDim.Dim)
  Dim.mean <- Dim.width <- numeric(length(dFixedDim.Dim))
  Dim.mean[1] <- abs(mean(matrFixedDim[1:dcumFixedDim.Dim[1], 3]))
  Dim.width[1] <- diff(range(matrFixedDim[1:dcumFixedDim.Dim[1], Dim.index]))
  if (length(dFixedDim.Dim) > 1L) {
    for (i in 2:length(dFixedDim.Dim)) {
      cluster.rows <- (dcumFixedDim.Dim[i - 1] + 1L):dcumFixedDim.Dim[i]
      Dim.mean[i] <- abs(mean(matrFixedDim[cluster.rows, 3]))
      Dim.width[i] <- diff(range(matrFixedDim[cluster.rows, Dim.index]))
    }
  }
  Dim.mean <- Dim.mean / max(Dim.mean)
  if (max(Dim.width) != 0) Dim.width <- Dim.width / max(Dim.width)
  Dim.weights <- beta * Dim.mean + (1 - beta) * Dim.width

  if (is.null(Dim.intknots)) Dim.intknots <- NA_real_
  xx <- findNewDimKnot(
    dcumFixedDim.Dim,
    Dim.weights,
    sort(c(Dim.intknots, range(matr[, Dim]))),
    matrFixedDim[, Dim.index],
    matrFixedDim[, 3]
  )

  list(
    Dim.newknot = as.numeric(xx$Dim.newknot),
    weightDim = as.numeric(xx$weightDim),
    flagDim = xx$flagDim
  )
}

