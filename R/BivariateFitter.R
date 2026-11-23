################################################################################
################################################################################
############################## Bivariate Fitters ###############################
################################################################################
################################################################################
fullRankBivariateLinear <- function(X, Y, Xintknots, Yintknots,
                                    Xextr, Yextr, weights) {
  positive <- weights > 0
  basis.size <- (length(Xintknots) + 2L) * (length(Yintknots) + 2L)
  if (basis.size > sum(positive)) return(FALSE)
  Xbasis <- splines::splineDesign(
    sort(c(Xintknots, rep(Xextr, 2L))), X[positive], ord = 2L,
    outer.ok = TRUE
  )
  Ybasis <- splines::splineDesign(
    sort(c(Yintknots, rep(Yextr, 2L))), Y[positive], ord = 2L,
    outer.ok = TRUE
  )
  basis <- tensorProd(Xbasis, Ybasis)
  qr(basis)$rank == NCOL(basis)
}

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
#' @param offset A numeric vector of length \eqn{N} added to the linear
#' predictor with coefficient fixed at one. Defaults to zero.
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
#' algorithm. Defaults to \code{300L} in these low-level fitters.
#' @param max.coef Optional upper bound on the number of linear Stage A
#' tensor-product coefficients plus parametric columns. Defaults to no
#' coefficient cap. An explicit limit does not cap a complete additive model
#' or its higher-order refits.
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
#' to be used. Choose \code{"SR"}, \code{"RD"}, \code{"LR"} or
#' \code{"SR-post"}. For SR-post, the \eqn{q}th root of the smoothed ratio is
#' compared with \code{phi}. Use the full name \code{"SR-post"}. See details
#' in \code{\link{NGeDS}} or \code{\link{GGeDS}}.
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
                            max.intknots = 300L, q = 2L, Xextr = range(X[weights > 0]),
                            Yextr = range(Y[weights > 0]), show.iters = TRUE,
                            tol = as.double(1e-12), stoptype = c("SR","RD","LR"),
                            higher_order = TRUE, Xintknots = NULL, Yintknots = NULL,
                            offset = rep(0, length(X)), max.coef = NULL)
  {
  # Capture the function call
  save <- match.call()
  nparam <- if (is.null(W)) 0L else NCOL(W)
  if (is.null(max.coef)) max.coef <- Inf
  if ((length(Xintknots) + 2L) * (length(Yintknots) + 2L) + nparam > max.coef) {
    stop("Initial spline exceeds 'max.coef'; supply a larger limit.", call. = FALSE)
  }
  coefficient_cap_reached <- FALSE
  if (length(offset) != length(X) || anyNA(offset) || any(!is.finite(offset))) {
    stop("'offset' must contain one finite value per observation.", call. = FALSE)
  }
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol, "offset" = offset,
               "max.coef" = max.coef)

  # Initialize rss and phis
  n_starting_intknots <- length(Xintknots) + length(Yintknots)
  rssnew <- numeric(n_starting_intknots)
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  # Stop type
  collect_path <- identical(stoptype, ".gamma_dynamic_path")
  force_refit <- identical(stoptype, ".gamma_dynamic_refit")
  phi_sr <- phi
  if (identical(stoptype, "SR-post")) {
    # Comparing the qth root of the smoothed ratio with phi is equivalent.
    phi_sr <- phi^q
    stoptype <- "SR"
  }
  if (!collect_path && !force_refit) stoptype <- match.arg(stoptype)

  # Initialize knots matrix
  previousX <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  previousY <- matrix(nrow = max.intknots + 1, ncol = max.intknots + 4)
  # Initialize coefficients matrix
  nw <- if(!is.null(W)) NCOL(W) else 0
  oldcoef   <- matrix(nrow = max.intknots + 1,
                      ncol = round((max.intknots/2 + 2)^2) + nw) # max number of coef; comes from maximizing f(x) = (x + 2)(max.intknots - x + 2)

  # Knot placement uses only positive-weight rows, as in the univariate fitter:
  # zero-weight rows must not split residual clusters or change the strips.
  place <- weights > 0
  Xplace <- X[place]; Yplace <- Y[place]; n.place <- length(Xplace)
  if (!all(place)) {
    indicator <- if (anyDuplicated(cbind(Xplace, Yplace))) table(Xplace, Yplace) else NULL
  }
  ordX <- order(Xplace, Yplace); ordY <- order(Yplace, Xplace)

  ##################################################################################
  ## STEP 1: Divide the sample space D into M_1/M_2 rectangular strips in X_1/X_2 ##
  ##################################################################################

  # Set the number of intervals for dividing the X_1 and X_2 dimensions (M_1 and M_2)
  nintX <- nintY <- as.integer(sqrt(n.place))

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
    if (collect_path && !fullRankBivariateLinear(
      X, Y, Xintknots, Yintknots, Xextr, Yextr, weights
    )) break
    first.deg <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                               InterKnotsX = Xintknots, InterKnotsY = Yintknots,
                               Xextr = Xextr, Yextr = Yextr, n = 2, fast = TRUE)

    # Store knots and coefficients
    previousX[j, 1:(length(Xintknots)+4)] <- sort(c(Xintknots, rep(Xextr,2)))
    previousY[j, 1:(length(Yintknots)+4)] <- sort(c(Yintknots, rep(Yextr,2)))
    lth <- length(first.deg$theta); ncoef <- c(ncoef, lth)
    oldcoef[j, 1:lth] <- first.deg$theta

    # Store weighted residuals
    matr <- cbind(X, Y, first.deg$residuals*weights)[place, , drop = FALSE]
    # Store rss
    rss.tmp <- first.deg$rss
    rssnew <- c(rssnew, rss.tmp)

    ###########################
    ## STEP 3: Stopping rule ##
    ###########################
    if (!collect_path && !force_refit && j > q + n_starting_intknots) {

      if (rssnew[j]/rssnew[j-q] > 1) break

      # Adding the current ratio of deviances to the 'phis' vector
      if (stoptype == "LR") {
        phis <- c(phis, rssnew[j-q]-rssnew[j])
        } else {
          phnew <- (rssnew[j]/rssnew[j-q])^(1/(ncoef[j]-ncoef[j-q]))
          if (stoptype == "SR" && (!is.finite(phnew) || phnew >= 1)) break
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
          if(phi_kappa >= phi_sr)  break

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

    if ((collect_path || force_refit) && j == max.intknots + 1L) break

    ###################################
    ## STEP 4. (i) X knot placement ##
    ##################################
    placeXKnot <- placeKnot(
      Dim = "X", Dim.intknots = Xintknots, matr = matr,
      indicator = indicator, FixedDim = Yplace, ordFixedDim = ordY,
      nintFixedDim = nintY, upperFixedDim = upperY, beta = beta
    )

    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim

    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(
      Dim = "Y", Dim.intknots = Yintknots, matr = matr,
      indicator = indicator, FixedDim = Xplace, ordFixedDim = ordX,
      nintFixedDim = nintX, upperFixedDim = upperX, beta = beta
    )

    Ynewknot = placeYKnot$Dim.newknot; weightY = placeYKnot$weightDim; flagY = placeYKnot$flagDim

    cap_flags <- .geds_bivariate_cap_flags(flagX, flagY, length(Xintknots),
                                          length(Yintknots), nw, max.coef)
    flagX <- cap_flags$flagX; flagY <- cap_flags$flagY
    coefficient_cap_reached <- cap_flags$capped && flagX && flagY
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
    if (show.iters && !collect_path && !force_refit) {
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
    if((length(Yintknots)+3)*(length(Xintknots)+3)>=n.place) {
      warning("Exiting stage A: Too many knots found")
      break # Exit the loop to avoid adding too many knots (prevent overfitting)
      }
  }

  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################

  if (collect_path) return(utils::tail(rssnew, length(rssnew) - n_starting_intknots))

  # Keep the non-NA columns from the "j"th row
  toBeSaved <- sum(!is.na(previousX[j,]))
  previousX <- previousX[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]
  toBeSaved <- sum(!is.na(previousY[j,]))
  previousY <- previousY[ ,-((toBeSaved + 1):max(max.intknots + 4, toBeSaved + 1)), drop = FALSE]

  # Retain both the tensor-spline and parametric coefficient history.
  oldcoef <- oldcoef[, seq_len((NCOL(previousX) - 2) * (NCOL(previousY) - 2) + nw), drop = FALSE]

  if (j == max.intknots + 1) {
    if (!force_refit) warning("Maximum number of iterations exceeded")
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    iter <- j
    } else {
      # Delete from the "j+1th" row until the "max.intknots+1th" row (i.e. keep the j first rows)
      previousX <- previousX[-((j+1):(max.intknots+1)), , drop = FALSE]
      previousY <- previousY[-((j+1):(max.intknots+1)), , drop = FALSE]
      oldcoef   <- oldcoef[-((j+1):(max.intknots+1)), , drop = FALSE]

      iter <- if (coefficient_cap_reached) j else j - q
      lastXknots <- sum(!is.na(previousX[iter, ]))
      lastYknots <- sum(!is.na(previousY[iter, ]))
    }

  selected.iter <- iter
  while (iter > 1L) {
    candidateX <- previousX[iter, !is.na(previousX[iter, ])]
    candidateY <- previousY[iter, !is.na(previousY[iter, ])]
    internalX <- if (length(candidateX) > 4L) {
      candidateX[3L:(length(candidateX) - 2L)]
    } else NULL
    internalY <- if (length(candidateY) > 4L) {
      candidateY[3L:(length(candidateY) - 2L)]
    } else NULL
    if (fullRankBivariateLinear(
      X, Y, internalX, internalY, Xextr, Yextr, weights
    )) break
    iter <- iter - 1L
  }
  if (iter < selected.iter) {
    warning("Selected linear tensor basis is rank deficient; using the latest full-rank Stage A fit instead.")
    lastXknots <- sum(!is.na(previousX[iter, ]))
    lastYknots <- sum(!is.na(previousY[iter, ]))
  }

  # 1. LINEAR
  if (iter < 2) {
    if (!force_refit) warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    llX <- llY <- NULL
    lin <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                         InterKnotsX = llX, InterKnotsY = llY,
                         Xextr = Xextr, Yextr = Yextr, n = 2)
    } else {
      ikX <- if (lastXknots > 4) previousX[iter, 3:(lastXknots - 2)] else NULL
      ikY <- if (lastYknots > 4) previousY[iter, 3:(lastYknots - 2)] else NULL
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                           InterKnotsX = llX, InterKnotsY = llY,
                           Xextr = Xextr, Yextr = Yextr, n = 2)
    }
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
  # 2. QUADRATIC
  if (iter < 3) {
    if (!force_refit) warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    qqX <- qqY <- NULL
    squ <- SplineReg_biv(
                         X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                         InterKnotsX = qqX, InterKnotsY = qqY,
                         Xextr = Xextr, Yextr = Yextr, n = 3)
    } else {
      # Stage B.1 (averaging knot location)
      qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
      qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
      # Stage B.2
      squ <- SplineReg_biv(
                           X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                           InterKnotsX = qqX, InterKnotsY = qqY,
                           Xextr = Xextr, Yextr = Yextr, n = 3)
      }
  # 3. CUBIC
  if(iter < 4) {
    if (!force_refit) warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ccX <- ccY <- NULL
    cub <- SplineReg_biv(
                         X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                         InterKnotsX = ccX, InterKnotsY = ccY,
                         Xextr = Xextr, Yextr = Yextr, n = 4)
    } else {
      # Stage B.1 (averaging knot location)
      ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
      ccY <- if (length(ikY) < 3) NULL else makenewknots(ikY, 4)
      # Stage B.2
      cub <- SplineReg_biv(
                           X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
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
                               Xextr = range(X[weights > 0]), Yextr = range(Y[weights > 0]),
                               show.iters=TRUE, tol = as.double(1e-12),
                               stoptype = c("SR","RD","LR"), higher_order = TRUE,
                               offset = rep(0, length(X)), max.coef = NULL)
{
  # Capture the function call
  save <- match.call()
  # Extract arguments
  if (length(offset) != length(X) || anyNA(offset) || any(!is.finite(offset))) {
    stop("'offset' must contain one finite value per observation.", call. = FALSE)
  }
  nparam <- if (is.null(W)) 0L else NCOL(W)
  if (is.null(max.coef)) max.coef <- Inf
  if (4L + nparam > max.coef) stop("'max.coef' is smaller than the initial spline.", call. = FALSE)
  coefficient_cap_reached <- FALSE
  args <- list("X" = X, "Y" = Y, "Z" = Z, "W" = W, "weights" = weights, "beta" = beta,
               "phi" = phi, "min.intknots" = min.intknots, "max.intknots" = max.intknots,
               "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol, family = family,
               "offset" = offset, "max.coef" = max.coef)

  # Initialize rss and phis
  rssnew <- numeric()
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL
  # Stop type
  phi_sr <- phi
  if (identical(stoptype, "SR-post")) {
    phi_sr <- phi^q
    stoptype <- "SR"
  }
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
  # Knot placement uses only positive-weight rows, as in the univariate fitter:
  # zero-weight rows must not split residual clusters or change the strips.
  place <- weights > 0
  Xplace <- X[place]; Yplace <- Y[place]; n.place <- length(Xplace)
  if (!all(place)) {
    indicator <- if (anyDuplicated(cbind(Xplace, Yplace))) table(Xplace, Yplace) else NULL
  }
  ordX <- order(Xplace, Yplace); ordY <- order(Yplace, Xplace)

  ##################################################################################
  ## STEP 1: Divide the sample space D into M_1/M_2 rectangular strips in X_1/X_2 ##
  ##################################################################################

  # Set the number of intervals for dividing the X_1 and X_2 dimensions (M_1 and M_2)
  nintX <- nintY <- as.integer(sqrt(n.place))

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
    first.deg <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
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
    matr <- cbind(X, Y, first.deg$residuals*working.weights*weights)[place, , drop = FALSE]

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
          if(phi_kappa >= phi_sr)  break
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
      indicator = indicator, FixedDim = Yplace, ordFixedDim = ordY,
      nintFixedDim = nintY, upperFixedDim = upperY, beta = beta
    )

    Xnewknot = placeXKnot$Dim.newknot; weightX = placeXKnot$weightDim; flagX = placeXKnot$flagDim

    ###################################
    ## STEP 4. (ii) Y knot placement ##
    ###################################
    placeYKnot <- placeKnot(
      Dim = "Y", Dim.intknots = Yintknots, matr = matr,
      indicator = indicator, FixedDim = Xplace, ordFixedDim = ordX,
      nintFixedDim = nintX, upperFixedDim = upperX, beta = beta
    )

    Ynewknot = placeYKnot$Dim.newknot; weightY = placeYKnot$weightDim; flagY = placeYKnot$flagDim


    cap_flags <- .geds_bivariate_cap_flags(flagX, flagY, length(Xintknots),
                                          length(Yintknots), nw, max.coef)
    flagX <- cap_flags$flagX; flagY <- cap_flags$flagY
    coefficient_cap_reached <- cap_flags$capped && flagX && flagY
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
    if((length(Yintknots)+3)*(length(Xintknots)+3)>=n.place) {
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

  # Retain both the tensor-spline and parametric coefficient history.
  oldcoef <- oldcoef[, seq_len((NCOL(previousX) - 2) * (NCOL(previousY) - 2) + nw), drop = FALSE]

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

    iter <- if (coefficient_cap_reached) j else j - q
    lastXknots <- sum(!is.na(previousX[iter, ]))
    lastYknots <- sum(!is.na(previousY[iter, ]))
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
    lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                             InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
                             n = 2, family = family, mustart = mustart)
    } else {
      ikX <- if (lastXknots > 4) previousX[iter, 3:(lastXknots-2)] else NULL
      ikY <- if (lastYknots > 4) previousY[iter, 3:(lastYknots-2)] else NULL
      # Stage B.1 (averaging knot location)
      llX <- if (length(ikX) < 1) NULL else makenewknots(ikX, 2)
      llY <- if (length(ikY) < 1) NULL else makenewknots(ikY, 2)
      # Stage B.2
      lin <- SplineReg_biv_GLM(X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                               InterKnotsX = llX, InterKnotsY = llY, Xextr = Xextr, Yextr = Yextr,
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
      squ <- SplineReg_biv_GLM(
                               X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                               InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                               n = 3, family = family, mustart = guess_lin)
      } else {
        # Stage B.1 (averaging knot location)
        qqX <- if (length(ikX) < 2) NULL else makenewknots(ikX, 3)
        qqY <- if (length(ikY) < 2) NULL else makenewknots(ikY, 3)
        # Stage B.2
        guess_lin <- lin$predicted
        squ <- SplineReg_biv_GLM(
                                 X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                                 InterKnotsX = qqX, InterKnotsY = qqY, Xextr = Xextr, Yextr = Yextr,
                                 n = 3, family = family, mustart = guess_lin)
        }
    # 3. CUBIC
    if (iter < 4) {
      warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      ccX <- ccY <- NULL
      guess_sq <- squ$predicted
      cub <- SplineReg_biv_GLM(
                               X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                               InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
                               n = 4, family = family, mustart = guess_sq)
      } else {
        # Stage B.1 (averaging knot location)
        ccX <- if (length(ikX) < 3) NULL else makenewknots(ikX, 4)
        ccY <- if (length(ikY) < 3) NULL else makenewknots(ikY, 4)
        # Stage B.2
        guess_sq <- squ$predicted
        cub <- SplineReg_biv_GLM(
                                 X = X, Y = Y, Z = Z, W = W, weights = weights, offset = offset,
                                 InterKnotsX = ccX, InterKnotsY = ccY, Xextr = Xextr, Yextr = Yextr,
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
## Shared tensor-product basis helper                  ##
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

  # This produces one strip ID per row.
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

