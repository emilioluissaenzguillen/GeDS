################################################################################
################################################################################
############################## Univariate Fitters ##############################
################################################################################
################################################################################
#' @title Functions Used to Fit GeDS Objects with a Univariate Spline Regression
#' @name UnivariateFitters
#' @aliases Fitters UnivariateFitter GenUnivariateFitter
#' @description
#' These are computing engines called by \code{\link{NGeDS}} and
#' \code{\link{GGeDS}}, needed for the underlying fitting procedures.
#' @param X A numeric vector containing \eqn{N} sample values of the covariate
#' chosen to enter the spline regression component of the predictor model.
#' @param Y A vector of size \eqn{N} containing the observed values of the
#' response variable \eqn{y}.
#' @param Z A design matrix with \eqn{N} rows containing other covariates
#' selected to enter the parametric component of the predictor model (see
#' \code{\link[=formula.GeDS]{formula}}). If no such covariates are selected, it
#' is set to \code{NULL} by default.
#' @param offset A vector of size \eqn{N} that can be used to specify a fixed
#' covariate to be included in the predictor model avoiding the estimation of
#' its corresponding regression coefficient. In case more than one covariate is
#' fixed, the user should sum the corresponding coordinates of the fixed
#' covariates to produce one common \eqn{N}-vector of coordinates. The
#' \code{offset} argument is particularly useful when using
#' \code{GenUnivariateFitter} if the link function used is not the identity.
#' @param weights An optional vector of size \eqn{N} of `prior weights' to be
#' put on the observations in the fitting process in case the user requires
#' weighted GeDS fitting. It is \code{NULL} by default.
#' @param beta Numeric parameter in the interval \eqn{[0,1]} tuning the knot
#' placement in stage A of GeDS. See the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param phi Numeric parameter in the interval \eqn{(0,1)} specifying the
#' threshold for the stopping rule  (model selector) in stage A of GeDS. See
#' also \code{stoptype} and Details in the description of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param min.intknots Optional parameter allowing the user to set a minimum
#' number of internal knots required. By default equal to zero.
#' @param max.intknots Optional parameter allowing the user to set a maximum
#' number of internal knots to be added by the GeDS estimation algorithm. By
#' default uses the per-spline sample-and-support Stage A coefficient budget described
#' for \code{\link{NGeDS}}. An explicit value overrides this safeguard.
#' @param q Numeric parameter which allows to fine-tune the stopping rule of
#' stage A of GeDS, by default equal to 2. See Details in the description of
#' \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param extr Numeric vector of 2 elements representing the left-most and
#' right-most limits of the interval embedding the sample values of \code{X}. By
#' default equal correspondingly to the smallest and largest values of \code{X}.
#' @param show.iters Logical variable indicating whether or not to print
#' information at each step. By default equal to \code{FALSE}.
#' @param tol Numeric value indicating the tolerance to be used in the knot
#' placement steps in stage A. By default equal to \code{1e-12}. See Details below.
#' @param stoptype A character string indicating the type of GeDS stopping rule
#' to be used. It should be either \code{"SR"}, \code{"RD"} or \code{"LR"},
#' partial match allowed. See Details of \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param higher_order A logical that defines whether to compute the higher
#' order fits (quadratic and cubic) after stage A is run. Default is
#' \code{TRUE}.
#' @param intknots_init Vector of initial internal knots from which to start the GeDS
#' Stage A iterations. See Section 3 of Kaishev et al. (2016). Default is \code{NULL}.
#' @param fit_init A list containing fitted values \code{pred}, along with
#' corresponding \code{intknots} and \code{coef}, representing the initial fit from
#' which to begin Stage A GeDS iteration (i.e. departing from step 2).
#' @param only_pred Logical, if \code{TRUE} only predictions are computed.
#' @param family A description of the error distribution and link function to be
#' used in the model. This can be a character string naming a family function
#' (e.g. \code{"gaussian"}), the family function itself (e.g.
#' \code{\link[stats]{gaussian}}) or the result of a call to a family function
#' (e.g. \code{gaussian()}). See \link[stats]{family} for details on family
#' functions.
#'
#' @return A \code{"GeDS"} class object, but without the \code{formula},
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
#' legend("topleft", c("Quadratic", "Cubic"),
#'        col = c("red", "blue"),
#'        lty = c(1,2))
#'
#' @seealso  \code{\link{NGeDS}} and \code{\link{GGeDS}}.
#' @importFrom stats .lm.fit qchisq pchisq
#' @export
#' @rdname UnivariateFitters

UnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0,NROW(Y)),
                             weights = rep(1,length(X)), beta=.5, phi = 0.5,
                             min.intknots = 0, max.intknots = 300, q = 2,
                             extr = range(X), show.iters = FALSE,
                             tol = as.double(1e-12), stoptype = c("SR","RD","LR"),
                             higher_order = TRUE, intknots_init = NULL, fit_init = NULL,
                             only_pred = FALSE)
  {
  # Capture the function call
  save <- match.call()
  if (missing(max.intknots)) {
    max.intknots <- .geds_spline_cap(X, weights, if (is.null(Z)) 0L else NCOL(Z))$max.intknots
  }
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "weights" = weights,
               "beta" = beta, "phi" = phi, "min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "extr" = extr, "tol" = tol)

  path_state <- if (inherits(stoptype, "geds_gd_path")) stoptype else NULL
  restart <- path_state$restart
  if (!is.null(restart)) fit_init <- restart
  # Initialize intknots, rss and phis
  intknots <- if (!is.null(fit_init)) fit_init$intknots else intknots_init
  n_starting_intknots <- length(intknots)
  if (n_starting_intknots > max.intknots) {
    stop("Initial knots exceed 'max.intknots'; supply a larger limit.", call. = FALSE)
  }
  rssnew <- if (is.null(restart)) numeric(n_starting_intknots) else restart$loss
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL

  # Stop type, Indicator, distinctX
  collect_path <- !is.null(path_state) || identical(stoptype, ".gamma_dynamic_path")
  force_refit <- identical(stoptype, ".gamma_dynamic_refit")
  if (!collect_path && !force_refit) stoptype <- match.arg(stoptype)
  if (length(weights) != length(X) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must provide one finite non-negative value per observation.",
         call. = FALSE)
  }
  positive.weight <- weights > 0
  if (!any(positive.weight)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  # NGeDS/GGeDS sort X before calling the fitter; rle() relies on this ordering.
  placement.X <- X[positive.weight]
  placement.weights <- weights[positive.weight]
  xr <- rle(placement.X)
  Indicator <- xr$lengths
  distinctX <- xr$values
  has_repeated_x <- any(Indicator > 1L)

  geometry_tol <- 32 * .Machine$double.eps * max(abs((Y - offset)[positive.weight]))
  site_weight_scale <- if (has_repeated_x) {
    makeNewRes(placement.weights, as.numeric(Indicator))
  } else {
    placement.weights
  }

  if (!collect_path) {
    # The GD calibration path only returns the losses.
    nz <- if (is.null(Z)) 0L else NCOL(Z)
    previous <- matrix(nrow = max.intknots + 1L, ncol = max.intknots + 4L)
    oldcoef <- matrix(nrow = max.intknots + 1L, ncol = max.intknots + 2L + nz)
  }
  plain_boundaries <- collect_path && is.numeric(extr) && length(extr) == 2L &&
    is.null(attributes(extr)) && all(is.finite(extr))
  boundary_tolerance <- sqrt(.Machine$double.eps)

  # GeDS iterations start by j = n_starting_intknots + 1
  init.iter <- if (is.null(intknots)) 1 else  n_starting_intknots + 1

  # GD path mode: run the Stage A loop below in C++ (gd_stage_a_path). It
  # returns NULL when its tridiagonal solver declines a fit, and the R loop
  # then recomputes the whole path.
  last.iter <- min(max.intknots + 1, length(Y) - 2)
  if (collect_path && is.null(Z) && plain_boundaries && !show.iters &&
      init.iter <= last.iter) {
    first_residuals <- if (is.null(fit_init)) NULL else if (is.null(restart)) {
      Y - fit_init$pred
    } else restart$residuals
    fast <- gd_stage_a_path(as.numeric(X), as.numeric(Y), as.numeric(offset),
                            as.numeric(weights), positive.weight, placement.weights,
                            Indicator, has_repeated_x, site_weight_scale, geometry_tol,
                            distinctX, beta, as.numeric(extr), tol, max.intknots,
                            as.numeric(intknots), as.numeric(rssnew), first_residuals,
                            is.null(restart), init.iter, last.iter)
    if (!is.null(fast)) {
      if (!is.null(path_state)) {
        knots <- if (!length(fast$intknots) && is.null(intknots)) NULL else fast$intknots
        path_state$restart <- if (fast$j == max.intknots + 1L) {
          list(intknots = knots, loss = fast$loss, residuals = fast$residuals)
        } else NULL
      }
      return(fast$loss)
    }
  }


  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  for(j in init.iter:min(max.intknots + 1, length(Y) - 2)) {


    if(is.null(fit_init) || j > (n_starting_intknots + 1) ) {
      #############################################################
      ## STEP 1/STEP 8: Find the least squares linear spline fit ##
      #############################################################
      # The GD refit uses the same linear solver as the GD path it reproduces.
      first.deg <- .spline_reg_lm_fast(X, Y, Z, offset, weights, intknots, 2L, extr,
                                      gd_path = collect_path || force_refit)
      # Store knots and coefficients
      if (!collect_path) {
        previous[j, 1:(j+3)] <- sort(c(intknots, rep(extr, 2)))
        oldcoef[j, 1:(j+1+nz)] <- first.deg$theta
      }

      #####################################
      ## STEP 9: Store residuals and rss ##
      ####################################
      res.tmp <- first.deg$residuals
      rssnew <- c(rssnew, first.deg$rss)

    } else {
      # Store knots and coefficients
      if (!collect_path) previous[j, 1:(j+3)] <- sort(c(intknots, rep(extr, 2)))
      # oldcoef[j, 1:(j+1+nz)] <- first.deg$theta # we don't need the coef for anything

      if (is.null(restart)) {
        res.tmp <- Y - fit_init$pred
        rssnew <- c(rssnew, .weighted_rss(res.tmp, weights))
      } else res.tmp <- restart$residuals
    }

    ############################
    ## STEP 10: Stopping Rule ##
    ############################
    if (collect_path || force_refit) {
      prnt <- ""
      if (j == max.intknots + 1L) break
    } else {
      res <- stopping_rule(j = j, q = q, n_starting_intknots = n_starting_intknots,
                           rssnew = rssnew, phis = phis, stoptype = stoptype,
                           intknots = intknots, min.intknots = min.intknots, phi = phi,
                           phis_star = phis_star, oldintc = oldintc, oldslp = oldslp)

      phis    <- res$phis; phis_star <- res$phis_star
      oldintc <- res$oldintc; oldslp <- res$oldslp
      prnt    <- res$prnt

      if (res$should_break) break
    }

    ##################
    ## STEPS 2 to 7 ##
    ##################
    # Weighted residual clusters, cluster weights and the new knot, in C++.
    # R reference: .geds_knot_step_R() in R/cpp_functions_R.R.
    support_order <- 2
    newknot <- geds_knot_step(res.tmp, positive.weight, placement.weights,
                              Indicator, has_repeated_x, site_weight_scale,
                              geometry_tol, distinctX, beta, as.numeric(extr),
                              as.numeric(intknots), tol, support_order)

    if (is.na(newknot)) break
    at_boundary <- if (plain_boundaries && is.finite(newknot)) {
      # Scalar all.equal() uses the target's magnitude as its relative scale.
      scale <- if (abs(newknot) > boundary_tolerance) abs(newknot) else 1
      any(abs(newknot - extr[1:2]) / scale <= boundary_tolerance)
    } else {
      isTRUE(all.equal(newknot, extr[1])) || isTRUE(all.equal(newknot, extr[2]))
    }
    if (at_boundary) break

    intknots <- c(intknots, newknot)

    # Print iteration
    if(show.iters) {
      indent <- rep(" ", nchar(options()$prompt))
      indent <- paste(indent, collapse="")
      if (j > q + n_starting_intknots) {
        toprint <- paste0(indent, "Iteration ",j,": New Knot = ", round(newknot, 3),
                          ", rss = " , round(rssnew[j],3), prnt,"\n")
      } else {
          toprint <- paste0(indent,"Iteration ",j,": New Knot = ", round(newknot, 3),
                            ", rss = " , round(rssnew[j],3), "\n")
      }
      cat(toprint)
    }

  }

  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  if (collect_path) {
    if (!is.null(path_state)) {
      path_state$restart <- if (j == max.intknots + 1L) {
        list(intknots = intknots, loss = rssnew, residuals = res.tmp)
      } else NULL
    }
    return(rssnew)
  }
  if (j == max.intknots + 1) {
    warning(warningCondition(
      paste0("Internal-knot cap reached (max.intknots = ", max.intknots, ")."),
      class = "geds_knot_cap_warning"
    ))
    iter <- j
    } else if (j <= max.intknots) {
      # Eliminate NAs in previous
      previous <- previous[-((j+1):(max.intknots+1)), , drop = FALSE] # eliminate all the NA rows
      previous <- previous[ ,-((j+4):(max.intknots+4)), drop = FALSE] #  k = j - 1 internal knots (i.e. k + 4 = j + 3 knots) in the last iteration
      # Eliminate NAs in oldcoef
      oldcoef  <- oldcoef[-((j+1):(max.intknots+1)), , drop = FALSE]  # eliminate all the NA rows
      oldcoef  <- oldcoef[ , 1:(j+1+nz), drop = FALSE]                # p = n + k = j + 1  B-splines

      # If stage A breaks due to NA, w take the jth fit as the final fit
      # Also, if intknots_init was set, the first length(intknots_init) rows of previous
      # will be empty
      if(force_refit || is.na(newknot) || !is.null(intknots_init) && length(intknots_init) >= (j-q) ) {
        iter <- j
        } else {
          iter <- j - q
        }
    }

  # 1. LINEAR
  if (iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, InterKnots = ll, n = 2, extr = extr,
                        only_pred = only_pred)
    } else {
      ik <- as.numeric(na.omit(previous[iter, -c(1, 2, iter + 2, iter + 3)])) # keep the internal knots in the "iter"th row
      # Stage B.1 (averaging knot location)
      ll <- makenewknots(ik, 2)
      # Stage B.2
      lin <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, InterKnots = ll, n = 2, extr = extr,
                          only_pred = only_pred)
    }
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
    # 2. QUADRATIC
    if (iter < 3) {
      warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      qq <- NULL
      squ <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                          extr = extr, InterKnots = qq, n = 3, only_pred = only_pred)
      } else {
        # Stage B.1 (averaging knot location)
        qq <- makenewknots(ik, 3)
        # Stage B.2
        squ <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                            extr = extr, InterKnots = qq, n = 3, only_pred = only_pred)
        }
    # 3. CUBIC
    if (iter < 4) {
      warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      cc <- NULL
      cub <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                          extr = extr, InterKnots = cc, n = 4, only_pred = only_pred)
      } else {
        # Stage B.1 (averaging knot location)
        cc <- makenewknots(ik, 4)
        # Stage B.2
        cub <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                            extr = extr, InterKnots = cc, n = 4, only_pred = only_pred)
      }
    } else {
      qq <- squ <- cc <- cub <- NULL
    }

  out <- list("type" = "LM - Univ", "linear.intknots" = ll, "quadratic.intknots" = qq, "cubic.intknots" = cc,
              "dev.linear" = lin$rss, "dev.quadratic" = squ$rss, "dev.cubic" = cub$rss,
              "rss" = rssnew, "linear.fit" = lin, "quadratic.fit" = squ, "cubic.fit" = cub, "stored" = previous,
              "args"= args, "call"= save, "Nintknots" = iter - 1, "iters" = j, "guesses" = NULL,
              "coefficients" = oldcoef, stopinfo = list("phis" = phis, "phis_star" = phis_star, "oldintc" = oldintc, "oldslp" = oldslp))

  class(out) <- "GeDS"
  return(out)
}


################################################################################
############################# GenUnivariateFitter ##############################
################################################################################
#' @rdname UnivariateFitters
#' @importFrom Matrix rankMatrix
#' @importFrom stats gaussian Gamma binomial poisson quasibinomial quasipoisson .lm.fit pchisq qchisq
#' @export
GenUnivariateFitter <- function(X, Y, Z = NULL, offset = rep(0, NROW(Y)),
                                weights = rep(1, length(X)), family = gaussian(),
                                beta = 0.5, phi = 0.5, min.intknots = 0,
                                max.intknots = 300, q = 2, extr = range(X),
                                show.iters = F, tol = as.double(1e-12),
                                stoptype = c("SR","RD","LR"), higher_order = TRUE)
{
  # Capture the function call
  save <- match.call()
  if (length(weights) != length(X) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must provide one finite non-negative value per observation.", call. = FALSE)
  }
  positive.weight <- weights > 0
  if (!any(positive.weight)) {
    stop("'weights' must contain at least one positive value.", call. = FALSE)
  }
  if (missing(extr)) extr <- range(X[positive.weight])
  if (missing(max.intknots)) {
    max.intknots <- .geds_spline_cap(X, weights, if (is.null(Z)) 0L else NCOL(Z))$max.intknots
  }
  # Extract arguments
  args <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "weights"=weights,
               "beta" = beta, "phi" = phi, "family"=family, "min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "extr" = extr)

  # Initialize rss and phis
  path_state <- if (inherits(stoptype, "geds_gd_path")) stoptype else NULL
  restart <- path_state$restart
  rssnew <- if (is.null(restart)) numeric() else restart$loss
  phis <- NULL
  # Initialize \hat{\phi}_\kappa, \hat{\gamma}_0 and \hat{\gamma}_\1 (stoptype = "SR"; see eq. 9 in Dimitrova et al. (2023))
  phis_star <- NULL; oldintc <- NULL; oldslp <- NULL

  # Stop type, Indicator, distinctX
  collect_path <- !is.null(path_state) || identical(stoptype, ".gamma_dynamic_path")
  force_refit <- identical(stoptype, ".gamma_dynamic_refit")
  if (!collect_path && !force_refit) stoptype <- match.arg(stoptype)
  # NGeDS/GGeDS sort X before calling the fitter; rle() relies on this ordering.
  placement.X <- X[positive.weight]
  xr <- rle(placement.X)
  Indicator <- xr$lengths
  distinctX <- xr$values
  has_repeated_x <- any(Indicator > 1L)

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

  # Accumulated number of IRLS iterations
  irlsAccumIterCount <- NULL
  # GD path collection keeps the fitted prefix when IRLS fails at a later knot.
  prefix_path <- collect_path && family$family %in% c("inverse.gaussian", "binomial")
  fit_linear <- function() {
    mustart <- if (j == 1L && isTRUE(path_state$inverse_gaussian_bootstrap)) {
      rep(mean(Y), length(Y))
    } else NULL
    fit <- function(start) {
      SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                    InterKnots = intknots, n = 2, extr = extr,
                    family = family, inits = guess, mustart = start, fast = TRUE)
    }
    if (j == 1L && !length(intknots) && family$family == "inverse.gaussian" &&
        (collect_path || force_refit || identical(stoptype, "SR"))) {
      .retry_ig_zero_start(fit, mustart, Y, weights)
    } else {
      fit(mustart)
    }
  }
  # Control basis matrix singularity
  flag <- FALSE
  insertion_failed <- FALSE
  last_iteration <- min(max.intknots + 1, length(Y) - 2)
  # GD path mode continues in compiled code after the first fit (see
  # gd_glm_path() in src/gd_glm_path.cpp); it hands any step it cannot
  # reproduce exactly back to the R loop below.
  glm_codes <- if (collect_path) .gd_glm_fast_codes(family) else NULL
  glm_fast <- !is.null(glm_codes) && is.null(Z) && !show.iters &&
    is.numeric(extr) && length(extr) == 2L && all(is.finite(extr)) &&
    !any(weights == 0 & (X < extr[1L] | X > extr[2L]))
  skip_to <- 0L

  ##############################################################################
  ################################## STAGE A ###################################
  ##############################################################################
  first_iteration <- if (is.null(restart)) 1L else restart$j
  for (j in seq.int(first_iteration, last_iteration)) {
    if (j < skip_to) next
    if (!is.null(restart) && j == first_iteration) {
      intknots <- restart$intknots
      guess <- restart$guess
      guess_z <- restart$guess_z
      res.tmp <- restart$residuals
      working.weights <- restart$weights
    } else {
      if (flag) j <- j - 2 # adjusting iteration count after removing a knot

      # Sort internal knots and update the oldguess matrix; create new guess vector
      if (j > 1) {
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
      first.deg <- if (prefix_path && length(rssnew)) {
        tryCatch(fit_linear(), error = function(e) {
          if (!.geds_irls_failure(e)) stop(e)
          e
        })
      } else fit_linear()
      if (inherits(first.deg, "error")) {
        if (!is.null(path_state)) path_state$restart <- NULL
        bad <- .gd_invalid_path_steps(rssnew)
        return(if (length(bad)) utils::head(rssnew, bad[1L] - 1L) else rssnew)
      }

      # 1. Check for NA values in the theta vector to handle potential singularities
      if(anyNA(first.deg$theta)) {
        design <- cbind(first.deg$basis, Z)
        working.weights <- first.deg$temporary$weights
        informative <- working.weights > 0
        weighted.design <- design[informative, , drop = FALSE] * sqrt(working.weights[informative])
        rank.basis <- qr(weighted.design)$rank
        cols <- NCOL(design)

        # (i) Handle the case when the basis matrix is singular
        if(rank.basis < cols) {
          # Knot deletion cannot repair aliased parametric terms or a deficiency
          # that exists only after weighting an otherwise full-rank spline basis.
          if (any(which(is.na(first.deg$theta)) > NCOL(first.deg$basis)) ||
              rankMatrix(first.deg$basis) == NCOL(first.deg$basis)) {
            stop("The weighted design matrix is rank deficient; check positive-weight ",
                 "spline support and parametric covariates.", call. = FALSE)
          }
          # (a) If basis was singular for second consecutive time, break loop
          if (flag) {
            warning("Matrix singular for the second time. Breaking the loop.")
            break
          }
          # (b) Identify and remove problematic knots and guesses based on NA positions
          check <- which(is.na(first.deg$theta))
          intknots <- intknots[-(check-1)]
          guess <- guess[1:length(first.deg$theta)][-check]
          toprint <- paste0("basis Matrix singular, deleting one knot")
          print(toprint)
          flag  <- TRUE
          # (c) Check if the number of knots equals the number of unique X values and issue a warning if true
          if(cols == length(distinctX)) {
            warning("Number of knots equal to number of unique Xs. Breaking the loop.")
            break
            } else {
              # Continue to the next iteration otherwise
              next
            }

        # (ii) NA values in the theta vector, but basis matrix is not singular (i.e. other issues)
          } else {
            stop("NA(s) in the coefficients")
          }

      # 2. If no NAs, update guess (coefficients initial value in the next iteration)
      } else {
        guess <- first.deg$theta[1:(j+1)]
      }

      # Accumulated number of IRLS iterations at each GeDS iteration
      if (flag) irlsAccumIterCount <- irlsAccumIterCount[1:(j-1)]
      irlsAccumIterCount <- c(irlsAccumIterCount, first.deg$temporary$iter)

      # Store knots and coefficients
      previous[j,1:(j+3)] <- sort(c(intknots,rep(extr,2)))
      oldcoef[j,1:(j+1+nz)] <- first.deg$theta
      guess_z <- if(nz > 0) first.deg$theta[(j+2):(j+1+nz)] else NULL

      # Store residuals and deviance
      res.tmp <- first.deg$residuals
      rss.tmp <- first.deg$temporary$deviance
      if (flag) rssnew <- rssnew[1:(j-1)]
      rssnew <- c(rssnew, rss.tmp)
      # Working weights (weights in the final iteration of the IRLS fit)
      working.weights <- first.deg$temporary$weights
    }

    ############################
    ## STEP 10: Stopping Rule ##
    ############################
    if (collect_path || force_refit) {
      prnt <- ""
      if (j == max.intknots + 1L) break
    } else {
      res <- stopping_rule(j = j, q = q, flag = flag,
                           rssnew = rssnew, phis = phis, stoptype = stoptype,
                           intknots = intknots, min.intknots = min.intknots, phi = phi,
                           phis_star = phis_star, oldintc = oldintc, oldslp = oldslp)

      phis    <- res$phis; phis_star <- res$phis_star
      oldintc <- res$oldintc; oldslp <- res$oldslp
      prnt    <- res$prnt

      if (res$should_break) break
    }

    if (glm_fast && !flag) {
      fast <- gd_glm_path(as.numeric(X), as.numeric(Y), as.numeric(offset),
                          as.numeric(weights), positive.weight, as.integer(Indicator),
                          has_repeated_x, as.numeric(distinctX), beta, as.numeric(extr),
                          tol, glm_codes[1L], glm_codes[2L], as.numeric(intknots),
                          as.numeric(guess), as.numeric(res.tmp),
                          as.numeric(working.weights), j = as.integer(j),
                          max_intknots = max.intknots, last_iter = last_iteration)
      rssnew <- c(rssnew, fast$loss)
      intknots <- fast$intknots
      guess <- fast$guess
      guess_z <- NULL
      res.tmp <- fast$residuals
      working.weights <- fast$weights
      if (fast$status == 3L) {
        # The R loop refits iteration fast$j from these knots and guess.
        skip_to <- fast$j
        next
      }
      j <- fast$j
      if (fast$status == 4L) {
        # The R code below places the next knot after the fit at iteration j.
        skip_to <- j + 1L
      } else {
        if (fast$status == 2L) insertion_failed <- TRUE
        break
      }
    }

    ###############################################
    ## STEPS 3-8 (i.e. steps 2-7 of Normal GeDS) ##
    ###############################################
    if (has_repeated_x) {
      # Average res.tmp*weights for repeated values of X
      res.weighted <- makeNewRes2(resold = res.tmp[positive.weight],
                                 weights = (weights * working.weights)[positive.weight],
                                 recurr = as.numeric(Indicator))
    } else {
      res.weighted <- (res.tmp * weights * working.weights)[positive.weight]
    }
    # Group consecutive residuals into same-sign clusters (run-length encoding)
    d <- rle(sign(res.weighted))$lengths
    ####################################################################
    ## STEP 4: within residual cluster means +  within-cluster ranges ##
    ####################################################################
    u    <- length(d)
    dcum <- cumsum(d)
    cs   <- cumsum(res.weighted)
    means    <- abs((cs[dcum] - c(0, cs[dcum[-u]])) / d)   # |cluster means|
    left_idx <- c(1, dcum[-u] + 1)
    wc.range <- distinctX[dcum] - distinctX[left_idx]      # cluster x-ranges
    ######################################################################
    ## STEP 5: calculate the normalized within-cluster means and ranges ##
    ######################################################################
    means <- means/max(means)
    # If the residual clusters are all singletons then all the wc.ranges will equal 0, and we cannot divide by 0
    if (max(wc.range) != 0) wc.range <- wc.range/max(wc.range)

    ###########################################
    ## STEP 6: calculate the cluster weights ##
    ###########################################
    w <- beta*means + (1 - beta)*wc.range

    ###############
    ## STEPS 7/8 ##
    ###############
    # i. Calculate new knot (R reference for Knotnew(): Knotnew_R_clean() in R/cpp_functions_R.R)
    support_order <- 2
    newknot <- Knotnew(weights = w, residuals = res.weighted, x = distinctX,
                       dcum = dcum, oldknots = c(rep(extr, support_order), intknots),
                       tol = tol, support_order = support_order)[1]

    if (!is.finite(newknot) || newknot <= extr[1L] || newknot >= extr[2L]) {
      insertion_failed <- TRUE
      break
    }

    # ii. Calculate guess-coefficient for newknot
    guess <- newknot.guess(intknots, extr, guess, newknot)

    # iii. Update internal knots vector
    intknots <- c(intknots, newknot)

    # Print iteration
    if (show.iters) {
      indent <- rep(" ", nchar(options()$prompt))
      indent <- paste(indent, collapse="")
      if(j > q) {
        toprint <- paste0(indent, "Iteration ",j,": New knot = ", round(newknot,3),
                          ", dev =  " , round(rssnew[j],3), prnt,"\n")
      } else {
        toprint <- paste0(indent, "Iteration ", j, ": New knot = ", round(newknot,3),
                          ", dev =  " , round(rssnew[j],3), "\n")
      }
      cat(toprint)
    }

  }

  ##############################################################################
  ################################## STAGE B ###################################
  ##############################################################################
  if (collect_path) {
    if (!is.null(path_state)) {
      path_state$restart <- if (j == max.intknots + 1L && !flag) {
        list(j = j, loss = rssnew, intknots = intknots, guess = guess, guess_z = guess_z,
             residuals = res.tmp, weights = working.weights)
      } else NULL
    }
    return(rssnew)
  }
  if (j == max.intknots + 1) {
    warning(warningCondition(
      paste0("Internal-knot cap reached (max.intknots = ", max.intknots, ")."),
      class = "geds_knot_cap_warning"
    ))
    iter <- j
  } else if (j <= max.intknots) {
    # Eliminate NAs in previous
    previous <- previous[-((j+1):(max.intknots+1)), , drop = FALSE] # eliminate all the NA rows
    previous <- previous[ ,-((j+4):(max.intknots+4)), drop = FALSE] #  k = j - 1 internal knots (i.e. k + 4 = j + 3 knots) in the last iteration
    # Eliminate NAs in oldcoef
    oldcoef  <- oldcoef[-((j+1):(max.intknots+1)), , drop = FALSE]  # eliminate all the NA rows
    oldcoef  <- oldcoef[ ,1:(j+1+nz), drop = FALSE]                 # p = n + k = j + 1  B-splines
    # An exhausted insertion search retains the latest successful fit.
    iter <- if (force_refit || insertion_failed) j else j - q
  }

  # If model selected is from first iteration
  if (iter == 1) {
    mustart <- NULL
    linear.inits <- NULL
  } else {
    mustart <- oldguess[iter, 1:(iter+1)]
    selected.guess.z <- if (nz > 0L) {
      oldcoef[iter, (iter + 2):(iter + 1 + nz)]
    } else {
      NULL
    }
    linear.inits <- c(mustart, selected.guess.z)
  }

  # 1. LINEAR
  if (iter < 2) {
    warning("Too few internal knots found: Linear spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
    ll <- NULL
    lin <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = ll, n = 2, family = family,
                         inits = linear.inits)
  } else {
    ik <- as.numeric(na.omit(previous[iter,-c(1,2,iter+2,iter+3)]))
    # Stage B.1 (averaging knot location)
    ll <- makenewknots(ik, 2)
    # Stage B.2
    lin <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                         extr = extr, InterKnots = ll, n = 2, family = family,
                         inits = linear.inits)
  }
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
    # 2. QUADRATIC
    if (iter < 3) {
      warning("Too few internal knots found: Quadratic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      qq <- NULL
      squ <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                           extr = extr, InterKnots = qq, n = 3, family = family,
                           mustart = lin$predicted)
      } else {
        # Stage B.1 (averaging knot location)
        qq <- makenewknots(ik, 3)
        # Stage B.2
        squ <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                             extr = extr, InterKnots = qq, n = 3, family = family,
                             mustart = lin$predicted)
        }
    # 3. CUBIC
    if (iter < 4) {
      warning("Too few internal knots found: Cubic spline will be computed with NULL internal knots. Try to set a different value for 'q' or a different treshold")
      cc <- NULL
      cub <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                           extr = extr, InterKnots = cc, n = 4, family = family,
                           mustart = squ$predicted)
      } else {
        # Stage B.1 (averaging knot location)
        cc <- makenewknots(ik, 4)
        # Stage B.2
        cub <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights,
                             extr = extr, InterKnots = cc, n = 4, family = family,
                             mustart = squ$predicted)
      }
    } else {
      qq <- squ <- cc <- cub <- NULL
    }

  out <- list("type" = "GLM - Univ", "linear.intknots" = ll, "quadratic.intknots" = qq, "cubic.intknots" = cc,
              "dev.linear" = lin$rss, "dev.quadratic" = squ$rss, "dev.cubic" = cub$rss,
              "rss" = rssnew, "linear.fit" = lin, "quadratic.fit" = squ, "cubic.fit" = cub, "stored" = previous,
              "args" = args, "call" = save, "Nintknots" = iter - 1, "iters" = j, "guesses" = oldguess,
              "coefficients" = oldcoef, "iterIrls" = irlsAccumIterCount,
              stopinfo = list("phis" = phis, "phis_star" = phis_star, "oldintc" = oldintc, "oldslp" = oldslp))

  class(out) <- "GeDS"
  return(out)
}

# Family and link codes for the compiled GD path, gd_glm_path(); NULL for any
# family or link it does not implement, or customised family callbacks, which
# then stay with the R loop.
.gd_glm_fast_codes <- function(family) {
  families <- c(gaussian = 0L, poisson = 1L, binomial = 2L, Gamma = 3L,
                inverse.gaussian = 4L)
  links <- c(identity = 0L, log = 1L, logit = 2L, inverse = 3L, "1/mu^2" = 4L)
  if (!is.list(family) || !is.character(family$family) || !is.character(family$link) ||
      !family$family %in% names(families) || !family$link %in% names(links)) {
    return(NULL)
  }
  standard <- get(family$family, envir = asNamespace("stats"))(link = family$link)
  callbacks <- c("variance", "dev.resids", "linkfun", "linkinv", "mu.eta",
                 "validmu", "valideta", "initialize")
  for (key in callbacks) {
    actual <- family[[key]]
    expected <- standard[[key]]
    if (!identical(actual, expected, ignore.environment = TRUE)) return(NULL)
    if (is.function(expected)) {
      # Built-in family closures have fresh environments on each call.
      # Compare the bindings their bodies use, rather than those environments'
      # identities; a copied body with changed captured values is still custom.
      symbols <- setdiff(all.names(body(expected)), names(formals(expected)))
      for (symbol in symbols) {
        actual_has <- exists(symbol, environment(actual), inherits = TRUE)
        expected_has <- exists(symbol, environment(expected), inherits = TRUE)
        if (!identical(actual_has, expected_has)) return(NULL)
        if (actual_has && !identical(get(symbol, environment(actual)),
                                    get(symbol, environment(expected)))) return(NULL)
      }
    }
  }
  c(families[[family$family]], links[[family$link]])
}
