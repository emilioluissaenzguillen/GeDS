#' @title NGeDSgam: Local Scoring Algorithm with GeD Splines in Backfitting
#' @name NGeDSgam
#' @description
#' Implements the Local Scoring Algorithm (Hastie and Tibshirani
#' (1986)), applying normal linear GeD splines (i.e., \code{\link{NGeDS}}
#' function) to fit the targets within each backfitting iteration. Higher order 
#' fits are computed by pursuing stage B of GeDS after the local-scoring algorithm
#' is run.
#' @param formula a description of the model structure to be fitted,
#' specifying both the dependent and independent variables. Unlike \code{\link{NGeDS}}
#' and \code{\link{GGeDS}}, this formula supports multiple additive (normal) GeD
#' spline regression components as well as linear components. For example, setting
#' \code{formula = Y ~ f(X1) + f(X2) + X3} implies using a normal linear GeD
#' spline as the smoother for \code{X1} and for \code{X2}, while for \code{X3} a
#' linear model would be used.
#' @param family a character string indicating the response variable distribution
#' and link function to be used. Default is \code{"gaussian"}. This should be a
#' character or a family object.
#' @param data a data frame containing the variables referenced in the formula.
#' @param weights an optional vector of `prior weights' to be put on the
#' observations during the fitting process. It should be \code{NULL} or a numeric
#' vector of the same length as the response variable defined in the formula.
#' @param offset a vector of size \eqn{N} that can be used to specify a fixed
#' component to be included in the linear predictor during fitting. In case
#' more than one covariate is fixed, the user should sum the corresponding
#' coordinates of the fixed covariates to produce one common \eqn{N}-vector of
#' coordinates.
#' @param normalize_data a logical that defines whether the data should be
#' normalized (standardized) before fitting the baseline linear model, i.e.,
#' before running the local-scoring algorithm. Normalizing the data involves
#' scaling the predictor variables to have a mean of 0 and a standard deviation
#' of 1. This process alters the scale and interpretation of the knots and
#' coefficients estimated. Default is equal to \code{FALSE}.
#' @param min_iterations optional parameter to manually set a minimum number of
#' boosting iterations to be run. If not specified, it defaults to 0L.
#' @param max_iterations optional parameter to manually set the maximum number
#' of boosting iterations to be run. If not specified, it defaults to 100L.
#' This setting serves as a fallback when the stopping rule, based on
#' consecutive deviances and tuned by \code{phi_gam_exit} and \code{q_gam},
#' does not trigger an earlier termination (see Dimitrova et al. (2025)).
#' Therefore, users can increase/decrease the number of boosting iterations,
#' by increasing/decreasing the value \code{phi_gam_exit} and/or \code{q_gam},
#' or directly specify \code{max_iterations}.
#' @param phi_gam_exit Convergence threshold for local-scoring and backfitting.
#' Both algorithms stop when the relative change in the deviance is below this
#' threshold. Default is \code{0.99}.
#' @param q_gam numeric parameter which allows to fine-tune the stopping rule of
#' the local-scoring and backfitting iterations. By default equal to \code{2L}.
#' @param beta numeric parameter in the interval \eqn{[0,1]}
#' tuning the knot placement in stage A of GeDS, for each of the GeD spline
#' components of the model. Default is equal to \code{0.5}.
#' See details in \code{\link{NGeDS}}.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the
#' threshold for the stopping rule  (model selector) in stage A of GeDS, for each
#' of the GeD spline components of the model. Default is equal to \code{0.99}.
#' See details in \code{\link{NGeDS}}.
#' @param internal_knots The maximum number of internal knots that can be added
#' by the GeDS base-learners in each boosting iteration, effectively setting the
#' value of \code{max.intknots} in \code{\link{NGeDS}} at each backfitting
#' iteration. Default is \code{500L}.
#' @param q numeric parameter which allows to fine-tune the stopping rule of
#' stage A of GeDS, for each of the GeD spline components of the model. By
#' default equal to \code{2L}. See details in \code{\link{NGeDS}}.
#' @param higher_order a logical that defines whether to compute the higher order
#' fits (quadratic and cubic) after the local-scoring algorithm is run. Default
#' is \code{TRUE}.
#'
#' @return \code{\link{GeDSgam-Class}} object, i.e. a list of items that
#' summarizes the main details of the fitted GAM-GeDS model. See
#' \code{\link{GeDSgam-Class}} for details. Some S3 methods are available in
#' order to make these objects tractable, such as
#' \code{\link[=coef.GeDSgam]{coef}}, \code{\link[=knots.GeDSgam]{knots}},
#' \code{\link[=print.GeDSgam]{print}} and \code{\link[=predict.GeDSgam]{predict}}.
#'   
#' @details
#' The  \code{NGeDSgam} function employs the local scoring algorithm to fit a
#' Generalized Additive Model (GAM). This algorithm iteratively fits weighted
#' additive models by backfitting. Normal linear GeD splines, as well as linear
#' learners, are supported as function smoothers within the backfitting
#' algorithm. The local-scoring algorithm ultimately produces a linear fit.
#' Higher order fits (quadratic and cubic) are then computed by calculating the
#' Schoenberg’s variation diminishing spline (VDS) approximation of the linear
#' fit.
#' 
#' On the one hand, \code{NGeDSgam} includes all the parameters of
#' \code{\link{NGeDS}}, which in this case tune the function smoother fit at each
#' backfitting iteration. On the other hand, \code{NGeDSgam} includes some
#' additional parameters proper to the local-scoring procedure. We describe
#' the main ones as follows. 
#' 
#' The \code{family} chosen determines the link function, adjusted dependent
#' variable and weights to be used in the local-scoring algorithm. The number of
#' local-scoring and backfitting iterations is controlled by a
#' \emph{Ratio of Deviances} stopping rule similar to the one presented for
#' \code{\link{GGeDS}}. In the same way \code{phi} and \code{q} tune the stopping
#' rule of \code{\link{GGeDS}}, \code{phi_boost_exit} and \code{q_boost} tune the
#' stopping rule of \code{NGeDSgam}. The user can also manually control the number
#' of local-scoring iterations through \code{min_iterations} and 
#' \code{max_iterations}.
#' 
#' @examples
#' 
#' # Load package
#' library(GeDS) 
#' 
#' data(airquality) 
#' data = na.omit(airquality)
#' data$Ozone <- data$Ozone^(1/3)
#' 
#' formula = Ozone ~ f(Solar.R) + f(Wind, Temp)
#' Gmodgam <- NGeDSgam(formula = formula, data = data,
#' phi = 0.8)
#' MSE_Gmodgam_linear <- mean((data$Ozone - Gmodgam$predictions$pred_linear)^2)
#' MSE_Gmodgam_quadratic <- mean((data$Ozone - Gmodgam$predictions$pred_quadratic)^2)
#' MSE_Gmodgam_cubic <- mean((data$Ozone - Gmodgam$predictions$pred_cubic)^2)
#' 
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#' "Linear NGeDSgam:", MSE_Gmodgam_linear, "\n",
#' "Quadratic NGeDSgam:", MSE_Gmodgam_quadratic, "\n",
#' "Cubic NGeDSgam:", MSE_Gmodgam_cubic, "\n")
#' 
#' ## S3 methods for class 'GeDSboost'
#' # Print 
#' print(Gmodgam)
#' # Knots
#' knots(Gmodgam, n = 2L)
#' knots(Gmodgam, n = 3L)
#' knots(Gmodgam, n = 4L)
#' # Coefficients
#' coef(Gmodgam, n = 2L)
#' coef(Gmodgam, n = 3L)
#' coef(Gmodgam, n = 4L)
#' # Deviances
#' deviance(Gmodgam, n = 2L)
#' deviance(Gmodgam, n = 3L)
#' deviance(Gmodgam, n = 4L)
#'
#' @seealso \code{\link{NGeDS}}; \code{\link{GGeDS}}; \code{\link{GeDSgam-Class}};
#' S3 methods such as \code{\link{knots.GeDSgam}}; \code{\link{coef.GeDSgam}};
#' \code{\link{deviance.GeDSgam}}; \code{\link{predict.GeDSgam}}
#' 
#' @export
#' @importFrom stats setNames
#' 
#' @references 
#' Hastie, T. and Tibshirani, R. (1986). Generalized Additive Models.
#' \emph{Statistical Science} \strong{1 (3)} 297 - 310. \cr
#' DOI: \doi{10.1214/ss/1177013604}
#' 
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
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
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. L. (2025).
#' \pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
#' Models and Functional Gradient Boosting, based on Geometrically Designed
#' (GeD) Splines. \emph{Manuscript submitted for publication.}

################################################################################
################################################################################
################# Local Scoring Algorithm for Additive Models ##################
################################################################################
################################################################################

#####################
### Local Scoring ###
#####################
NGeDSgam <- function(formula, family = "gaussian", data, weights = NULL, offset = NULL,
                     normalize_data = FALSE, min_iterations, max_iterations,
                     phi_gam_exit = 0.99, q_gam = 2,
                     beta = 0.5, phi = 0.99, internal_knots = 500, q = 2,
                     higher_order = TRUE)
{
  
  # Capture the function call
  extcall <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  
  # Convert integer variables to numeric
  data <- data.frame(lapply(data, function(x) if(is.integer(x)) as.numeric(x) else x))
  
  # Formula
  read.formula <- read.formula.gam(formula, data)
  terms <-  read.formula$terms
  response <- read.formula$response
  predictors <- read.formula$predictors
  base_learners <- read.formula$base_learners
  
  # Weights and offset
  nobs = length(data[[response]])
  if (is.null(weights)) weights <- rep.int(1, nobs)
  else weights <- rescale_weights(weights)
  offset_index <- match("offset", names(extcall), 0L)
  offset <- if (offset_index != 0) eval(extcall[[offset_index]], data) else NULL
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  
  # Eliminate indexes and keep only relevant variables
  rownames(data) <- NULL
  all_vars <- c(response, predictors)
  if(all(all_vars %in% names(data))) {
    data <- data[, all_vars]
  } else {
    stop("Error: Some of the variables are not found in the data.")
  }
  # Check for NA in each column and print a warning message if found
  na_vars <- apply(data, 2, function(x) any(is.na(x)))
  for(var in names(data)[na_vars]){
    warning(paste("variable", var, "contains missing values.\n",
                  "These observations will be excluded from the fitting process."))
  }
  # Eliminate rows containing NA values
  data <- na.omit(data)
  # Check if there are any character variables in the data
  if(any(sapply(data, is.character))) {
    # Convert all character variables to factors
    data[] <- lapply(data, function(x) if(is.character(x)) as.factor(x) else x)
    warning("Character variables are converted into factors.")
  }
  
  # Family arguments
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  
  # Ensure the response variable is a factor when using a binomial family
  if (family$family == "binomial" && !is.factor(data[[response]])) {
    data[[response]] <- as.factor(data[[response]])
    if (nlevels(data[[response]]) != 2) 
      stop("response is not a factor at two levels but ", sQuote("family = binomial"))
  }
  
  variance <- family$variance
  dev.resids <- family$dev.resids
  link <- family$linkfun # g
  linkinv <- family$linkinv # g^{-1}
  mu.eta <- family$mu.eta   # dg^{-1}/d\eta
  
  # Min/max iterations
  min_iterations <- validate_iterations(min_iterations, 0L, "min_iterations")
  max_iterations <- validate_iterations(max_iterations, 100L, "max_iterations")
  
  # Save arguments
  args <- list(
    "predictors" = data[predictors], 
    "base_learners"= base_learners,
    "family" = family,
    "normalize_data" = normalize_data,
    "offset" = offset
    )
  
  # Add 'response' as the first element of the list 'args'
  if (args$family$family != "binomial") {
    args <- c(list(response = data[response]), args)
  } else {
    args <- c(list(response = data.frame(as.numeric(data[[response]])-1)), args) # encode to 0/1
    names(args$response) <- response
  }
  
  # Normalize data if necessary
  if (normalize_data == TRUE) {
    if (args$family$family != "binomial") {
      # Mean and SD of the original response and predictor(s) variables (to de-normalize afterwards)
      args$Y_mean <- mean(data[[response]])
      args$Y_sd <- sd(data[[response]])
      data[response] <- as.vector(scale(data[response]))
      
      numeric_predictors <- names(data[predictors])[sapply(data[predictors], is.numeric)]
      args$X_mean <- colMeans(data[numeric_predictors])
      args$X_sd <- sapply(data[numeric_predictors], sd)
      # Scale only numeric predictors
      data[numeric_predictors] <- scale(data[numeric_predictors])
      
    } else {
      # If the family is "binomial" only normalize predictors
      args$X_mean <- colMeans(data[predictors])
      args$X_sd <- sapply(data[predictors], sd)
      data[predictors] <- scale(data[predictors])
    }
  }
  
  # Boundary knots
  args$extr <- lapply(data[ , sapply(data, is.numeric)], range)
 
  # Initialize the knots, intervals, and coefficients for each base-learner
  base_learners_list <- list()
  
  for (bl_name in names(base_learners)) {
    # Extract predictor variable(s) of base-learner
    pred_vars <- base_learners[[bl_name]]$variables
    ## (A) GeDS base-learners
    if (base_learners[[bl_name]]$type=="GeDS") {
      # Get min and max values for each predictor
      min_vals <- sapply(pred_vars, function(p) min(data[[p]]))
      max_vals <- sapply(pred_vars, function(p) max(data[[p]]))
      # (A.1) UNIVARIATE BASE-LEARNERS
      if (length(pred_vars) == 1) {
        knots <- c(min_vals[[1]], max_vals[[1]])
        coefficients <- list("b0" = 0, "b1" = 0)
        base_learners_list[[bl_name]] <- list("knots" = knots, "coefficients" = coefficients)
        # (A.2) BIVARIATE BASE-LEARNERS
      } else {
        knots <- list(Xk = c(min_vals[1], max_vals[1]), Yk = c(min_vals[2], max_vals[2]))
        base_learners_list[[bl_name]] <- list("knots" = knots, "coefficients" = list())
      }
      ## (B) Linear base-learners (initialize only coefficients)
    } else if (base_learners[[bl_name]]$type=="linear"){
      # Factor
      if (is.factor(data[[pred_vars]])){
        n <- nlevels(data[[pred_vars]]); coefficients <- vector("list", length = n)
        names(coefficients) <- paste0("b", 0:(n-1)); coefficients[] <- 0
        # Continuous
      } else {
        coefficients <- list("b0" = 0, "b1" = 0)
      }
      base_learners_list[[bl_name]] <- list("coefficients" = coefficients)
    }
  }
  
  # 1. INITIALIZE
  Y <- data[[response]]
  temp_env <- list2env(list(y = Y, nobs = nobs,
                            etastart = NULL, start = NULL, mustart = NULL ))
  eval(family$initialize, envir = temp_env)
  Y <- temp_env$y # e.g. to encode factor variables to 0-1
  mustart <- temp_env$mustart
  # eta and mu
  eta <- link(mustart); mu <- linkinv(eta)
  # old.dev
  old.dev <- sum(dev.resids(Y, mu, weights))
  # models
  models = list(); n_iters <- NULL
  
  # 2. ITERATE
  for (iter in 1:max_iterations) {
    good <- weights > 0
    varmu <- variance(mu)
    if (any(is.na(varmu[good])))
      stop("NAs in V(mu)")
    if (any(varmu[good] == 0))
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good])))
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    z <- eta - offset
    z[good] <- z[good] + (Y - mu)[good]/mu.eta.val[good]
    wz <- weights
    wz[!good] <- 0
    wz[good] <- wz[good] * mu.eta.val[good]^2/varmu[good]
    
    # 2.3 Fit an additive model to z using the backfitting algorithm and weights wz
    bf <- backfitting(z = z, base_learners = base_learners, base_learners_list = base_learners_list,
                      data = data, wz = wz, phi_gam_exit = phi_gam_exit, q_gam = q_gam, iter = iter,
                      internal_knots = internal_knots, beta = beta, phi = phi, q = q)
    
    eta <- bf$z_hat + offset
    mu <- linkinv(eta)
    base_learners_list <- bf$base_learners_list
    n_iters <- c(n_iters, bf$n_iters)
    
    models[[paste0("iter", iter)]] <- list(Y_hat = list(eta = eta, mu = mu, z = z),
                                           base_learners = base_learners_list)
    
    # 3. UNTIL
    # If family = "gaussian", then localscoring = backfitting, hence no need of further iterations
    if(args$family$family == "gaussian") break
    
    new.dev <- sum(dev.resids(Y, mu, weights))
    if (is.infinite(new.dev)) stop("Infinite deviance: please decrease shrinkage parameter")
    if (new.dev < 1e-9) {
      cat("Deviance is almost zero: stopping boosting iterations\n\n")
      break
    }
    old.dev <- c(old.dev, new.dev)
    if (is.na(new.dev)) {
      break
      warning("iterations terminated prematurely because of singularities")
    }
    else if (!is.null(phi_gam_exit) && !is.null(q_gam) && iter >= q_gam && iter >= min_iterations) {
      # Check if the change in deviance is less than the tolerance level
      phi_gam = old.dev[iter+1]/old.dev[iter+1-q_gam]
      if (phi_gam >= phi_gam_exit) {
        cat("Stopping iterations due to small improvement in deviance\n\n")
        break
      }
    }
    # original stopping rule: else if ( abs(new.dev-old.dev)/(old.dev + 0.1) < stop ) break
    # old.dev <- new.dev
  }
  
  ## 4. Set the "final model" to be the one with lower deviance
  # 4.1. De-normalize predictions if necessary
  if (normalize_data == TRUE && args$family$family != "binomial") {
    models <- lapply(models, function(model) {
      model$Y_hat$mu = model$Y_hat$mu * args$Y_sd + args$Y_mean
      model$Y_hat$eta = link(model$Y_hat$mu)
      return(model)
    })
  }
  # 4.2. Calculate each model's deviance
  deviances <- sapply(models, function(model) {
    mu <- model$Y_hat$mu
    sum(dev.resids(args$response[[response]], mu, weights))
  })
  # 4.3. Find the model with the smallest deviance
  model_min_deviance <- names(models)[which.min(deviances)]
  final_model <- list(
    model_name = model_min_deviance,
    DEV = min(deviances),
    Y_hat = models[[model_min_deviance]]$Y_hat,      
    base_learners = models[[model_min_deviance]]$base_learners
  )
  # Save linear internal knots for each base-learner
  for (bl_name in names(base_learners)) {
    final_model$base_learners[[bl_name]]$linear.int.knots <- get_internal_knots(final_model$base_learners[[bl_name]]$knots)
  }
  
  # Extract variables of GeDS/linear base-learners
  GeDS_variables <- unname(unlist(lapply(base_learners, function(x) if (x$type == "GeDS") x$variables)))
  linear_variables <- unname(unlist(lapply(base_learners, function(x) if (x$type == "linear") x$variables)))
  
  #############################
  ## B-Spline Representation ##
  #############################
  # (i) Knots
  ll_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                   args$X_sd, args$X_mean, normalize_data, n = 2)
  
  # (ii) Coefficients
  # Filter GeDS base learners based on type and number of variables
  univariate_GeDS_learners <- Filter(function(bl) bl$type == "GeDS" && length(bl$variables) == 1, base_learners)
  bivariate_GeDS_learners <- Filter(function(bl) bl$type == "GeDS" && length(bl$variables) == 2, base_learners)
  
  # Extract B-spline coefficients
  univariate_GeDS_theta <- unlist(bSpline.coef(final_model, univariate_learners = univariate_GeDS_learners))
  bivariate_GeDS_theta <- lapply(names(bivariate_GeDS_learners),
                                 function(bl) final_model$base_learners[[bl]]$coefficients)
  for (i in seq_along(bivariate_GeDS_theta)) names(bivariate_GeDS_theta[[i]]) <- NULL
  bivariate_GeDS_theta <- unlist(setNames(bivariate_GeDS_theta, names(bivariate_GeDS_learners)))
    
  # Handle linear and factor bl/variables
  linear_coef <- NULL
  if (length(linear_variables) > 0) {
    # Initialize the intercept and slopes
    intercept <- 0; slopes <- NULL
    # Loop through each variable in linear_variables
    for (var in linear_variables) {
      coef <- final_model$base_learners[[var]]$coefficients
      int <- coef$b0
      slp <- unlist(coef[names(coef) != "b0"]); names(slp) <- paste0(var, names(slp))
      # Sum up the intercepts
      intercept <- intercept + int 
      # Vector of slopes
      slopes <- c(slopes, slp)
    }
    names(intercept) <- "b0"
    linear_coef <- c(intercept, slopes)
  } else {
    Z <- NULL  # Set Z to NULL if there are no factor variables
  }
  
  theta <- c(univariate_GeDS_theta, bivariate_GeDS_theta, linear_coef)
  
  linear_fit <- tryCatch({
    suppressMessages(
      SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$response[[response]],
                         Z = args$predictors[linear_variables], offset = args$offset + mean(final_model$Y_hat$z),
                         base_learners = args$base_learners, weights = weights,
                         InterKnotsList = ll_list, n = 2, family = args$family,
                         coefficients = theta, linear_intercept = TRUE, de_mean = TRUE)
      )
    }, error = function(e) {
      cat(paste0("Error computing linear fit:", e))
      return(NULL)
      })
  
  # De-normalize if necessary
  if (normalize_data == TRUE && family$family != "binomial") {
    linear_fit$Predicted <- as.numeric(linear_fit$Predicted) * args$Y_sd + args$Y_mean
  }
  
  final_model$Linear.Fit <- linear_fit
  
  pred_linear <- as.numeric(final_model$Y_hat$mu)
  
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
    
    # Quadratic fit
    qq_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                   args$X_sd, args$X_mean, normalize_data, n = 3)
    quadratic_fit <- tryCatch({
      suppressMessages(
        SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$response[[response]],
                           Z = args$predictors[linear_variables], offset = args$offset,
                           base_learners = args$base_learners, weights = weights,
                           InterKnotsList = qq_list, n = 3, family = family)
        )
      }, error = function(e) {
        cat(paste0("Error computing quadratic fit:", e))
        return(NULL)
        })
    final_model$Quadratic.Fit <- quadratic_fit
    pred_quadratic <- as.numeric(quadratic_fit$Predicted)
    
    # Cubic fit
    cc_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                     args$X_sd, args$X_mean, normalize_data, n = 4)
    cubic_fit <- tryCatch({
      suppressMessages(
        SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$response[[response]],
                           Z = args$predictors[linear_variables], offset = args$offset,
                           base_learners = args$base_learners, weights = weights,
                           InterKnotsList = cc_list, n = 4, family = family)
        )
      }, error = function(e) {
        cat(paste0("Error computing cubic fit:", e))
        return(NULL)
        })
    final_model$Cubic.Fit <- cubic_fit
    pred_cubic <- as.numeric(cubic_fit$Predicted)
    
    # Save quadratic and cubic knots for each base-learner
    for (bl_name in names(base_learners)){
      final_model$base_learners[[bl_name]]$quadratic.int.knots <- qq_list[[bl_name]]
      final_model$base_learners[[bl_name]]$cubic.int.knots <- cc_list[[bl_name]]
    }
  } else {
    pred_quadratic <-  pred_cubic <- NULL
    final_model$base_learners[[bl_name]]$quadratic.int.knots <- NULL
    final_model$base_learners[[bl_name]]$cubic.int.knots <- NULL
  }
  
  # Simplify final_model structure
  for (bl_name in names(base_learners)){
    final_model$base_learners[[bl_name]]$knots <- NULL
  }
  
  # Store predictions
  preds <- list(pred_linear = pred_linear, pred_quadratic = pred_quadratic, pred_cubic = pred_cubic)
  
  # Store internal knots
  linear.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  quadratic.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  cubic.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  # Loop through each base learner and extract the int.knots
  for(bl_name in names(base_learners)){
    # Extract and store linear internal knots
    linear_int.knt <- final_model$base_learners[[bl_name]]$linear.int.knots
    linear.int.knots[bl_name] <- list(linear_int.knt)
    # Extract and store quadratic internal knots
    quadratic_int.knt <- final_model$base_learners[[bl_name]]$quadratic.int.knots
    quadratic.int.knots[bl_name] <- list(quadratic_int.knt)
    # Extract and store cubic internal knots
    cubic_int.knt <- final_model$base_learners[[bl_name]]$cubic.int.knots
    cubic.int.knots[bl_name] <- list(cubic_int.knt)
  }
  # Combine the extracted knots into a single list
  internal_knots <- list(linear.int.knots = linear.int.knots, quadratic.int.knots = quadratic.int.knots,
                         cubic.int.knots = cubic.int.knots)
  
  output <- list(extcall = extcall, formula = formula, args = args, final_model = final_model, predictions = preds,
                 internal_knots = internal_knots, iters = list(local_scoring = iter, backfitting = n_iters))
  
  class(output) <- "GeDSgam"
  
  return(output)
}

###################
### backfitting ###
###################

#' @importFrom stats formula

backfitting <- function(z, base_learners, base_learners_list, data, wz, phi_gam_exit, q_gam,
                        iter, internal_knots, beta, phi, q)
{
  # 1. Initialize
  # (I) Intercept
  alpha <- rep(mean(z), length(z))
  # (II) Base-learners
  f <- matrix(0, nrow = length(z), ncol = length(base_learners))
  colnames(f) <- names(base_learners)
  
  ok <- TRUE; rss0 <- sum((z-alpha)^2); model_formula_template <- "partial_resid ~ "
  
  # 2. Cycle
  n_iters <- 0
  while (ok) { # backfitting loop
    n_iters <- n_iters + 1
    for (bl_name in names(base_learners)) { # loop through the smooth terms
      
      # 2.1. Fit bl_name to partial residuals
      partial_resid <- z - alpha - rowSums(f[, !colnames(f) %in% c(bl_name), drop = FALSE])
      
      pred_vars <- base_learners[[bl_name]]$variables
      data_loop <- cbind(partial_resid = partial_resid, data[pred_vars])
      
      # (A) GeDS base-learners
      if (base_learners[[bl_name]]$type == "GeDS") {
        
        if (length(pred_vars) == 2 && internal_knots == 0) {
          stop("internal_knots must be > 0 for bivariate learners")
        } else {
          max.intknots <- internal_knots
        }
        
        model_formula <- formula(paste0(model_formula_template, bl_name))
        error <- FALSE
        suppressWarnings({
          fit <- tryCatch(
            NGeDS(model_formula, data = data_loop, weights = wz, beta = beta, phi = phi,
                  min.intknots = 0, max.intknots = max.intknots, q = q, Xextr = NULL, Yextr = NULL,
                  show.iters = FALSE, stoptype = "RD", higher_order = FALSE, only_pred = TRUE),
            error = function(e) {
              message(paste0("Error occurred in NGeDS() for base learner ", bl_name, ": ", e))
              error <<- TRUE
              }
            )
          })
        # Continue to next bl if error
        if (error) {
          message(paste0("Skipped iteration ", iter," for base_learner ", bl_name, "."))
          next
          }
        # Compute predictions
        pred <- if (length(pred_vars) == 1) {
          predict_GeDS_linear(fit, data_loop[[pred_vars]])
          } else if (length(pred_vars) == 2) {
            predict_GeDS_linear(fit, X = data[pred_vars[1]], Y = data[pred_vars[2]], Z = data_loop[["partial_resid"]])
          }
        
        ## Update knots and coefficients ##
        # UNIVARIATE BASE-LEARNERS
        if (length(base_learners[[bl_name]]$variables) == 1) {
          base_learners_list[[bl_name]]$knots <- pred$knt
          base_learners_list[[bl_name]]$coefficients$b0 <- pred$b0
          base_learners_list[[bl_name]]$coefficients$b1 <- pred$b1
        # BIVARIATE BASE-LEARNERS
          } else if(length(base_learners[[bl_name]]$variables) == 2) {
            base_learners_list[[bl_name]]$knots <- pred$knt
            base_learners_list[[bl_name]]$coefficients <- pred$theta
          }
        
      # (B) Linear base-learners
      } else if (base_learners[[bl_name]]$type=="linear") {
        model_formula <- formula(paste0(model_formula_template, bl_name))
        error <- FALSE
        suppressWarnings({
          fit <- tryCatch(
            lm(model_formula, data = data_loop, weights = wz),
            error = function(e) {
              message(paste0("Error occurred in lm() for base learner ", bl_name, ": ", e))
              error <<- TRUE
            }
          )
        })
        # Continue to next bl if error
        if (error) {
          message(paste0("Skipped iteration ", iter," for base_learner ", bl_name, "."))
          next
        }
        
        # Compute predictions
        pred <- list(Y_hat = fit$fitted.values)
        ## Update coefficients ##
        for (i in 1:length(fit$coefficients)) {
          coef_name <- paste0("b", i - 1)
          base_learners_list[[bl_name]]$coefficients[[coef_name]] <- fit$coefficients[[i]]
        }
      }
      
      ## 2.2. Center predicted values
      f[, bl_name] <- pred$Y_hat - mean(pred$Y_hat)
    }
    
  # 3. Stopping rule
    rss <- sum((z-alpha-rowSums(f))^2)
    rss0 <- c(rss0, rss)
    if (!is.null(phi_gam_exit) && !is.null(q_gam) && length(rss0) > q_gam) {
      # Check if the change in SSR is less than the tolerance level
      phi_gam = rss0[length(rss0)]/rss0[length(rss0)-q_gam]
      if (phi_gam >= phi_gam_exit) ok <- FALSE
    }
    # original stopping rule: if (abs(rss-rss0)<stop*rss0) ok <- FALSE
    # rss0 <- rss
  }
  
  out <- list(z_hat = alpha+rowSums(f), base_learners_list = base_learners_list, n_iters = n_iters)
  return(out)
}
