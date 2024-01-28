#' NGeDSgam: Local Scoring Algorithm with GeD Splines in Backfitting
#'
#' @description Implements the Local Scoring Algorithm (Hastie and Tibshirani (1986)),
#' applying GeD splines to fit the targets within the backfitting iterations.
#' 
#' @importFrom stats setNames
#'
#' @param formula a formula specifying the model structure. It is expected to be in the
#'                format `Y ~ f(x1) + f(x2) + ...`.
#' @param data a data frame containing the variables referenced in the formula.
#' @param weights optional vector of weights for each observation.
#' @param offset optional offset to be added to the linear predictor during the backfitting.
#' @param normalize_data a logical that defines whether the data is normalized before fitting the model or not.
#' @param family a character string indicating the response distribution and link function to be used. Default is "gaussian".
#'               This should be a character or a family object. Accepted characters are "gaussian" and "binomial".
#' @param min_iterations Minimum number of iterations in the algorithm. Default is 0.
#' @param max_iterations Maximum number of iterations in the algorithm. Default is 100.
#' @param phi_gam_exit Convergence threshold for local-scoring. The algorithm stops when the relative change in the deviance is below this threshold.
#'             Default is 1e-3.
#' @param q_gam numeric parameter which allows to fine-tune the stopping rule of stage A of GeDS, by default equal to 2.
#' @param beta numeric parameter in the interval \eqn{[0,1]}
#' tuning the knot placement in stage A of GeDS. See details.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the threshold for
#'  the stopping rule  (model selector) in stage A of GeDS. See also \code{stoptype} and details below.
#' @param internal_knots maximum number of internal knots to be added by the GeDS learners within the backfitting iterations. Default is 500.
#' @param q numeric parameter which allows to fine-tune the stopping rule of stage A of GeDS, by default equal to 2.
#' @param higher_order a logical that defines whether to compute the higher order fits (quadratic and cubic).
#'
#' @return A list of class "GeDSgam" containing:
#'   * formula: The formula provided as input.
#'   * args: List of arguments and relevant data used internally.
#'   * final_model: The final fitted model after local scoring.
#'   * predictions: Predictions for linear, quadratic, and cubic fits.
#'
#' @examples
#' # Load package
#' library(GeDS) 
#' 
#' data(airquality) 
#' data = na.omit(airquality)
#' data$Ozone <- data$Ozone^(1/3)
#' 
#' formula = Ozone ~ f(Solar.R) + f(Wind, Temp)
#' Gmod_gam <- NGeDSgam(formula = formula, data = data,
#' phi_gam_exit = 0.995, phi = 0.995, q = 2)
#' mse_ngedsgam1 <- mean((data$Ozone - Gmod_gam$predictions$pred_linear)^2)
#' mse_ngedsgam2 <- mean((data$Ozone - Gmod_gam$predictions$pred_quadratic)^2)
#' mse_ngedsgam3 <- mean((data$Ozone - Gmod_gam$predictions$pred_cubic)^2)
#' 
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#' "Linear NGeDSgam:", mse_ngedsgam1, "\n",
#' "Quadratic NGeDSgam:", mse_ngedsgam2, "\n",
#' "Cubic NGeDSgam:", mse_ngedsgam3, "\n")
#'
#' @export
#' @seealso \link{gam}, \link{glm}
#' @references
#' Hastie, T., & Tibshirani, R. (1986). Generalized Additive Models. Statistical Science, 1(3), 297–310. http://www.jstor.org/stable/2245459

################################################################################
################################################################################
################# Local Scoring Algorithm for Additive Models ##################
################################################################################
################################################################################

#####################
### Local Scoring ###
#####################
NGeDSgam <- function(formula, data, weights = NULL, offset = NULL, normalize_data = FALSE, family = "gaussian",
                     min_iterations = 0, max_iterations = 100,
                     phi_gam_exit = 0.995, q_gam = 2,
                     beta = 0.5, phi = 0.99, internal_knots = 500, q = 2,
                     higher_order = TRUE)
  {
  
  # Capture the function call
  extcall <- match.call()
  
  # Formula
  read.formula <- read.formula.gam(formula, data)
  terms <-  read.formula$terms
  outcome <- read.formula$outcome
  predictors <- read.formula$predictors
  base_learners <- read.formula$base_learners
  
  # Eliminate indexes and keep only relevant variables
  rownames(data) <- NULL
  all_vars <- c(outcome, predictors)
  if(all(all_vars %in% names(data))) {
    data <- data[, all_vars]
  } else {
    stop("Error: Some of the variables are not found in the data.")
  }
  # Check for NA in each column and print a warning message if found
  na_vars <- apply(data, 2, function(x) any(is.na(x)))
  for(var in names(data)[na_vars]){
    warning(paste("variable", var, "contains missing values;\n",
                  "these observations are not used for fitting"))
  }
  # Eliminate rows containing NA values
  data <- na.omit(data)
  
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
  
  variance <- family$variance
  dev.resids <- family$dev.resids
  link <- family$linkfun # g
  linkinv <- family$linkinv # g^{-1}
  mu.eta <- family$mu.eta   # dg^{-1}/d\eta
  
  
  # Save arguments
  args <- list(
    "predictors" = data[predictors], 
    "base_learners"= base_learners,
    "normalize_data" = normalize_data,
    "family" = family
  )
  
  if (args$family$family != "binomial") {
    args$outcome <- data[outcome]
  } else {
    args$outcome <- data.frame(as.numeric(data[[outcome]])-1)
    names(args$outcome) <- outcome
  }
  
  # Normalize data if necessary
  if (normalize_data == TRUE) {
    if (args$family$family != "binomial") {
      # Mean and SD of the original outcome and predictor(s) variables (to de-normalize afterwards)
      args$Y_mean <- mean(data[[outcome]])
      args$Y_sd <- sd(data[[outcome]])
      data[outcome] <- as.vector(scale(data[outcome]))
      
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
  
  # Data matrices
  Y <- data[[outcome]]
  X <- as.matrix(data[predictors])
  
  # Weights and offset
  nobs = length(Y)
  if (is.null(weights)) weights <- rep.int(1, nobs)
  else weights <- rescale_weights(weights)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  
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
  models = list()
  
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
    
    models[[paste0("iter", iter)]] <- list(Y_hat = list(eta = eta, mu = mu, z = z),
                                           base_learners = base_learners_list)
    
    # 3. UNTIL
    # If family = "gaussian", then localscoring = backfitting, hence no need of further iterations
    if(args$family$family == "gaussian") break
    
    new.dev <- sum(dev.resids(Y, mu, weights))
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
  # Use sapply to calculate deviances for each model
  deviances <- sapply(models, function(model) {
    mu <- model$Y_hat$mu
    sum(dev.resids(args$outcome[[outcome]], mu, weights))
  })
  # Find the model with the smallest deviance
  model_min_deviance <- names(models)[which.min(deviances)]
  # Set the final model to the model with the smallest MSE
  final_model <- models[[model_min_deviance]]
  final_model$model_name <- model_min_deviance
  final_model$RSS <- min(deviances)
  final_model$best_bl <- NULL # to simplify the final_model output
  # Save predictions
  if (normalize_data == TRUE && args$family$family != "binomial"){
    final_model$Y_hat$mu = final_model$Y_hat$mu * args$Y_sd + args$Y_mean
    final_model$Y_hat$eta = link(final_model$Y_hat$mu)
  } 
  
  #######################
  ## Higher order fits ##
  #######################
  if (higher_order) {
    
    # Extract variables of GeDS/linear base-learners
    GeDS_variables <- lapply(base_learners, function(x) {if(x$type == "GeDS") return(x$variables) else return(NULL)})
    GeDS_variables <- unname(unlist(GeDS_variables))
    linear_variables <- lapply(base_learners, function(x) {if(x$type == "linear") return(x$variables) else return(NULL)})
    linear_variables <- unname(unlist(linear_variables))
    
    # Quadratic fit
    qq_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                   args$X_sd, args$X_mean, normalize_data, n = 3)
    quadratic_fit <- tryCatch({
      SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$outcome[[outcome]],
                         Z = args$predictors[linear_variables],
                         base_learners = args$base_learners, InterKnotsList = qq_list,
                         n = 3, family = family)}, error = function(e) {
                           cat(paste0("Error computing quadratic fit:", e))
                           return(NULL)
                           })
    final_model$Quadratic.Fit <- quadratic_fit
    pred_quadratic <- as.numeric(quadratic_fit$Predicted)
    
    # Cubic fit
    cc_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                     args$X_sd, args$X_mean, normalize_data, n = 4)
    cubic_fit <- tryCatch({
      SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$outcome[[outcome]],
                         Z = args$predictors[linear_variables],
                         base_learners = args$base_learners, InterKnotsList = cc_list,
                         n = 4, family = family)}, error = function(e) {
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
  
  preds <- list(pred_linear = as.numeric(final_model$Y_hat$mu), pred_quadratic = pred_quadratic, pred_cubic = pred_cubic)
  
  # Store internal knots
  linear.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  quadratic.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  cubic.int.knots <- setNames(vector("list", length(base_learners)), names(base_learners))
  # Loop through each base learner and extract the int.knots
  for(bl_name in names(base_learners)){
    # Extract and store linear knots
    linear_knots <- get_internal_knots(final_model$base_learners[[bl_name]]$knots)
    linear.int.knots[bl_name] <- list(linear_knots)
    # Extract and store quadratic knots
    quadratic_knots <- final_model$base_learners[[bl_name]]$quadratic.int.knots
    quadratic.int.knots[bl_name] <- list(quadratic_knots)
    # Extract and store cubic knots
    cubic_knots <- final_model$base_learners[[bl_name]]$cubic.int.knots
    cubic.int.knots[bl_name] <- list(cubic_knots)
  }
  # Combine the extracted knots into a single list
  internal_knots <- list(linear.int.knots = linear.int.knots, quadratic.int.knots = quadratic.int.knots,
                         cubic.int.knots = cubic.int.knots)
  
  output <- list(extcall = extcall, formula = formula, args = args, final_model = final_model, predictions = preds,
                 internal_knots = internal_knots)
  
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
  while (ok) { # backfitting loop
    for (bl_name in names(base_learners)) { # loop through the smooth terms
      
      # 2.1. Fit bl_name to partial residuals
      partial_resid <- z - alpha - rowSums(f[, !colnames(f) %in% c(bl_name), drop = FALSE])
      
      pred_vars <- base_learners[[bl_name]]$variables
      data_loop <- cbind(partial_resid = partial_resid, data[pred_vars])
      # (A) GeDS base-learners
      if (base_learners[[bl_name]]$type == "GeDS") {
        max.intknots <- max.intknots <- if (length(pred_vars) == 1) {internal_knots + q
        } else if (length(pred_vars) == 2 && internal_knots == 0) {stop("internal_knots must be > 0 for bivariate learners")
            } else {internal_knots}
        
        model_formula <- formula(paste0(model_formula_template, bl_name))
        error <- FALSE
        suppressWarnings({
          fit <- tryCatch(
            NGeDS(model_formula, data = data_loop, weights = wz, beta = beta, phi = phi,
                  min.intknots = 0, max.intknots = max.intknots, q = q, Xextr = NULL, Yextr = NULL,
                  show.iters = FALSE, stoptype = "RD"),
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
        pred <- if(length(pred_vars) == 1){
          predict_GeDS_linear(fit, data_loop[[pred_vars]])
          } else if (length(pred_vars) == 2) {
            predict_GeDS_linear(fit, X = data[pred_vars[1]], Y = data[pred_vars[2]], Z = data_loop[["partial_resid"]])
          }
        
        ## Update knots and coefficients ##
        # UNIVARIATE BASE-LEARNERS
        if(length(base_learners[[bl_name]]$variables) == 1) {
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
  
  out <- list(z_hat = alpha+rowSums(f), base_learners_list = base_learners_list)
  return(out)
}
