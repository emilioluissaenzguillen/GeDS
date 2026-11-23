#########################################################################################
#########################################################################################
#################################### Cross-validation ###################################
#########################################################################################
#########################################################################################
# Declare loop variables
globalVariables(c("i", "j", "k", "l", "m"))

# Define a function to perform k-fold cross-validation
#' @title K-Fold Cross-Validation
#' @name crossv_GeDS
#' @description
#' \code{crossv_GeDS} performs k-fold cross-validation for tuning the relevant
#' parameters of \code{NGeDS}, \code{GGeDS}, \code{NGeDSgam}, and
#' \code{NGeDSboost} functions.
#'
#' @param formula A description of the structure of the model structure,
#' including the dependent and independent variables.
#' @param data A \code{data.frame} containing the variables referenced in the formula.
#' @param model_fun The GeDS model to cross-validate, that is, \code{NGeDS},
#' \code{GGeDS}, \code{NGeDSgam} or \code{NGeDSboost}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the spline fit evaluated during cross-validation. The default
#' is \code{2L}.
#' @param n_folds Integer from 2 to the number of observations specifying the
#' number of cross-validation folds. The default is \code{5L}. For \code{GGeDS} and \code{NGeDSgam},
#' binomial responses are stratified across folds.
#' @param parameters A named list of numeric grids to be tuned. The supported
#' grids are \code{beta_grid}, \code{phi_grid}, and \code{q_grid} for
#' \code{NGeDS}, \code{GGeDS}, and \code{NGeDSgam}; \code{NGeDSboost} also
#' supports \code{int.knots_init_grid} and \code{shrinkage_grid}. Omitted
#' parameters retain the fitting function's defaults and are not shown in the
#' results. The default is an empty list, which evaluates the default fit.
#' Select \code{stoptype = "RD"} in
#' \code{...} to tune its threshold with \code{phi_grid}. If a threshold grid
#' is supplied without a rule, cross-validation uses the previous default
#' (univariate \code{NGeDS}, \code{NGeDSgam} smoothers and \code{NGeDSboost} use RD, and joint
#' \code{NGeDS} and \code{GGeDS} use SR). Use the same rule when refitting
#' with the selected parameters. An explicit \code{stoptype = "GD"} cannot
#' be combined with \code{phi_grid}, since GD does not use \code{phi}.
#' @param n_cores Integer specifying the number of cores to use for parallel
#' cross-validation. If \code{NULL}, the function uses one fewer than the number
#' of detected cores, subject to CRAN check limits.
#' @param ... Additional arguments forwarded to the fitting function in each
#' fold, such as \code{family = binomial()} for \code{NGeDSgam} or \code{GGeDS}.
#'
#' @return Two data frames, \code{best_params} and \code{results}.
#' \code{best_params} contains the best parameter combination and its scores.
#' Only supplied parameter grids appear as parameter columns. For
#' \code{NGeDS}, \code{GGeDS}, and \code{NGeDSgam}, \code{results} reports
#' cross-validated mean squared error (\code{cv_mse}) for Gaussian fits or
#' mean deviance (\code{cv_deviance}) for non-Gaussian fits. For
#' \code{NGeDSboost}, it reports \code{cv_mse} and the cross-validated number
#' of boosting iterations. Average internal-knot counts are included where
#' available.
#'
#' @examples
#' ###################################################
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#' data = data.frame(X = X, Y = Y)
#'
#' \dontrun{
#' ## NGeDS
#' # Define different combinations of parameters to cross-validate
#' param = list(beta_grid = c(0.5),
#'              phi_grid = c(0.9, 0.95),
#'              q_grid = c(2))
#'
#' cv_NGeDS <- crossv_GeDS(Y ~ f(X), data = data, NGeDS, n = 3,
#'                         parameters = param)
#'
#' print(cv_NGeDS$best_params)
#' View(cv_NGeDS$results)
#'
#' ## NGeDSboost
#' param = list(int.knots_init_grid = c(1, 2),
#'              shrinkage_grid = 1,
#'              beta_grid = c(0.3, 0.5),
#'              phi_grid = c(0.95, 0.99),
#'              q_grid = 2)
#'
#' cv_NGeDSboost <- crossv_GeDS(Y ~ f(X), data = data, NGeDSboost, n = 2L,
#'                              n_folds = 2L, parameters = param)
#'
#' print(cv_NGeDSboost$best_params)
#' View(cv_NGeDSboost$results)
#' }
#'
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom stats predict
#' @importFrom utils globalVariables
#' @export
crossv_GeDS <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                        parameters = list(), n_cores = NULL, ...)
  {
  # Ensure model_fun is a supported function
  model_fun_name <- deparse(substitute(model_fun))
  supported_models <- c("NGeDS", "GGeDS", "NGeDSgam", "NGeDSboost")
  if (!model_fun_name %in% supported_models) {
    stop("Unsupported model function. Please use one of the following: ", paste(supported_models, collapse=", "))
  }

  # Cross-validation based on the model function
  if (model_fun_name %in% c("NGeDS", "GGeDS", "NGeDSgam")) {
    cv <- cross_validate.GeDS(formula = formula, data = data, model_fun = model_fun,
                              n = n, n_folds = n_folds, parameters = parameters,
                              n_cores = n_cores, ...)
    } else if (model_fun_name %in% c("NGeDSboost")) {
      cv <- cross_validate.GeDSboost(formula = formula, data = data, model_fun = model_fun,
                                     n = n, n_folds = n_folds, parameters = parameters,
                                     n_cores = n_cores, ...)
    }

  return(list(best_params = cv$best_params, results = as.data.frame(cv$results)))
}

################################################################################
############################ NGeDS, GGeDS, NGeDSgam ############################
################################################################################
cross_validate.GeDS <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                                parameters = list(),
                                n_cores = NULL, ...)
{
  # Parse formula
  terms <- all.vars(formula)
  outcome <- terms[1]
  predictors <- terms[-1]
  # Response and predictor data
  Y <- data[[outcome]]
  # Get the number of observations
  n_obs <- nrow(data)

  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }

  # If cv model is linear, avoid computing higher order fits
  if (n == 2) higher_order <- FALSE else higher_order <- TRUE

  # Validate n_folds
  if (is.null(n_folds) || !is.numeric(n_folds) || n_folds <= 0) {
    n_folds <- 5L
    warning("n_folds was NULL, non-numeric, or less than or equal to 0; set to default: 5.")
  } else {
    n_folds <- as.integer(n_folds)
  }
  if (n_folds < 2L || n_folds > n_obs) {
    stop("'n_folds' must be between 2 and the number of observations.")
  }
  fit_args <- list(...)
  family_obj <- fit_args$family
  if (is.null(family_obj)) family_obj <- gaussian()
  if (is.character(family_obj)) {
    family_obj <- get(family_obj, mode = "function", envir = parent.frame())
  }
  if (is.function(family_obj)) family_obj <- family_obj()
  if (!inherits(family_obj, "family")) stop("'family' must specify an R family.")
  is_binomial <- family_obj$family %in% c("binomial", "quasibinomial")
  score_name <- if (family_obj$family == "gaussian") "cv_mse" else "cv_deviance"

  # Balance the response classes when forming binomial folds.
  if (is_binomial) {
    shuffled <- unlist(lapply(split(seq_len(n_obs), Y), function(idx) {
      idx[sample.int(length(idx))]
    }), use.names = FALSE)
    fold_id <- integer(n_obs)
    fold_id[shuffled] <- rep(sample.int(n_folds), length.out = n_obs)
    folds <- split(seq_len(n_obs), fold_id)
  } else {
    folds <- split(sample(seq_len(n_obs)), rep(seq_len(n_folds), length.out = n_obs))
  }

  grids <- prepare_cv_grids(parameters, c("beta_grid", "phi_grid", "q_grid"))
  beta_grid <- grids$values$beta_grid
  phi_grid <- grids$values$phi_grid
  q_grid <- grids$values$q_grid
  tuned <- grids$tuned

  if ("phi" %in% tuned && is.null(fit_args$stoptype)) {
    if (identical(model_fun, NGeDS)) {
      fit_args$stoptype <- if (ncol(read.formula(formula, data)$X) == 1L) "RD" else "SR"
    } else if (identical(model_fun, NGeDSgam)) {
      fit_args$stoptype <- "RD"
    } else if (identical(model_fun, GGeDS)) {
      fit_args$stoptype <- "SR"
    }
  }
  if ("phi" %in% tuned && is.character(fit_args$stoptype) &&
      identical(pmatch(fit_args$stoptype, c("GD", "RD", "SR", "LR")), 1L)) {
    stop("Choose a fixed-threshold stopping rule to tune 'phi'.", call. = FALSE)
  }



  # Register parallel backend
  available_cores <- detectCores()
  if (is.na(available_cores)) available_cores <- 1L

  if (is.null(n_cores)) {
    n_cores <- max(1L, available_cores - 1L)
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }

  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    n_cores <- min(n_cores, 2L)
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Parallel processing
  results_list <- foreach(i = 1:length(beta_grid), .combine = rbind, .packages = "GeDS") %:%
    foreach(j = 1:length(phi_grid), .combine = rbind) %:%
    foreach(k = 1:length(q_grid), .combine = rbind) %dopar% {

      # A failed fold makes the candidate ineligible rather than improving its score.
      score_fold <- rep(Inf, n_folds)
      intknots_fold <- list(univbl = rep(NA_real_, n_folds),
                             bivblX = rep(NA_real_, n_folds),
                             bivblY = rep(NA_real_, n_folds))
      time_fold <- rep(NA_real_, n_folds)
      last_fit <- NULL

      for (fold_idx in 1:n_folds) {
        # Create the training and test sets for the current fold
        train_data <- data[-folds[[fold_idx]], ]
        test_data <- data[folds[[fold_idx]], ]
        rownames(train_data) <- rownames(test_data) <- NULL
        # Fit the model
        start_time <- Sys.time()
        candidate <- list(beta = beta_grid[i], phi = phi_grid[j], q = q_grid[k])
        fit <- tryCatch(
          do.call(model_fun, c(list(formula = formula, data = train_data,
                                    higher_order = higher_order),
                               candidate[tuned], fit_args)),
          error = function(e) {
            message(paste0("Error occurred for fold ", fold_idx, ": ", e))
            NULL
          }
        )
        if (is.null(fit)) next
        last_fit <- fit
        end_time <- Sys.time()
        time_fold[fold_idx] <- difftime(end_time, start_time, units = "secs")
        # Evaluate the model on the test set
        test_pred <- tryCatch(predict(fit, newdata = test_data, n = n),
                              error = function(e) NULL)
        if (is.null(test_pred) || any(!is.finite(test_pred))) next

        response <- test_data[[outcome]]
        if (is_binomial && is.factor(response)) response <- as.numeric(response) - 1L
        fold_score <- tryCatch({
          if (score_name == "cv_mse") {
            mean((response - test_pred)^2)
          } else {
            mean(family_obj$dev.resids(
              response, test_pred, rep(1, length(response))))
          }
        }, error = function(e) Inf)
        if (is.finite(fold_score)) score_fold[fold_idx] <- fold_score

        # Number of internal knots
        if (inherits(fit, "GeDS")) {
          if (fit$type == "LM - Univ" || fit$type == "GLM - Univ") {
            intknots_fold$univbl[fold_idx] <- fit$Nintknots
          } else if (fit$type == "LM - Biv" || fit$type == "GLM - Biv") {
            intknots_fold$bivblX[fold_idx] <- fit$Nintknots$X
            intknots_fold$bivblY[fold_idx] <- fit$Nintknots$Y
          }
        }

      }

      if (is.null(last_fit) && !identical(model_fun, NGeDSgam)) {
        stop("All folds failed for a GeDS parameter combination.")
      }
      cv_score <- if (score_name == "cv_mse") {
        mean(score_fold)
      } else {
        weighted.mean(score_fold, lengths(folds))
      }

      result <- data.frame(score = cv_score)
      if (length(tuned)) result <- cbind(as.data.frame(candidate[tuned]), result)
      if (inherits(last_fit, "GeDS")) {
        if (last_fit$type == "LM - Univ" || last_fit$type == "GLM - Univ") {
          result$avg_intknots <- mean(intknots_fold$univbl, na.rm = TRUE)
        } else if (last_fit$type == "LM - Biv" || last_fit$type == "GLM - Biv") {
          result$avg_intknots_X <- mean(intknots_fold$bivblX, na.rm = TRUE)
          result$avg_intknots_Y <- mean(intknots_fold$bivblY, na.rm = TRUE)
        }
      }
      result$avg_fitting_time <- mean(time_fold, na.rm = TRUE)
      result

    }

  # Stop the parallel backend
  stopCluster(cl)

  results_df <- as.data.frame(results_list)
  names(results_df)[names(results_df) == "score"] <- score_name
  if (!any(is.finite(results_df[[score_name]]))) {
    stop("All cross-validation candidates failed.")
  }

  # Find the optimal parameters that minimize the cross-validated score
  best_params <- results_df[
    which.min(results_df[[score_name]]),
    ,
    drop = FALSE
  ]

  # Return the optimal parameters
  return(list(best_params = best_params, results = results_df))
}

################################################################################
################################## NGeDSboost ##################################
################################################################################
cross_validate.GeDSboost <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                                     parameters = list(),
                                     n_cores = NULL, ...)
{
  fit_args <- list(...)
  # Parse formula
  terms <- all.vars(formula)
  outcome <- terms[1]
  predictors <- terms[-1]
  # Data matrices
  Y <- data[[outcome]]
  X <- as.matrix(data[predictors])
  # Get the number of observations
  n_obs <- nrow(data)

  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }

  # If cv GeDSboost/GeDSgam model is linear, avoid computing higher order fits
  if (n == 2) higher_order <- FALSE else higher_order <- TRUE

  # Validate n_folds
  if (is.null(n_folds) || !is.numeric(n_folds) || n_folds <= 0) {
    n_folds <- 5L
    warning("n_folds was NULL, non-numeric, or less than or equal to 0; set to default: 5.")
  } else {
    n_folds <- as.integer(n_folds)
  }
  # Generate a list of fold indices
  folds <- split(sample(1:n_obs), rep(1:n_folds, len = n_obs))

  grids <- prepare_cv_grids(parameters, c("int.knots_init_grid", "shrinkage_grid",
                                         "beta_grid", "phi_grid", "q_grid"))
  int.knots_init_grid <- grids$values$int.knots_init_grid
  shrinkage_grid <- grids$values$shrinkage_grid
  beta_grid <- grids$values$beta_grid
  phi_grid <- grids$values$phi_grid
  q_grid <- grids$values$q_grid
  tuned <- grids$tuned

  if ("phi" %in% tuned) {
    if (is.null(fit_args$stoptype)) fit_args$stoptype <- "RD"
    if (is.character(fit_args$stoptype) &&
        identical(pmatch(fit_args$stoptype, c("GD", "RD")), 1L)) {
      stop("Choose a fixed-threshold stopping rule to tune 'phi'.", call. = FALSE)
    }
  }

  # Register parallel backend
  available_cores <- detectCores()
  if (is.na(available_cores)) available_cores <- 1L

  if (is.null(n_cores)) {
    n_cores <- max(1L, available_cores - 1L)
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }

  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    n_cores <- min(n_cores, 2L)
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Parallel processing
  results_list <- foreach(i = 1:length(int.knots_init_grid), .combine = rbind, .packages = "GeDS",
                          .export = c("n_intknots")) %:%
    foreach(j = 1:length(shrinkage_grid), .combine = rbind) %:%
    foreach(k = 1:length(beta_grid), .combine = rbind) %:%
    foreach(l = 1:length(phi_grid), .combine = rbind) %:%
    foreach(m = 1:length(q_grid), .combine = rbind) %dopar% {

      # A failed fold makes the candidate ineligible.
      mse_fold <- rep(Inf, n_folds)
      n_iterations_fold <- rep(NA_real_, n_folds)
      intknots_fold <- list(univbl = rep(NA_real_, n_folds),
                            bivblX = rep(NA_real_, n_folds),
                            bivblY = rep(NA_real_, n_folds))
      time_fold <- rep(NA_real_, n_folds)

      for (fold_idx in 1:n_folds) {
        # Create the training and test sets for the current fold
        train_data <- data[-folds[[fold_idx]], ]
        test_data <- data[folds[[fold_idx]], ]
        rownames(train_data) <- rownames(test_data) <- NULL
        # Fit the model
        start_time <- Sys.time()
        error <- FALSE
        candidate <- list(int.knots_init = int.knots_init_grid[i],
                          shrinkage = shrinkage_grid[j], beta = beta_grid[k],
                          phi = phi_grid[l], q = q_grid[m])
        fit <- tryCatch(
          do.call(model_fun, c(list(formula = formula, data = train_data,
                                    higher_order = higher_order),
                               candidate[tuned], fit_args)),
          error = function(e) {
            message(paste0("Error occurred for fold ", fold_idx, ": ", e))
            error <<- TRUE
          }
        )
        if (error) next
        end_time <- Sys.time()
        time_fold[fold_idx] <- difftime(end_time, start_time, units = "secs")
        # Evaluate the model on the test set
        test_pred <- predict(fit, newdata = test_data, n = n)

        # Calculate the mean squared error for the test set
        mse_fold[fold_idx] <- mean((test_data[[outcome]] - test_pred)^2)

        # Number of boosting iterations and total number of internal knots
        aux <- n_intknots(fit, avg = FALSE)
        n_iterations_fold[fold_idx] <- aux$n_iterations
        intknots_fold$univbl[fold_idx] <- aux$n_intknots$univbl
        intknots_fold$bivblX[fold_idx] <- aux$n_intknots$bivblX
        intknots_fold$bivblY[fold_idx] <- aux$n_intknots$bivblY
      }

      result <- data.frame(cv_mse = mean(mse_fold),
                           n_iterations = mean(n_iterations_fold),
                           avg_intknots_univbl = mean(intknots_fold$univbl),
                           avg_intknots_bivblX = mean(intknots_fold$bivblX),
                           avg_intknots_bivblY = mean(intknots_fold$bivblY),
                           avg_fitting_time = mean(time_fold))
      if (length(tuned)) result <- cbind(as.data.frame(candidate[tuned]), result)
      result
    }

  # Stop the parallel backend
  stopCluster(cl)

  results_df <- as.data.frame(results_list)
  if (!any(is.finite(results_df$cv_mse))) {
    stop("All cross-validation candidates failed.")
  }

  best_params <- results_df[
    which.min(results_df$cv_mse),
    ,
    drop = FALSE
  ]

  return(list(best_params = best_params, results = results_df))
}



################################################################################
############################ NGeDS, GGeDS, NGeDSgam ############################
################################################################################
# cross_validate.GeDS <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
#                                 parameters = list(beta_grid = c(0.3, 0.5, 0.7),
#                                                   phi_grid = c(0.9, 0.95, 0.99),
#                                                   q_grid = c(2, 3)),...)
# {
#   # Parse formula
#   terms <- all.vars(formula)
#   outcome <- terms[1]
#   predictors <- terms[-1]
#   # Data matrices
#   Y <- data[[outcome]]
#   X <- as.matrix(data[predictors])
#   # Get the number of observations
#   n_obs <- nrow(data)
#
#   # Check if order is correctly set
#   n <- as.integer(n)
#   if(!(n %in% 2L:4L)) {
#     n <- 3L
#     warning("'n' incorrectly specified. Set to 3.")
#   }
#
#   # If cv model is linear, avoid computing higher order fits
#   if (n == 2) higher_order <- FALSE else higher_order <- TRUE
#
#   # Validate n_folds
#   if (is.null(n_folds) || !is.numeric(n_folds) || n_folds <= 0) {
#     n_folds <- 5L
#     warning("n_folds was NULL, non-numeric, or less than or equal to 0; set to default: 5.")
#   } else {
#     n_folds <- as.integer(n_folds)
#   }
#   # Generate a list of fold indices
#   folds <- split(sample(1:n_obs), rep(1:n_folds, len = n_obs))
#
#   # Extract parameters
#   beta_grid <- param_gridCheck(parameters, "beta_grid", c(0.3, 0.5, 0.7))
#   phi_grid <- param_gridCheck(parameters, "phi_grid", c(0.9, 0.95, 0.99))
#   q_grid <- param_gridCheck(parameters, "q_grid", c(2, 3))
#
#   # Initialize a dataframe to store hyperparameters and corresponding MSE values
#   results_df <- data.frame()
#
#   # Loop through each combination
#   for (i in 1:length(beta_grid)) {
#     for (j in 1:length(phi_grid)) {
#       for (k in 1:length(q_grid)) {
#         # Use tryCatch to handle errors
#         tryCatch({
#           # Initialize a vector to store MSE values for each fold
#           mse_fold <- numeric(n_folds)
#           n_iterations_fold <- numeric(n_folds)
#           avg_intknots_fold <- list(avg_intknots_univbl = numeric(n_folds),
#                                     avg_intknots_bivblX = numeric(n_folds),
#                                     avg_intknots_bivblY = numeric(n_folds))
#           time_fold <- numeric(n_folds)
#           for (fold_idx in 1:n_folds) {
#             # Create the training and test sets for the current fold
#             train_data <- data[-folds[[fold_idx]], ]
#             test_data <- data[folds[[fold_idx]], ]
#             rownames(train_data) <- rownames(test_data) <- NULL
#             # Fit the model
#             start_time <- Sys.time()
#             error <- FALSE
#             fit <- tryCatch(
#               model_fun(formula = formula, data = train_data, beta = beta_grid[i], phi = phi_grid[j],
#                         q = q_grid[k], higher_order = higher_order, ...),
#               error = function(e) {
#                 message(paste0("Error occurred for fold ", fold_idx, ": ", e))
#                 error <<- TRUE
#               }
#             )
#             if (error) next
#             end_time <- Sys.time()
#             time_fold[fold_idx] <- difftime(end_time, start_time, units = "secs")
#             # Evaluate the model on the test set
#             test_pred <- predict(fit, newdata = test_data, n = n)
#
#             # Calculate the mean squared error for the test set
#             mse_fold[fold_idx] <- mean((test_data[[outcome]] - test_pred)^2)
#
#             # Number of boosting iterations and average number of internal knots per boosting iteration
#             if (inherits(fit, "GeDS")) {
#               n_iterations_fold[fold_idx] <- fit$iters
#               if (fit$type == "LM - Univ" || fit$type == "GLM - Univ") {
#                 avg_intknots_fold$avg_intknots_univbl[fold_idx] <- fit$Nintknots
#                 avg_intknots_fold$avg_intknots_bivblX[fold_idx] <- NA
#                 avg_intknots_fold$avg_intknots_bivblY[fold_idx] <- NA
#                 } else if (fit$type == "LM - Biv" || fit$type == "GLM - Biv") {
#                   avg_intknots_fold$avg_intknots_univbl[fold_idx] <- NA
#                   avg_intknots_fold$avg_intknots_bivblX[fold_idx] <- fit$Nintknots$X
#                   avg_intknots_fold$avg_intknots_bivblY[fold_idx] <- fit$Nintknots$Y
#                 }
#             } else if (inherits(fit, "GeDSgam")) {
#               aux <- n_intknots(fit)
#               n_iterations_fold[fold_idx] <- mean(fit$iters$backfitting)
#             }
#
#           }
#           # Store hyperparameters and corresponding MSE value in results_df
#           results_df <- rbind(results_df, data.frame(beta = beta_grid[i],
#                                                      phi = phi_grid[j],
#                                                      q = q_grid[k],
#                                                      mse = mean(mse_fold),
#                                                      n_iterations = mean(n_iterations_fold),
#                                                      avg_intknots_univbl = mean(avg_intknots_fold$avg_intknots_univbl),
#                                                      avg_intknots_bivblX = mean(avg_intknots_fold$avg_intknots_bivblX),
#                                                      avg_intknots_bivblY = mean(avg_intknots_fold$avg_intknots_bivblY),
#                                                      avg_fitting_time = mean(time_fold)
#           ))
#         },
#         # Specify what to do if an error occurs
#         error = function(e) {
#           warning(paste("Error for beta=", beta_grid[i], ", phi=", phi_grid[j],
#                         ", q=", q_grid[k], ". Iteration was skipped."))
#           return()
#         })
#       }
#     }
#   }
#
#   # Find the optimal internal_knots, max_iterations and shrinkage values that minimize the cross-validated mse
#
#    results_df <- as.data.frame(results_matrix)
#    best_params <- results_df[
#      which.min(results_df$mse),
#      ,
#      drop = FALSE
#      ]
#
#   # Return the optimal internal_knots, max_iterations and shrinkage values
#   return(list(best_params = best_params, results = results_df))
# }

################################################################################
################################## NGeDSboost ##################################
################################################################################
# cross_validate.GeDSboost <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
#                                  parameters = list(int.knots_init_grid = c(0, 1, 2),
#                                                    shrinkage_grid = c(0.1, 0.5, 1),
#                                                    beta_grid = c(0.3, 0.5, 0.7),
#                                                    phi_grid = c(0.9, 0.95, 0.99),
#                                                    q_grid = c(2, 3)),...)
# {
#   # Parse formula
#   terms <- all.vars(formula)
#   outcome <- terms[1]
#   predictors <- terms[-1]
#   # Data matrices
#   Y <- data[[outcome]]
#   X <- as.matrix(data[predictors])
#   # Get the number of observations
#   n_obs <- nrow(data)
#
#   # Check if order is correctly set
#   n <- as.integer(n)
#   if(!(n %in% 2L:4L)) {
#     n <- 3L
#     warning("'n' incorrectly specified. Set to 3.")
#   }
#
#   # If cv GeDSboost/GeDSgam model is linear, avoid computing higher order fits
#   if (n == 2) higher_order <- FALSE else higher_order <- TRUE
#
#   # Validate n_folds
#   if (is.null(n_folds) || !is.numeric(n_folds) || n_folds <= 0) {
#     n_folds <- 5L
#     warning("n_folds was NULL, non-numeric, or less than or equal to 0; set to default: 5.")
#   } else {
#     n_folds <- as.integer(n_folds)
#   }
#   # Generate a list of fold indices
#   folds <- split(sample(1:n_obs), rep(1:n_folds, len = n_obs))
#
#   # Extract parameters
#   int.knots_init_grid <- param_gridCheck(parameters, "int.knots_init_grid", c(0, 1, 2))
#   shrinkage_grid <- param_gridCheck(parameters, "shrinkage_grid", c(0.8, 0.9, 1))
#   beta_grid <- param_gridCheck(parameters, "beta_grid", c(0.3, 0.5, 0.7))
#   phi_grid <- param_gridCheck(parameters, "phi_grid", c(0.9, 0.95, 0.99))
#   q_grid <- param_gridCheck(parameters, "q_grid", c(2, 3))
#
#   # Initialize a dataframe to store hyperparameters and corresponding MSE values
#   results_df <- data.frame()
#
#   # Loop through each combination
#   for (i in 1:length(int.knots_init_grid)) {
#     for (j in 1:length(shrinkage_grid)) {
#       for (k in 1:length(beta_grid)) {
#         for (l in 1:length(phi_grid)) {
#           for (m in 1:length(q_grid)) {
#             # Use tryCatch to handle errors
#             tryCatch({
#               # Initialize a vector to store MSE values for each fold
#               mse_fold <- numeric(n_folds)
#               n_iterations_fold <- numeric(n_folds)
#               avg_intknots_fold <- list(avg_intknots_univbl = numeric(n_folds),
#                                         avg_intknots_bivblX = numeric(n_folds),
#                                         avg_intknots_bivblY = numeric(n_folds))
#               time_fold <- numeric(n_folds)
#               for (fold_idx in 1:n_folds) {
#                 # Create the training and test sets for the current fold
#                 train_data <- data[-folds[[fold_idx]], ]
#                 test_data <- data[folds[[fold_idx]], ]
#                 rownames(train_data) <- rownames(test_data) <- NULL
#                 # Fit the model
#                 start_time <- Sys.time()
#                 error <- FALSE
#                 fit <- tryCatch(
#                   model_fun(formula = formula, data = train_data, int.knots_init = int.knots_init_grid[i],
#                             shrinkage = shrinkage_grid[j], beta = beta_grid[k], phi = phi_grid[l],
#                             q = q_grid[m], higher_order = higher_order),
#                   error = function(e) {
#                     message(paste0("Error occurred for fold ", fold_idx, ": ", e))
#                     error <<- TRUE
#                   }
#                 )
#                 if (error) next
#                 end_time <- Sys.time()
#                 time_fold[fold_idx] <- difftime(end_time, start_time, units = "secs")
#                 # Evaluate the model on the test set
#                 test_pred <- predict.GeDSboost(fit, newdata = test_data, n = n)
#
#                 # Calculate the mean squared error for the test set
#                 mse_fold[fold_idx] <- mean((test_data[[outcome]] - test_pred)^2)
#
#                 # Number of boosting iterations and average number of internal knots per boosting iteration
#                 aux <- n_intknots(fit)
#                 n_iterations_fold[fold_idx] <- aux$n_iterations
#                 avg_intknots_fold$avg_intknots_univbl[fold_idx] <- aux$avg_intknots$avg_intknots_univbl
#                 avg_intknots_fold$avg_intknots_bivblX[fold_idx] <- aux$avg_intknots$avg_intknots_bivblX
#                 avg_intknots_fold$avg_intknots_bivblY[fold_idx] <- aux$avg_intknots$avg_intknots_bivblY
#               }
#               # Store hyperparameters and corresponding MSE value in results_df
#               results_df <- rbind(results_df, data.frame(int.knots_init = int.knots_init_grid[i],
#                                                          shrinkage = shrinkage_grid[j],
#                                                          beta = beta_grid[k],
#                                                          phi = phi_grid[l],
#                                                          q = q_grid[m],
#                                                          mse = mean(mse_fold),
#                                                          n_iterations = mean(n_iterations_fold),
#                                                          avg_intknots_univbl = mean(avg_intknots_fold$avg_intknots_univbl),
#                                                          avg_intknots_bivblX = mean(avg_intknots_fold$avg_intknots_bivblX),
#                                                          avg_intknots_bivblY = mean(avg_intknots_fold$avg_intknots_bivblY),
#                                                          avg_fitting_time = mean(time_fold)
#               ))
#             },
#             # Specify what to do if an error occurs
#             error = function(e) {
#               warning(paste("Error for int.knots_init=", int.knots_init_grid[i],
#                             ", shrinkage=", shrinkage_grid[j], ", beta=", beta_grid[k],
#                             ", phi=", phi_grid[l], ", q=", q_grid[m], ". Iteration was skipped."))
#               return()
#             })
#           }
#         }
#       }
#     }
#   }
#
#   # Find the optimal internal_knots, max_iterations and shrinkage values that minimize the cross-validated mse
#    results_df <- as.data.frame(results_matrix)
#    best_params <- results_df[
#      which.min(results_df$mse),
#      ,
#      drop = FALSE
#      ]
#
#   # Return the optimal internal_knots, max_iterations and shrinkage values
#   return(list(best_params = best_params, results = results_df))
# }


##############
# n_intknots #
##############
n_intknots <- function(object, avg = TRUE){
  # Compute the # of boosting iterations
  if (inherits(object, "GeDSboost")) {
    n_iterations <- length(object$models) - if (object$args$initial_learner) 0L else 1L
    } else if (inherits(object, "GeDSgam")) {
      n_iterations <- sum(object$iters$backfitting)
    }
  n_intknots <- n_intknotsX <- n_intknotsY <- 0

  if (n_iterations == 0) {
    n_intknots_univbl <- n_intknots_bivblX <- n_intknots_bivblY <- 0
  } else {
    # Compute the average # of int. knots per boosting iteration
    base_learners <- object$final_model$base_learners
    for (bl_name in names(base_learners)) {
      n_vars <- length(object$args$base_learners[[bl_name]]$variables)
      if (n_vars == 1) {
        n_intknots <- length(base_learners[[bl_name]]$linear.int.knots)
      } else if (n_vars == 2) {
        n_intknotsX <- length(base_learners[[bl_name]]$linear.int.knots$ikX)
        n_intknotsY <- length(base_learners[[bl_name]]$linear.int.knots$ikY)
      }
    }

    if(avg) {
      n_intknots_univbl <- n_intknots / n_iterations
      n_intknots_bivblX <- n_intknotsX / n_iterations
      n_intknots_bivblY <- n_intknotsY / n_iterations
      } else {
        n_intknots_univbl <- n_intknots
        n_intknots_bivblX <- n_intknotsX
        n_intknots_bivblY <- n_intknotsY
      }

  }

  n_intknots <- list(univbl = n_intknots_univbl,
                     bivblX = n_intknots_bivblX,
                     bivblY = n_intknots_bivblY)



  return(list(n_intknots = n_intknots, n_iterations = n_iterations))
}

prepare_cv_grids <- function(parameters, allowed) {
  if (!is.list(parameters) ||
      (length(parameters) > 0L &&
       (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!nzchar(names(parameters))) || anyDuplicated(names(parameters))))) {
    stop("'parameters' must be a named list of numeric grids.")
  }
  unknown <- setdiff(names(parameters), allowed)
  if (length(unknown)) stop("Unknown parameter grid: ", paste(unknown, collapse = ", "))
  parameters <- Filter(Negate(is.null), parameters)
  for (name in names(parameters)) {
    values <- parameters[[name]]
    if (!is.numeric(values) || !length(values) || any(!is.finite(values))) {
      stop("'", name, "' must contain finite numeric values.")
    }
  }
  values <- setNames(lapply(allowed, function(name) {
    if (name %in% names(parameters)) parameters[[name]] else NA_real_
  }), allowed)
  list(values = values, tuned = sub("_grid$", "", names(parameters)))
}


