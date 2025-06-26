#########################################################################################
#########################################################################################
#################################### Cross-validation ###################################
#########################################################################################
#########################################################################################
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
#' @param parameters A set of parameters to be tuned via cross-validation.
#' These are: \code{beta}, \code{phi} and \code{q} in the case of \code{NGeDS},
#' \code{GGeDS} and \code{NGeDSgam}. In addition, for \code{NGeDSboost},
#' \code{int.knots_init} and \code{shrinkage} can also be tuned. Default values
#' are:
#' \itemize{
#' \item \code{int.knots_init_grid = c(0, 1, 2)},
#' \item \code{shrinkage_grid = c(0.1, 0.5, 1)},
#' \item \code{beta_grid = c(0.3, 0.5, 0.7)},
#' \item \code{phi_grid = c(0.9, 0.95, 0.99)} and
#' \item \code{q_grid = c(2, 3))}.
#' }
#'
#' @return Two data frames, \code{best_params} and \code{results}.
#' \code{best_params} contains the best combination of parameters according to
#' the cross-validated MSE. \code{results} presents the cross-validated MSE and 
#' the average number of internal knots across the folds for each possible 
#' combination of parameters, given the input \code{parameters}. In the case of
#' \code{model_fun = NGeDSboost}, it also provides the cross-validated number of
#' boosting iterations.
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

# Declare loop variables
globalVariables(c("i", "j", "k", "l", "m"))

# Define a function to perform k-fold cross-validation
crossv_GeDS <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                        parameters, ...)
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
                              n = n, n_folds = n_folds, parameters = parameters)
    } else if (model_fun_name %in% c("NGeDSboost")) {
      cv <- cross_validate.GeDSboost(formula = formula, data = data, model_fun = model_fun,
                                     n = n, n_folds = n_folds, parameters = parameters)
    }
  
  return(list(best_params = cv$best_params, results = as.data.frame(cv$results)))
}

################################################################################
############################ NGeDS, GGeDS, NGeDSgam ############################
################################################################################
cross_validate.GeDS <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                                parameters = list(beta_grid = c(0.3, 0.5, 0.7),
                                                  phi_grid = c(0.9, 0.95, 0.99),
                                                  q_grid = c(2, 3)), ...)
{
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

  # If cv model is linear, avoid computing higher order fits
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

  # Extract parameters
  beta_grid <- param_gridCheck(parameters, "beta_grid", c(0.3, 0.5, 0.7))
  phi_grid <- param_gridCheck(parameters, "phi_grid", c(0.9, 0.95, 0.99))
  q_grid <- param_gridCheck(parameters, "q_grid", c(2, 3))



  # Initialize a matrix to store parameters and corresponding MSE values
  results_matrix <- matrix(nrow = length(beta_grid) * length(phi_grid) * length(q_grid),
                           ncol = 11)

  # Register parallel backend
  n_cores <- max(1, detectCores() - 1) # Avoid using all cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Parallel processing
  results_list <- foreach(i = 1:length(beta_grid), .combine = rbind, .packages = "GeDS") %:%
    foreach(j = 1:length(phi_grid), .combine = rbind) %:%
    foreach(k = 1:length(q_grid), .combine = rbind) %dopar% {

      # Initialize a vector to store MSE values for each fold
      mse_fold <- numeric(n_folds)
      intknots_fold <- list(univbl = numeric(n_folds),
                            bivblX = numeric(n_folds),
                            bivblY = numeric(n_folds))
      time_fold <- numeric(n_folds)

      for (fold_idx in 1:n_folds) {
        # Create the training and test sets for the current fold
        train_data <- data[-folds[[fold_idx]], ]
        test_data <- data[folds[[fold_idx]], ]
        rownames(train_data) <- rownames(test_data) <- NULL
        # Fit the model
        start_time <- Sys.time()
        error <- FALSE
        fit <- tryCatch(
          model_fun(formula = formula, data = train_data, beta = beta_grid[i], phi = phi_grid[j],
                    q = q_grid[k], higher_order = higher_order),
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

        # Number of internal knots
        if (inherits(fit, "GeDS")) {
          if (fit$Type == "LM - Univ" || fit$Type == "GLM - Univ") {
            intknots_fold$univbl[fold_idx] <- fit$Nintknots
          } else if (fit$Type == "LM - Biv" || fit$Type == "GLM - Biv") {
            intknots_fold$bivblX[fold_idx] <- fit$Nintknots$X
            intknots_fold$bivblY[fold_idx] <- fit$Nintknots$Y
          }
        }

      }
      
      if (inherits(fit, "GeDS")) {
        if (fit$Type == "LM - Univ" || fit$Type == "GLM - Univ") {
          data.frame(beta = beta_grid[i],
                     phi = phi_grid[j],
                     q = q_grid[k],
                     cv_mse = mean(mse_fold),
                     avg_intknots = mean(intknots_fold$univbl),
                     avg_fitting_time = mean(time_fold))
          
          } else if (fit$Type == "LM - Biv" || fit$Type == "GLM - Biv")
            data.frame(beta = beta_grid[i],
                       phi = phi_grid[j],
                       q = q_grid[k],
                       cv_mse = mean(mse_fold),
                       avg_intknots_X = mean(intknots_fold$bivblX),
                       avg_intknots_Y = mean(intknots_fold$bivblY),
                       avg_fitting_time = mean(time_fold))
        } else {
          # NGeDSgam
          data.frame(beta = beta_grid[i],
                     phi = phi_grid[j],
                     q = q_grid[k],
                     cv_mse = mean(mse_fold),
                     avg_fitting_time = mean(time_fold))
        }
      
    }

  # Stop the parallel backend
  stopCluster(cl)

  # Convert results to a dataframe
  results_matrix <- do.call(cbind, results_list)

  # Find the optimal parameters that minimize the cross-validated mse
  best_params <- results_matrix[which.min(results_matrix[,"cv_mse"]), ]

  # Return the optimal parameters
  return(list(best_params = best_params, results = results_matrix))
}

################################################################################
################################## NGeDSboost ##################################
################################################################################
cross_validate.GeDSboost <- function(formula, data, model_fun, n = 2L, n_folds = 5L,
                                     parameters = list(int.knots_init_grid = c(0, 1, 2),
                                                       shrinkage_grid = c(0.1, 0.5, 1),
                                                       beta_grid = c(0.3, 0.5, 0.7),
                                                       phi_grid = c(0.9, 0.95, 0.99),
                                                       q_grid = c(2, 3)), ...)
{
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

  # Extract parameters
  int.knots_init_grid <- param_gridCheck(parameters, "int.knots_init_grid", c(0, 1, 2))
  shrinkage_grid <- param_gridCheck(parameters, "shrinkage_grid", c(0.8, 0.9, 1))
  beta_grid <- param_gridCheck(parameters, "beta_grid", c(0.3, 0.5, 0.7))
  phi_grid <- param_gridCheck(parameters, "phi_grid", c(0.9, 0.95, 0.99))
  q_grid <- param_gridCheck(parameters, "q_grid", c(2, 3))

  # Initialize a matrix to store parameters and corresponding MSE values
  results_matrix <- matrix(nrow = length(int.knots_init_grid) * length(shrinkage_grid) * length(beta_grid) * length(phi_grid) * length(q_grid),
                           ncol = 11)

  # Register parallel backend
  n_cores <- detectCores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Parallel processing
  results_list <- foreach(i = 1:length(int.knots_init_grid), .combine = rbind, .packages = "GeDS",
                          .export = c("n_intknots")) %:%
    foreach(j = 1:length(shrinkage_grid), .combine = rbind) %:%
    foreach(k = 1:length(beta_grid), .combine = rbind) %:%
    foreach(l = 1:length(phi_grid), .combine = rbind) %:%
    foreach(m = 1:length(q_grid), .combine = rbind) %dopar% {

      # Initialize a vector to store MSE values for each fold
      mse_fold <- numeric(n_folds)
      n_iterations_fold <- numeric(n_folds)
      intknots_fold <- list(univbl = numeric(n_folds),
                            bivblX = numeric(n_folds),
                            bivblY = numeric(n_folds))
      time_fold <- numeric(n_folds)

      for (fold_idx in 1:n_folds) {
        # Create the training and test sets for the current fold
        train_data <- data[-folds[[fold_idx]], ]
        test_data <- data[folds[[fold_idx]], ]
        rownames(train_data) <- rownames(test_data) <- NULL
        # Fit the model
        start_time <- Sys.time()
        error <- FALSE
        fit <- tryCatch(
          model_fun(formula = formula, data = train_data, int.knots_init = int.knots_init_grid[i],
                    shrinkage = shrinkage_grid[j], beta = beta_grid[k], phi = phi_grid[l],
                    q = q_grid[m], higher_order = higher_order, ...),
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

      data.frame(int.knots_init = int.knots_init_grid[i],
                 shrinkage = shrinkage_grid[j],
                 beta = beta_grid[k],
                 phi = phi_grid[l],
                 q = q_grid[m],
                 cv_mse = mean(mse_fold),
                 n_iterations = mean(n_iterations_fold),
                 avg_intknots_univbl = mean(intknots_fold$univbl),
                 avg_intknots_bivblX = mean(intknots_fold$bivblX),
                 avg_intknots_bivblY = mean(intknots_fold$bivblY),
                 avg_fitting_time = mean(time_fold))
    }

  # Stop the parallel backend
  stopCluster(cl)

  # Convert results to a dataframe
  results_matrix <- do.call(cbind, results_list)

  # Find the optimal internal_knots, max_iterations and shrinkage values that minimize the cross-validated mse
  best_params <- results_matrix[which.min(results_matrix[,"cv_mse"]), ]

  # Return the optimal internal_knots, max_iterations and shrinkage values
  return(list(best_params = best_params, results = results_matrix))
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
#               if (fit$Type == "LM - Univ" || fit$Type == "GLM - Univ") {
#                 avg_intknots_fold$avg_intknots_univbl[fold_idx] <- fit$Nintknots
#                 avg_intknots_fold$avg_intknots_bivblX[fold_idx] <- NA
#                 avg_intknots_fold$avg_intknots_bivblY[fold_idx] <- NA
#                 } else if (fit$Type == "LM - Biv" || fit$Type == "GLM - Biv") {
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
#   best_params <- results_df[which.min(results_df$mse), ]
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
#   best_params <- results_df[which.min(results_df$mse), ]
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
    n_iterations <- if (object$args$initial_learner) object$iters + 1
    else object$iters
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

###################
# param_gridCheck #
###################
param_gridCheck <- function(param_grid, param_name, default_value) {
  # Check if the parameter is missing or null
  if (is.null(param_grid[[param_name]]) || !param_name %in% names(param_grid)) {
    warning(paste0(param_name, " was NULL or missing; set to default."))
    param_grid[[param_name]] <- default_value
  } else {
    # Check if the parameter is numeric
    if (!is.numeric(param_grid[[param_name]])) {
      warning(paste0(param_name, " is not numeric; set to default."))
      param_grid[[param_name]] <- default_value
    }
  }
  return(param_grid[[param_name]])
}


