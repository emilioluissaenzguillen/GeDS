################################################################################
################################################################################
################################## NGeDSboost ##################################
################################################################################
################################################################################
#' @title Component-wise gradient boosting with NGeDS base-learners
#' @name NGeDSboost
#' @description
#' \code{NGeDSboost} performs component-wise gradient boosting (Bühlmann and Yu
#' (2003), Bühlmann and Hothorn (2007)) using normal GeD splines (i.e., fitted
#' with \code{\link{NGeDS}} function) as base-learners (see Dimitrova et al. (2024)).
#' @param formula a description of the structure of the model to be fitted,
#' including the dependent and independent variables. Unlike \code{\link{NGeDS}}
#' and \code{\link{GGeDS}}, the formula specified allows for multiple additive
#' GeD spline regression components (as well as linear components) to be
#' included (e.g., \code{Y ~ f(X1) + f(X2) + X3}).
#' See \code{\link[=formula.GeDS]{formula}} for further details.
#' @param data a data frame containing the variables referenced in the formula.
#' @param weights an optional vector of `prior weights' to be put on the
#' observations during the fitting process. It should be \code{NULL} or a
#' numeric vector of the same length as the response variable defined in the
#' formula.
#' @param normalize_data a logical that defines whether the data should be
#' normalized (standardized) before fitting the baseline linear model, i.e.,
#' before running the FGB algorithm. Normalizing the data involves scaling the
#' predictor variables to have a mean of 0 and a standard deviation of 1. This
#' process alters the scale and interpretation of the knots and coefficients
#' estimated. Default is equal to \code{FALSE}.
#' @param family determines the loss function to be optimized by the boosting
#' algorithm. In case \code{initial_learner = FALSE} it also determines the
#' corresponding empirical risk minimizer to be used as offset initial learner.
#' By default, it is set to \code{mboost::Gaussian()}. Users can specify any
#' \code{\link[mboost]{Family}} object from the \pkg{mboost} package.
#' @param initial_learner a logical value. If set to \code{TRUE}, the model's
#' initial learner will be a normal GeD spline. If set to FALSE, then the
#' initial predictor will consist of the empirical risk minimizer corresponding
#' to the specified family. Note that if \code{initial_learner = TRUE},
#' \code{family} must be \code{mboost::Gaussian()}.
#' @param int.knots_init optional parameter allowing the user to set a
#' maximum number of internal knots to be added by the initial GeDS learner in
#' case \code{initial_learner = TRUE}. Default is equal to \code{2L}.
#' @param min_iterations optional parameter to manually set a minimum number of
#' boosting iterations to be run. If not specified, it defaults to 0L.
#' @param max_iterations optional parameter to manually set the maximum number
#' of boosting iterations to be run. If not specified, it defaults to 100L.
#' This setting serves as a fallback when the stopping rule, based on
#' consecutive deviances and tuned by \code{phi_boost_exit} and \code{q_boost},
#' does not trigger an earlier termination (see Dimitrova et al. (2024)).
#' Therefore, users can increase/decrease the number of boosting iterations,
#' by increasing/decreasing the value \code{phi_boost_exit} and/or#
#' \code{q_boost}, or directly specify \code{max_iterations}.
#' @param shrinkage numeric parameter in the interval \eqn{[0,1]} defining the
#' step size or shrinkage parameter. This controls the size of the steps taken
#' in the direction of the gradient of the loss function. In other words, the
#' magnitude of the update each new iteration contributes to the final model.
#' Default is equal to \code{1}.
#' @param phi_boost_exit numeric parameter in the interval \eqn{[0,1]}
#' specifying the threshold for the boosting iterations stopping rule. Default
#' is equal to \code{0.995}.
#' @param q_boost numeric parameter which allows to fine-tune the boosting
#' iterations stopping rule, by default equal to \code{2L}.
#' @param beta numeric parameter in the interval \eqn{[0,1]} tuning the knot
#' placement in stage A of GeDS. Default is equal to \code{0.5}. See details in
#' \code{\link{NGeDS}}.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the
#' threshold for the stopping rule  (model selector) in stage A of GeDS.
#' Default is equal to \code{0.99}. See details in \code{\link{NGeDS}}.
#' @param internal_knots The maximum number of internal knots that can be added
#' by the GeDS base-learners in each boosting iteration, effectively setting the
#' value of \code{max.intknots} in \code{\link{NGeDS}} at each boosting
#' iteration. Default is \code{500L}.
#' @param q numeric parameter which allows to fine-tune the stopping rule of
#' stage A of GeDS, by default equal to \code{2L}. See details in
#' \code{\link{NGeDS}}.
#' @param higher_order a logical that defines whether to compute the higher
#' order fits (quadratic and cubic) after the FGB algorithm is run. Default is
#' \code{TRUE}.
#' 
#' @return \code{\link{GeDSboost-Class}} object, i.e. a list of items that
#' summarizes the main details of the fitted FGB-GeDS model. See
#' \code{\link{GeDSboost-Class}} for details. Some S3 methods are available in
#' order to make these objects tractable, such as
#' \code{\link[=coef.GeDSboost]{coef}}, \code{\link[=knots.GeDSboost]{knots}},
#' \code{\link[=print.GeDSboost]{print}} and
#' \code{\link[=predict.GeDSboost]{predict}}. Also variable importance measures
#' (\code{\link[=bl_imp.GeDSboost]{bl_imp}}) and improved plotting facilities
#' (\code{\link[=visualize_boosting.GeDSboost]{visualize_boosting}}).
#' 
#' @details
#' The  \code{NGeDSboost} function implements functional gradient boosting
#' algorithm for some pre-defined loss function, using linear GeD splines as
#' base learners. At each boosting iteration, the negative gradient vector is
#' fitted through the base procedure encapsulated within the \code{\link{NGeDS}}
#' function. The latter constructs a Geometrically Designed variable knots
#' spline regression model for a response having a Normal distribution. The FGB
#' algorithm yields a final linear fit. Higher order fits (quadratic and cubic)
#' are then computed by calculating the Schoenberg’s variation diminishing
#' spline (VDS) approximation of the linear fit.
#' 
#' On the one hand, \code{NGeDSboost} includes all the parameters of
#' \code{\link{NGeDS}}, which in this case tune the base-learner fit at each
#' boosting iteration. On the other hand, \code{NGeDSboost} includes some
#' additional parameters proper to the FGB procedure. We describe the main ones
#' as follows. 
#' 
#' First, \code{family} allows to specify the loss function and corresponding
#' risk function to be optimized by the boosting algorithm. If
#' \code{initial_learner = FALSE}, the initial learner employed will be the
#' empirical risk minimizer corresponding to the family chosen. If
#' \code{initial_learner = TRUE} then the initial learner will be an
#' \code{\link{NGeDS}} fit with maximum number of internal knots equal to
#' \code{int.knots_init}.
#' 
#' \code{shrinkage} tunes the step length/shrinkage parameter which helps to 
#' control the learning rate of the model. In other words, when a new base
#' learner is added to the ensemble, its contribution to the final prediction is
#' multiplied by the shrinkage parameter. The smaller \code{shrinkage} is, the
#' slower/more gradual the learning process will be, and viceversa.
#' 
#' The number of boosting iterations is controlled by a
#' \emph{Ratio of Deviances} stopping rule similar to the one presented for
#' \code{\link{GGeDS}}. In the same way \code{phi} and \code{q} tune the
#' stopping rule of \code{\link{GGeDS}}, \code{phi_boost_exit} and
#' \code{q_boost} tune the stopping rule of \code{NGeDSboost}. The user can also
#' manually control the number of boosting iterations through
#' \code{min_iterations} and \code{max_iterations}.
#' 
#' @examples
#' 
#' ################################# Example 1 #################################
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
#' Y <- rnorm(N, means, sd = 0.2)
#' data = data.frame(X, Y)
#'
#' # Fit a Normal FGB-GeDS regression using NGeDSboost
#'
#' Gmodboost <- NGeDSboost(Y ~ f(X), data = data)
#' MSE_Gmodboost_linear <- mean((sapply(X, f_1) - Gmodboost$predictions$pred_linear)^2)
#' MSE_Gmodboost_quadratic <- mean((sapply(X, f_1) - Gmodboost$predictions$pred_quadratic)^2)
#' MSE_Gmodboost_cubic <- mean((sapply(X, f_1) - Gmodboost$predictions$pred_cubic)^2)
#'
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#'     "Linear NGeDSboost:", MSE_Gmodboost_linear, "\n",
#'     "Quadratic NGeDSboost:", MSE_Gmodboost_quadratic, "\n",
#'     "Cubic NGeDSboost:", MSE_Gmodboost_cubic, "\n")
#' 
#' # Compute predictions on new randomly generated data
#' X <- sort(runif(100, min = -2, max = 2))
#' 
#' pred_linear <- predict(Gmodboost, newdata = data.frame(X), n = 2L)
#' pred_quadratic <- predict(Gmodboost, newdata = data.frame(X), n = 3L)
#' pred_cubic <- predict(Gmodboost, newdata = data.frame(X), n = 4L)
#' 
#' MSE_Gmodboost_linear <- mean((sapply(X, f_1) - pred_linear)^2)
#' MSE_Gmodboost_quadratic <- mean((sapply(X, f_1) - pred_quadratic)^2)
#' MSE_Gmodboost_cubic <- mean((sapply(X, f_1) - pred_cubic)^2)
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#'     "Linear NGeDSboost:", MSE_Gmodboost_linear, "\n",
#'     "Quadratic NGeDSboost:", MSE_Gmodboost_quadratic, "\n",
#'     "Cubic NGeDSboost:", MSE_Gmodboost_cubic, "\n")
#' 
#' ## S3 methods for class 'GeDSboost'
#' # Print 
#' print(Gmodboost)
#' # Knots
#' knots(Gmodboost, n = 2L)
#' knots(Gmodboost, n = 3L)
#' knots(Gmodboost, n = 4L)
#' # Coefficients
#' coef(Gmodboost, n = 2L)
#' coef(Gmodboost, n = 3L)
#' coef(Gmodboost, n = 4L)
#' # Deviances
#' deviance(Gmodboost, n = 2L)
#' deviance(Gmodboost, n = 3L)
#' deviance(Gmodboost, n = 4L)
#' 
#' ############################ Example 2 - Bodyfat ############################
#' library(TH.data)
#' data("bodyfat", package = "TH.data")
#' 
#' Gmodboost <- NGeDSboost(formula = DEXfat ~ age + f(hipcirc, waistcirc) + f(kneebreadth),
#' data = bodyfat)
#' 
#' MSE_Gmodboost_linear <- mean((bodyfat$DEXfat - Gmodboost$predictions$pred_linear)^2)
#' MSE_Gmodboost_quadratic <- mean((bodyfat$DEXfat - Gmodboost$predictions$pred_quadratic)^2)
#' MSE_Gmodboost_cubic <- mean((bodyfat$DEXfat - Gmodboost$predictions$pred_cubic)^2)
#' # Comparison
#' cat("\n", "MSE", "\n",
#'     "Linear NGeDSboost:", MSE_Gmodboost_linear, "\n",
#'     "Quadratic NGeDSboost:", MSE_Gmodboost_quadratic, "\n",
#'     "Cubic NGeDSboost:", MSE_Gmodboost_cubic, "\n")
#'
#' @seealso \code{\link{NGeDS}}; \code{\link{GGeDS}}; \code{\link{GeDSboost-Class}};
#' S3 methods such as \code{\link{knots.GeDSboost}}; \code{\link{coef.GeDSboost}};
#' \code{\link{deviance.GeDSboost}}; \code{\link{predict.GeDSboost}}
#'      
#' @export
#' @import foreach
#' @import doParallel
#' @import doFuture
#' @import future
#' @import doRNG
#' @import TH.data
#' @importFrom parallel detectCores
#' @importFrom stats setNames
#' 
#' @references 
#' Friedman, J.H. (2001).
#' Greedy function approximation: A gradient boosting machine.
#' \emph{The Annals of Statistics}, \strong{29 (5)}, 1189--1232. \cr
#' DOI: \doi{10.1214/aos/1013203451}
#' 
#' Bühlmann P., Yu B. (2003).
#' Boosting With the L2 Loss.
#' \emph{Journal of the American Statistical Association},
#' \strong{98(462)}, 324–339.
#' \doi{10.1198/016214503000125}
#' 
#' Bühlmann P., Hothorn T. (2007).
#' Boosting Algorithms: Regularization, Prediction and Model Fitting.
#' \emph{Statistical Science}, \strong{22(4)}, 477 – 505. \cr
#' DOI: \doi{10.1214/07-STS242}
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
#' Dimitrova, D. S., Guillen, E. S. and Kaishev, V. K.  (2024).
#' \pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
#' Models and Functional Gradient Boosting, based on Geometrically Designed
#' (GeD) Splines. \emph{Manuscript submitted for publication.}

##################
### NGeDSboost ###
##################
NGeDSboost <- function(formula, data, weights = NULL, normalize_data = FALSE,
                       family = mboost::Gaussian(), initial_learner = TRUE,
                       int.knots_init = 2L, min_iterations,
                       max_iterations, shrinkage = 1,
                       phi_boost_exit = 0.995, q_boost = 2L, beta = 0.5,
                       phi = 0.99, internal_knots = 500L, q = 2L,
                       higher_order = TRUE)
  {
  # Capture the function call
  extcall <- match.call()
  
  # Models list
  models <- list()
  
  # Formula
  read.formula <- read.formula.boost(formula, data)
  response <- read.formula$response
  predictors <- read.formula$predictors
  base_learners <- read.formula$base_learners
  
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
    warning(paste("variable", var, "contains missing values;\n",
                  "these observations are not used for fitting"))
  }
  # Eliminate rows containing NA values
  data <- na.omit(data)
  
  # Family arguments
  ngradient <- family@ngradient
  risk <- family@risk
  offset <- family@offset
  check_y_family <- family@check_y
  
  # Min/max iterations
  min_iterations <- validate_iterations(min_iterations, 0L, "min_iterations")
  max_iterations <- validate_iterations(max_iterations, 100L, "max_iterations")
  
  # Save arguments
  args <- list(
    "predictors" = data[predictors], 
    "base_learners"= base_learners,
    "family" = family,
    "initial_learner" = initial_learner, 
    "int.knots_init" = if(!is.null(initial_learner)){int.knots_init} else {NULL},
    "shrinkage" = shrinkage, 
    "normalize_data" = normalize_data
  )
  
  # Response variable check
  data[response] <- check_y_family(data[[response]])
  # Add 'response' as the first element of the list 'args'
  args <- c(list(response = data.frame(data[[response]])), args)
  names(args$response) <- response
  
  # Normalize data if necessary
  if (normalize_data == TRUE) {
    if (family@name != "Negative Binomial Likelihood (logit link)") {
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
        # If the family is "Negative Binomial Likelihood (logit link)" only normalize predictors
        args$X_mean <- colMeans(data[predictors])
        args$X_sd <- sapply(data[predictors], sd)
        data[predictors] <- scale(data[predictors])
    }
  }
  
  # Weights
  nobs = length(data[[response]])
  if (is.null(weights)) weights <- rep.int(1, nobs)
  else weights <- rescale_weights(weights)
  if (!family@weights(weights))
    stop(sQuote("family"), " is not able to deal with weights")

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
        intervals <- matrix(c(min_vals[[1]], max_vals[[1]]), ncol = 2)
        colnames(intervals) <- c("start", "end")
        coefficients <- list("b0" = 0, "b1" = 0)
        base_learners_list[[bl_name]] <- list("knots" = knots, "intervals" = intervals, "coefficients" = coefficients)
      # (A.2) BIVARIATE BASE-LEARNERS
        } else {
          knots <- list(Xk = c(min_vals[1], max_vals[1]), Yk = c(min_vals[2], max_vals[2]))
          base_learners_list[[bl_name]] <- list("knots" = knots, "iterations" = list())
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
  
  #####################
  ## Initial Learner ##
  #####################
  
  ## (I) GeDS (or Linear) Initial Learner ##
  if (initial_learner) {
    if (family@name != "Squared Error (Regression)") {
      stop("Initial learner only allowed for family = Gaussian()")
    }
    
    # Pre-allocate memory for results
    results <- vector("list", length = length(base_learners))
    names(results) <- names(base_learners)
    best_pred <- best_bl <- NULL; best_resid <- best_ssr <- Inf
    # Template for model formula
    model_formula_template <- paste0(response, " ~ ")
    
    ## 1. Loop through predictors and fit an initial base learner to data (NGeDS with # int.knots_init)
    model_results <- lapply(names(base_learners), function(bl_name) {
      componentwise_fit(bl_name, data = data, response = response, model_formula_template, weights, base_learners, m = 0,
                        internal_knots = int.knots_init, beta, phi, q)
      })
    
    ## 2. Calculate SSR and find model that fits best U according to SSR
    ssr_values <- numeric(length(model_results))
    for (i in seq_along(model_results)) {
      if (!model_results[[i]]$error) {
        ssr_values[i] <- model_results[[i]]$ssr
        } else {
          message(model_results[[i]]$message)
          ssr_values[i] <- Inf
        }
    }
    # Finding the best model
    min_ssr_index <- which.min(ssr_values)
    best_resid    <- model_results[[min_ssr_index]]$resid 
    best_ssr      <- model_results[[min_ssr_index]]$ssr
    best_pred     <- model_results[[min_ssr_index]]$pred
    best_bl       <- model_results[[min_ssr_index]]$bl_name
    
    ## (A) GeDS best base-learner
    if (base_learners[[best_bl]]$type == "GeDS") {
      
      # (A.1) UNIVARIATE BASE-LEARNERS
      if (length(base_learners[[best_bl]]$variables) == 1) {
        ## 3. GeDS fit into linear piecewise polynomial form
        pred_linear <- best_pred$Y_hat
        int.knt <- best_pred$int.knt
        b0 <- best_pred$b0
        b1 <- best_pred$b1
        coef <- list("b0" = b0, "b1" = b1, "mat" = best_pred$mat)
        ## 4. Knots, intervals and coefficients
        # (i) Knots
        knots = sort(c(base_learners_list[[best_bl]]$knots, int.knt))
        int.knots = knots[-c(1, length(knots))]
        # (ii) Intervals
        # Combine the 1st value from the "start" col of intervals, with the values from the "end" col of intervals and int.knt
        combined <- sort(c(unname(base_learners_list[[best_bl]]$intervals[1, "start"]),
                           unname(base_learners_list[[best_bl]]$intervals[, "end"]), int.knots))
        intervals <- matrix(c(combined[-length(combined)], combined[-1]), ncol = 2)
        colnames(intervals) <- c("start", "end")
        
        # Save the coefficients, knots and intervals for the best predictor
        base_learners_list[[best_bl]]$knots <- knots
        base_learners_list[[best_bl]]$intervals <- intervals
        base_learners_list[[best_bl]]$coefficients <- list("b0" = b0, "b1" = b1)
        
        # Number of int. knots placed on initial learner
        n_intknots <- length(knots) - 2
        
      # (A.2) BIVARIATE BASE-LEARNERS
        } else if (length(base_learners[[best_bl]]$variables) == 2){
          ## 3/4. Save knots and coefficients (bivariate base-learners are not converted into PP form)
          pred_linear <- best_pred$Y_hat
          int.knt <- best_pred$int.knt
          knots = list(Xk = sort(c(base_learners_list[[best_bl]]$knots$Xk, int.knt$Xint.knt)),
                       Yk = sort(c(base_learners_list[[best_bl]]$knots$Yk, int.knt$Yint.knt)))
          base_learners_list[[best_bl]]$knots <- knots
          # Number of int. knots placed on initial learner
          n_intknots <- list(X = length(knots$Xk) - 2, Y = length(knots$Yk) - 2)
          
          coef <-  best_pred$theta
          base_learners_list[[best_bl]]$iterations[["model0"]] <- list("int.knt" = int.knt, "coef" = coef)
        }
    
    ## (B) Linear best base-learner
      } else if (base_learners[[best_bl]]$type == "linear") {
        pred_linear <- best_pred$Y_hat
        n <- length(best_pred) - 1; coef <- vector("list", length = n)
        names(coef) <- paste0("b", 0:(n-1))
        coef[] <- best_pred[-1]
        
        int.knt <- n_intknots <- NULL
        
        base_learners_list[[best_bl]]$coefficients <-  coef
      }
    
    ## 3. Evaluate model and compute residuals
    Y_hat <- best_pred$Y_hat
    U     <- best_resid
    # Initialize old.dev for stopping rule
    old.dev <- best_ssr
    
    # Save model
    models[["model0"]] <- list("best_bl" = list("name" = best_bl, "n_intknots" = n_intknots, "int.knt" = int.knt,
                                                "coef" = coef, "pred_linear" = pred_linear),
                               "Y_hat" = Y_hat, "base_learners" = base_learners_list)
    
    
  ## (II) Offset Initial Learner ##
  } else {
    # 1. Initialize Y_hat with an offset value
    Y_hat <- offset(data[[response]], weights)
    if (length(Y_hat) == 1){Y_hat <- rep(Y_hat, NROW(data[[response]]))}
    old.dev <- risk(data[[response]], Y_hat, weights)
    # 2. Compute the negative gradient of the loss function
    U <- ngradient(data[[response]], Y_hat, weights)
    
    # Save model
    models[["model0"]] <- list("Y_hat" = Y_hat, "base_learners" = base_learners_list)
  }
  
  ## Initialize parallel processing if required (i.e. if # base-learners > 1000)
  pprocessing_threshold <- 1000
  if (length(base_learners) >= pprocessing_threshold) {
    # Number of cores
    n_cores <- detectCores() - 3 # Leave 3 cores free
    plan(multisession, workers = n_cores)
    registerDoFuture()
    # Number of predictors
    n_bl <- length(base_learners)
    # Number of predictors per batch
    bl_per_batch <- ceiling(n_bl / n_cores)
    # Split predictors into batches
    bl_batches <- split(names(base_learners), ceiling(seq_along(names(base_learners)) / bl_per_batch))
    # Call necessary internal functions
    componentwise_fit   <- componentwise_fit
    predict_GeDS_linear <- predict_GeDS_linear
    last_row_with_value <- last_row_with_value
  }
  
  #########################
  ## Boosting iterations ##
  #########################
  for (m in 1:max_iterations){
    
    previous_model <- models[[paste0("model", m-1)]]
    data_loop <- cbind(U = U, data[predictors])
    
    # Pre-allocate memory for results
    results <- vector("list", length = length(base_learners))
    names(results) <- names(base_learners)
    best_pred <- best_bl <- NULL; best_ssr <- Inf
    # Template for model formula
    model_formula_template <- "U ~ "
    
    if (length(base_learners) < pprocessing_threshold) {
      ## 3. Loop through base learners and fit to negative gradient (NGeDS with # internal_knots or linear model)
      model_results <- lapply(names(base_learners), function(bl_name) {
        componentwise_fit(bl_name, data = data_loop, response = "U", model_formula_template, weights, base_learners, m,
                          internal_knots = internal_knots, beta, phi, q)
      })
      
      ## 4. Calculate SSR and find model that fits best U according to SSR
      ssr_values <- numeric(length(model_results))
      for (i in seq_along(model_results)) {
        if (!model_results[[i]]$error) {
          ssr_values[i] <- model_results[[i]]$ssr
        } else {
          message(model_results[[i]]$message)
          ssr_values[i] <- Inf
        }
      }
      
      # Finding the best model
      min_ssr_index <- which.min(ssr_values)
      best_ssr      <- model_results[[min_ssr_index]]$ssr
      best_pred     <- model_results[[min_ssr_index]]$pred
      best_bl       <- model_results[[min_ssr_index]]$bl_name
      
    #########################
    ## PARALLEL PROCESSING ##
    #########################
      } else {
       ## 3. Loop through base-learners and fit to negative gradient (NGeDS with # internal_knots)
        results <- foreach(bl_name = names(base_learners), .combine = 'rbind', .packages = "GeDS",
                           .export = c("predict_GeDS_linear", "last_row_with_value")) %dorng% {
          componentwise_fit(bl_name, data = data_loop, response = "U", model_formula_template, weights, base_learners, m,
                            internal_knots = internal_knots, beta, phi, q)
        }
      ## 4. Calculate SSR and find model that fits best U according to SSR
        best_result_index <- which.min(sapply(1:nrow(results), function(i) results[i, "ssr"]))
        best_ssr <- unlist(results[best_result_index, "ssr"])
        best_pred <- results[best_result_index, "pred"][[1]]
        best_bl <- unlist(results[best_result_index, "bl_name"])
        
      }
    
    # If all base-learners were skipped, stop the iterations
    if (is.infinite(best_ssr)) {
      message(paste0("All predictors were skipped in iteration ", m, ". Stopping boosting iterations."))
      break
    }
    
    ## (A) GeDS best base-learner
    if (base_learners[[best_bl]]$type == "GeDS") {
      
      # (A.1) UNIVARIATE BASE-LEARNERS
      if (length(base_learners[[best_bl]]$variables) == 1){
        ## 4.1. GeDS fit into linear piecewise polynomial form
        pred_linear <- best_pred$Y_hat
        int.knt <- best_pred$int.knt
        b0 <- best_pred$b0
        b1 <- best_pred$b1
        coef <- list("b0" = b0, "b1" = b1, "mat" = best_pred$mat)
        
        ## 4.2. Knots, intervals and coefficients
        # (i) Update knots vector
        knots = sort(c(previous_model$base_learners[[best_bl]]$knots, int.knt))
        int.knots = knots[-c(1, length(knots))]
        
        # (ii) Update intervals matrix
        # Combine the 1st value from the "start" col of intervals, with the values from the "end" col of intervals and int.knt
        combined <- sort(c(unname(previous_model$base_learners[[best_bl]]$intervals[1, "start"]),
                           unname(previous_model$base_learners[[best_bl]]$intervals[, "end"]), int.knt))
        intervals <- matrix(c(combined[-length(combined)], combined[-1]), ncol = 2)
        colnames(intervals) <- c("start", "end")
        
        # (iii) Coefficients
        # Combine old intervals with old coefficients
        old_aux <- cbind(previous_model$base_learners[[best_bl]]$intervals,
                         data.frame(b0 = previous_model$base_learners[[best_bl]]$coefficients$b0),
                         data.frame(b1 = previous_model$base_learners[[best_bl]]$coefficients$b1))
        # Combine new intervals with coefficients from the piecewise polynomial
        new_aux <- data.frame(intervals)
        new_aux$interval <- findInterval(new_aux$end, c(-Inf, int.knt, Inf))
        new_aux$interval[new_aux$end %in% int.knt] <- new_aux$interval[new_aux$end %in% int.knt] - 1
        
        new_aux$b0 <- b0[new_aux$interval]
        new_aux$b1 <- b1[new_aux$interval]
        new_aux$interval <- NULL
        
        # Initialize new coefficients
        b0_new <- numeric(nrow(new_aux))
        b1_new <- numeric(nrow(new_aux))
        
        # Sum old coefficients and coefficients of piecewise polynomial
        # Using a matrix approach for faster computation
        mat_start <- matrix(old_aux$start, nrow = nrow(new_aux), ncol = nrow(old_aux), byrow = TRUE)
        mat_end <- matrix(old_aux$end, nrow = nrow(new_aux), ncol = nrow(old_aux), byrow = TRUE)
        
        # Find intervals that satisfy the conditions
        indices <- which(mat_start <= new_aux$start & new_aux$end <= mat_end, arr.ind = TRUE)
        
        # Update coefficients based on the identified intervals
        b0_new[indices[, 1]] <- old_aux$b0[indices[, 2]] + shrinkage * new_aux$b0[indices[, 1]]
        b1_new[indices[, 1]] <- old_aux$b1[indices[, 2]] + shrinkage * new_aux$b1[indices[, 1]]
        
        # Update the coefficients, knots, and intervals for the best predictor
        base_learners_list[[best_bl]]$knots <- knots
        base_learners_list[[best_bl]]$intervals <- intervals
        base_learners_list[[best_bl]]$coefficients <- list("b0" = b0_new, "b1" = b1_new)
        
        # Calculate number of internal knots placed on current iteration
        n_intknots <- length(knots) - length(previous_model$base_learners[[best_bl]]$knots)
        
      # (A.2) BIVARIATE BASE-LEARNERS
        } else if (length(base_learners[[best_bl]]$variables) == 2){
          ## 3/4. Save knots and coefficients (bivariate base-learners are not converted into PP form)
          pred_linear <- best_pred$Y_hat
          int.knt <- best_pred$int.knt
          knots = list(Xk = sort(c(base_learners_list[[best_bl]]$knots$Xk, int.knt$Xint.knt)),
                       Yk = sort(c(base_learners_list[[best_bl]]$knots$Yk, int.knt$Yint.knt)))
          base_learners_list[[best_bl]]$knots <- knots
          # Number of int. knots placed on current iteration
          n_intknots <- list(X = length(knots$Xk) - length(previous_model$base_learners[[best_bl]]$knots$Xk),
                             Y = length(knots$Yk) - length(previous_model$base_learners[[best_bl]]$knots$Yk))
          coef <-  best_pred$theta
          base_learners_list[[best_bl]]$iterations[[paste0("model", m)]] <- list("int.knt" = int.knt, "coef" = coef)
        }
      
    ## (B) Linear best base-learner
      } else if (base_learners[[best_bl]]$type == "linear") {
        pred_linear <- best_pred$Y_hat
        n <- length(best_pred) - 1; coef <- vector("list", length = n)
        names(coef) <- paste0("b", 0:(n-1))
        coef[] <- best_pred[-1]
        
        int.knt <- n_intknots <- NULL
        
        base_learners_list[[best_bl]]$coefficients <- mapply(function(old_coef, new_coef) old_coef + shrinkage * new_coef,
                                                             previous_model$base_learners[[best_bl]]$coefficients,
                                                             coef, SIMPLIFY=FALSE)
      }
    
    
    ## 5. Evaluate model and recompute gradient vector and residuals
    Y_hat <- models[[paste0("model", m-1)]]$Y_hat + shrinkage*best_pred$Y_hat
    U <- ngradient(data[[response]], Y_hat, weights)
    new.dev <- risk(data[[response]], Y_hat, weights)
    
    # Save model
    models[[paste0("model", m)]] <- list("best_bl" = list("name" = best_bl, "n_intknots" = n_intknots, 
                                                          "int.knt"= int.knt, "coef"= coef, "pred_linear" = pred_linear),
                                         "Y_hat" = Y_hat, "base_learners" = base_learners_list)
    
    
    ## 6. Stopping Rule (Ratio of Deviances)
    # Update the previous dev with current dev
    old.dev <- c(old.dev, new.dev)
   
    if (!is.null(phi_boost_exit) && !is.null(q_boost) && m >= q_boost && m >= min_iterations) {
      # Check if the change in Deviance is less than the tolerance level
      phi_boost =  old.dev[m+1]/ old.dev[m+1-q_boost]
      if (phi_boost >= phi_boost_exit) {
        cat("Stopping iterations due to small improvement in Deviance\n\n")
        break
      }
    }
  }
  
  
  ## 7. Set the "final model" to be the one with lower deviance
  # 7.1. De-normalize predictions if necessary
  if (normalize_data == TRUE && family@name != "Negative Binomial Likelihood (logit link)") {
    models <- lapply(models, function(model) {
      model$Y_hat <- model$Y_hat * args$Y_sd + args$Y_mean
      return(model)
      })
  }
  # 7.2. Calculate each model's deviance
  family_stats <- get_mboost_family(family@name)
  deviances <- sapply(models, function(model) {
    mu <- family@response(model$Y_hat)
    if(family_stats$family == "binomial" && all(args$response[[response]] %in% c(-1, 1))) {
      # since for NGeDSboost(family = "binomial") the encoding is -1/1 and for stats::binomial() the encoding is 0/1
      sum(family_stats$dev.resids((args$response[[response]] + 1) / 2, mu, weights))
      } else {
        sum(family_stats$dev.resids(args$response[[response]], mu, weights))
        }
    })
  # 7.3. Find the model with the smallest deviance
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
      SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$response[[response]],
                         Z = args$predictors[linear_variables],
                         base_learners = args$base_learners, InterKnotsList = qq_list,
                         n = 3, family = get_mboost_family(family@name))}, error = function(e) {
                           cat(paste0("Error computing quadratic fit:", e))
                           return(NULL)
                           })
    final_model$Quadratic.Fit <- quadratic_fit
    pred_quadratic <- as.numeric(quadratic_fit$Predicted)
    
    # Cubic fit
    cc_list <- compute_avg_int.knots(final_model, base_learners = base_learners,
                                   args$X_sd, args$X_mean, normalize_data, n = 4)
    cubic_fit <- tryCatch({
      SplineReg_Multivar(X = args$predictors[GeDS_variables], Y = args$response[[response]],
                         Z = args$predictors[linear_variables],
                         base_learners = args$base_learners, InterKnotsList = cc_list,
                         n = 4, family = get_mboost_family(family@name))}, error = function(e) {
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
  
  # Response (for example, if family is binomial, define the probability estimates)
  f <- as.numeric(final_model$Y_hat); pred_linear <- family@response(f)
  
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
  
  output <- list(extcall = extcall, formula = formula, args = args, models = models, final_model = final_model, predictions = preds,
                 internal_knots = internal_knots, iters = m)
  class(output) <- "GeDSboost"
  return(output)
}


#######################################
### Component-wise fitting function ###
#######################################
componentwise_fit <- function(bl_name, data, response, model_formula_template, weights, base_learners, m, 
                              internal_knots, beta, phi, q) {
  
  # Initialize a results list with default values
  results <- list(
    bl_name = bl_name,
    resid = Inf,
    ssr = Inf,
    pred = NULL,
    error = FALSE,
    message = NULL
  )
  
  pred_vars <- base_learners[[bl_name]]$variables
  
  ## (A) GeDS base-learners
  if (base_learners[[bl_name]]$type == "GeDS") {
    max.intknots <- if (length(pred_vars) == 1) internal_knots
    else if (length(pred_vars) == 2 && internal_knots == 0) stop("internal_knots must be > 0 for bivariate learners")
    else internal_knots
    
    model_formula <- formula(paste0(model_formula_template, bl_name))
    error <- FALSE
    suppressWarnings({
      fit <- tryCatch(
        NGeDS(model_formula, data = data, weights = weights, beta = beta, phi = phi,
              min.intknots = 0, max.intknots = max.intknots, q = q, Xextr = NULL, Yextr = NULL,
              show.iters = FALSE, stoptype = "RD", higher_order = FALSE),
        error = function(e) {
          message(paste0("Error occurred in NGeDS() for base-learner ", bl_name, ": ", e))
          error <<- TRUE
        }
      )
    })
    if (error) {
      results$error <- TRUE
      results$message <- paste0("Skipped iteration ", m, " for base-learner ", bl_name, ".")
      return(results)
    }
    
    # Compute predictions
    pred <- if(length(pred_vars) == 1) {
      predict_GeDS_linear(fit, data[[pred_vars]])
    } else if(length(pred_vars) == 2) {
      predict_GeDS_linear(fit, X = data[pred_vars[1]], Y = data[pred_vars[2]], Z = data[[response]])
    }
    
    ## (B) Linear base-learners
  } else if (base_learners[[bl_name]]$type == "linear") {
    model_formula <- formula(paste0(model_formula_template, bl_name))
    error <- FALSE
    suppressWarnings({
      fit <- tryCatch(
        lm(model_formula, data = data),
        error = function(e) {
          message(paste0("Error occurred in lm() for base-learner ", bl_name, ": ", e))
          error <<- TRUE
        }
      )
    })
    if (error) {
      results$error <- TRUE
      results$message <- paste0("Skipped iteration ", m, " for base-learner ", bl_name, ".")
      return(results)
    }
    
    # Compute predictions
    pred <- list(Y_hat = fit$fitted.values)
    # Create coefficient entries
    for (i in 1:length(fit$coefficients)) {
      coef_name <- paste0("b", i - 1)
      pred[[coef_name]] <- fit$coefficients[[i]]
    }
  }
  
  # Calculate SSR
  resid <- data[[response]] - pred$Y_hat
  ssr <- sum((resid)^2)
  
  # Save results
  results$resid <- resid
  results$ssr <- ssr
  results$pred <- pred
  
  return(results)
}


