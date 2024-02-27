################################################################################
################################################################################
############################# Auxiliar functions ###############################
################################################################################
################################################################################

########################
## 1. Rescale weights ##
########################
rescale_weights<- function (weights) {
  if (max(abs(weights - floor(weights))) < sqrt(.Machine$double.eps))
    return(weights)
  return(weights/sum(weights) * sum(weights > 0))
}

##########################
## 2. get_mboost_family ##
#########################
#' @import mboost
#' @import stats
get_mboost_family <- function(family) {
  if (!is.character(family)) {
    stop("Family must be a character string.")
  }
  
  switch(family,
         # Gaussian
         "Squared Error (Regression)" = stats::gaussian(link = "identity"),
         # Gamma
         "Negative Gamma Likelihood" = stats::Gamma(link = "inverse"),
         # Binomial 
         "Negative Binomial Likelihood (logit link)" = stats::binomial(link = "logit"),
         "Negative Binomial Likelihood -- probit link" = stats::binomial(link = "probit"),
         "Negative Binomial Likelihood -- cauchit link" = stats::binomial(link = "cauchit"),
         "Negative Binomial Likelihood -- log link" = stats::binomial(link = "log"),
         "Negative Binomial Likelihood -- cloglog link" = stats::binomial(link = "cloglog"),
         # Poisson
         "Poisson Likelihood" = stats::poisson(link = "log"),
         stop("Invalid distribution type")
  )
}


#######################################################################################################
## 3. Function to get which is the last row with at least one non-NA in Coefficients/Stored matrices ##
#######################################################################################################
last_row_with_value <- function(matrix) {
  non_na_rows <- which(!apply(is.na(matrix), 1, all))
  if(length(non_na_rows) == 0) {
    return(NULL)
  } else {
    return(max(non_na_rows))
  }
}

###########################################
## 4. Function to get the internal knots ##
###########################################
get_internal_knots <- function(knots, depth = 1) {
  # Helper function to extract internal knots
  extract_knots <- function(k) {
    if (length(k) > 2) {
      return(k[-c(1, length(k))])
    } else {
      return(NULL)
    }
  }
  # Recursive call based on depth
  if (depth > 1) {
    knots <- get_internal_knots(knots, depth - 1)
  }
  # Bivariate
  if (is.list(knots)) {
    # Get internal knots for each component and store them in a list
    ikX <- extract_knots(knots$Xk)
    ikY <- extract_knots(knots$Yk)
    return(list(ikX = ikX, ikY = ikY))
  } else { 
    # Univariate
    return(extract_knots(knots))
  }
}

###################################################################################
## 5. Function for computing GeDS-class object linear predictions (i.e. stage A) ##
###################################################################################
#' @importFrom stats na.omit
predict_GeDS_linear <- function(Gmod, X, Y, Z){
  
  terms <- all.vars(Gmod$Formula)
  num_predictors <- length(terms[-1])
  
  q <- Gmod$Args$q
  
  ## Univariate GeDS
  if(num_predictors == 1){
    
    X <- as.matrix(X)
    
    # Knots
    if (last_row_with_value(Gmod$Stored) > q){
      knt <- na.omit(Gmod$Stored[last_row_with_value(Gmod$Stored) - q,]) # an exit from stage A is performed with the spline fit l = k + 1 - q
    } else {
      knt <- na.omit(Gmod$Stored[last_row_with_value(Gmod$Stored),])
    }
    knt <- knt[-c(1, length(knt))]
    int.knt <- knt[-c(1, length(knt))]
    if(length(knt) == 2) {int.knt <- NULL}
    # Coefficients
    if (last_row_with_value(Gmod$Coefficients) > q){
      coeff <- na.omit(Gmod$Coefficients[last_row_with_value(Gmod$Coefficients) - q,])
    } else {
      coeff <- na.omit(Gmod$Coefficients[last_row_with_value(Gmod$Coefficients),])
    }
    
    # Initialize output vector
    Y_hat <- numeric(nrow(X))
    # Number of intervals
    n_intervals <- length(knt) - 1
    # Add first and last intervals to cover all possible X values
    intervals <- c(-Inf, int.knt, Inf)
    
    # Initialize lists for coefficients
    b0_list <- numeric(n_intervals)
    b1_list <- numeric(n_intervals)
    
    # Calculate predictions for the current predictor
    for (int in 1:n_intervals) {
      # Find X values within current interval
      idx <- X >= intervals[int] & X < intervals[int+1]
      # Calculate predicted Y values for current interval
      epsilon <- 1e-10 # to avoid b1 = Inf in the rare occasions where knt[int+1] = knt[int];
      # check 'cross_validation\HPC\GeDS_initial_learner\ex5\additive\scripts\cv_ex5_1160.R' (seed=6967) for an example.
      b1 <- (coeff[int+1] - coeff[int]) / max(knt[int+1] - knt[int], epsilon)
      b0 <- coeff[int] - b1 * knt[int]
      Y_hat[idx] <- b0 + b1 * X[idx]
      # Store coefficients
      b0_list[int] <- b0
      b1_list[int] <- b1
    }
    return(list(Y_hat = Y_hat, knt = knt, int.knt = int.knt, b0 = b0_list, b1 = b1_list, theta = coeff))
    
  ## Bivariate GeDS
  } else if(num_predictors == 2){
    
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    
    # Knots and Coefficients
    last_row_with_value <- last_row_with_value(Gmod$Stored$previousX)
    
    if (last_row_with_value > q) {
      qplus1_last_rows_X <- Gmod$Stored$previousX[(nrow(Gmod$Stored$previousX) - q):nrow(Gmod$Stored$previousX), ]
      qplus1_last_rows_Y <- Gmod$Stored$previousY[(nrow(Gmod$Stored$previousY) - q):nrow(Gmod$Stored$previousY), ]
      no_NAs <- !any(is.na(qplus1_last_rows_X)) || !any(is.na(qplus1_last_rows_Y))
    } else {
      no_NAs <- FALSE
    }
    
    if (no_NAs || is.null(Gmod$Linear)){
      index <-  last_row_with_value - q
      Xknt <- na.omit(Gmod$Stored$previousX[index,]) # an exit from stage A is performed with the spline fit l = k + 1 - q
      Yknt <- na.omit(Gmod$Stored$previousY[index,])
      theta <- as.numeric(na.omit(Gmod$Coefficients[index,]))
    } else {
      index <-  last_row_with_value
      Xknt <- na.omit(Gmod$Stored$previousX[index,])
      Yknt <- na.omit(Gmod$Stored$previousY[index,])
      theta <- as.numeric(na.omit(Gmod$Coefficients[index,]))
    }
    
    Xknt <- Xknt[-c(1, length(Xknt))]
    Yknt <- Yknt[-c(1, length(Yknt))]
    Xint.knt <- Xknt[-c(1, length(Xknt))]
    Yint.knt <- Yknt[-c(1, length(Yknt))]
    
    if(length(Xknt) == 2) {Xint.knt <- NULL}
    if(length(Yknt) == 2) {Yint.knt <- NULL}
    
    # SplineReg_biv
    lin <- SplineReg_biv(X, Y, Z, InterKnotsX=Xint.knt, InterKnotsY=Yint.knt,
                         Xextr=Gmod$Args$Xextr, Yextr=Gmod$Args$Yextr, n=2,
                         coefficients = theta)
    Y_hat <- lin$Predicted
    
    return(list(Y_hat = Y_hat, knt = list("Xknt" = Xknt, "Yknt" = Yknt),
                int.knt = list("Xint.knt" = Xint.knt, "Yint.knt" = Yint.knt),
                "theta" = lin$Theta))
    
  } else {
    print("Only 1 (i.e. Y ~ f(X)) or 2 (i.e. Z ~ f(X, Y)) predictors are allowed for NGeDS models.")
  }
}

############################################################################
## 6. Function for computing piecewise multivariate additive linear model ##
############################################################################
piecewise_multivar_linear_model <- function(X, model, base_learners = NULL) {
  
  X <- data.frame(X)
  
  # Initialize output vector
  Y_hat <- numeric(nrow(X))
  
  # To compute the boosted prediction made by specific base learners:
  if(!is.null(base_learners)){
    model$base_learners <- model$base_learners[names(model$base_learners) %in% names(base_learners)]
  }
  
  # For each base learner (and its corresponding model)
  for (base_learner in names(model$base_learners)) {
    # Get the model for the current base learner
    model_bl <-  model$base_learners[[base_learner]]
    # Get the predictor variable
    predictor <- base_learners[[base_learner]]$variables
    # Add first and last intervals to cover all possible X values
    intervals <- c(-Inf, model_bl$linear.int.knots, Inf)
    # Number of intervals
    n_intervals <- length(intervals)-1
    
    # Calculate predictions for the current predictor
    for (int in 1:n_intervals) {
      # Find X values within current interval
      idx <- X[,predictor] >= intervals[int] & X[,predictor] < intervals[int+1]
      # Add predicted Y values for current interval to Y_hat
      Y_hat[idx] <- Y_hat[idx] + model_bl$coefficients$b0[int] + model_bl$coefficients$b1[int] * X[idx,predictor]
    }
  }
  return(Y_hat)
}

#####################################################
## 7. Function for computing additive linear model ##
#####################################################
# 7.1. Univariate base-learners
univariate_bl_linear_model <- function(pred_vars, model, shrinkage, base_learners = NULL) {
  
  # Initialize output vector
  Y_hat <- numeric(nrow(pred_vars))
  
  # To compute the boosted prediction made by specific base learners:
  if(!is.null(base_learners)){
    model$base_learners <- model$base_learners[names(model$base_learners) %in% names(base_learners)]
  }
  
  # For each base learner (and its corresponding model)
  for (base_learner in names(model$base_learners)) {
    
    X <- pred_vars[, base_learners[[base_learner]][1]]
    
    bl      <- model$base_learners[[base_learner]]
    int.knt <- bl$linear.int.knots
    theta   <- bl$coefficients
    
    lin <- SplineReg_LM(X, InterKnots=int.knt, extr=range(X), n=2,
                        coefficients = theta)
    
    Y_hat <- Y_hat + lin$Predicted
    
  }
  return(Y_hat)
}

# 7.2. Bivariate base-learners
bivariate_bl_linear_model <- function(pred_vars, model, shrinkage, base_learners = NULL, type = "boost") {
  
  # Initialize output vector
  Y_hat <- numeric(nrow(pred_vars))
  
  # To compute the boosted prediction made by specific base learners:
  if(!is.null(base_learners)){
    model$base_learners <- model$base_learners[names(model$base_learners) %in% names(base_learners)]
  }
  
  # For each base learner (and its corresponding model)
  for (base_learner in names(model$base_learners)) {
    
    X <- pred_vars[, base_learners[[base_learner]]$variables[1]]
    Y <- pred_vars[, base_learners[[base_learner]]$variables[2]]
    
    bl <- model$base_learners[[base_learner]]
    # (I) Bivariate boosted base learners
    if(type == "boost"){
      for (mod in names(bl$iterations)){
        Xint.knt <- bl$iterations[[mod]]$int.knt$Xint.knt
        Yint.knt <- bl$iterations[[mod]]$int.knt$Yint.knt
        theta <- bl$iterations[[mod]]$coef
        lin <- SplineReg_biv(X, Y, InterKnotsX=Xint.knt, InterKnotsY=Yint.knt,
                             Xextr=range(X), Yextr=range(Y), n=2,
                             coefficients = theta)
        if (mod=="model0") {Y_hat <- Y_hat + lin$Predicted}
        else {Y_hat <- Y_hat + shrinkage*lin$Predicted}
      }
    # (II) Bivariate gam
    } else if (type == "gam"){
      Xint.knt <- bl$linear.int.knots$ikX
      Yint.knt <- bl$linear.int.knots$ikY
      theta <- bl$coefficients
      lin <- SplineReg_biv(X, Y, InterKnotsX=Xint.knt, InterKnotsY=Yint.knt,
                           Xextr=range(X), Yextr=range(Y), n=2,
                           coefficients = theta)
      Y_hat <- Y_hat + lin$Predicted
    }
  }
  return(Y_hat)
}


################################
### 8. compute_avg_int.knots ###
################################
# Function for computing the internal knots of base learners (averaging knot location, stage B.1.)
compute_avg_int.knots <- function(final_model, base_learners = base_learners, X_sd, X_mean, normalize_data, n) {
  warning_displayed <- FALSE  # Initialize variable to track warning display
  # Select GeDS base-learners
  base_learners =  base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
  # Check if base_learners is empty
  if (length(base_learners) == 0) {
    cat(paste0("No GeDS base-learners found. Returning NULL when computing averaging knot location for n = ", n))
    return(kk_list = NULL)
  }
  
  kk_list <- lapply(names(base_learners), function(bl) {
    
    pred_vars <- base_learners[[bl]]$variables
    
    # Univariate
    if (length(pred_vars) == 1) {
      intknt <- get_internal_knots(final_model$base_learners[[bl]]$knots)
      if (normalize_data){
        if (!is.null(intknt)) intknt <- intknt * X_sd[[pred_vars[1]]] + X_mean[[pred_vars[1]]]
      }
      if (length(intknt) >= n - 1) return(makenewknots(intknt, n))
      cat(paste0(bl, " has less than ", n - 1, " linear internal knots. "))
      warning_displayed <<- TRUE  # Update the warning flag
      return(intknt)
      # Bivariate
    } else if (length(pred_vars) == 2) {
      intknt <- list(X = get_internal_knots(final_model$base_learners[[bl]]$knots$Xk), 
                     Y = get_internal_knots(final_model$base_learners[[bl]]$knots$Yk))
      if (normalize_data){
        if (!is.null(intknt$X)) {intknt$X <- intknt$X * X_sd[[pred_vars[1]]] + X_mean[[pred_vars[1]]]}
        if (!is.null(intknt$Y)) {intknt$Y <- intknt$Y * X_sd[[pred_vars[2]]] + X_mean[[pred_vars[2]]]}
      }
      if (length(intknt$X) >= n - 1) {
        ikX <- makenewknots(intknt$X, n) 
      } else {
        cat(paste0(bl, " has less than ", n - 1, " linear internal knots for ", pred_vars[1], ". "))
        warning_displayed <<- TRUE  # Update the warning flag
        ikX <- intknt$X
      }
      
      # Number of linear knots warning
      if (length(intknt$Y) >= n - 1) {
        ikY <- makenewknots(intknt$Y, n) 
      } else {
        cat(paste0(bl, " has less than ", n - 1, " linear internal knots for ", pred_vars[2], ". "))
        warning_displayed <<- TRUE  # Update the warning flag
        ikY <- intknt$Y
      }
      return(list(ikX = ikX, ikY = ikY))
    }
  })
  # Check if any warning was displayed and show the appropriate follow-up message
  if (warning_displayed) {
    if (n == 3) {
      cat("Quadratic averaging knot location is not computed for base learners with less than 2 linear internal knots.\n")
    } else if (n == 4) {
      cat("Cubic averaging knot location is not computed for base learners with less than 3 linear internal knots.\n")
    }
  }
  # Naming the elements of the list with the base-learner names
  names(kk_list) <- names(base_learners)
  return(kk_list = kk_list)
}
