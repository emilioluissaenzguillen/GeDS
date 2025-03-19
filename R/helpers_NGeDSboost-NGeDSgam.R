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

############################
## 2. validate_iterations ##
############################
validate_iterations <- function(iterations, default, name) {
  if (missing(iterations) || is.null(iterations)) {
    iterations <- default
  } else if (!is.numeric(iterations) || iterations %% 1 != 0 || iterations < 0) {
    warning(sprintf("%s must be a non-negative integer; was set to %dL", name, default))
    iterations <- default
  } else {
    iterations <- as.integer(iterations)
  }
  return(iterations)
}

##########################
## 3. get_mboost_family ##
#########################
#' @import mboost
#' @import stats
#' @import mi

# Input should be "family`at`name", with family being a boost_family_glm object
get_mboost_family <- function(family) {
  if (!is.character(family)) {
    stop("Family must be a character string.")
  }
  
  switch(family,
         # Gaussian
         "Squared Error (Regression)" = stats::gaussian(link = "identity"),
         # Gamma; GammaReg():The implemented loss function is the negative Gamma log-likelihood with logarithmic link function
         "Negative Gamma Likelihood" = stats::Gamma(link = "log"),
         # Binomial 
         "Negative Binomial Likelihood (logit link)" = stats::binomial(link = "logit"),
         "Negative Binomial Likelihood -- probit link" = stats::binomial(link = "probit"),
         "Negative Binomial Likelihood -- cauchit link" = stats::binomial(link = "cauchit"),
         "Negative Binomial Likelihood -- log link" = stats::binomial(link = "log"),
         "Negative Binomial Likelihood -- cloglog link" = stats::binomial(link = "cloglog"),
         # Poisson
         "Poisson Likelihood" = stats::poisson(link = "log"),
         "Negative Multinomial Likelihood" = mi::multinomial(link = "logit"),
         stop("Invalid distribution type")
  )
}


#######################################################################################################
## 4. Function to get which is the last row with at least one non-NA in Coefficients/Stored matrices ##
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
## 5. Function to get the internal knots ##
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
## 6. Function for computing GeDS-class object linear predictions (i.e. stage A) ##
###################################################################################
#' @importFrom stats na.omit
predict_GeDS_linear <- function(Gmod, X, Y, Z){
  
  terms <- all.vars(Gmod$Formula)
  num_predictors <- length(terms[-1])
  
  max.intknots <- Gmod$Args$max.intknots
  q <- Gmod$Args$q
  
  ## Univariate GeDS
  if (num_predictors == 1) {
    
    X <- as.matrix(X)
    j <- NROW(Gmod$Stored) # last_row_with_value(Gmod$Stored), last_row_with_value(Gmod$Coefficients)
    
    # Knots
    if (j > q && j != max.intknots + 1) {
      knt <- na.omit(Gmod$Stored[j - q,]) # an exit from stage A is performed with the spline fit l = k + 1 - q
      } else {
        knt <- na.omit(Gmod$Stored[j,])
      }
    
    knt <- knt[-c(1, length(knt))]
    int.knt <- knt[-c(1, length(knt))]
    if (length(knt) == 2) int.knt <- NULL
    
    # Coefficients
    if (j > q && j != max.intknots + 1) {
      coeff <- na.omit(Gmod$Coefficients[j - q,])
      } else {
        coeff <- na.omit(Gmod$Coefficients[j,])
      }
    
    # Initialize output vector
    F_hat <- numeric(nrow(X))
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
      F_hat[idx] <- b0 + b1 * X[idx]
      # Store coefficients
      b0_list[int] <- b0
      b1_list[int] <- b1
    }
    
    if (is.null(Gmod$Args$family)){
      Y_hat <- F_hat
    } else{
      Y_hat <- Gmod$Args$family$linkinv(F_hat)
    }
    
    return(list(Y_hat = Y_hat, knt = knt, int.knt = int.knt, b0 = b0_list, b1 = b1_list, theta = coeff))
    
  ## Bivariate GeDS
  } else if (num_predictors == 2){
    
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    j <- NROW(Gmod$Stored$previousX) #last_row_with_value(Gmod$Stored$previousX)
    
    # Knots and Coefficients
    if (j > q && j != max.intknots + 1) {
      Xknt <- na.omit(Gmod$Stored$previousX[j - q,]) # an exit from stage A is performed with the spline fit l = k + 1 - q
      Yknt <- na.omit(Gmod$Stored$previousY[j - q,])
      theta <- as.numeric(na.omit(Gmod$Coefficients[j - q,]))
      } else {
        Xknt <- na.omit(Gmod$Stored$previousX[j,]) # an exit from stage A is performed with the spline fit l = k + 1 - q
        Yknt <- na.omit(Gmod$Stored$previousY[j,])
        theta <- as.numeric(na.omit(Gmod$Coefficients[j,]))
      }
    
    Xknt <- Xknt[-c(1, length(Xknt))]
    Yknt <- Yknt[-c(1, length(Yknt))]
    Xint.knt <- Xknt[-c(1, length(Xknt))]
    Yint.knt <- Yknt[-c(1, length(Yknt))]
    
    if(length(Xknt) == 2) {Xint.knt <- NULL}
    if(length(Yknt) == 2) {Yint.knt <- NULL}
    
    # SplineReg_biv
    if (Gmod$Type == "LM - Biv") {
      lin <- SplineReg_biv(X, Y, Z, InterKnotsX = Xint.knt, InterKnotsY = Yint.knt,
                           Xextr = Gmod$Args$Xextr, Yextr = Gmod$Args$Yextr, n = 2,
                           coefficients = theta)
      } else if (Gmod$Type == "GLM - Biv") {
        lin <- SplineReg_biv_GLM(X, Y, Z, InterKnotsX = Xint.knt, InterKnotsY = Yint.knt,
                                 Xextr = Gmod$Args$Xextr, Yextr = Gmod$Args$Yextr, n = 2,
                                 family = Gmod$Args$family, coefficients = theta)
        }
    
    Y_hat <-lin$Predicted
    
    return(list(Y_hat = Y_hat, knt = list("Xknt" = Xknt, "Yknt" = Yknt),
                int.knt = list("Xint.knt" = Xint.knt, "Yint.knt" = Yint.knt),
                "theta" = lin$Theta))
    } else {
      print("Only 1 (i.e. Y ~ f(X)) or 2 (i.e. Z ~ f(X, Y)) predictors are allowed for NGeDS models.")
    }
}

############################################################################
## 7. Function for computing piecewise multivariate additive linear model ##
############################################################################
# 7.1
lin_model <- function(pred_vars, model, lin_bl, nobs) {
  pred_lin <- numeric(nobs)  # Initialize prediction vector
  
  # Loop over base learners
  for (bl in names(lin_bl)) {
    # Extract coefficients
    coeffs <- model$base_learners[[bl]]$coefficients
    
    # (i) Categorical
    if (is.factor(pred_vars[[bl]])) {
      # Compute fitted values manually
      baseline <- levels(pred_vars[[bl]])[1]
      pred_bl <- numeric(nobs)
      
      # For baseline level
      pred_bl[pred_vars[[bl]] == baseline] <- coeffs[1]
      
      # Loop through all levels except the baseline
      for (i in 2:length(levels(pred_vars[[bl]]))) {
        level <- levels(pred_vars[[bl]])[i]
        mask <- pred_vars[[bl]] == level
        pred_bl[mask] <- coeffs[[1]] + coeffs[[i]]
      }
      pred_bl <- as.numeric(pred_bl)
    } else {  # (ii) Continuous
      pred_bl <- coeffs$b0 + coeffs$b1 * pred_vars[[bl]]
    }
    
    # Add to overall prediction
    pred_lin <- pred_lin + pred_bl
  }
  
  return(pred_lin)
}

# 7.2
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
## 8. Function for computing additive linear model ##
#####################################################
# 8.1. Univariate base-learners
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

# 8.2. Bivariate base-learners
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
    if(type == "boost") {
      
      for (mod in names(bl$iterations)){
        Xint.knt <- bl$iterations[[mod]]$int.knt$Xint.knt
        Yint.knt <- bl$iterations[[mod]]$int.knt$Yint.knt
        theta <- bl$iterations[[mod]]$coef
        lin <- SplineReg_biv(X, Y, InterKnotsX=Xint.knt, InterKnotsY=Yint.knt,
                             Xextr=range(X), Yextr=range(Y), n=2,
                             coefficients = theta)
        
        if (mod=="model0") {
          Y_hat <- Y_hat + lin$Predicted # initial learner is added with no shrinking
        } else {
            Y_hat <- Y_hat + shrinkage*lin$Predicted
        }
        
      }
    # (II) Bivariate gam
    } else if (type == "gam") {
      
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
### 9. compute_avg_int.knots ###
################################
# Function for computing the internal knots of base learners (averaging knot location, stage B.1.)
compute_avg_int.knots <- function(final_model, base_learners = base_learners, X_sd, X_mean, normalize_data, n) {
  # Select GeDS base-learners
  base_learners =  base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
  # Check if base_learners is empty
  if (length(base_learners) == 0) {
    cat(paste0("No GeDS base-learners found. Returning NULL when computing averaging knot location for n = ", n,". "))
    return(kk_list = NULL)
  }
  
  warning_messages <- c()
  
  kk_list <- lapply(names(base_learners), function(bl) {
    pred_vars <- base_learners[[bl]]$variables
    
    # Univariate case
    if (length(pred_vars) == 1) {
      intknt <- get_internal_knots(final_model$base_learners[[bl]]$knots)
      if (normalize_data) {
        if (!is.null(intknt)) intknt <- intknt * X_sd[[pred_vars[1]]] + X_mean[[pred_vars[1]]]
      }
      if (length(intknt) >= n - 1) return(makenewknots(intknt, n))
      
      warning_messages <<- c(warning_messages, paste0(bl, " has less than ", n - 1, " linear internal knots."))
      return(intknt)
      
      # Bivariate case
    } else if (length(pred_vars) == 2) {
      intknt <- list(
        X = get_internal_knots(final_model$base_learners[[bl]]$knots$Xk), 
        Y = get_internal_knots(final_model$base_learners[[bl]]$knots$Yk)
      )
      if (normalize_data) {
        if (!is.null(intknt$X)) intknt$X <- intknt$X * X_sd[[pred_vars[1]]] + X_mean[[pred_vars[1]]]
        if (!is.null(intknt$Y)) intknt$Y <- intknt$Y * X_sd[[pred_vars[2]]] + X_mean[[pred_vars[2]]]
      }
      
      ikX <- if (length(intknt$X) >= n - 1) {
        makenewknots(intknt$X, n)
      } else {
        warning_messages <<- c(warning_messages, paste0(bl, " has less than ", n - 1, " linear internal knots for ", pred_vars[1], "."))
        intknt$X
      }
      
      ikY <- if (length(intknt$Y) >= n - 1) {
        makenewknots(intknt$Y, n)
      } else {
        warning_messages <<- c(warning_messages, paste0(bl, " has less than ", n - 1, " linear internal knots for ", pred_vars[2], "."))
        intknt$Y
      }
      
      return(list(ikX = ikX, ikY = ikY))
    }
  })
  
  # Print warnings if they are below the threshold
  warning_threshold <- 10
  if (n >= 3 && length(warning_messages) > 0 && length(warning_messages) <= warning_threshold) {
    cat(paste(warning_messages, collapse = "\n"), sep = "")
  }
  # Check if any warning was displayed and show the appropriate follow-up message
  if (length(warning_messages) > 0) {
    if (n == 3) {
      cat(" Quadratic averaging knot location is not computed for base learners with less than 2 internal knots.\n")
    } else if (n == 4) {
      cat(" Cubic averaging knot location is not computed for base learners with less than 3 internal knots.\n")
    }
  }
  # Naming the elements of the list with the base-learner names
  names(kk_list) <- names(base_learners)
  return(kk_list = kk_list)
}

################################################################################

bSpline.coef <- function(final_model, univariate_learners)
  {
  learner_names <- names(univariate_learners)
  knots <- lapply(learner_names, function(bl) final_model$base_learners[[bl]]$knots)
  names(knots) <- learner_names
  poly_coef <- lapply(learner_names, function(bl) final_model$base_learners[[bl]]$coefficients)
  names(poly_coef) <- learner_names
  
  epsilon <- 1e-10
  coeff <- setNames(vector("list", length(learner_names)), learner_names) 
  
  # Loop through each learner
  for (bl in learner_names) {
    knt <- knots[[bl]]
    b0 <- poly_coef[[bl]]$b0
    b1 <- poly_coef[[bl]]$b1
    
    # Calculate coefficient
    coeff[[bl]] <- b0 + b1 * knt[-length(knt)]
    # Calculate the last element to append
    last_knot_diff <- max(knt[length(knt)] - knt[length(knt) - 1], epsilon)
    last_coeff <- coeff[[bl]][length(coeff[[bl]])] + b1[length(b1)] * last_knot_diff 
    
    # Append the last value to the coefficient vector
    coeff[[bl]] <- c(coeff[[bl]], last_coeff)
    
    # # Check
    # round(b1 - diff(coeff[[bl]]) / pmax(diff(knt), epsilon), 10)
    # round(b0 - ( coeff[[bl]][-length(coeff[[bl]])] - b1 * knt[-length(knt)] ), 10)
  }
  
  
  return(coeff)
}


biv_bSpline.coef <- function(final_model, bivariate_learners, shrinkage)
{
  bivbl_names <- names(bivariate_learners)
  biv_models <- lapply(bivbl_names, function(bl) final_model$base_learners[[bl]]$iterations)
  names(biv_models) <- bivbl_names
  
  theta_biv <- setNames(vector("list", length(bivbl_names)), bivbl_names) 
  for (bl in bivbl_names) {
    bl_models <- biv_models[[bl]]
    theta <- list()
    
    # Loop through each model in the base-learner
    for (model_name in names(bl_models)) {
      # Extract model number
      model_num <- as.numeric(gsub("model", "", model_name))
      #Extract coefs
      coefs <- biv_models[[bl]][[model_name]]$coef
      # Apply shrinkage
      if (model_name != "model0") coefs <- shrinkage * coefs
      # Store coefficients
      theta[[paste0("theta", model_num)]] <- coefs
      }
    
    # Store results for the base learner
    theta_biv[[bl]] <- theta
    
  }
  return(theta_biv)
}


### Validate Y_hat ##
# At each boosting iteration, you need to transform Y_hat into F_hat
# in order to recompute the gradient vector. Sometimes Y_hat might violate the requirements
# of the corresponding distribution. This function ensures Y_hat meets those requirements
validate_Y_hat <- function(Y_hat, family_stats) {
  valid <- family_stats$validmu(Y_hat)
  
  if (all(valid)) {
    return(Y_hat)
  }
  
  if (family_stats$family == "binomial") {
    # Ensure values are between 0 and 1
    Y_hat <- pmax(pmin(Y_hat, 1), 0)
  } else if (family_stats$family == "poisson" || family_stats$family == "Gamma" || family_stats$family == "inverse.gaussian") {
    # Ensure values are positive
    Y_hat <- pmax(Y_hat, .Machine$double.eps)
  }
  # Add more conditions for other families if necessary
  
  valid <- family_stats$validmu(Y_hat)
  
  if (!all(valid)) {
    warning("Some Y_hat values are still invalid after correction.")
  }
  
  return(Y_hat)
}

##########################
## count_initial_zeroes ##
##########################
# To know the rounding when defining the shrinkage rate as 2/max(Y)
count_initial_zeroes <- function(x) {
  # Convert the number to a string
  num_str <- as.character(x)
  
  # Split the string at the decimal point
  parts <- strsplit(num_str, split = "\\.")[[1]]
  
  # Check if there's a decimal part
  if (length(parts) < 2) {
    return(0)  # No decimal part means no zeroes after the decimal
  }
  
  # Get the decimal part
  decimal_part <- parts[2]
  
  # Count the number of leading zeroes in the decimal part
  leading_zeroes <- nchar(gsub("^0*", "", decimal_part, perl = TRUE)) - nchar(decimal_part)
  return(-leading_zeroes + 1)
} 

