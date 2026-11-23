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
#' @importFrom stats gaussian Gamma binomial poisson

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
         #"Negative Multinomial Likelihood" = mi::multinomial(link = "logit"),
         stop("Invalid distribution type")
  )
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
    # Preserve the historical bivariate names and named coordinate knot lists.
    if (any(c("Xk", "Xknt") %in% names(knots))) {
      Xk <- if (!is.null(knots$Xk)) knots$Xk else knots$Xknt
      Yk <- if (!is.null(knots$Yk)) knots$Yk else knots$Yknt
      return(list(ikX = extract_knots(Xk), ikY = extract_knots(Yk)))
    }
    return(lapply(knots, extract_knots))
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

  terms <- all.vars(Gmod$formula)
  num_predictors <- length(terms[-1])

  ## Univariate GeDS
  if (num_predictors == 1) {

    X <- as.matrix(X)
    # Use the selected refit, including exits without an admissible next knot.
    int.knt <- Gmod$linear.intknots
    knt <- sort(c(Gmod$args$extr, int.knt))
    coeff <- as.numeric(Gmod$linear.fit$theta)

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

    if (is.null(Gmod$args$family)){
      Y_hat <- F_hat
    } else{
      Y_hat <- Gmod$args$family$linkinv(F_hat)
    }

    return(list(Y_hat = Y_hat, knt = knt, int.knt = int.knt, b0 = b0_list, b1 = b1_list, theta = coeff))

  ## Bivariate GeDS
  } else if (num_predictors == 2){

    X <- as.matrix(X)
    Y <- as.matrix(Y)
    # The fitter may have rejected a later rank-deficient tensor basis.
    # Keep its selected knots and coefficients together, as above for one dimension.
    Xint.knt <- Gmod$linear.intknots$Xk
    Yint.knt <- Gmod$linear.intknots$Yk
    theta <- as.numeric(Gmod$linear.fit$theta)
    Xknt <- sort(c(Gmod$args$Xextr, Xint.knt))
    Yknt <- sort(c(Gmod$args$Yextr, Yint.knt))

    # SplineReg_biv
    if (Gmod$type == "LM - Biv") {
      lin <- SplineReg_biv(X, Y, Z, InterKnotsX = Xint.knt, InterKnotsY = Yint.knt,
                           Xextr = Gmod$args$Xextr, Yextr = Gmod$args$Yextr, n = 2,
                           coefficients = theta)
      } else if (Gmod$type == "GLM - Biv") {
        lin <- SplineReg_biv_GLM(X, Y, Z, InterKnotsX = Xint.knt, InterKnotsY = Yint.knt,
                                 Xextr = Gmod$args$Xextr, Yextr = Gmod$args$Yextr, n = 2,
                                 family = Gmod$args$family, coefficients = theta)
        }

    Y_hat <-lin$predicted

    return(list(Y_hat = Y_hat, knt = list("Xk" = Xknt, "Yk" = Yknt),
                int.knt = list("Xint.knt" = Xint.knt, "Yint.knt" = Yint.knt),
                "theta" = lin$theta))
    } else {
      print("Only 1 (i.e. Y ~ f(X)) or 2 (i.e. Z ~ f(X, Y)) predictors are allowed for NGeDS models.")
    }
}

# Fit the linear (Stage A) tensor GeDS learner used inside GAM backfitting and
# component-wise boosting.  The surrounding algorithms always fit a Gaussian
# working response here, irrespective of the response family's final link.

############################################################################
## 7. Function for computing piecewise multivariate additive linear model ##
############################################################################
# 7.1
lin_model <- function(pred_vars, model, lin_bl) {

  pred_vars <- as.data.frame(pred_vars)
  nobs <- nrow(pred_vars)
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
piecewise_multivar_linear_model <- function(X, model, base_learners = NULL,
                                             clamp_linear = FALSE, extr = NULL) {

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
    if (clamp_linear) {
      boundary <- extr[[predictor]]
      X[[predictor]] <- pmax(boundary[1], pmin(boundary[2], X[[predictor]]))
    }
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


# 8.2. Bivariate base-learners
bivariate_bl_linear_model <- function(pred_vars, model, shrinkage, base_learners = NULL,
                                      extr, type = "boost", clamp_linear = FALSE) {

  # Initialize output vector
  Y_hat <- numeric(nrow(pred_vars))

  # To compute the boosted prediction made by specific base learners:
  if(!is.null(base_learners)){
    model$base_learners <- model$base_learners[names(model$base_learners) %in% names(base_learners)]
  }

  # For each base learner (and its corresponding model)
  for (base_learner in names(model$base_learners)) {

    vars <- base_learners[[base_learner]]$variables
    bl <- model$base_learners[[base_learner]]

    X <- pred_vars[, vars[1]]
    Y <- pred_vars[, vars[2]]

    # Boundary knots
    extrX <- extr[[base_learners[[base_learner]]$variables[1]]]
    extrY <- extr[[base_learners[[base_learner]]$variables[2]]]

    predict_bivariate <- function(Xint.knt, Yint.knt, theta) {
      basis <- tensorProd(
        .prediction_spline_basis(sort(c(Xint.knt, rep(extrX, 2))), X, 2L, clamp_linear),
        .prediction_spline_basis(sort(c(Yint.knt, rep(extrY, 2))), Y, 2L, clamp_linear)
      )
      if (length(theta) != ncol(basis)) {
        stop("Stored bivariate coefficients are not conformable with the fitted spline basis.",
             call. = FALSE)
      }
      theta[is.na(theta)] <- 0
      as.numeric(basis %*% theta)
    }

    # (I) Bivariate boosted base learners
    if(type == "boost") {

      for (mod in names(bl$iterations)){
        Xint.knt <- bl$iterations[[mod]]$int.knt$Xint.knt
        Yint.knt <- bl$iterations[[mod]]$int.knt$Yint.knt
        theta <- bl$iterations[[mod]]$coef
        predicted <- predict_bivariate(Xint.knt, Yint.knt, theta)

        if (mod=="model0") {
          Y_hat <- Y_hat + predicted # initial learner is added with no shrinking
        } else {
            Y_hat <- Y_hat + shrinkage*predicted
        }

      }
    # (II) Bivariate gam
    } else if (type == "gam") {

      Xint.knt <- bl$linear.int.knots$ikX
      Yint.knt <- bl$linear.int.knots$ikY
      theta <- bl$coefficients
      Y_hat <- Y_hat + predict_bivariate(Xint.knt, Yint.knt, theta)
    }
  }
  return(as.numeric(Y_hat))
}


# Intercept of the linear (order-2) GeDSgam predictor on the internal (offset-free,
# possibly normalized) scale. Backfitting centres the intercept and each learner with
# the IRLS weights, whereas prediction re-centres each learner with ordinary training
# means. Using the ordinary mean of the fitted additive predictor z_hat keeps both
# representations equal on the training data. Objects fitted before z_hat was stored
# fall back to the previous mean of the working response.
.gam_linear_intercept <- function(Y_hat) {
  if (!is.null(Y_hat$z_hat)) mean(Y_hat$z_hat) else mean(Y_hat$z)
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
      return(NULL)

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
        NULL
      }

      ikY <- if (length(intknt$Y) >= n - 1) {
        makenewknots(intknt$Y, n)
      } else {
        warning_messages <<- c(warning_messages, paste0(bl, " has less than ", n - 1, " linear internal knots for ", pred_vars[2], "."))
        NULL
      }

      return(list(ikX = ikX, ikY = ikY))
    } else {
      stop("Each GeDS spline component supports at most two covariates.", call. = FALSE)
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
      cat(" Quadratic fits use zero internal knots in dimensions with fewer than 2 linear internal knots.\n")
    } else if (n == 4) {
      cat(" Cubic fits use zero internal knots in dimensions with fewer than 3 linear internal knots.\n")
    }
  }
  # Naming the elements of the list with the base-learner names
  names(kk_list) <- names(base_learners)
  return(kk_list = kk_list)
}

################################################################################

lin.coef <- function(final_model, linear_variables) {
  # Handle linear and factor bl/variables
  linear_coef <- NULL
  if (length(linear_variables) > 0) {
    # Initialize the intercept and slopes
    intercept <- 0; slopes <- NULL
    # Loop through each variable in linear_variables
    for (var in linear_variables) {
      coef <- final_model$base_learners[[var]]$coefficients
      int <- coef$b0
      slp <- unlist(coef[names(coef) != "b0"])
      names(slp) <- paste0(var, names(slp))
      # Sum up the intercepts
      intercept <- intercept + int
      # Concatenate the slopes
      slopes <- c(slopes, slp)
    }
    names(intercept) <- "b0"
    linear_coef <- c(intercept, slopes)
  }
  return(linear_coef)
}

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


########################
### 10. Stage B fits ###
########################
stageB_fit <- function(n, args, GeDS_variables, linear_variables,
                       response, final_model, normalize_data,
                       weights = FALSE, offset = FALSE, link = FALSE,
                       only_pred = TRUE, penalty = "none", stabilise = TRUE) {

  X_mat <- args$predictors[, GeDS_variables, drop = FALSE]
  Z_mat <- args$predictors[, linear_variables, drop = FALSE]
  Y_vec <- args$response[[response]]
  n_obs <- length(Y_vec)

  # Turn boolean flags into actual values (or NULL)
  wts <- if (isTRUE(weights)) args$weights else rep(1, n_obs)
  off <- if (isTRUE(offset))  args$offset  else rep(0, n_obs)
  lnk <- if (isTRUE(link))    args$link    else NULL
  extr <- lapply(
    X_mat[wts > 0, , drop = FALSE],
    range
  )

  n_informative <- sum(wts > 0)

  # 1) knots
  iklist <- compute_avg_int.knots(
    final_model,
    base_learners = args$base_learners[names(final_model$base_learners)],
    args$X_sd, args$X_mean, normalize_data,
    n = n
  )

  spline_learners <- args$base_learners[names(iklist)]
  basis_sizes <- vapply(names(iklist), function(bl) {
    dims <- length(spline_learners[[bl]]$variables)
    knots <- iklist[[bl]]
    if (dims == 1L) length(knots) + n else prod(lengths(knots) + n)
  }, numeric(1))
  total_coef <- sum(basis_sizes) + NCOL(Z_mat) + as.integer(NCOL(Z_mat) > 0L)
  max.coef <- if (is.null(args$max.coef)) 100000L else args$max.coef
  if (total_coef > max.coef) {
    if (penalty != "none")
      stop("The penalised Stage B basis exceeds 'max.coef'.", call. = FALSE)
    warning("Omitting the order-", n, " fit because its additive tensor basis requires ",
            total_coef, " coefficients, exceeding 'max.coef' = ", max.coef, ".",
            call. = FALSE)
    return(list(fit = NULL, preds = rep(NA_real_, n_obs), int.knots = iklist))
  }

  # A joint tensor with at least as many basis coefficients as informative
  # observations is not identifiable by the unpenalized Stage B regression.
  # In that case lm/glm may return rank-deficient coefficients which later
  # turn otherwise valid predictions into NA/NaN values. Additive collections
  # of univariate learners retain their existing behaviour.
  joint_not_estimable <- any(vapply(names(iklist), function(bl) {
    length(spline_learners[[bl]]$variables) > 1L &&
      basis_sizes[[bl]] >= n_informative
  }, logical(1)))
  if (joint_not_estimable && penalty == "none") {
    warning(
      "Omitting the order-", n,
      " fit because a joint tensor basis is not estimable from the available observations.",
      call. = FALSE
    )
    return(list(fit = NULL, preds = rep(NA_real_, n_obs), int.knots = iklist))
  }

  if (penalty != "none") {
    fit <- .stageB_penalised_fit(
      X_mat, Y_vec, Z_mat, off, wts, spline_learners, iklist, n,
      args$family, lnk, penalty
    )
    return(list(fit = fit, preds = fit$predicted, int.knots = iklist))
  }

  guarded <- if (stabilise) .stageB_stabilise(
    X_mat, Z_mat, wts, spline_learners, iklist, n,
    boundaries = lapply(X_mat, range)
  ) else list(internal = iklist, stability = NULL)
  iklist <- guarded$internal

  # (IR)LS fit
  fit <- tryCatch(
    suppressMessages(
      SplineReg_Multivar(X = X_mat, Y = Y_vec,
                         Z = Z_mat, offset = off,
                         weights = wts,
                         base_learners = args$base_learners[names(final_model$base_learners)],
                         InterKnotsList = iklist, n = n, extrList = extr,
                         family = args$family,
                         link = lnk, only_pred = only_pred)
    ),
    error = function(e) { warning(sprintf("Error computing n=%d fit: %s", n, conditionMessage(e))); NULL }
  )

  preds <- if (!is.null(fit) && !is.null(fit$predicted)) as.numeric(fit$predicted) else rep(NA_real_, n_obs)
  if (!is.null(fit)) fit$stability <- guarded$stability
  list(fit = fit, preds = preds, int.knots = iklist)
}

