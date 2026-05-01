predict_all_orders <- function(object, newdata = NULL, type = c("response", "link")) {
  
  type <- match.arg(type)
  model <- object$final_model
  
  if (is.null(newdata)) {
    if (type == "response") {
      return(object$predictions)
    }
    
    return(list(
      pred_linear    = extract_link_pred(model$linear.fit$temporary),
      pred_quadratic = extract_link_pred(model$quadratic.fit$temporary),
      pred_cubic     = extract_link_pred(model$cubic.fit$temporary)
    ))
  }
  
  orders <- 2L:4L
  names_preds <- c("pred_linear", "pred_quadratic", "pred_cubic")
  
  preds <- lapply(
    orders,
    function(k) predict(object, newdata = newdata, n = k, type = type)
  )
  
  setNames(preds, names_preds)
}

predict_training_data <- function(object, model, n, type) {
  
  if (type == "response") {
    key <- .prediction_key(n)
    pred <- object$predictions[[key]]
    
    if (is.null(pred)) {
      stop(sprintf("Prediction for '%s' not found.", key))
    }
    
    return(pred)
  }
  
  key <- .fit_key(n)
  Fit <- model[[key]]$temporary
  
  if (inherits(Fit, "glm")) {
    return(as.numeric(Fit$linear.predictors))
  }
  
  if (inherits(Fit, "lm")) {
    return(as.numeric(Fit$fitted.values))
  }
  
  stop(sprintf("No valid fitted object found for '%s'.", key))
}

.prediction_key <- function(n) {
  switch(
    as.character(n),
    "2" = "pred_linear",
    "3" = "pred_quadratic",
    "4" = "pred_cubic",
    stop("'n' must be 2, 3, or 4.")
  )
}

.fit_key <- function(n) {
  switch(
    as.character(n),
    "2" = "linear.fit",
    "3" = "quadratic.fit",
    "4" = "cubic.fit",
    stop("'n' must be 2, 3, or 4.")
  )
}


predict_newdata_full_model <- function(object, model, newdata, n, type) {
  
  if (!all(names(object$args$predictors) %in% colnames(newdata))) {
    missing_vars <- setdiff(names(object$args$predictors), colnames(newdata))
    stop(
      paste(
        "The following predictors are missing in newdata:",
        paste(missing_vars, collapse = ", ")
      )
    )
  }
  
  nobs <- nrow(newdata)
  
  selected_bl <- intersect(
    names(model$base_learners),
    names(object$args$base_learners)
  )
  
  base_learners <- object$args$base_learners[selected_bl]
  
  offset <- rep(0, nobs)
  
  if (n == 2L) {
    
    # Extract predictor variables from newdata as a data.frame
    X_df <- newdata[, names(object$args$predictors), drop = FALSE]
    # 1.1 Extract linear base learners
    lin_bl <- base_learners[sapply(base_learners, function(x) x$type == "linear")]
    # 1.2 Extract univariate base learners
    univariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 1)]
    # 1.3 Extract bivariate base learners
    bivariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 2)]
    
    # Family
    if (inherits(object, "GeDSboost")) {
      family_name <- get_mboost_family(object$args$family@name)$family 
    } else {
      family_name <- object$args$family$family
    }
    
    # Normalized model
    if (object$args$normalize_data) {
      # Identify the numeric predictors
      numeric_predictors <- names(X_df)[sapply(X_df, is.numeric)]
      # Scale only the numeric predictors
      if (length(numeric_predictors) > 0) {
        X_df[numeric_predictors] <- scale(
          X_df[numeric_predictors],
          center = object$args$X_mean[numeric_predictors],
          scale  = object$args$X_sd[numeric_predictors]
        )
      }
      if (family_name != "binomial") {
        Y_mean <- object$args$Y_mean; Y_sd <- object$args$Y_sd # Normalized non-binary response
      }
    }
    
    # Initialize prediction vectors
    pred0 <- pred_lin <- pred_univ <- pred_biv <- numeric(nobs)
    
    ####################
    ## 1.1. GeDSboost ##
    ####################
    if (inherits(object, "GeDSboost")) {
      
      # shrinkage
      shrinkage <- object$args$shrinkage
      
      # 1.0 Offset initial learner
      if (!object$args$initial_learner) {
        pred0 <- object$args$family@offset(object$args$response[[1]],  object$args$weights)
        pred0 <- rep(pred0, nrow(newdata)); 
        if(object$args$normalize_data && family_name != "binomial") pred0 <- (pred0-Y_mean)/Y_sd
      }
      # 1.1 Linear base learners
      if (length(lin_bl)!=0) {
        pred_lin <- lin_model(X_df, model, lin_bl)
      }
      # 1.2 GeDS univariate base learners
      if (length(univariate_bl)!=0) {
        pred_univ <- piecewise_multivar_linear_model(X = X_df, model, base_learners = univariate_bl)
      }
      # 1.3 GeDS bivariate base learners
      if (length(bivariate_bl) != 0) {
        pred_biv <- bivariate_bl_linear_model(pred_vars = X_df, model,
                                              shrinkage, base_learners = bivariate_bl,
                                              extr = object$args$extr)
      }
      
      if (object$args$normalize_data && family_name != "binomial") {
        pred <- (pred0 + pred_lin + pred_univ + pred_biv)*Y_sd + Y_mean
      } else {
        pred <- pred0 + pred_lin + pred_univ + pred_biv
      }
      
      pred <- if (type == "response") object$args$family@response(pred) else if (type == "link") pred
      
      return(as.numeric(pred))
      
      ##################
      ## 1.2. GeDSgam ##
      ##################
    } else if (inherits(object, "GeDSgam")) {
      
      # 1.0 Initial learner
      pred0 <- rep(mean(model$Y_hat$z), nrow(newdata))
      # 1.1 Linear base learners
      if (length(lin_bl)!=0) {
        pred_lin <- lin_model(X_df, model, lin_bl)
        if (object$args$normalize_data) {
          # Scale numeric predictors
          scaled <- as.data.frame(lapply(object$args$predictors, function(x) {
            if (is.numeric(x)) scale(x) else x
          }))
          alpha_lin <- mean(lin_model(scaled, model, lin_bl))
        } else {
          alpha_lin <- mean(lin_model(object$args$predictors, model, lin_bl))
        }
        # alpha_lin <- mean(pred_lin)
        pred_lin <- pred_lin - alpha_lin
      }
      # 1.2 GeDS univariate base learners
      if (length(univariate_bl)!= 0) {
        pred_univ <- piecewise_multivar_linear_model(X_df, model, base_learners = univariate_bl)
        if (object$args$normalize_data) {
          alpha_univ <- mean(piecewise_multivar_linear_model(scale(object$args$predictors[numeric_predictors]),
                                                             model, base_learners = univariate_bl))
        } else {
          alpha_univ <- mean(piecewise_multivar_linear_model(object$args$predictors,
                                                             model, base_learners = univariate_bl))
        }
        # alpha_univ <- mean(pred_univ)
        pred_univ <- pred_univ - alpha_univ
      }
      # 1.3 GeDS bivariate base learners
      if (length(bivariate_bl) != 0) {
        pred_biv <- bivariate_bl_linear_model(X_df, model, base_learners = bivariate_bl,
                                              extr = object$args$extr, type = "gam")
        if (object$args$normalize_data) {
          alpha_biv <- mean(bivariate_bl_linear_model(scale(object$args$predictors[numeric_predictors]),
                                                      model, base_learners = bivariate_bl,
                                                      extr = object$args$extr, type = "gam"))
        } else {
          alpha_biv <- mean(bivariate_bl_linear_model(object$args$predictors, model,
                                                      base_learners = bivariate_bl,
                                                      extr = object$args$extr, type = "gam"))
        }
        # alpha_biv <- mean(pred_biv)
        pred_biv <- pred_biv - alpha_biv
      }
      
      if(object$args$normalize_data && family_name != "binomial") {
        pred <- (pred0 + pred_lin + pred_univ + pred_biv)*Y_sd + Y_mean
      } else {
        pred <- pred0 + pred_lin + pred_univ + pred_biv
      }
      
      pred <- if (type == "response") object$args$family$linkinv(pred) else if (type == "link") pred
      
      return(as.numeric(pred))
      
    }
    
  } else if (n == 3L || n == 4L) {
    
    # Extract family from GeDSboost/GeDSgam object
    if (inherits(object, "GeDSboost")) {
      family <- get_mboost_family(object$args$family@name)
      if (!is.null(object$args$link)) family <- get(family$family)(link = object$args$link)
    } else {
      family <- object$args$family
    }
    
    GeDS_variables <- lapply(base_learners, function(x) {if(x$type == "GeDS") return(x$variables) else return(NULL)})
    GeDS_variables <- unname(unlist(GeDS_variables))
    linear_variables <- lapply(base_learners, function(x) {if(x$type == "linear") return(x$variables) else return(NULL)})
    linear_variables <- unname(unlist(linear_variables))
    
    X <- newdata[GeDS_variables]
    Z <- newdata[linear_variables]
    
    if (n == 3) {
      if (is.null(model$quadratic.fit)) {
        cat("No Quadratic Fit to compute predictions.\n")
        return(NULL)
      }
      int.knots <- "quadratic.int.knots"
      Fit <- "quadratic.fit"
    } else if (n == 4) {
      if (is.null(model$cubic.fit)) {
        cat("No Cubic Fit to compute predictions.\n")
        return(NULL)
      }
      int.knots <- "cubic.int.knots"
      Fit <- "cubic.fit"
    }
    
    # Internal Knots
    InterKnotsList <- if (length(object$internal_knots[[int.knots]]) == 0) {
      # if averaging knot location was not computed use linear internal knots
      object$internal_knots$linear.int.knots
    } else {
      object$internal_knots[[int.knots]]
    }
    # Exclude linear learners
    InterKnotsList <- InterKnotsList[setdiff(names(InterKnotsList), names(Z))]
    # GeDS learners
    InterKnotsList_univ <- list()
    for (bl in names(InterKnotsList)) {
      # Check if the length of the variables is equal to 1
      if (length(base_learners[[bl]]$variables) == 1) {
        InterKnotsList_univ[bl] <- InterKnotsList[bl]
      }
    }
    InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
    
    # Select GeDS base-learners
    base_learners <- base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
    # Univariate
    if (length(InterKnotsList_univ) != 0){
      # Create a list to store individual design matrices
      univariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 1]
      univariate_vars <- sapply(univariate_learners, function(bl) bl$variables)
      X_univ <- X[, univariate_vars, drop = FALSE]
      
      # If new data exceeds boundary knots limits, redefine boundary knots
      extrListfit <- lapply(univariate_vars, function(var) range(object$args$predictors[[var]]))
      extrListnew <- lapply(X_univ, range)
      extrList <- mapply(function(var, range1, range2) {
        if(range2[1] < range1[1] || range2[2] > range1[2])
          warning(sprintf("Input values for variable '%s' exceed original boundary knots; extending boundary knots to cover new data range.", var))
        c(min(range1[1], range2[1]), max(range1[2], range2[2]))
      }, var = univariate_vars, range1 = extrListfit, range2 = extrListnew, SIMPLIFY = FALSE)
      
      matrices_univ_list <- vector("list", length = ncol(X_univ))
      # Generate design matrices for each predictor
      for (j in 1:ncol(X_univ)) {
        matrices_univ_list[[j]] <- splineDesign(knots = sort(c(InterKnotsList_univ[[j]], rep(extrList[[j]], n))),
                                                derivs = rep(0, length(X_univ[,j])), x = X_univ[,j], ord = n, outer.ok = TRUE)
      }
      names(matrices_univ_list) <- names(univariate_learners)
    } else {
      matrices_univ_list <- NULL
    }
    # Bivariate
    if (length(InterKnotsList_biv) != 0) {
      bivariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 2]
      matrices_biv_list <- list()
      
      for (learner_name in names(bivariate_learners)) {
        vars <- bivariate_learners[[learner_name]]$variables
        X_biv <- X[, vars, drop = FALSE]
        Xextr <- range(object$args$predictors[vars[1]])
        Yextr <- range(object$args$predictors[vars[2]])
        
        # If new data exceeds boundary knots limits, redefine boundary knots
        if(min(X_biv[,1]) < Xextr[1] || max(X_biv[,1]) > Xextr[2] || min(X_biv[,2]) < Yextr[1] || max(X_biv[,2]) > Yextr[2]) {
          Xextr <- range(c(Xextr, X_biv[,1]))
          Yextr <- range(c(Yextr, X_biv[,2]))
          warning("Input values exceed original boundary knots; extending boundary knots to cover new data range.")
        }
        
        knots <- InterKnotsList_biv[[learner_name]]
        basisMatrixX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                                     x=X_biv[,1],ord=n,outer.ok = TRUE)
        basisMatrixY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                                     x=X_biv[,2],ord=n,outer.ok = TRUE)
        matrices_biv_list[[learner_name]] <- tensorProd(basisMatrixX, basisMatrixY)
      }
    } else {
      matrices_biv_list <- NULL
    }
    
    # Combine all matrices side-by-side
    matrices_list <- c(matrices_univ_list, matrices_biv_list)
    if (!is.null(matrices_list) && length(matrices_list) > 0) {
      full_matrix <- do.call(cbind, matrices_list)
    } else {
      full_matrix <- matrix(ncol = 0, nrow = nrow(Z))
    }
    # Convert any factor columns in Z to dummy variables
    if (NCOL(Z) != 0) {
      Z <- model.matrix(~ ., data = Z)
      Z <-  Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    }
    
    basisMatrix2 <- cbind(full_matrix, as.matrix(Z))
    
    # # Alternative 1 to compute predictions (required @importFrom stats predict)
    # tmp <- model[[Fit]]$temporary
    # # Set the environment of the model's terms to the current environment for variable access
    # environment(tmp$terms) <- environment()
    # pred <- predict(tmp, newdata=data.frame(basisMatrix2), type = "response")
    
    # Alternative 2 to compute predictions
    coefs <- model[[Fit]]$theta
    coefs[is.na(coefs)] <- 0
    pred <- basisMatrix2 %*% coefs + offset
    
    pred <- if (type == "response") family$linkinv(pred) else if (type == "link") pred
    
    return(as.numeric(pred))
    
  }
}

predict_newdata_base_learner <- function(object, model, newdata, n,
                                         base_learner, type) {
  
  # Single base-learner prediction
  bl_name <-  gsub(",\\s*", ", ", as.character(base_learner))
  bl <- object$args$base_learners[[bl_name]]
  if(is.null(bl)) stop(paste0(bl_name, " not found in the model."))
  
  # If GeDSboost and the requested base-learner was never selected by the
  # boosting algorithm, its contribution to the fitted predictor is zero.
  if (inherits(object, "GeDSboost") && !bl_name %in% names(model$base_learners)) {
    return(rep(0, nrow(newdata)))
  }
  
  Y <- object$args$response[[1]]
  pred_vars <- object$args$predictors[bl$variables]
  X_df <- newdata[, intersect(bl$variables, colnames(newdata)), drop = FALSE]
  
  if (object$args$normalize_data && n == 2) {
    # Family
    if (inherits(object, "GeDSboost")) {
      family_name <- get_mboost_family(object$args$family@name)$family 
    } else {
      family_name <- object$args$family$family
    }
    
    # Normalized model
    # Identify the numeric predictors
    numeric_predictors <- names(X_df)[sapply(X_df, is.numeric)]
    # Scale only the numeric predictors
    if (length(numeric_predictors) > 0) {
      X_df[numeric_predictors] <- scale(
        X_df[numeric_predictors],
        center = object$args$X_mean[bl$variables],
        scale  = object$args$X_sd[bl$variables]
      )
    }
    if (family_name != "binomial") {
      Y_mean <- object$args$Y_mean; Y_sd <- object$args$Y_sd # Normalized non-binary response
    }
  }
  
  # Check whether newdata includes all the necessary predictors
  if (!all(bl$variables %in% colnames(newdata))) {
    missing_vars <- setdiff(bl$variables, colnames(newdata))
    stop(paste("The following predictors are missing in newdata:", paste(missing_vars, collapse = ", ")))
  }
  
  # Extract estimated knots and coefficients
  if (n == 2) {
    
    if (inherits(object, "GeDSboost") && is.character(model$linear.fit$theta) ) {
      # model$linear.fit == "When using bivariate base-learners, the 'single spline representation' (in pp form or B-spline form) of the boosted fit is not available.") {
      object$args$base_learners <- object$args$base_learners[bl_name]
      object$args$predictors <- pred_vars
      object$args$extr <- object$args$extr[names(pred_vars)]
      if(object$args$normalize_data) {
        object$args$X_mean <- object$args$X_mean[names(pred_vars)]
        object$args$X_sd <- object$args$X_sd[names(pred_vars)]
      }
      object$args$family <- mboost::Gaussian() # to guarantee pred is returned at the linear predictor level
      
      # Substract the offset initial learner
      if (!object$args$initial_learner) {
        pred0 <- object$args$family@offset(object$args$response[[1]],  object$args$weights)
        if(object$args$normalize_data && family_name != "binomial") pred0 <- (pred0-Y_mean)/Y_sd
      } else {
        pred0 <- 0
      }
      pred0 <- rep(pred0, nrow(newdata))
      
      # Return only the corresponding (normalized) bl-prediction
      pred <- predict(object, newdata = newdata, n = n) - pred0
      if (object$args$normalize_data) {
        pred <- (pred - object$args$Y_mean)/object$args$Y_sd
      }
      
      return(pred)
    }
    
    theta <- model$linear.fit$theta
    int.knots <- object$internal_knots$linear.int.knots
  } else if (n == 3) {
    theta <- model$quadratic.fit$theta
    int.knots <- object$internal_knots$quadratic.int.knots
  } else if (n == 4) {
    theta <- model$cubic.fit$theta
    int.knots <- object$internal_knots$cubic.int.knots
  }
  
  int.knt <- int.knots[[bl_name]]
  
  pattern <- paste0("^", gsub("([()])", "\\\\\\1", bl_name))
  theta <- theta[grep(pattern, names(theta))]
  # Replace NA values with 0
  theta[is.na(theta)] <- 0
  
  # 1. Univariate learners
  if (NCOL(X_df) == 1) {
    
    if (bl$type == "GeDS") {
      
      if (n != 2 && object$args$normalize_data) {
        extr <- range(pred_vars)
      } else {
        extr <- object$args$extr[[bl$variables]]
      }
      
      # If new data exceeds boundary knots limits, redefine boundary knots
      if(min(X_df) < extr[1] || max(X_df) > extr[2]) {
        extr <- range(c(extr, X_df))
        warning("Input values exceed original boundary knots; extending boundary knots to cover new data range.")
      }
      
      # Create spline basis matrix using specified knots, evaluation points and order
      basisMatrix <- splineDesign(knots = sort(c(int.knt,rep(extr,n))),
                                  x = X_df[,1], ord = n, derivs = rep(0,length(X_df[,1])),
                                  outer.ok = T)
      # To recover backfitting predictions need de_mean
      pred <- if (n == 2) basisMatrix %*% theta - mean(basisMatrix %*% theta) else basisMatrix %*% theta 
      
    } else if (bl$type == "linear") {
      # Linear
      if (!is.factor(X_df[,1])) {
        pred <- theta * X_df[,1]
        # Factor
      } else {
        names(theta) <- levels(X_df[,1])[-1]
        theta[levels(X_df[,1])[1]] <- 0 # set baseline coef to 0
        pred <- theta[as.character(X_df[,1])]
      }
    }
    
    return(as.numeric(pred))
    
    # 2. Bivariate learners
  } else if (NCOL(X_df) == 2) {
    
    if (n != 2 && object$args$normalize_data) {
      Xextr <- range(pred_vars[,1])
      Yextr <- range(pred_vars[,2])
    } else {
      Xextr <- object$args$extr[bl$variables][[1]]
      Yextr <- object$args$extr[bl$variables][[2]]
    }
    # If new data exceeds boundary knots limits, redefine boundary knots
    if(min(X_df[,1]) < Xextr[1] || max(X_df[,1]) > Xextr[2] || min(X_df[,2]) < Yextr[1] || max(X_df[,2]) > Yextr[2]) {
      Xextr <- range(c(Xextr, X_df[,1]))
      Yextr <- range(c(Yextr, X_df[,2]))
      warning("Input values exceed original boundary knots; extending boundary knots to cover new data range.")
    }
    
    # Generate spline basis matrix for X and Y dimensions using object knots and given order
    basisMatrixX <- splineDesign(knots = sort(c(int.knt$ikX,rep(Xextr,n))), derivs = rep(0,length(X_df[,1])),
                                 x = X_df[,1], ord = n, outer.ok = T)
    basisMatrixY <- splineDesign(knots = sort(c(int.knt$ikY,rep(Yextr,n))), derivs = rep(0,length(X_df[,2])),
                                 x = X_df[,2], ord = n, outer.ok = T)
    # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
    basisMatrixbiv <- tensorProd(basisMatrixX, basisMatrixY)
    # Multiply the bivariate spline basis by model coefficients to get fitted values
    f_hat_XY_val <- basisMatrixbiv %*% theta[1:dim(basisMatrixbiv)[2]]
    
    return(as.numeric(f_hat_XY_val))
    
  }
  
}