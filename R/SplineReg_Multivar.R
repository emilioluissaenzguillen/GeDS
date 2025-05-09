
########################
## SplineReg_Multivar ##
########################
SplineReg_Multivar <- function(X, Y, Z = NULL, offset = rep(0, NROW(Y)), base_learners,
                               weights = rep(1, NROW(Y)), InterKnotsList, n, 
                               extrList = lapply(X, range), prob = 0.95, family, link,
                               mustart = NULL, inits = NULL, coefficients = NULL,
                               linear_intercept = FALSE, de_mean = FALSE) {
  
  if (inherits(family, "boost_family_glm") || inherits(family, "boost_family")) {
    family_name <- get_mboost_family(family@name)$family
    if (!is.null(link)) family <- get(family_name)(link = link)
    
    # Since for NGeDSboost(family = "binomial") the encoding is -1/1
    # and for stats::binomial() the encoding is 0/1
    if (family_name == "binomial") {
      Y <- (Y + 1) / 2
      }
    
    } else {
      family_name <- family$family
    }
  
  # Check if family is "gaussian"
  if (family_name == "gaussian") {
    # Call Spline_LM_Multivar function
    result <- SplineReg_LM_Multivar(X, Y, Z, offset, base_learners, weights, InterKnotsList,
                                    n, extrList, prob, coefficients, linear_intercept, de_mean)
    } else {
      # Call Spline_GLM_Multivar function
      result <- SplineReg_GLM_Multivar(X, Y, Z, offset, base_learners, weights, InterKnotsList,
                                       n, extrList, family, mustart, inits, coefficients, linear_intercept, de_mean)
    }
  
  return(result)
}

###########################
## SplineReg_LM_Multivar ##
###########################
#'
#' @importFrom MASS ginv
#' 
SplineReg_LM_Multivar <- function(X, Y, Z = NULL, offset = rep(0, NROW(Y)), base_learners,
                                  weights = rep(1, NROW(Y)), InterKnotsList, n, extrList = lapply(X, range),
                                  prob = 0.95, coefficients, linear_intercept = FALSE, de_mean = FALSE)
{
  n <- as.integer(n)
  
  InterKnotsList_univ <- list()
  for (bl in names(InterKnotsList)) {
    # Check if the length of the variables is equal to 1
    if (length(base_learners[[bl]]$variables) == 1) {
      InterKnotsList_univ[bl] <- InterKnotsList[bl]
    }
  }
  InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
  
  # Select GeDS base-learners
  base_learners <- if (length(base_learners) > 0) base_learners[sapply(base_learners, function(x) x$type == "GeDS")] else NULL
  univariate_learners <- bivariate_learners <- NULL
  # Univariate
  if (length(InterKnotsList_univ) != 0) {
    # Create a list to store individual design matrices
    univariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 1]
    univariate_vars <- sapply(univariate_learners, function(bl) bl$variables)
    X_univ <- X[, univariate_vars, drop = FALSE]
    extrList = lapply(X_univ, range)
    matrices_univ_list <- vector("list", length = ncol(X_univ))
    # Generate design matrices for each predictor
    for (j in 1:ncol(X_univ)) {
      matrices_univ_list[[j]] <- splineDesign(knots = sort(c(InterKnotsList_univ[[j]], rep(extrList[[j]], n))),
                                              derivs = rep(0, length(X_univ[,j])), x = X_univ[,j], ord = n, outer.ok = TRUE)
    }
    names(matrices_univ_list) <- names(univariate_learners)
    
    # Assign base-learner name to the columns of each of the matrices
    for (matrix_name in names(matrices_univ_list)) {
      num_cols <- ncol(matrices_univ_list[[matrix_name]])
      col_names <- paste(matrix_name, 1:num_cols, sep = "_")
      colnames(matrices_univ_list[[matrix_name]]) <- col_names
    }
    
  } else {
    matrices_univ_list <- NULL
  }
  # Bivariate
  if (length(InterKnotsList_biv) != 0){
    bivariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 2]
    matrices_biv_list <- list()
    matrices_biv_list_aux <- list()
    
    for(learner_name in names(bivariate_learners)){
      vars <- bivariate_learners[[learner_name]]$variables
      X_biv <- X[, vars, drop = FALSE]
      Xextr = range(X_biv[,1])
      Yextr = range(X_biv[,2])
      knots <- InterKnotsList_biv[[learner_name]]
      basisMatrixX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                               x=X_biv[,1],ord=n,outer.ok = TRUE)
      basisMatrixY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                               x=X_biv[,2],ord=n,outer.ok = TRUE)
      
      # To help saving control polygon knots afterwards
      matrices_biv_list_aux[[learner_name]] <- list(basisMatrixX, basisMatrixY)
      names(matrices_biv_list_aux[[learner_name]]) <- vars
      
      matrices_biv_list[[learner_name]] <- tensorProd(basisMatrixX, basisMatrixY)
      
      # Assign base-learner name to the columns of each of the matrices
      for (matrix_name in names(matrices_biv_list)) {
        num_cols <- ncol(matrices_biv_list[[matrix_name]])
        col_names <- paste(matrix_name, 1:num_cols, sep = "_")
        colnames(matrices_biv_list[[matrix_name]]) <- col_names
      }
      
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
  if (!is.null(Z) && NCOL(Z) > 0) {
    Z <- model.matrix(~ ., data = Z)
    if (!linear_intercept) Z <-  Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    } else {
      Z <- NULL
    }
  
  basisMatrix2 <- cbind(full_matrix, Z)
  
  # 1) If coefficients are NOT provided estimate the corresponding regression model
  Y0 <- Y - offset
  if (is.null(coefficients)) {
    tmp <- lm(Y0 ~ -1 + basisMatrix2, weights = as.numeric(weights))
    # the ‘-1’ serving to suppress the redundant extra intercept that would be added by default
    # 'splineDesign' already includes a basis that accounts for the intercept
    theta <- coef(tmp)
    # Check if any coefficient is NA, which indicates a rank deficiency and recover theta
    if (any(is.na(theta))) {
      # # Compute the minimal-norm solution for theta using the Moore-Penrose generalized inverse.
      # theta <- as.numeric(ginv(basisMatrix2) %*% Y0)
      # # Now theta contains the computed coefficients that reproduce lm()'s fitted values.
      
      # Compute t(basisMatrix2) %*% basisMatrix2 using crossprod (more efficient)
      matcb <- crossprod(basisMatrix2)
      matcbinv <- MASS::ginv(matcb)
      theta <- as.numeric(matcbinv %*% crossprod(basisMatrix2, tmp$fitted.values))
      
    }
    names(theta) <- sub("basisMatrix2", "", names(coef(tmp)))
    predicted <- tmp$fitted.values + offset
    
    # Reset environment of lm object
    f <- tmp$terms
    environment(f) <- .GlobalEnv
    tmp$terms <- f
    terms_tmp <- attr(tmp$model, "terms")
    environment(terms_tmp) <- .GlobalEnv
    attr(tmp$model, "terms") <- terms_tmp
  
  # 2) If coefficients are provided, use them to compute predicted values directly
  } else {
    theta <- coefficients
    fitted.values <- if (de_mean) basisMatrix2 %*% theta - mean(basisMatrix2 %*% theta) else basisMatrix2 %*% theta # to recover backfitting predictions need de_mean
    predicted <- fitted.values + offset
    # Manually reconstruct lm object to facilitate S3 methods application
    form <- Y0 ~ -1 + basisMatrix2
    mf   <- model.frame(form, weights = as.numeric(weights))
    mm    <- model.matrix(form, data = mf)
    tmp <- list(
      coefficients = setNames(theta, colnames(mm)),
      residuals = setNames(as.numeric(Y0 - fitted.values), rownames(mf)),
      fitted.values = setNames(as.numeric(fitted.values), rownames(mf)),
      weights = as.numeric(weights),
      rank = qr(mm)$rank,  # or qr(basisMatrix2)$rank or rankMatrix(basisMatrix2) or length(coefficients)
      qr = qr(mm), # or qr(basisMatrix2)
      df.residual = as.numeric(nrow(basisMatrix2) - rankMatrix(basisMatrix2)), # residual degrees of freedom
      y = Y0,
      terms = terms(mf),
      model = mf
    )
    class(tmp) <- "lm"
    
  }
  
  # Control polygon knots
  polyknots_list <- list()
  # Univariate 
  if (length(InterKnotsList_univ) != 0) {
    for (learner_name in names(InterKnotsList_univ)) {
      if (!is.null(InterKnotsList_univ[[learner_name]]) &&
          length(InterKnotsList_univ[[learner_name]]) >= n - 1) {
        polyknots_list[[learner_name]] <-
          makenewknots(
            sort(c(InterKnotsList_univ[[learner_name]],
                   rep(extrList[[univariate_learners[[learner_name]]$variables]], n)))[-c(1, NCOL(matrices_univ_list[[learner_name]]) + 1)],
            degree = n
            )
        } else {
          polyknots_list[[learner_name]] <- NULL
      }
    }
  }
  # Bivariate
  if (length(InterKnotsList_biv) != 0) {
    for (learner_name in names(InterKnotsList_biv)) {
      learner_vars <- bivariate_learners[[learner_name]]$variables
      learner_knots <- list()
      for (i in seq_along(learner_vars)) {
        var_name <- learner_vars[i]
        if (i == 1 && !is.null(InterKnotsList_biv[[learner_name]]$ikX) &&
            length(InterKnotsList_biv[[learner_name]]$ikX) >= n - 1) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikX,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else if (i == 2 && !is.null(InterKnotsList_biv[[learner_name]]$ikY) &&
                   length(InterKnotsList_biv[[learner_name]]$ikY) >= n - 1) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikY,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else {
          learner_knots[[var_name]] <- NULL
        }
      }
      polyknots_list[[learner_name]] <- learner_knots
    }
  }
  
  resid <- Y - predicted
  df <- tmp$df.residual
  sigma_hat <- sqrt(sum(resid^2)/df)
  prob <- 1 - 0.5 * (1 - prob)
  # CI_j =\hat{y_j} ± t_{α/2,df}*\hat{σ}*\sqrt{H_{jj}}; H = X(X'X)^{−1}X'
  H_diag <- stats::hat(basisMatrix2, intercept = FALSE) # or influence(tmp)$hat
  band <- qt(prob,df) * sigma_hat * H_diag^.5
  
  # Huang (2003) method for confidence band width (see Theorem 6.1)
  n_obs <- length(Y)
  dim_threshold <- 1500
  if (NCOL(full_matrix) != 0 && n_obs < dim_threshold) {
    
    # i. E_n[B(X)B^t(X)] = (1/n)*\sum_{i=1}^nB(X_i)B^t(X_i)
    matcb <- crossprod(full_matrix) / n_obs
    matcbinv <- tryCatch({
      chol2inv(chol(matcb))  # Primary method: Fastest for SPD matrices
    }, error = function(e1) {
      message("SplineReg_LM_Multivar: Variance matrix for computing asymptotic confidence intervals is not SPD; falling back to solve().")
      tryCatch({
        solve(matcb)  # Secondary method: Standard inverse
      }, error = function(e2) {
        message("SplineReg_LM_Multivar: Variance matrix for computing asymptotic confidence intervals is singular; falling back to ginv().")
        MASS::ginv(matcb)  # Final fallback
      })
    })
    
    # ii. Var(\hat{f} | X) = (1/n)*B^t(x) * E_n[B(X)B^t(X)]^-1 * B(x) * \hat{σ}^2
    S <- full_matrix %*% matcbinv
    conditionalVariance <- (sigma_hat^2 / n_obs) * rowSums(S * full_matrix)
    
    # iii. ± z_{1-α/2} * Var(\hat{f} | X)
    band_width_huang <- qnorm(prob) * sqrt(conditionalVariance)
    
    Upp = predicted + band_width_huang
    Low = predicted - band_width_huang
    
  } else {
    Upp = NULL
    Low = NULL
  }
  
  out <- list(Theta = theta, Predicted = predicted, Residuals = resid, 
              RSS = as.numeric(crossprod(resid)), NCI = list(Upp = predicted + 
                                                     band, Low = predicted - band), Basis = full_matrix, 
              Polygon = list(Kn = polyknots_list, Thetas = theta[1:NCOL(full_matrix)]), 
              temporary = tmp, ACI = list(Upp = Upp, 
                                          Low = Low))
  return(out)
  
}

############################
## SplineReg_GLM_Multivar ##
############################
SplineReg_GLM_Multivar <- function(X, Y, Z = NULL, offset = rep(0, NROW(Y)), base_learners,
                                   weights = rep(1, NROW(Y)), InterKnotsList, n, extrList = lapply(X, range),
                                   family, mustart = NULL, inits = NULL, coefficients, linear_intercept = FALSE,
                                   de_mean = FALSE)
{
  n <- as.integer(n)
  # If boosting family save some relevant functions before converting to glm familt
  if (inherits(family, "boost_family_glm") || inherits(family, "boost_family")) {
    family_name <- family@name
    family_linkinv <- family@response
    family <- get_mboost_family(family_name)
  } else {
    family_name <- family$family
    family_linkinv <- family$linkinv
  }
  
  InterKnotsList_univ <- list()
  for (bl in names(InterKnotsList)) {
    # Check if the length of the variables is equal to 1
    if (length(base_learners[[bl]]$variables) == 1) {
      InterKnotsList_univ[bl] <- InterKnotsList[bl]
    }
  }
  InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
  
  # Select GeDS base-learners
  if (length(base_learners) != 0) {
    base_learners =  base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
    univariate_learners <- bivariate_learners <- NULL
    # Univariate
    if (length(InterKnotsList_univ) != 0){
      # Create a list to store individual design matrices
      univariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 1]
      univariate_vars <- sapply(univariate_learners, function(bl) bl$variables)
      X_univ <- X[, univariate_vars, drop = FALSE]
      extrList = lapply(X_univ, range)
      matrices_univ_list <- vector("list", length = ncol(X_univ))
      # Generate design matrices for each predictor
      for (j in 1:ncol(X_univ)) {
        matrices_univ_list[[j]] <- splineDesign(knots = sort(c(InterKnotsList_univ[[j]], rep(extrList[[j]], n))),
                                                derivs = rep(0, length(X_univ[,j])), x = X_univ[,j], ord = n, outer.ok = TRUE)
      }
      names(matrices_univ_list) <- names(univariate_learners)
      
      # Assign base-learner name to the columns of each of the matrices
      for (matrix_name in names(matrices_univ_list)) {
        num_cols <- ncol(matrices_univ_list[[matrix_name]])
        col_names <- paste(matrix_name, 1:num_cols, sep = "_")
        colnames(matrices_univ_list[[matrix_name]]) <- col_names
      }
      
    } else {
      matrices_univ_list <- NULL
    }
    # Bivariate
    if (length(InterKnotsList_biv) != 0) {
      bivariate_learners <- base_learners[sapply(base_learners, function(bl) length(bl$variables)) == 2]
      matrices_biv_list <- list()
      matrices_biv_list_aux <- list()
      
      for(learner_name in names(bivariate_learners)){
        vars <- bivariate_learners[[learner_name]]$variables
        X_biv <- X[, vars, drop = FALSE]
        Xextr = range(X_biv[,1])
        Yextr = range(X_biv[,2])
        knots <- InterKnotsList_biv[[learner_name]]
        basisMatrixX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                                     x=X_biv[,1],ord=n,outer.ok = TRUE)
        basisMatrixY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                                     x=X_biv[,2],ord=n,outer.ok = TRUE)
        
        # To help saving control polygon knots afterwards
        matrices_biv_list_aux[[learner_name]] <- list(basisMatrixX, basisMatrixY)
        names(matrices_biv_list_aux[[learner_name]]) <- vars
        
        matrices_biv_list[[learner_name]] <- tensorProd(basisMatrixX, basisMatrixY)
        
        # Assign base-learner name to the columns of each of the matrices
        for (matrix_name in names(matrices_biv_list)) {
          num_cols <- ncol(matrices_biv_list[[matrix_name]])
          col_names <- paste(matrix_name, 1:num_cols, sep = "_")
          colnames(matrices_biv_list[[matrix_name]]) <- col_names
        }
        
      }
    } else {
      matrices_biv_list <- NULL
    }
  } else {
    matrices_univ_list <- matrices_biv_list <- NULL
  }
  
  # Combine all matrices side-by-side
  matrices_list <- c(matrices_univ_list, matrices_biv_list)
  if (!is.null(matrices_list) && length(matrices_list) > 0) {
    full_matrix <- do.call(cbind, matrices_list)
  } else {
    full_matrix <- matrix(ncol = 0, nrow = nrow(Z))
  }
  
  # Convert any factor columns in Z to dummy variables
  if (!is.null(Z) && NCOL(Z) > 0) {
    Z <- model.matrix(~ ., data = Z)
    if (!linear_intercept) Z <-  Z[, colnames(Z) != "(Intercept)"]
    } else {
      Z <- NULL
    }
  
  basisMatrix2 <- cbind(full_matrix, Z)
  
  # 1) If coefficients are NOT provided estimate the corresponding regression model
  if (is.null(coefficients)) {
    
    # Initialization
    if (missing(mustart)||is.null(mustart)) {
      env <- parent.frame()
      if (is.null(inits)) {
        temp_env <- list2env(list(y = Y, nobs = length(Y),
                                  etastart = NULL, start = NULL, mustart = NULL ))
        eval(family$initialize, envir = temp_env)
        mustart <- temp_env$mustart
      } else {
        if(length(inits)!= NCOL(basisMatrix2)) stop("'inits' must be of length length(InterKnots) + n + NCOL(Z)")
        mustart <- family$linkinv(basisMatrix2 %*% inits)
      }
    }
    
    # tmp <- IRLSfit(basisMatrix2, Y, offset = offset,
    #                family=family, mustart = mustart, weights = weights)
    # tmp <- glm.fit(basisMatrix2, Y, family = family,
    #                weights = weights, mustart = mustart)
    tmp <- tryCatch({
      glm(Y ~ -1 + basisMatrix2, family = family, weights = weights, mustart = mustart)
      }, error = function(e) {
        # If glm throws an error, fallback to IRLSfit
        IRLSfit(basisMatrix2, Y, offset = offset, family = family, mustart = mustart, weights = weights)
        })
    
    # Extract fitted coefficients
    theta <- coef(tmp)
    names(theta) <- sub("basisMatrix2", "", names(theta))
    # Compute predicted mean values of the response variable
    # predicted <- family$linkinv(basisMatrix2 %*% theta + offset)
    predicted <- family$linkinv(tmp$linear.predictors + offset)
    
    # 2) If coefficients are provided, use them to compute predicted values directly
  } else {
    theta <- coefficients
    linear.predictors <- if (de_mean) basisMatrix2 %*% theta - mean(basisMatrix2 %*% theta) else basisMatrix2 %*% theta # to recover backfitting predictions need de_mean
    fitted.values <- family_linkinv(linear.predictors)
    predicted <- family_linkinv(linear.predictors + offset)
    
    # Manually reconstruct glm object to facilitate S3 methods application
    form <- Y ~ -1 + basisMatrix2
    mf   <- model.frame(form, weights = as.numeric(weights))
    mm    <- model.matrix(form, data = mf)
    tmp <- list(
      coefficients = setNames(theta, colnames(mm)),
      residuals = NULL, # working IWLS residuals
      fitted.values = fitted.values,
      null.deviance = sum(family$dev.resids(Y, family$linkinv(offset), weights)), # Y ~ -1: no intercept
      weights = as.numeric(weights),
      rank = qr(mm)$rank,  # or qr(basisMatrix2)$rank or qr(mm)$rank or rankMatrix(basisMatrix2) or length(coefficients)
      qr = qr(mm), # or qr(basisMatrix2)
      family = family,
      linear.predictors = linear.predictors,
      deviance =  sum(family$dev.resids(Y, fitted.values, weights)),
      prior.weights = weights, # the weights initially supplied
      df.residual = as.numeric(nrow(basisMatrix2) - rankMatrix(basisMatrix2)), # residual degrees of freedom
      y = Y,
      terms = terms(mf),
      model = mf
    )
    logLik_val <- -0.5 * tmp$deviance
    tmp$aic <- -2 * logLik_val + 2 * length(coefficients)
    class(tmp) <- "glm"
  }
  
  # Control polygon knots
  polyknots_list <- list()
  # Univariate 
  if (length(InterKnotsList_univ) != 0) {
    for (learner_name in names(InterKnotsList_univ)) {
      if (!is.null(InterKnotsList_univ[[learner_name]]) &&
          length(InterKnotsList_univ[[learner_name]]) >= n - 1) {
        polyknots_list[[learner_name]] <- list(
          learner_name = makenewknots(
            sort(c(InterKnotsList_univ[[learner_name]],
                   rep(extrList[[univariate_learners[[learner_name]]$variables]], n)))[-c(1, NCOL(matrices_univ_list[[learner_name]]) + 1)],
            degree = n
          )
        )
      } else {
        polyknots_list[[learner_name]] <- NULL
      }
    }
  }
  # Bivariate
  if (length(InterKnotsList_biv) != 0) {
    for (learner_name in names(InterKnotsList_biv)) {
      learner_vars <- bivariate_learners[[learner_name]]$variables
      learner_knots <- list()
      for (i in seq_along(learner_vars)) {
        var_name <- learner_vars[i]
        if (i == 1 && !is.null(InterKnotsList_biv[[learner_name]]$ikX) &&
            length(InterKnotsList_biv[[learner_name]]$ikX) >= n - 1) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikX,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else if (i == 2 && !is.null(InterKnotsList_biv[[learner_name]]$ikY) &&
                   length(InterKnotsList_biv[[learner_name]]$ikY) >= n - 1) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikY,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else {
          learner_knots[[var_name]] <- NULL
        }
      }
      polyknots_list[[learner_name]] <- learner_knots
    }
  }
  
  resid <- tmp$residuals
  deviance <- tmp$deviance
  
  out <- list(Theta = theta, Predicted = predicted,
              Residuals = resid, RSS = deviance,
              Basis = full_matrix,
              Polygon = list(Kn = polyknots_list, Thetas = theta[1:NCOL(full_matrix)]),
              temporary = tmp, deviance = deviance)
  return(out)
}





