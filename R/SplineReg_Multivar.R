
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
  
  InterKnotsList_univ <- InterKnotsList[sapply(InterKnotsList, is.atomic)]
  InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
  
  # Select GeDS base-learners
  base_learners =  if (length(base_learners) > 0) base_learners[sapply(base_learners, function(x) x$type == "GeDS")] else NULL
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
      matriceX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                               x=X_biv[,1],ord=n,outer.ok = TRUE)
      matriceY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                               x=X_biv[,2],ord=n,outer.ok = TRUE)
      
      # To help saving control polygon knots afterwards
      matrices_biv_list_aux[[learner_name]] <- list(matriceX, matriceY)
      names(matrices_biv_list_aux[[learner_name]]) <- vars
      
      matriceY_noint <- cut_int(matriceY)
      matrices_biv_list[[learner_name]] <- tensorProd(matriceX,matriceY_noint)
      
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
    if (linear_intercept) {
      Z <- model.matrix(~ ., data = Z)
      } else {
        Z <- model.matrix(~ . -1, data = Z)
        }
    } else {
      Z <- NULL
    }
  
  matrice2 <- cbind(full_matrix, Z)
  
  # 1) If coefficients are NOT provided estimate the corresponding regression model
  if (is.null(coefficients)) {
    Y0 <- Y - offset
    tmp <- lm(Y0 ~ -1 + matrice2, weights=as.numeric(weights))
    # the ‘-1’ serving to suppress the redundant extra intercept that would be added by default
    # 'splineDesign' already includes a basis that accounts for the intercept
    theta <- coef(tmp)
    names(theta) <- sub("matrice2", "", names(theta))
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
    tmp <- NULL
    theta <- coefficients
    f <- if (de_mean) matrice2 %*% theta - mean(matrice2 %*% theta) else matrice2 %*% theta # to recover backfitting predictions need de_mean
    predicted <- f + offset
  }
  
  # Control polygon knots
  polyknots_list <- list()
  # Univariate 
  if (length(InterKnotsList_univ) != 0) {
    for (learner_name in names(InterKnotsList_univ)) {
      if (!is.null(InterKnotsList_univ[[learner_name]]) && length(InterKnotsList_univ[[learner_name]]) >= n - 1) {
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
        if (i == 1 && !is.null(InterKnotsList_biv[[learner_name]]$ikX) && length(InterKnotsList_biv[[learner_name]]$ikX) >= n - 1) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikX,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else if (i == 2 && !is.null(InterKnotsList_biv[[learner_name]]$ikY) && length(InterKnotsList_biv[[learner_name]]$ikY) >= n - 1) {
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
  df <- if(!is.null(tmp)) tmp$df.residual else as.numeric(nrow(matrice2) - rankMatrix(matrice2)) # residual degrees of freedom
  sigma_hat <- sqrt(sum(resid^2)/df)
  prob <- 1 - 0.5 * (1 - prob)
  # CI_j =\hat{y_j} ± t_{α/2,df}*\hat{σ}*\sqrt{H_{jj}}; H = X(X'X)^{−1}X'
  H_diag <- if(!is.null(tmp)) influence(tmp)$hat else stats::hat(matrice2, intercept = FALSE)
  band <- qt(prob,df) * sigma_hat * H_diag^.5
  
  n <- length(Y)
  N <- NCOL(full_matrix)
  if (N != 0 && n < 1500) {
    matcb <- matrix(0, N, N)
    for (i in 1:n) {
      matcb <- matcb + full_matrix[i, ] %*% t(full_matrix[i, ])
      }
    matcb <- matcb / n
    matcbinv <- tryCatch({
      solve(matcb)
      }, error = function(e) {
        # If there's an error with solve(), use ginv() as a fallback
        # Moore-Penrose pseudo-inverse to skip multicolinearity issues that make matcb singular
        message("Warning message in SplineReg_LM_Multivar: Matrix is singular for computing confidence intervals; using ginv() as a fallback.")
        MASS::ginv(matcb)
        })
    
    band_width_huang <- qnorm(prob) * n^(-0.5) * diag(sigma_hat * full_matrix %*% matcbinv %*% t(full_matrix))
    
    Upp = predicted + band_width_huang
    Low = predicted - band_width_huang
    
    } else {
      Upp = NULL
      Low = NULL
    }
  
  out <- list(Fit = tmp, Theta = theta, Predicted = predicted, Residuals = resid, 
              RSS = t(resid) %*% resid, NCI = list(Upp = predicted + 
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
  
  InterKnotsList_univ <- InterKnotsList[sapply(InterKnotsList, is.atomic)]
  InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
  
  # Select GeDS base-learners
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
      matriceX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                               x=X_biv[,1],ord=n,outer.ok = TRUE)
      matriceY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                               x=X_biv[,2],ord=n,outer.ok = TRUE)
      
      # To help saving control polygon knots afterwards
      matrices_biv_list_aux[[learner_name]] <- list(matriceX, matriceY)
      names(matrices_biv_list_aux[[learner_name]]) <- vars
      
      matriceY_noint <- cut_int(matriceY)
      matrices_biv_list[[learner_name]] <- tensorProd(matriceX,matriceY_noint)
      
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
    if (linear_intercept) {
      Z <- model.matrix(~ ., data = Z)
    } else {
      Z <- model.matrix(~ . -1, data = Z)
    }
  } else {
    Z <- NULL
  }
  
  matrice2 <- cbind(full_matrix, Z)
  
  # 1) If coefficients are NOT provided estimate the corresponding regression model
  if (is.null(coefficients)) {
    
    # Since for NGeDSboost(family = "binomial") the encoding is -1/1
    # and for stats::binomial() the encoding is 0/1
    if (family_name == "Negative Binomial Likelihood (logit link)") {
      Y <- (Y + 1) / 2
    }
    
    # Initialization
    if (missing(mustart)||is.null(mustart)) {
      env <- parent.frame()
      if (is.null(inits)) {
        temp_env <- list2env(list(y = Y, nobs = length(Y),
                                  etastart = NULL, start = NULL, mustart = NULL ))
        eval(family$initialize, envir = temp_env)
        mustart <- temp_env$mustart
      } else {
        if(length(inits)!= NCOL(matrice2)) stop("'inits' must be of length length(InterKnots) + n + NCOL(Z)")
        mustart <- family$linkinv(matrice2%*%inits)
      }
    }
    
    # tmp <- IRLSfit(matrice2, Y, offset = offset,
    #                family=family, mustart = mustart, weights = weights)
    # tmp <- glm.fit(matrice2, Y, family = family,
    #                weights = as.numeric(weights), mustart = mustart)
    
    tmp <- glm(Y ~ -1 + matrice2, family = family, weights = weights)
    
    theta <- coef(tmp)
    names(theta) <- sub("matrice2", "", names(theta))
    # predicted <- family$linkinv(matrice2%*%theta + offset)
    predicted <- family$linkinv(tmp$linear.predictors + offset)
    
    # 2) If coefficients are provided, use them to compute predicted values directly
  } else {
    tmp <- NULL
    theta <- coefficients
    f <- if (de_mean) matrice2 %*% theta - mean(matrice2 %*% theta) else matrice2 %*% theta # to recover backfitting predictions need de_mean
    predicted <- family_linkinv(f + offset)
  }
  
  # Control polygon knots
  polyknots_list <- list()
  # Univariate 
  if (length(InterKnotsList_univ) != 0) {
    for (learner_name in names(InterKnotsList_univ)) {
      if (!is.null(InterKnotsList_univ[[learner_name]])){
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
        if (i == 1 && !is.null(InterKnotsList_biv[[learner_name]]$ikX)) {
          learner_knots[[var_name]] <- makenewknots(
            sort(c(InterKnotsList_biv[[learner_name]]$ikX,
                   rep(extrList[[var_name]], n)))[-c(1, NCOL(matrices_biv_list_aux[[learner_name]][[var_name]]) + 1)],
            degree = n)
        } else if (i == 2 && !is.null(InterKnotsList_biv[[learner_name]]$ikY)) {
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
  
  resid <- if(!is.null(tmp)) tmp$residuals else NULL
  deviance <- if(!is.null(tmp)) tmp$deviance else sum(family$dev.resids(Y, predicted, weights))
  
  out <- list(Fit = tmp, Theta = theta, Predicted = predicted,
              Residuals = resid, RSS = deviance,
              Basis = full_matrix,
              Polygon = list(Kn = polyknots_list, Thetas = theta[1:NCOL(full_matrix)]),
              temporary = tmp, deviance = deviance)
  return(out)
}





