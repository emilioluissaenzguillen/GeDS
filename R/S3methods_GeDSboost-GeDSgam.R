################################################################################
################################# COEFFICIENTS #################################
################################################################################
#' @noRd
coef.GeDSboost_GeDSgam <- function(object, n = 3L, ...)
  {
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")
  
  # Check if object is of class "GeDSboost" or "GeDSgam"
  if (!inherits(object, "GeDSboost") && !inherits(object, "GeDSgam")) {
    stop("The input 'object' must be of class 'GeDSboost' or 'GeDSgam'")
  }
  
  # Check if order was wrongly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  model <- object$final_model
  base_learners <- object$args$base_learners
  
  if(n == 2L){
    
    # Piecewise Polynomial coefficients of univariate learners
    univariate_bl <- Filter(function(bl) length(bl$variables) == 1, base_learners)
    univ_coeff <- lapply(model$base_learners[names(univariate_bl)],
                         function(bl) bl$coefficients)
    names(univ_coeff) <- names(univariate_bl)
    
    
    # B-spline coefficients of bivariate learners
    bivariate_bl <- Filter(function(bl) length(bl$variables) == 2, base_learners)
    if (inherits(object, "GeDSboost")) {
      biv_coeff <- lapply(model$base_learners[names(bivariate_bl)],
                          function(bl) lapply(bl$iterations, function(x) x$coef))
      names(biv_coeff) <- names(bivariate_bl)
      } else if (inherits(object, "GeDSgam")) {
        biv_coeff <- lapply(model$base_learners[names(bivariate_bl)],
                            function(bl) bl$coefficients)
      }
    
    coefficients <- c(univ_coeff, biv_coeff)
  }
  if(n == 3L){
    coefficients <- model$Quadratic.Fit$Theta
  }
  if(n == 4L){
    coefficients <- model$Cubic.Fit$Theta
  }
  
  return(coefficients)
}

#' @title Coef method for GeDSboost, GeDSgam
#' @name coef.GeDSboost,gam
#' @description
#' Methods for the functions \code{\link[stats]{coef}} and
#' \code{\link[stats]{coefficients}} that allow to extract the estimated
#' coefficients of \code{\link{GeDSboost-Class}} or \code{\link{GeDSgam-Class}}
#' object.
#' @param object the  \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-Class}} object from which the coefficients should be
#' extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the FGB-GeDS/GAM-GeDS fit whose coefficients should be
#' extracted. If \code{n = 2L} piecewise polynomial coefficients of the
#' univariate GeDS base-learners are provided. In the case of bivariate GeDS
#' base learners and \code{class(object) == "GeDSboost"}, the B-spline
#' coefficients for each iteration where a particular base-learner was selected
#' are provided. In the case of bivariate base learners and
#' \code{class(object) == "GeDSgam"}, the final local-scoring B-spline
#' coefficients for each base-learner are provided. If \code{n = 3L} or
#' \code{n = 4L} B-spline coefficients are provided. By default \code{n} is
#' equal to \code{3L}. Non-integer values will be passed to the function
#' \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the
#' generic function). They will be ignored, but with a warning.
#' 
#' @return
#' A named vector containing the required coefficients of the fitted
#' multivariate predictor model.
#' 
#' @aliases coef.GeDSboost, coef.GeDSgam
#' @seealso \code{\link[stats]{coef}} for the standard definition;
#' \code{\link{NGeDSboost}} and \code{\link{NGeDSgam}} for examples.

#' @export
#' @rdname coef.GeDSboost_GeDSgam
coef.GeDSboost <- coef.GeDSboost_GeDSgam

#' @export
#' @rdname coef.GeDSboost_GeDSgam
coefficients.GeDSboost <- coef.GeDSboost_GeDSgam

#' @export
#' @rdname coef.GeDSboost_GeDSgam
coef.GeDSgam <- coef.GeDSboost_GeDSgam

#' @export
#' @rdname coef.GeDSboost_GeDSgam
coefficients.GeDSgam <- coef.GeDSboost_GeDSgam


################################################################################
################################### DEVIANCE ###################################
################################################################################
#' @noRd
deviance.GeDSboost_GeDSgam <- function(object, n = 3L, ...)
  {
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object','newdata', and 'n' arguments will be considered")
  
  # Check if object is of class "GeDSboost" or "GeDSgam"
  if (!inherits(object, "GeDSboost") && !inherits(object, "GeDSgam")) {
    stop("The input 'object' must be of class 'GeDSboost' or 'GeDSgam'")
  }
  
  # Check if order was wrongly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  model <- object$final_model
  
  if(n == 2L){
    dev <- model$DEV
  }
  if(n == 3L){
    dev <- model$Quadratic.Fit$RSS
  }
  if(n == 4L){
    dev <- model$Cubic.Fit$RSS
  }
  dev <- as.numeric(dev)
  return(dev)
}

#' @export
#' @rdname deviance.GeDS
deviance.GeDSboost <- deviance.GeDSboost_GeDSgam

#' @export
#' @rdname deviance.GeDS
deviance.GeDSgam <- deviance.GeDSboost_GeDSgam


################################################################################
##################################### KNOTS ####################################
################################################################################
#' @noRd
knots.GeDSboost_GeDSgam <-  function(Fn, n = 3L,
                                     options = c("all","internal"), ...)
  {
  # Handle additional arguments
  if(!missing(...)) warning("Only 'Fn','n', and 'options' arguments will be considered")
  
  # Check if Fn is of class "GeDSboost" or "GeDSgam"
  if (!inherits(Fn, "GeDSboost") && !inherits(Fn, "GeDSgam")) {
    stop("The input 'Fn' must be of class 'GeDSboost' or 'GeDSgam'")
  }
  
  # Check if order was wrongly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # Ensure that the options argument is one of the allowed choices ("all" or "internal")
  options <- match.arg(options)
  
  if(n == 2L){
    kn <- Fn$internal_knots$linear.int.knots
  }
  if(n == 3L){
    kn <- Fn$internal_knots$quadratic.int.knots
  }
  if(n == 4L){
    kn <- Fn$internal_knots$cubic.int.knots
  }
  
  # Add the n left and right most knots (assumed to be coalescent)
  if(options=="all") {
    base_learners <- Fn$args$base_learners
    GeDS_learners <- base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
    GeDS_vars <- unlist(sapply(GeDS_learners, function(bl) bl$variables))
    X_GeDS <- Fn$args$predictors[, GeDS_vars, drop = FALSE]
    
    extr <- rbind(apply(X_GeDS, 2, min), apply(X_GeDS, 2, max))
    
    # Univariate GeDS
    univariate_GeDS_learners <- GeDS_learners[sapply(GeDS_learners, function(bl) length(bl$variables)) == 1]
    for (bl_name in names(univariate_GeDS_learners)) {
      bl <- univariate_GeDS_learners[[bl_name]]
      variables <- bl$variables 
      kn[[bl_name]] <- sort(c(rep(extr[,variables], n), kn[[bl_name]]))
    }
    # Bivariate GeDS
    bivariate_GeDS_learners  <- GeDS_learners[sapply(GeDS_learners, function(bl) length(bl$variables)) == 2]
    for (bl_name in names(bivariate_GeDS_learners)) {
      bl <- bivariate_GeDS_learners[[bl_name]]
      variables <- bl$variables 
      kn[[bl_name]]$ikX <- sort(c(rep(extr[,variables[1]], n), kn[[bl_name]]$ikX))
      kn[[bl_name]]$ikY <- sort(c(rep(extr[,variables[2]], n), kn[[bl_name]]$ikY))
    }
  }
  
  return(kn)
}

#' @export
#' @rdname knots
knots.GeDSboost <- knots.GeDSboost_GeDSgam

#' @export
#' @rdname knots
knots.GeDSgam <- knots.GeDSboost_GeDSgam


################################################################################
#################################### PREDICT ###################################
################################################################################
#' @noRd
predict.GeDSboost_GeDSgam <- function(object, newdata, n = 2L, ...)
{
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object','newdata', and 'n' arguments will be considered")
  
  # Check if object is of class "GeDSboost" or "GeDSgam"
  if (!inherits(object, "GeDSboost") && !inherits(object, "GeDSgam")) {
    stop("The input 'object' must be of class 'GeDSboost' or 'GeDSgam'")
  }
  
  # Check if order was wrongly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  newdata <- as.data.frame(newdata)
  # Check whether newdata includes all the necessary predictors
  if (!all(names(object$args$predictors) %in% colnames(newdata))) {
    missing_vars <- setdiff(names(object$args$predictors), colnames(newdata))
    stop(paste("The following predictors are missing in newdata:", paste(missing_vars, collapse = ", ")))
  }
  
  nobs <- nrow(newdata)
  model <- object$final_model
  base_learners <- object$args$base_learners
  offset <- rep(0, nobs)
  
  ###############
  ## 1. LINEAR ##
  ###############
  if (n == 2) {
    ####################
    ## 1.1. GeDSboost ##
    ####################
    if (inherits(object, "GeDSboost")) {
      # Extract predictor variables from newdata
      pred_vars <- newdata[, names(object$args$predictors)]
      # family and shrinkage
      family_name <- object$args$family@name; shrinkage <- object$args$shrinkage
      
      # 1.1 Extract linear base learners
      lin_bl <- base_learners[sapply(base_learners, function(x) x$type == "linear")]
      # 1.2 Extract univariate base learners
      univariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 1)]
      # 1.3 Extract bivariate base learners
      bivariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 2)]
      
      # Normalized model
      if(object$args$normalize_data) {
        # (i) Normalized non-binary response
        if(family_name != "Negative Binomial Likelihood (logit link)") {
          pred_vars <- data.frame(scale(pred_vars)); Y_mean <- object$args$Y_mean; Y_sd <- object$args$Y_sd
          # (ii) Normalized binary response
        } else if (family_name == "Negative Binomial Likelihood (logit link)") {
          pred_vars <- data.frame(scale(pred_vars))
        }
      }
      
      # Initialize prediction vectors
      pred0 <- pred_lin <- pred_univ <- pred_biv <- numeric(nobs)
      
      # 1.0 Offset initial learner
      if (!object$args$initial_learner) {
        pred0 <- rep(mean(object$models$model0$Y_hat), nrow(newdata))
        if(object$args$normalize_data && family_name != "Negative Binomial Likelihood (logit link)") {pred0 <- (pred0-Y_mean)/Y_sd}
      }
      # 1.1 Linear base learners
      if (length(lin_bl)!=0) {
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
            # (ii) Continuous
          } else {
            pred_bl <- coeffs$b0 + coeffs$b1 * pred_vars[[bl]]
          }
          # Add the result to pred_lin
          pred_lin <- pred_lin + pred_bl
        }
      }
      # 1.2 GeDS univariate base learners
      if (length(univariate_bl)!=0) {
        pred_univ <- piecewise_multivar_linear_model(X = pred_vars, model, base_learners = univariate_bl)
      }
      # 1.3 GeDS bivariate base learners
      if (length(bivariate_bl) != 0) {
        pred_biv <- bivariate_bl_linear_model(pred_vars = pred_vars, model, shrinkage, base_learners = bivariate_bl)
      }
      
      if (object$args$normalize_data && family_name != "Negative Binomial Likelihood (logit link)") {
        pred <- (pred0 + pred_lin + pred_univ + pred_biv)*Y_sd + Y_mean
      } else {
        pred <- pred0 + pred_lin + pred_univ + pred_biv
      }
      
      pred <- object$args$family@response(as.numeric(pred))
      return(pred)
      
      ##################
      ## 1.2. GeDSgam ##
      ##################
    } else if (inherits(object, "GeDSgam")) {
      # Extract predictor variables from newdata
      pred_vars <- newdata[, names(object$args$predictors)]
      # Base learners and family
      base_learners <- object$args$base_learners; family_name <- object$args$family$family
      
      # 1.1 Extract linear base learners
      lin_bl <- base_learners[sapply(base_learners, function(x) x$type == "linear")]
      # 1.2 Extract univariate base learners
      univariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 1)]
      # 1.3 Extract bivariate base learners
      bivariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 2)]
      
      # Normalized model
      if (object$args$normalize_data) {
        # (i) Normalized non-binary response
        if (family_name != "binomial") {
          numeric_predictors <- names(pred_vars)[sapply(pred_vars, is.numeric)]
          pred_vars[numeric_predictors] <- data.frame(scale(pred_vars[numeric_predictors]))
          Y_mean <- object$args$Y_mean; Y_sd <- object$args$Y_sd
          # (ii) Normalized binary response
        } else if (family_name == "Negative Binomial Likelihood (logit link)") {
          numeric_predictors <- names(pred_vars)[sapply(pred_vars, is.numeric)]
          pred_vars[numeric_predictors] <- data.frame(scale(pred_vars[numeric_predictors]))
        }
      }
      
      # Initialize prediction vectors
      pred0 <- pred_lin <- pred_univ <- pred_biv <- numeric(nobs)
      
      # 1.0 Initial learner
      pred0 <- rep(mean(model$Y_hat$z), nrow(newdata))
      # 1.1 Linear base learners
      if (length(lin_bl)!=0) {
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
            # (ii) Continuous
          } else {
            pred_bl <- coeffs$b0 + coeffs$b1 * pred_vars[[bl]]
          }
          # Add the result to pred_lin
          pred_lin <- pred_lin + pred_bl
        }
        pred_lin <- pred_lin - mean(pred_lin)
      }
      # 1.2 GeDS univariate base learners
      if (length(univariate_bl)!= 0) {
        pred_univ <- piecewise_multivar_linear_model(pred_vars, model, base_learners = univariate_bl)
        pred_univ <- pred_univ - mean(pred_univ)
      }
      # 1.3 GeDS bivariate base learners
      if (length(bivariate_bl) != 0){
        pred_biv <- bivariate_bl_linear_model(pred_vars, model, base_learners = bivariate_bl, type = "gam")
        pred_biv <- pred_biv - mean(pred_biv)
      }
      
      if(object$args$normalize_data && family_name != "binomial") {
        pred <- (pred0 + pred_lin + pred_univ + pred_biv)*Y_sd + Y_mean
      } else {
        pred <- pred0 + pred_lin + pred_univ + pred_biv
      }
      pred <- object$args$family$linkinv(pred)
      
      return(as.numeric(pred))
    }
    
    ##################
    ## 2. Higher Order ##
    ##################
  } else if (n==3 || n==4) {
    
    # Extract family from GeDSboost/GeDSgam object
    if (inherits(object, "GeDSboost")) {
      family <- get_mboost_family(object$args$family@name)
    } else {
      family <- object$args$family
    }
    
    GeDS_variables <- lapply(base_learners, function(x) {if(x$type == "GeDS") return(x$variables) else return(NULL)})
    GeDS_variables <- unname(unlist(GeDS_variables))
    linear_variables <- lapply(base_learners, function(x) {if(x$type == "linear") return(x$variables) else return(NULL)})
    linear_variables <- unname(unlist(linear_variables))
    
    X = newdata[GeDS_variables]
    Z = newdata[linear_variables]
    
    if (n == 3) {
      if (is.null(model$Quadratic.Fit)) {
        cat("No Quadratic Fit to compute predictions.\n")
        return(NULL)
      }
      int.knots <- "quadratic.int.knots"
      Fit <- "Quadratic.Fit"
    } else if (n == 4) {
      if (is.null(model$Cubic.Fit)) {
        cat("No Cubic Fit to compute predictions.\n")
        return(NULL)
      }
      int.knots <- "cubic.int.knots"
      Fit <- "Cubic.Fit"
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
    InterKnotsList_univ <- InterKnotsList[sapply(InterKnotsList, is.atomic)]
    InterKnotsList_biv <- InterKnotsList[!names(InterKnotsList) %in% names(InterKnotsList_univ)]
    
    # Select GeDS base-learners
    base_learners =  base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
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
        Xextr = range(X_biv[,1])
        Yextr = range(X_biv[,2])
        knots <- InterKnotsList_biv[[learner_name]]
        matriceX <- splineDesign(knots=sort(c(knots$ikX,rep(Xextr,n))), derivs=rep(0,length(X_biv[,1])),
                                 x=X_biv[,1],ord=n,outer.ok = TRUE)
        matriceY <- splineDesign(knots=sort(c(knots$ikY,rep(Yextr,n))),derivs=rep(0,length(X_biv[,2])),
                                 x=X_biv[,2],ord=n,outer.ok = TRUE)
        matriceY_noint <- cut_int(matriceY)
        matrices_biv_list[[learner_name]] <- tensorProd(matriceX,matriceY_noint)
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
    if (NCOL(Z) != 0) Z <- model.matrix(~ . - 1, data = Z)
    
    matrice2 <- cbind(full_matrix, as.matrix(Z))
    
    # # Alternative 1 to compute predictions (required @importFrom stats predict)
    # tmp <- model[[Fit]]$Fit
    # # Set the environment of the model's terms to the current environment for variable access
    # environment(tmp$terms) <- environment()
    # pred <- predict(tmp, newdata=data.frame(matrice2), type = "response")
    
    # Alternative 2 to compute predictions
    coefs <- model[[Fit]]$Fit$coefficients
    coefs[is.na(coefs)] <- 0
    pred <- matrice2%*%coefs+offset
    pred <- family$linkinv(as.numeric(pred))
    return(pred)
  }
}

#' @title Predict method for GeDSboost, GeDSgam
#' @name predict.GeDSboost,gam
#' @description 
#' This method computes predictions from GeDSboost and GeDSgam objects. 
#' It is designed to be user-friendly and accommodate different orders of the
#' GeDSboost or GeDSgam fits.
#' @param object The \code{\link{GeDSboost-class}} or
#' \code{\link{GeDSgam-class}} object.
#' @param newdata An optional data frame for prediction.
#' @param n The order of the GeDS fit (2 for linear, 3 for quadratic, and 4 for
#' cubic). 
#' Default is 2.
#' @param ... potentially further arguments.
#' 
#' @return Numeric vector of predictions (vector of means).
#' @aliases predict.GeDSboost, predict.GeDSgam

#' @export
#' @rdname predict.GeDSboost_GeDSgam
predict.GeDSboost <- predict.GeDSboost_GeDSgam

#' @export
#' @rdname predict.GeDSboost_GeDSgam
predict.GeDSgam <- predict.GeDSboost_GeDSgam


################################################################################
##################################### PRINT ####################################
################################################################################
#' @noRd 
print.GeDSboost_GeDSgam <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    ...)
  {
  cat("\nCall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  kn <- knots(x,n=2,options="int")
  DEVs <- numeric(3)
  names(DEVs) <- c("Order 2","Order 3","Order 4")
  
  # Iterate over each base learner
  for (bl_name in names(x$internal_knots$linear.int.knots)) {
    # Get the linear internal knots
    int.knt <- x$internal_knots$linear.int.knots[[bl_name]]
    
    # 1) No knots or linear base-learner
    if (is.null(int.knt)) {
      cat(paste0("Number of internal knots of the second order spline for '",
                 bl_name, "': 0\n"))
    } else {
      # 2) Bivariate GeDS
      if (is.list(int.knt)) {
        for (component_name in names(int.knt)) {
          component <- int.knt[[component_name]]
          cat(paste0("Number of internal knots of the second order spline for '",
                     bl_name, "', ", component_name, ": ", length(component), "\n"))
        }
        # 3) Univariate GeDS
      } else { 
        cat(paste0("Number of internal knots of the second order spline for '",
                   bl_name, "': ", length(int.knt), "\n"))
      }
    }
  }
  cat("\n")
  DEVs[1] <- deviance(x, n=2L)
  DEVs[2] <- if(!is.null(deviance(x, n=3L))) deviance(x, n=3L) else NA
  DEVs[3] <- if(!is.null(deviance(x, n=4L))) deviance(x, n=4L) else NA
  cat("Deviances:\n")
  print.default(format(DEVs, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  print <- list("Nknots" = length(kn), "Deviances" = DEVs, "Call" = x$extcall)
  
  x$print <- print
  invisible(x)
}

#' @rdname print.GeDS
#' @export
print.GeDSboost <- print.GeDSboost_GeDSgam

#' @rdname print.GeDS
#' @export
print.GeDSgam <- print.GeDSboost_GeDSgam


################################################################################
################################################################################
######################## Visualization/Analysis functions ######################
################################################################################
################################################################################

######################################
##### Fitting process plotting #######
######################################

# Helper function to split the text into lines
split_into_lines <- function(text, max_length)
  {
  words <- strsplit(text, split = ", ")[[1]]
  lines <- character(0)
  current_line <- character(0)
  for (word in words) {
    if (nchar(paste(paste(current_line, collapse = ", "), word, sep = ", ")) > max_length) {
      lines <- c(lines, paste(current_line, collapse = ", "))
      current_line <- word
    } else {
      current_line <- c(current_line, word)
    }
  }
  lines <- c(lines, paste(current_line, collapse = ", "))
  # Add braces
  lines[1] <- paste("{", lines[1], sep = "")
  lines[length(lines)] <- paste(lines[length(lines)], "}", sep = "")
  lines
}

#' @title Visualize Boosting Iterations
#' @name visualize_boosting
#' @description
#' This function plots the \code{\link{NGeDSboost}} fit to the data at the
#' beginning of a given boosting iteration and then plots the subsequent
#' \code{\link{NGeDS}} fit on the corresponding residual (negative gradient).
#' Note: Applicable only for \code{\link{NGeDSboost}} models with one covariate
#' and \code{family = mboost::Gaussian()}.
#' @param M Numeric, specifies the iteration number.
#' @param object A \code{\link{GeDSboost-Class}} object.
#' 
#' @method visualize_boosting GeDSboost
#' @examples
#' # Load packages
#' library(GeDS)
#' 
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
#' object <- NGeDSboost(Y ~ f(X), data = data, normalize_data = TRUE)
#' 
#' # Plot
#' plot(X, Y, pch=20, col=c("darkgrey"))
#' lines(X, sapply(X, f_1), col = "black", lwd = 2)
#' lines(X, object$predictions$pred_linear, col = "green4", lwd = 2)
#' lines(X, object$predictions$pred_quadratic, col="red", lwd=2)
#' lines(X, object$predictions$pred_cubic, col="purple", lwd=2)
#' legend("topright",
#' legend = c("Order 2 (degree=1)", "Order 3 (degree=2)", "Order 4 (degree=3)"),
#' col = c("green4", "red", "purple"),
#' lty = c(1, 1),
#' lwd = c(2, 2, 2),
#' cex = 0.75,
#' bty="n",
#' bg = "white")
#' # Visualize boosting iterations
#' par(mfrow=c(2,2))
#' visualize_boosting(0, object)
#' visualize_boosting(1, object)
#' par(mfrow=c(1,1))
#' 
#' @export
#' @aliases visualize_boosting visualize_boosting.GeDSboost
#' @rdname visualize_boosting
#' @importFrom graphics mtext

visualize_boosting.GeDSboost <- function(M, object)
  {
  # Check if object is of class "GeDSboost"
  if(!inherits(object, "GeDSboost")) {
    stop("The input 'object' must be of class 'GeDSboost'")
  }
  # Check if object has only one predictor
  if(length(object$args$predictors) > 1) {
    stop("Visualization only available for models with a single predictor")
  }
  # Check if family = Gaussian
  if(object$args$family@name != "Squared Error (Regression)") {
    stop("Visualization only available for family = 'Gaussian'")
  }
  
  Y <- object$args$response[[1]]; X <- object$args$predictors[[1]]
  
  model <- object$models[[paste0("model", M)]]
  Y_hat <- model$Y_hat
  next_model <- object$models[[paste0("model", M+1)]]
  
  # 1. Data plot
  # Plot range
  x_range <- range(X)
  y_range <- range(Y, Y_hat)
  x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
  y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
  
  # Split int knots into lines
  int.knots <- round(as.numeric(model$base_learners[[1]]$linear.int.knots), 2) 
  int.knots_lines <- split_into_lines(paste(int.knots, collapse=", "), 75)
  # Create the title
  title_text_1 <- bquote(atop(plain(Delta[.(M) * ",2"] == .(int.knots_lines[1]))))
  
  # Plot
  plot(X, Y, pch = 20, col = "darkgrey", tck = 0.02, main = "",
       xlim = x_range, ylim = y_range, cex.axis = 2.5, cex.lab = 2.5)
  lines(X, Y_hat, lwd = 2)
  legend("topleft",
         legend = c(2*M),
         bty = "n",
         text.font = 2,
         cex=2.5)
  
  # Add the title
  if(length(int.knots_lines)==1){
    mtext(title_text_1, side = 3, line=-0.5625, cex = 1.65)
  } else if(length(int.knots_lines)==2){
    mtext(title_text_1, side = 3, cex = 1.65)
    mtext(int.knots_lines[2], side = 3, line=0.5, cex = 1.65)
  } else if(length(int.knots_lines)==3){
    mtext(title_text_1, side = 3, line=0.75, cex = 0.85)
    mtext(int.knots_lines[2], side = 3, line = 1.5, cex = 0.85)
    mtext(int.knots_lines[3], side = 3, line = 0.25, cex = 0.85)
  } else if(length(int.knots_lines)==4){
    mtext(title_text_1, side = 3, line=1, cex = 0.85)
    mtext(int.knots_lines[2], side = 3, line = 1.5, cex = 0.85)
    mtext(int.knots_lines[3], side = 3, line = 0.65, cex = 0.85)
    mtext(int.knots_lines[4], side = 3, line = 0, cex = 0.85)
  }
  
  # 2. Residuals plot
  if (M < length(object$models) - 1){
    
    residuals <- Y - model$Y_hat
    
    # Plot range
    x_range <- range(X)
    y_range <- range(residuals, next_model$best_bl$pred_linear)
    x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
    y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
    
    # Split int knots into lines
    int.knots <- round(as.numeric(next_model$best_bl$int.knt), 2)
    int.knots_lines <- split_into_lines(paste(int.knots, collapse=", "), 75)
    # Create the title
    title_text_2 <- bquote(atop(plain(delta[.(M) * ",2"] == .(int.knots_lines[1]))))
    
    plot(X, residuals, pch=20, col=c("darkgrey"), tck = 0.02, main = "",
         xlim = x_range, ylim = y_range, cex.axis = 2.5, cex.lab = 2.5)
    lines(X, next_model$best_bl$pred_linear, col="blue", lty = 2, lwd = 1)
    points(next_model$best_bl$coef$mat[,1], next_model$best_bl$coef$mat[,2], col="blue", pch=21)
    legend("topleft",
           legend = c(2*M+1),
           bty = "n",
           text.font = 2,
           cex=2.5)
    # Add the title
    if(length(int.knots_lines)==1){
      mtext(title_text_2, side = 3, line=-0.5625, cex = 1.65)
    } else if(length(int.knots_lines)==2){
      mtext(title_text_2, side = 3, cex = 1.65)
      mtext(int.knots_lines[2], side = 3, line=0.5, cex = 1.65)
    } else if(length(int.knots_lines)==3){
      mtext(title_text_2, side = 3, line=0.75, cex = 0.85)
      mtext(int.knots_lines[2], side = 3, line = 1.5, cex = 0.85)
      mtext(int.knots_lines[3], side = 3, line = 0.25, cex = 0.85)
    } else if(length(int.knots_lines)==4){
      mtext(title_text_2, side = 3, line=1, cex = 0.85)
      mtext(int.knots_lines[2], side = 3, line = 1.5, cex = 0.85)
      mtext(int.knots_lines[3], side = 3, line = 0.65, cex = 0.85)
      mtext(int.knots_lines[4], side = 3, line = 0, cex = 0.85)
    }
  }
}

#' @export
visualize_boosting <- visualize_boosting.GeDSboost


################################################################################
############################ BASE LEARNER IMPORTANCE ###########################
################################################################################
#' @title Base Learner Importance for GeDSboost objects
#' @name bl_imp.GeDSboost
#' @description
#' This function calculates the in-bag mean squared error (MSE) reduction
#' ascribable to each of the base-learners with regards to the final prediction
#' of the component-wise gradient boosted model encapsulated in a
#' \code{\link{GeDSboost-Class}} object. Essentially, it measures the decrease
#' in MSE attributable to each base-learner for every time it is selected across
#' the boosting iterations, and aggregates them. This provides a measure on how
#' much each base-learner contributes to the overall improvement in the model's
#' accuracy, as reflected by the decrease in MSE. This function is adapted from
#' \code{\link[mboost]{varimp}} and is compatible with the available
#' \code{\link[mboost]{mboost-package}} methods for \code{\link[mboost]{varimp}},
#' including \code{plot}, \code{print} and \code{as.data.frame}.
#' @param object An object of class \code{\link{GeDSboost-Class}}.
#' @param ... potentially further arguments.
#'
#' @return An object of class \code{varimp} with available \code{plot},
#' \code{print} and \code{as.data.frame} methods.
#' @details
#' See \code{\link[mboost]{varimp}} for details.
#' @method bl_imp GeDSboost
#' @examples
#' library(GeDS)
#' library(TH.data)
#' set.seed(290875)
#' data("bodyfat", package = "TH.data")
#' data = bodyfat
#' Gmodboost <- NGeDSboost(formula = DEXfat ~ f(hipcirc) + f(kneebreadth) + f(anthro3a),
#'                         data = data, initial_learner = FALSE)
#' MSE_Gmodboost_linear <- mean((data$DEXfat - Gmodboost$predictions$pred_linear)^2)
#' MSE_Gmodboost_quadratic <- mean((data$DEXfat - Gmodboost$predictions$pred_quadratic)^2)
#' MSE_Gmodboost_cubic <- mean((data$DEXfat - Gmodboost$predictions$pred_cubic)^2)
#'
#' # Print MSE
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#'     "Linear NGeDSboost:", MSE_Gmodboost_linear, "\n",
#'     "Quadratic NGeDSboost:", MSE_Gmodboost_quadratic, "\n",
#'     "Cubic NGeDSboost:", MSE_Gmodboost_cubic, "\n")
#'
#' # Base Learner Importance
#' bl_imp <- bl_imp(Gmodboost)
#' print(bl_imp)
#' plot(bl_imp)
#' 
#' @export
#' @aliases bl_imp bl_imp.GeDSboost
#' @references
#' Hothorn T., Buehlmann P., Kneib T., Schmid M. and Hofner B. (2022).
#' mboost: Model-Based Boosting. R package version 2.9-7, \url{https://CRAN.R-project.org/package=mboost}.
#' @rdname bl_imp
bl_imp.GeDSboost <- function(object, ...)
  {
  # Check if object is of class "GeDSboost"
  if(!inherits(object, "GeDSboost")) {
    stop("The input 'object' must be of class 'GeDSboost'")
  }
  
  # 1. Variables
  ## Response variable
  response <- names(object$args$response)
  ## Base-learners
  bl_names <- names(object$args$base_learners)
  bl_selected <- sapply(object$models, function(model) model$best_bl$name)
  bl_indices <- match(bl_selected, bl_names)
  
  # 2. In-bag risk
  # Calculate initial risk
  initial_risk <- mean((object$args$response[[response]] - 0)^2)
  # Get the riskdiff of model0
  model0_risk <- mean((object$args$response[[response]] - object$models[[paste0("model", 0)]]$Y_hat)^2)
  first_riskdiff <- initial_risk - model0_risk
  
  # Get the number of models
  n_models <- length(object$models)
  if (!object$args$initial_learner) n_models <- n_models - 1
  # Initialize an empty vector to store the inbag risks
  inbag_risks <- vector(mode = "numeric", length = n_models)
  # Loop over each model from model1 onwards
  for (i in seq_len(n_models)) {
    # Extract the selected predictor for each model
    Y_hat <- object$models[[paste0("model", i)]]$Y_hat
    # Calculate the inbag risk for this iteration
    inbag_risks[i] <- mean((object$args$response[[response]] - Y_hat)^2)
  }
  inbag_risks <- c(model0_risk, inbag_risks)
  
  # 3. Calculate differences in inbag risk between consecutive boosting iterations
  riskdiff <- -1 * diff(inbag_risks)
  # Scale these differences by the number of observations
  riskdiff <- riskdiff / length(object$args$response[[response]])
  # Prepend the first risk difference to the riskdiff vector
  riskdiff <- c(first_riskdiff, riskdiff)
  
  # For offset initial-learner bl_imp is just measured within boosting iterations
  if(!object$args$initial_learner){
    riskdiff <- riskdiff[-1]
    bl_selected$model0 <- NULL
    bl_indices <- bl_indices[-1]
  }
  
  # 4. Sum up riskdiffs
  explained <- sapply(seq_along(bl_names), FUN = function(i) {
    sum(riskdiff[which(bl_indices == i)])
  })
  names(explained) <- bl_names
  class(explained) <- "varimp"
  # Calculate the proportion of models where the i-th learner was selected
  attr(explained, "selfreqs") <- sapply(seq_along(bl_names), 
                                        function(i) {
                                          mean(bl_indices == i)
                                        })
  # To accomodate structure requirements:
  var_names <- bl_names
  var_names <- sapply(strsplit(var_names, ", "), function(x) {
    do.call(function(...) paste(..., sep = ", "), as.list(x[order(x)]))
  })
  var_order <- order(sapply(var_names, function(i) {
    sum(explained[var_names == i])
  }))
  attr(explained, "variable_names") <- ordered(var_names, levels = unique(var_names[var_order]))
  
  return(explained)
}

#' @export
bl_imp <- bl_imp.GeDSboost
