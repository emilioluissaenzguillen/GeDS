################################################################################
################################# COEFFICIENTS #################################
################################################################################
#' @noRd
coef.GeDSboost_GeDSgam <- function(object, n = 3L, ...){
  
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
  
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")
  
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
    biv_coeff <- lapply(model$base_learners[names(bivariate_bl)], 
                        function(bl) lapply(bl$iterations, function(x) x$coef))
    names(biv_coeff) <- names(bivariate_bl)
    
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

#' Coef method for GeDSboost and GeDSgam objects
#' 
#' @description Methods for the functions \code{\link[stats]{coef}} and \code{\link[stats]{coefficients}}
#'  that allow to extract the estimated coefficients of \code{\link{GeDSboost-Class}} 
#'  or \code{\link{GeDSgam-Class}} object.
#' @param object the  \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-Class}} object from which
#' the coefficients should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' whose coefficients should be extracted. If \code{n = 2L} piecewise polynomial coefficients of univariate
#' base-learners are provided. In the case of bivariate base learners,
#' the B-spline coefficients for each iteration on which such base-learner was selected are provided.
#' If \code{n = 3L} or \code{n = 4L} B-spline coefficients are provided.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the generic function).
#' They will be ignored, but with a warning.
#'
#' @details
#' See \code{\link{coef.GeDS}} for details.
#'
#' @return
#' A named vector containing the required coefficients of the fitted multivariate predictor model.
#' The coefficient of linear base-learners is named as the base learner itself. The coefficients of GeDS
#' base-learners are named as the base learners themselves, followed by the index of the corresponding B-spline.
#'
#' @seealso \code{\link[stats]{coef}} for the standard definition; \code{\link{NGeDSboost}}
#' and \code{\link{NGeDSgam}} for examples.

#' @rdname coef.GeDSboost_GeDSgam
#' @export
coef.GeDSboost <- coef.GeDSboost_GeDSgam

#' @rdname coef.GeDSboost_GeDSgam
#' @export
coefficients.GeDSboost <- coef.GeDSboost_GeDSgam

#' @rdname coef.GeDSboost_GeDSgam
#' @export
coef.GeDSgam <- coef.GeDSboost_GeDSgam

#' @rdname coef.GeDSboost_GeDSgam
#' @export
coefficients.GeDSgam <- coef.GeDSboost_GeDSgam


################################################################################
################################### DEVIANCE ###################################
################################################################################
#' @noRd
deviance.GeDSboost_GeDSgam <- function(object, n = 3L, ...){
  
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
  
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object','newdata', and 'n' arguments will be considered")
  
  model <- object$final_model
  
  if(n == 2L){
    dev <- model$RSS
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

#' Deviance method for GeDSboost and GeDSgam objects
#' @description Method for the function \code{\link[stats]{deviance}} that allows the user to extract the value
#' of the deviance corresponding to a selected GeDSboost or GeDSgam fit from a
#' \code{\link{GeDSboost-Class}} or \code{\link{GeDSgam-Class}} object.
#' @param object the \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object from which
#' the deviance should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the
#' GeDSboost/GeDSgam fit whose deviance should be extracted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the generic function).
#' They will be ignored, but with a warning.
#'
#' @details This is a method for the function \code{\link[stats]{deviance}}.
#' As \code{\link{GeDSboost-class}} and \code{\link{GeDSgam-class}} objects contain three different
#' fits (linear, quadratic and cubic), it is possible
#' to specify the order of the GeDS fit for which the deviance is required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{deviance}} for the standard definition;
#' @return A numeric value corresponding to the  deviance of the selected GeDSboost/GeDSgam fit.

#' @rdname deviance.GeDSboost_GeDSgam
#' @export
deviance.GeDSboost <- deviance.GeDSboost_GeDSgam

#' @rdname deviance.GeDSboost_GeDSgam
#' @export
deviance.GeDSgam <- deviance.GeDSboost_GeDSgam


################################################################################
##################################### KNOTS ####################################
################################################################################
#' @noRd
knots.GeDSboost_GeDSgam <-  function(Fn, n = 3L, options = c("all","internal"), ...) {
  
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
  
  # Handle additional arguments
  if(!missing(...)) warning("Only 'Fn','n', and 'options' arguments will be considered")
  
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

#' Knots method for GeDSboost and GeDSgam objects
#' @description Method for the generic function \code{\link[stats]{knots}} that allows the user
#' to extract vector of knots of a GeDSboost or GeDSgam fit of a specified order
#' contained in a \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object, respectively.
#'
#' @param Fn the \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object from which the vector of knots for the
#' specified FGB-GeDS or GAM-GeDS fit should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the FGB-GeDS or GAM-GeDS fit
#' whose knots should be extracted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param options a character string specifying whether "\code{all}" knots, including
#' the left-most and the right-most limits of the interval embedding the observations (the default) or
#' only the "\code{internal}" knots should be extracted.
#' @param ... potentially further arguments (required for compatibility with the definition of
#' the generic function). Currently ignored, but with a warning.
#' @details This is a method for the function \code{\link[stats]{knots}} in the \pkg{stats} package.
#'
#' As \code{\link{GeDSboost-class}} and \code{\link{GeDSgam-class}} objects contain three different
#' fits (linear, quadratic and cubic), it is possible
#' to specify the order of the GeDS fit for which the deviance is required via the input argument \code{n}.
#'
#'
#' @seealso \code{\link[stats]{knots}} for the definition of the generic function; \code{\link{NGeDS}} and \code{\link{GGeDS}} for examples.
#'
#' @return A vector in which each element represents a knot of the FGB-GeDS/GAM-GeDS fit of the required order.

#' @rdname knots.GeDSboost_GeDSgam
#' @export
knots.GeDSboost <- knots.GeDSboost_GeDSgam

#' @rdname knots.GeDSboost_GeDSgam
#' @export
knots.GeDSgam <- knots.GeDSboost_GeDSgam


################################################################################
#################################### PREDICT ###################################
################################################################################
#' @noRd
predict.GeDSboost_GeDSgam <- function(object, newdata, n = 2L, ...)
{
  
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
  
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object','newdata', and 'n' arguments will be considered")
  
  
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
  if (n==2) {
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
    ## Higher Order ##
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

#' Predict method for GeDSboost and GeDSgam objects
#'
#' @description This method computes predictions from GeDSboost and GeDSgam objects. 
#' It is designed to be user-friendly and accommodate different orders of the GeDS fit.
#' @param object The \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object.
#' @param newdata An optional data frame for prediction.
#' @param n The order of the GeDS fit (2 for linear, 3 for quadratic, and 4 for cubic). 
#' Default is 2.
#' @param ... potentially further arguments.
#' @return Numeric vector of predictions.

#' @rdname predict.GeDSboost_GeDSgam
#' @export
predict.GeDSboost <- predict.GeDSboost_GeDSgam

#' @rdname predict.GeDSboost_GeDSgam
#' @export
predict.GeDSgam <- predict.GeDSboost_GeDSgam


################################################################################
##################################### PRINT ####################################
################################################################################
#' @noRd 
print.GeDSboost_GeDSgam <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
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
      cat(paste0("Number of internal knots of the second order spline for '", bl_name, "': 0\n"))
    } else {
      # 2) Bivariate GeDS
      if (is.list(int.knt)) {
        for (component_name in names(int.knt)) {
          component <- int.knt[[component_name]]
          cat(paste0("Number of internal knots of the second order spline for '", bl_name, "', ", component_name, ": ", length(component), "\n"))
        }
        # 3) Univariate GeDS
      } else { 
        cat(paste0("Number of internal knots of the second order spline for '", bl_name, "': ", length(int.knt), "\n"))
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

#' Print method for GeDS objects
#' 
#' @description Method for the generic function \code{\link[base]{print}} that allows to
#' print on screen the main information related to the fitted predictor model that can be extracted
#' from a \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object.
#' @param x the \code{\link{GeDSboost-class}} or \code{\link{GeDSgam-class}} object for which the main information should be printed on screen.
#' @param digits number of digits to be printed.
#' @param ... potentially further arguments (required by the definition of the generic function).
#'
#' @details This method allows to print on screen basic information related to the fitted predictor model such as the
#' function \code{call}, the number of internal knots for the linear FGB-GeDS/GAM-GeDS fit and the deviances
#' for the three (linear, quadratic and cubic) fitted predictor models embedded in the \code{\link{GeDSboost-class}}
#' or \code{\link{GeDSgam-class}} object.
#'
#' @seealso \code{\link[base]{print}} for the standard definition.
#' @return This function returns (invisibly) the same input object, but adding the slot \code{Print}
#' that contains the three sub-slots:
#' \item{Nknots}{ the number of internal knots for each base-learner of the linear FGB-GeDS/GAM-GeDS fit}
#' \item{Deviances}{ the deviances of the three (linear, quadratic and cubic) FGB-GeDS/GAM-GeDS fits}
#' \item{Call}{ the \code{call} to the function that produced the \code{x} object}

#' @rdname print.GeDSboost_GeDSgam
#' @export
print.GeDSboost <- print.GeDSboost_GeDSgam

#' @rdname print.GeDSboost_GeDSgam
#' @export
print.GeDSgam <- print.GeDSboost_GeDSgam



#########################################################################################
#########################################################################################
############################ Visualization/Analysis functions ###########################
#########################################################################################
#########################################################################################

######################################
##### Fitting process plotting #######
######################################

# Helper function to split the text into lines
split_into_lines <- function(text, max_length) {
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

#' Visualize Boosting Process
#'
#' @description This function presents  the fit over the data at the start of each boosting iteration and the subsequent fit on the residuals.
#'
#' @param M Numeric, specifies the iteration number.
#' @param Gmodboost An object of class  \code{\link{GeDSboost-Class}}, containing the results of the GeDS boosting model.
#' @importFrom graphics mtext
#' @examples
#' \dontrun{
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
#' Gmodboost <- NGeDSboost(Y ~ f(X), data = data, initial_learner = "GeDS", int.knots_init = 2,
#' internal_knots = 500, q_boost = 2, phi_boost_exit = 0.995, shrinkage = 1,
#' phi = 0.99, normalize_data = TRUE, family = mboost::Gaussian())
#' 
#' # Plot
#' plot(X, Y, pch=20, col=c("darkgrey"))
#' lines(X, sapply(X, f_1), col = "black", lwd = 2)
#' lines(X, Gmodboost$predictions$pred_linear, col = "green4", lwd = 2)
#' lines(X, Gmodboost$predictions$pred_quadratic, col="red", lwd=2)
#' lines(X, Gmodboost$predictions$pred_cubic, col="purple", lwd=2)
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
#' visualize_boosting(0, Gmodboost)
#' visualize_boosting(1, Gmodboost)
#' par(mfrow=c(1,1))
#' }
#'
#' @export

visualize_boosting = function(M, Gmodboost) {
  
  # Check if Gmodboost is of class "GeDSboost"
  if(!inherits(Gmodboost, "GeDSboost")) {
    stop("The input 'Gmodboost' must be of class 'GeDSboost'")
  }
  # Check if Gmodboost has only one predictor
  if(length(Gmodboost$args$predictors) > 1) {
    stop("Visualization only available for models with a single predictor")
  }
  # Check if family = Gaussian
  if(Gmodboost$args$family@name != "Squared Error (Regression)") {
    stop("Visualization only available for family = 'Gaussian'")
  }
  
  Y <- Gmodboost$args$outcome[[1]]; X <- Gmodboost$args$predictors[[1]]
  
  model <- Gmodboost$models[[paste0("model", M)]]
  Y_hat <- model$Y_hat
  next_model <- Gmodboost$models[[paste0("model", M+1)]]
  
  # 1. Data plot
  # Plot range
  x_range <- range(X)
  y_range <- range(Y, Y_hat)
  x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
  y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
  
  # Split int knots into lines
  int.knots <- round(as.numeric(get_internal_knots(model$base_learners[[1]]$knots)), 2)
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
  if (M < length(Gmodboost$models) - 1){
    
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


#' Variable Importance for GeDS Boosting Model
#'
#' This function calculates the importance of each base-learner in the final prediction of the 
#' Functional Gradient Boosting model that uses Geometrically Designed splines (GeDS) as base-learners.
#'
#' @param object An object of class \code{\link{GeDSboost-Class}}, which contains the results of the GeDS boosting model.
#' @param initial_model Logical, indicates whether to include the initial model in the importance calculation in the case of GeDS initial learner. 
#'        If \code{TRUE}, it includes the initial model's risk difference; otherwise, it is set to 0. Default is \code{TRUE}.
#'
#' @return A numeric vector of variable importances for each base-learner, with class "varimp". 
#'         Attributes include "selfreqs", which gives the proportion of models where the i-th learner was selected, 
#'         and "variable_names", which provides the ordered names of the variables.
#'
#' @examples
#' \dontrun{
#' # Load packages
#' library(GeDS)
#' library(TH.data)
#' set.seed(290875)
#' data("bodyfat", package = "TH.data")
#' data = bodyfat
#' Gmodboost <- NGeDSboost(formula = DEXfat ~ f(hipcirc) + f(kneebreadth) + f(anthro3a),
#' data = data, initial_learner = NULL, internal_knots = 500,
#' max_iterations = 100, q_boost = 2, phi_boost_exit = 0.995,
#' shrinkage=1, phi=0.99, normalize_data = FALSE)
#' mse_ngedsboost1 <- mean((data$DEXfat-Gmodboost$predictions$pred_linear)^2)
#' mse_ngedsboost2 <- mean((data$DEXfat-Gmodboost$predictions$pred_quadratic)^2)
#' mse_ngedsboost3 <- mean((data$DEXfat-Gmodboost$predictions$pred_cubic)^2)
#' 
#' # Print
#' cat("\n", "MEAN SQUARED ERROR", "\n",
#' "Linear NGeDSboost:", mse_ngedsboost1, "\n",
#' "Quadratic NGeDSboost:", mse_ngedsboost2, "\n",
#' "Cubic NGeDSboost:", mse_ngedsboost3, "\n")
#' # Variable importance
#' varimp <- varimp_GeDS(Gmodboost)
#' print(varimp)
#' plot(varimp)
#' }
#'
#' @export

######################################
############### varimp ###############
######################################

varimp_GeDS <- function (object, initial_model = TRUE)
{
  
  # Check if Gmodboost is of class "GeDSboost"
  if(!inherits(object, "GeDSboost")) {
    stop("The input 'Gmodboost' must be of class 'NGeDS_boosted'")
  }
  
  # 1. Variables
  ## Response variable
  response <- names(object$args$outcome)
  ## Covariates
  learner_names <- names(object$args$base_learners)
  learner_selected <- sapply(object$models, function(model) model$best_bl$name)
  learner_indices <- match(learner_selected, learner_names)
  
  # 2. In-bag risk
  # Calculate initial risk
  initial_risk <- mean((object$args$outcome[[response]] - 0)^2)
  # Get the riskdiff of model0
  initial_model_risk <- mean((object$args$outcome[[response]] - object$models[[paste0("model", 0)]]$Y_hat)^2)
  first_riskdiff <- initial_risk - initial_model_risk
  
  # Get the number of models
  n_models <- length(object$models)
  # Initialize an empty vector to store the inbag risks
  inbag_risks <- vector(mode = "numeric", length = n_models)
  
  # Loop over each model
  for (i in seq_len(n_models)) {
    # Extract the selected predictor for each model
    Y_hat <- object$models[[paste0("model", i - 1)]]$Y_hat
    # Calculate the inbag risk for this iteration
    inbag_risks[i] <- mean((object$args$outcome[[response]] - Y_hat)^2)
  }
  
  # 3. Calculate differences in inbag risk between consecutive boosting iterations
  riskdiff <- -1 * diff(inbag_risks)
  # Scale these differences by the number of observations
  riskdiff <- riskdiff / length(object$args$outcome[[response]])
  # Prepend the first risk difference to the riskdiff vector
  if (initial_model==FALSE){first_riskdiff = 0}
  riskdiff <- c(first_riskdiff, riskdiff)
  
  # For offset initial-learner varimp is just measured within boosting iterations
  if(is.null(object$args$initial_learner)){
    riskdiff <- riskdiff[-1]
    learner_indices <- learner_indices[-1]
  }
  
  # 4. Sum up riskdiffs
  explained <- sapply(seq_along(learner_names), FUN = function(i) {
    sum(riskdiff[which(learner_indices == i)])
  })
  names(explained) <- learner_names
  class(explained) <- "varimp"
  # Calculate the proportion of models where the i-th learner was selected
  attr(explained, "selfreqs") <- sapply(seq_along(learner_names), 
                                        function(i) {
                                          mean(learner_indices == i)
                                        })
  
  var_names <- learner_names
  var_names <- sapply(strsplit(var_names, ", "), function(x) {
    do.call(function(...) paste(..., sep = ", "), as.list(x[order(x)]))
  })
  var_order <- order(sapply(var_names, function(i) {
    sum(explained[var_names == i])
  }))
  attr(explained, "variable_names") <- ordered(var_names, levels = unique(var_names[var_order]))
  
  return(explained)
}
