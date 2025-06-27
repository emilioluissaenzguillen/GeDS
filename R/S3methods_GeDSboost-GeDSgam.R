################################################################################
##################################### COEF #####################################
################################################################################
#' @noRd
coef.GeDSgam_GeDSboost <- function(object, n = 3L, ...)
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
    coefficients <- model$quadratic.fit$theta
  }
  if(n == 4L){
    coefficients <- model$cubic.fit$theta
  }
  
  return(coefficients)
}

#' @title Coef Method for GeDSgam, GeDSboost
#' @name coef.GeDSgam,boost
#' @description
#' Method for the function \code{\link[stats]{coef}} that allows to extract the
#' estimated coefficients of a \code{"GeDSgam"} class or \code{"GeDSboost"} class
#' object.
#' @param object The \code{"GeDSgam"} class or
#' \code{"GeDSboost"} class object from which the coefficients should be
#' extracted.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{ + 1}) of the FGB-GeDS/GAM-GeDS fit whose coefficients should be
#' extracted. 
#' \itemize{
#'  \item If \code{n = 2L}: the piecewise polynomial coefficients of the univariate
#'  GeDS base-learners and linear base-learners are provided (the B-spline coefficients
#'  are accessible through \code{x$final_model$linear.fit$theta}); for bivariate GeDS
#'  base-learners and \code{class(object) == "GeDSboost"}, the B-spline coefficients
#'  for each iteration where the corresponding bivariate base-learner was selected are
#'  provided; for bivariate base-learners and \code{class(object) == "GeDSgam"},
#'  the final local-scoring B-spline coefficients are provided.
#'   \item If \code{n = 3L} or \code{n = 4L}: B-spline coefficients are provided
#'  for both univariate and bivariate GeDS learners.
#' }
#' By default \code{n} is equal to \code{3L}. Non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#' 
#' @return
#' A named vector containing the required coefficients of the fitted
#' multivariate predictor model.
#' 
#' @seealso \code{\link[stats]{coef}} for the standard definition;
#' \code{\link{NGeDSgam}} and \code{\link{NGeDSboost}} for examples.
#' 
#' @examples
#' data(mtcars)
#' # Convert specified variables to factors
#' categorical_vars <- c("cyl", "vs", "am", "gear", "carb")
#' mtcars[categorical_vars] <- lapply(mtcars[categorical_vars], factor)
#' 
#' Gmodgam <- NGeDSgam(mpg ~ f(disp, hp) + f(qsec) + carb,
#'                     data = mtcars, family = gaussian, phi = 0.7)
#' 
#' # Piecewise polynomial coefficients of the (univariate) GeDS and linear learners
#' coef(Gmodgam, n = 2)$`f(qsec)`
#' coef(Gmodgam, n = 2)$carb
#' # B-spline coefficients of the bivariate learner
#' coef(Gmodgam, n = 2)$`f(disp, hp)`
#' 
#' Gmodboost <- NGeDSboost(mpg ~ f(disp, hp) + f(qsec) + carb,
#'                         data = mtcars, family = mboost::Gaussian(), shrinkage = 0.7)
#' 
#' # Piecewise polynomial coefficients of the (univariate) GeDS and linear learners
#' coef(Gmodboost, n = 2)$`f(qsec)`
#' coef(Gmodboost, n = 2)$carb
#' # B-spline coefficients of the bivariate learner at each boosting iteration
#' # this was selected
#' coef(Gmodboost, n = 2)$`f(disp, hp)`
#' 
#' @aliases coef.GeDSgam, coef.GeDSboost

#' @export
#' @rdname coef.GeDSgam_GeDSboost
coef.GeDSgam <- coef.GeDSgam_GeDSboost

#' @export
#' @rdname coef.GeDSgam_GeDSboost
coef.GeDSboost <- coef.GeDSgam_GeDSboost

################################################################################
################################### CONFINT ####################################
################################################################################
#' @export
#' @rdname confint.GeDS
confint.GeDSgam <- confint.GeDS

#' @export
#' @rdname confint.GeDS
confint.GeDSboost <- confint.GeDS

################################################################################
################################### devIANCE ###################################
################################################################################
#' @noRd
deviance.GeDSgam_GeDSboost <- function(object, n = 3L, ...)
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
    dev <- model$dev
  }
  if(n == 3L){
    dev <- model$quadratic.fit$rss
  }
  if(n == 4L){
    dev <- model$cubic.fit$rss
  }
  dev <- as.numeric(dev)
  return(dev)
}

#' @export
#' @rdname deviance.GeDS
deviance.GeDSgam <- deviance.GeDSgam_GeDSboost

#' @export
#' @rdname deviance.GeDS
deviance.GeDSboost <- deviance.GeDSgam_GeDSboost

################################################################################
#################################### FAMILY ####################################
################################################################################
#' @noRd
family.GeDSgam_GeDSboost <- function(object, ...) object$args$family

#' @export
#' @rdname family.GeDS
family.GeDSgam <- family.GeDSgam_GeDSboost

#' @export
#' @rdname family.GeDS
family.GeDSboost <- family.GeDSgam_GeDSboost

################################################################################
################################### FORMULA ####################################
################################################################################
#' @export
#' @rdname formula.GeDS
formula.GeDSgam <- formula.GeDS

#' @export
#' @rdname formula.GeDS
formula.GeDSboost <- formula.GeDS

################################################################################
##################################### KNOTS ####################################
################################################################################
#' @noRd
knots.GeDSgam_GeDSboost <-  function(Fn, n = 3L,
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
knots.GeDSgam <- knots.GeDSgam_GeDSboost

#' @export
#' @rdname knots
knots.GeDSboost <- knots.GeDSgam_GeDSboost

################################################################################
#################################### LOGLIK ####################################
################################################################################
#' @export
#' @rdname logLik.GeDS
logLik.GeDSgam <- logLik.GeDS

#' @export
#' @rdname logLik.GeDS
logLik.GeDSboost <- logLik.GeDS

################################################################################
################################# N.boost.iter #################################
################################################################################
#' @title Extract Number of Boosting Iterations from a GeDSboost Object
#' @name N.boost.iter
#' @description
#' Method for \code{N.boost.iter} that returns the number of boosting iterations 
#' used in fitting a \code{"GeDSboost"} class object.
#'
#' @param object A \code{"GeDSboost"} class object.
#' @param ... Further arguments (ignored).
#'
#' @return An integer indicating the total number of boosting iterations.
#' 
#' @method N.boost.iter GeDSboost
#' 
#' @rdname N.boost.iter
#' @aliases N.boost.iter N.boost.iter.GeDSboost
#' @importFrom graphics mtext abline
#' @export

N.boost.iter.GeDSboost <- function(object, ...) {
  object$iters
}
#' @export
N.boost.iter <- N.boost.iter.GeDSboost

################################################################################
#################################### PREDICT ###################################
################################################################################
# Helper function to extract predictions
extract_link_pred <- function(fit) {
  if (inherits(fit, "glm")) {
    return(as.numeric(fit$linear.predictors))
  } else if (inherits(fit, "lm")) {
    return(as.numeric(fit$fitted.values))
  } else {
    stop("Unsupported model type in extract_link_pred().")
  }
}

#' @importFrom splines splineDesign
#' @noRd
predict.GeDSgam_GeDSboost <- function(object, newdata, n = 3L,
                                      base_learner = NULL, type = c("response", "link"), ...)
{
  # Handle additional arguments
  if(!missing(...)) warning("Only 'object','newdata', and 'n' arguments will be considered")
  
  # Check if object is of class "GeDSboost" or "GeDSgam"
  if (!inherits(object, "GeDSboost") && !inherits(object, "GeDSgam")) {
    stop("The input 'object' must be of class 'GeDSboost' or 'GeDSgam'")
  }
  # Validate and match the 'type' argument 
  type <- match.arg(type, c("response", "link"))
  # Select final model
  model <- object$final_model
  
  # ─────────── n = "all" ─────────── #
  if (identical(n, "all")) {
    if (missing(newdata)) {
      if (type == "response") {
        return(object$predictions)
      } else if (type == "link") {
        return(list(
          pred_linear    = extract_link_pred(model$linear.fit$temporary),
          pred_quadratic = extract_link_pred(model$quadratic.fit$temporary),
          pred_cubic     = extract_link_pred(model$cubic.fit$temporary)
        ))
      }
      
    } else {
      orders <- 2L:4L
      names <- c("pred_linear", "pred_quadratic", "pred_cubic")
      preds <- lapply(orders, function(k) predict(object, newdata = newdata, n = k, type = type))
      return(setNames(preds, names))
    }
  }
  
  # ─────────── Specific n ─────────── #
  n <- as.integer(n)
  if (!(n %in% 2L:4L)) {
    warning(sprintf("'n' incorrectly specified (%s). Set to 3.", n))
    n <- 3L
  }
  
  if (missing(newdata)) {
    if (type == "response") {
      key <- switch(as.character(n),
                    "2" = "pred_linear",
                    "3" = "pred_quadratic",
                    "4" = "pred_cubic")
      pred <- object$predictions[[key]]
      if (is.null(pred)) stop(sprintf("Prediction for '%s' not found.", key))
      return(pred)
    } else if (type == "link") {
      key <- switch(as.character(n),
                    "2" = "linear.fit",
                    "3" = "quadratic.fit",
                    "4" = "cubic.fit")
      Fit <- object$final_model[[key]]$temporary
      
      if (inherits(Fit, "glm")) {
        return(as.numeric(Fit$linear.predictors))
        } else if (inherits(Fit, "lm")) {
          return(as.numeric(Fit$fitted.values))
        }
    }
  }
  
  # Ensure newdata is a data.frame
  newdata <- as.data.frame(newdata)
  
  ###############################################
  # I. Compute predictions for all base_leaners #
  ###############################################
  if (is.null(base_learner)) {
    
    # Check whether newdata includes all the necessary predictors
    if (!all(names(object$args$predictors) %in% colnames(newdata))) {
      missing_vars <- setdiff(names(object$args$predictors), colnames(newdata))
      stop(paste("The following predictors are missing in newdata:",
                 paste(missing_vars, collapse = ", ")))
    }
    
    nobs <- nrow(newdata)
    base_learners <- object$args$base_learners
    offset <- rep(0, nobs)
    
    ###############
    ## 1. LINEAR ##
    ###############
    if (n == 2) {
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
      
      #####################
      ## 2. Higher Order ##
      #####################
    } else if (n==3 || n==4) {
      
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
      
      X = newdata[GeDS_variables]
      Z = newdata[linear_variables]
      
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
      base_learners =  base_learners[sapply(base_learners, function(x) x$type == "GeDS")]
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
        Z <-  Z[, colnames(Z) != "(Intercept)"]
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
    
    ####################################################
    # II. Compute predictions for a single base_leaner #
    # (at the linear predictor level)                  #
    ####################################################
    # This may be simplified but it is useful for checking that the linear B-spline
    # representation of the models is correct
    } else if (!is.null(base_learner)) {
      
      # Single base-learner prediction
      bl_name <-  gsub(",\\s*", ", ", as.character(base_learner))
      bl <- object$args$base_learners[[bl_name]]
      if(is.null(bl)) stop(paste0(bl_name, " not found in the model."))
      
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
}

#' @title Predict Method for GeDSgam, GeDSboost
#' @name predict.GeDSgam,boost
#' @description 
#' This method computes predictions from GeDSboost and GeDSgam objects. 
#' It is designed to be user-friendly and accommodate different orders of the
#' GeDSboost or GeDSgam fits.
#' @param object The \code{"GeDSgam"} class or
#' \code{"GeDSboost"} class object.
#' @param newdata An optional \code{data.frame} containing values of the independent
#' variables at which predictions are to be computed. If omitted, the fitted 
#' values are extracted from the \code{object} itself.
#' @param type Character string indicating the type of prediction required. By
#' default it is equal to \code{"response"}, i.e., the result is on the scale of
#' the response variable. Alternatively if one wants the predictions to be on the
#' predictor scale, it is necessary to set \code{type = "link"}.
#' @param n The order of the GeDS fit (\code{2L} for linear, \code{3L} for
#' quadratic, and \code{4L} for cubic). Alternatively, \code{"all"} can be used
#' to return a list with the predictions from all three fits. Default is \code{3L}.
#' @param base_learner Either \code{NULL} or a \code{character} string specifying
#' the base-learner of the model for which predictions should be computed. Note
#' that single base-learner predictions are provided on the linear predictor scale.
#' @param ... Potentially further arguments.
#' 
#' @return Numeric vector (or list) of predictions.
#' 
#' @references
#' Gu, C. and Wahba, G. (1991).
#' Minimizing GCV/GML Scores with Multiple Smoothing Parameters via the Newton Method.
#' \emph{SIAM J. Sci. Comput.}, \strong{12}, 383--398.
#' 
#' @examples
#' ## Gu and Wahba 4 univariate term example ##
#' # Generate a data sample for the response variable
#' # y and the covariates x0, x1 and x2; include a noise predictor x3
#' set.seed(123)
#' N <- 400
#' f_x0x1x2 <- function(x0,x1,x2) {
#'   f0 <- function(x0) 2 * sin(pi * x0)
#'   f1 <- function(x1) exp(2 * x1)
#'   f2 <- function(x2) 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
#'   f <- f0(x0) + f1(x1) + f2(x2)
#'   return(f)
#'}
#' x0 <- runif(N, 0, 1)
#' x1 <- runif(N, 0, 1)
#' x2 <- runif(N, 0, 1)
#' x3 <- runif(N, 0, 1)
#' # Specify a model for the mean of y
#' f <- f_x0x1x2(x0 = x0, x1 = x1, x2 = x2)
#' # Add (Normal) noise to the mean of y
#' y <- rnorm(N, mean = f, sd = 0.2)
#' data <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2, x3 = x3)
#' 
#' # Fit GAM-GeDS model
#' Gmodgam <- NGeDSgam(y ~ f(x0) + f(x1) + f(x2) + f(x3), data = data)
# 
#' # Check that the sum of the individual base-learner predictions equals the final
#' # model prediction
#' 
#' pred0 <- predict(Gmodgam, n = 2, newdata = data, base_learner = "f(x0)")
#' pred1 <- predict(Gmodgam, n = 2, newdata = data, base_learner = "f(x2)")
#' pred2 <- predict(Gmodgam, n = 2, newdata = data, base_learner = "f(x1)")
#' pred3 <- predict(Gmodgam, n = 2, newdata = data, base_learner = "f(x3)")
# 
#' round(predict(Gmodgam, n = 2, newdata = data) -
#' (mean(predict(Gmodgam, n = 2, newdata = data)) + pred0 + pred1 + pred2 + pred3), 12)
#' 
#' pred0 <- predict(Gmodgam, n = 3, newdata = data, base_learner = "f(x0)")
#' pred1 <- predict(Gmodgam, n = 3, newdata = data, base_learner = "f(x2)")
#' pred2 <- predict(Gmodgam, n = 3, newdata = data, base_learner = "f(x1)")
#' pred3 <- predict(Gmodgam, n = 3, newdata = data, base_learner = "f(x3)")
#' 
#' round(predict(Gmodgam, n = 3, newdata = data) - (pred0 + pred1 + pred2 + pred3), 12)
#' 
#' pred0 <- predict(Gmodgam, n = 4, newdata = data, base_learner = "f(x0)")
#' pred1 <- predict(Gmodgam, n = 4, newdata = data, base_learner = "f(x2)")
#' pred2 <- predict(Gmodgam, n = 4, newdata = data, base_learner = "f(x1)")
#' pred3 <- predict(Gmodgam, n = 4, newdata = data, base_learner = "f(x3)")
#' 
#' round(predict(Gmodgam, n = 4, newdata = data) - (pred0 + pred1 + pred2 + pred3), 12)
#' 
#' # Plot GeDSgam partial fits to f(x0), f(x1), f(x2)
#' par(mfrow = c(1,3))
#' for (i in 1:3) {
#'   # Plot the base learner
#'   plot(Gmodgam, n = 3, base_learners = paste0("f(x", i-1, ")"), col = "seagreen",
#'        cex.lab = 1.5, cex.axis = 1.5)
#'   # Add legend
#'   if (i == 2) {
#'     position <- "topleft"
#'     } else if (i == 3) {
#'       position <- "topright"
#'       } else {
#'         position <- "bottom"
#'       }
#'   legend(position, legend = c("GAM-GeDS Quadratic", paste0("f(x", i-1, ")")),
#'          col = c("seagreen", "darkgray"),
#'          lwd = c(2, 2),
#'          bty = "n",
#'          cex = 1.5)
#' }
#' 
#' @aliases predict.GeDSboost, predict.GeDSgam
#' @export
#' @rdname predict.GeDSgam_GeDSboost
predict.GeDSgam <- predict.GeDSgam_GeDSboost

#' @export
#' @rdname predict.GeDSgam_GeDSboost
predict.GeDSboost <- predict.GeDSgam_GeDSboost

################################################################################
##################################### PRINT ####################################
################################################################################
#' @noRd 
print.GeDSgam_GeDSboost <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    ...)
  {
  cat("\ncall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  kn <- knots(x,n=2,options="int")
  devs <- numeric(3)
  names(devs) <- c("Order 2","Order 3","Order 4")
  
  # Iterate over each base learner
  for (bl_name in names(x$internal_knots$linear.int.knots)) {
    # Get the linear internal knots
    int.knt <- x$internal_knots$linear.int.knots[[bl_name]]
    
    # 1) No knots or linear base-learner
    if (is.null(int.knt)) {
      cat(paste0("Number of internal knots of the second order (linear) spline for '",
                 bl_name, "': 0\n"))
    } else {
      # 2) Bivariate GeDS
      if (is.list(int.knt)) {
        for (component_name in names(int.knt)) {
          component <- int.knt[[component_name]]
          cat(paste0("Number of internal knots of the second order (linear) spline for '",
                     bl_name, "', ", component_name, ": ", length(component), "\n"))
        }
        # 3) Univariate GeDS
      } else { 
        cat(paste0("Number of internal knots of the second order (linear) spline for '",
                   bl_name, "': ", length(int.knt), "\n"))
      }
    }
  }
  cat("\n")
  devs[1] <- deviance(x, n = 2L)
  devs[2] <- if(!is.null(deviance(x, n = 3L))) deviance(x, n = 3L) else NA
  devs[3] <- if(!is.null(deviance(x, n = 4L))) deviance(x, n = 4L) else NA
  cat("deviances:\n")
  print.default(format(devs, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  print <- list("n.knots" = length(kn), "deviances" = devs, "call" = x$extcall)
  
  x$print <- print
  invisible(x)
}

#' @rdname print.GeDS
#' @export
print.GeDSgam <- print.GeDSgam_GeDSboost

#' @rdname print.GeDS
#' @export
print.GeDSboost <- print.GeDSgam_GeDSboost

################################################################################
################################### SUMMARY ####################################
################################################################################
#' @noRd 
summary.GeDSgam_GeDSboost <- function(object, ...)
{
  
  # 1) Header
  if (inherits(object, "GeDSboost")) {
    cat("\nFGB-GeDS Model\n\n")
  } else if (inherits(object, "GeDSgam")) {
    cat("\nGAM-GeDS Model\n\n")
  }
  
  # 2) Core print
  print.GeDSgam_GeDSboost(object, ...)
  
  # 3) Normalize
  if (object$args$normalize_data) cat("Data was normalized (standardized) before fitting.\n")

  # 4) Boosting/GAM details
  
  if (inherits(object, "GeDSboost")) {
    cat("\nBoosting details:\n")
    print(family.GeDSboost(object,...),...)
    cat("Number of boosting iterations:", object$iters, "\n")
    cat("Final model corresponds to boosting iteration", sub("^model", "", object$final_model$model_name), "\n")
    if (object$args$initial_learner) {
      if (object$args$family@name == "Squared Error (Regression)") {
        cat("A NGeDS() initial learner was used; ")
      } else {
        cat("A GGeDS() initial learner was used; ")
      }
      
    } else {
      cat("An offset initial learner was used; ")
    }
    cat("number of base learners:", length(names(object$final_model$base_learners)), "\n\n")
  } else if (inherits(object, "GeDSgam")) {
    cat("\nLocal-scoring details:\n")
    print(family.GeDSgam(object,...),...)
    cat("Number of local-scoring iterations:", object$iters$local_scoring, "\n")
  }
  
  invisible(object)
}

#' @rdname summary.GeDS
#' @export
summary.GeDSgam <- summary.GeDSgam_GeDSboost

#' @rdname summary.GeDS
#' @export
summary.GeDSboost <- summary.GeDSgam_GeDSboost

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
  if (length(lines) == 1) {
    lines[1] <- paste0("{", lines[1])
    } else if (length(lines) == 2) {
      lines[1] <- paste0("{", lines[1], ",")
      } else if (length(lines) == 3) {
        lines[1] <- paste0("{", lines[1], ",")
        lines[2] <- paste0(lines[2], ",")
        } else if (length(lines) == 4) {
          lines[1] <- paste0("{", lines[1], ",")
          lines[2] <- paste0(lines[2], ",")
          lines[3] <- paste0(lines[3], ",")
          }
  lines[length(lines)] <- paste0(lines[length(lines)], "}")
  
  lines
}

#' @title Visualize Boosting Iterations
#' @name visualize_boosting
#' @description
#' This function plots the \code{\link{NGeDSboost}} fit to the data at the
#' beginning of a given boosting iteration and then plots the subsequent
#' \code{\link{NGeDS}} fit on the corresponding negative gradient.
#' Note that this is only applicable to \code{\link{NGeDSboost}} models with a
#' single univariate base-learner.
#' @param iters Numeric, specifies the iteration(s) number.
#' @param object A \code{"GeDSboost"} class object.
#' @param final_fits Logical indicating whether the final linear, quadratic and
#' cubic fits should be plotted.
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
#' Gmodboost <- NGeDSboost(Y ~ f(X), data = data, normalize_data = TRUE)
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
#' # Visualize boosting iterations + final fits
#' par(mfrow=c(4,2))
#' visualize_boosting(Gmodboost, iters = 0:3, final_fits = TRUE)
#' par(mfrow=c(1,1))
#' 
#' @rdname visualize_boosting
#' @aliases visualize_boosting visualize_boosting.GeDSboost
#'  
#' @importFrom graphics plot lines legend rug points abline mtext
#' @export

visualize_boosting.GeDSboost <- function(object, iters = NULL, final_fits = FALSE)
{
  # Check if object is of class "GeDSboost"
  if(!inherits(object, "GeDSboost")) {
    stop("The input 'object' must be of class 'GeDSboost'")
  }
  # Check if object has only one predictor
  if(length(object$args$predictors) > 1) {
    stop("Visualization only available for models with a single predictor")
  }
  # If M is NULL, set it to be the initial model
  if (is.null(iters)) iters <- 0
  # Check if M > than # of boosting iterations
  if(any(iters > object$iters)) {
    stop(paste0("iters = ", iters, " but only ", object$iters, " boosting iterations were run."))
  }
  
  Y <- object$args$response[[1]]; X <- object$args$predictors[[1]]
  
  for (M in iters) {
    
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
    knots <- model$base_learners[[1]]$knots
    int.knots_round <- round(as.numeric(get_internal_knots(knots)), 2) 
    int.knots_lines <- split_into_lines(paste(int.knots_round, collapse=", "), 65)
    # Create the title
    title_text_1 <- bquote(atop(bold(Delta)[.(M) * ",2"] == .(int.knots_lines[1])))
    
    # Plot
    plot(X, Y, pch = 20, col = "darkgrey", tck = 0.02, main = "",
         xlim = x_range, ylim = y_range, cex.axis = 1.5, cex.lab = 1.5)
    lines(X, Y_hat, lwd = 2)
    legend("topleft",
           legend = c(2*M),
           bty = "n",
           text.font = 2,
           cex = 1.5)
    
    # Trace knots w/vertical lines
    if (length(knots) < 20) {
      for (knot in knots) {
        abline(v = knot, col = "gray", lty = 2)
      }
    } else {
      rug(knots)
    }
    
    # Add the title
    if(length(int.knots_lines) == 1) {
      mtext(title_text_1, side = 3, line = -0.5, cex = 1.15)
    } else if (length(int.knots_lines) == 2) {
      mtext(title_text_1, side = 3, line = 0, cex = 1.15)
      mtext(int.knots_lines[2], side = 3, line = 0.5, cex = 1.15)
    } else if (length(int.knots_lines) == 3) {
      mtext(title_text_1, side = 3, line=0.75, cex = 0.9)
      mtext(int.knots_lines[2], side = 3, line = 1.4, cex = 0.9)
      mtext(int.knots_lines[3], side = 3, line = 0.5, cex = 0.9)
    } else if (length(int.knots_lines) == 4) {
      mtext(title_text_1, side = 3, line = 1, cex = 0.9)
      mtext(int.knots_lines[2], side = 3, line = 1.8, cex = 0.9)
      mtext(int.knots_lines[3], side = 3, line = 0.9, cex = 0.9)
      mtext(int.knots_lines[4], side = 3, line = 0, cex =0.9)
    }
    
    # 2. Negative gradient plot
    if (M < length(object$models) - 1) {
      
      family <- object$args$family
      family_stats <- get_mboost_family(family@name)
      
      negative_gradient <- family@ngradient(y = Y, f = family_stats$linkfun(Y_hat), w = object$args$weights)
      
      # Plot range
      x_range <- range(X)
      y_range <- range(negative_gradient, next_model$best_bl$pred_linear)
      x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
      y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
      
      # Split int knots into lines
      int.knots <- next_model$best_bl$int.knt
      int.knots_round <- round(as.numeric(next_model$best_bl$int.knt), 2)
      int.knots_lines <- split_into_lines(paste(int.knots_round, collapse=", "), 65)    
      # Create the title
      title_text_2 <- bquote(atop(bold(delta)[.(M) * ",2"] == .(int.knots_lines[1])))
      
      # Plot
      plot(X, negative_gradient, pch=20, col=c("darkgrey"), tck = 0.02, main = "",
           xlim = x_range, ylim = y_range, cex.axis = 1.5, cex.lab = 1.5)
      lines(X, next_model$best_bl$pred_linear, col="blue", lty = 2, lwd = 1)
      points(next_model$best_bl$coef$mat[,1], next_model$best_bl$coef$mat[,2], col="blue", pch=21)
      legend("topleft",
             legend = c(2*M+1),
             bty = "n",
             text.font = 2,
             cex = 1.5)
      
      # Trace knots w/vertical lines
      if (length(int.knots) < 20) {
        for(int.knot in int.knots) {
          abline(v = int.knot, col = "gray", lty = 2)
        }
      } else {
        rug(int.knots)
      }
      
      # Add the title
      if (length(int.knots_lines) == 1) {
        mtext(title_text_2, side = 3, line = -0.5, cex = 1.15)
      } else if (length(int.knots_lines) == 2) {
        mtext(title_text_2, side = 3, line = 0, cex = 1.15)
        mtext(int.knots_lines[2], side = 3, line = 0.5, cex = 1.15)
      } else if (length(int.knots_lines) == 3) {
        mtext(title_text_2, side = 3, line = 0.75, cex = 0.9)
        mtext(int.knots_lines[2], side = 3, line = 1.4, cex = 0.9)
        mtext(int.knots_lines[3], side = 3, line = 0.5, cex = 0.9)
      } else if (length(int.knots_lines) == 4) {
        mtext(title_text_2, side = 3, line = 1, cex = 0.9)
        mtext(int.knots_lines[2], side = 3, line = 1.8, cex = 0.9)
        mtext(int.knots_lines[3], side = 3, line = 0.9, cex = 0.9)
        mtext(int.knots_lines[4], side = 3, line = 0, cex = 0.9)
      }
    }
  }
  
  # Add a final plot with the final fits
  if (final_fits) {
    
    x_range <- range(X)
    y_range <- range(Y, object$predictions$pred_linear, object$predictions$pred_quadratic, object$predictions$pred_cubic)
    x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
    y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
    plot(X, Y, pch = 20, col = "darkgrey", tck = 0.02, main = "Final fits",
         xlim = x_range, ylim = y_range, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
    lines(X, object$predictions$pred_linear, col = "green4", lwd = 2, lty = 1)
    lines(X, object$predictions$pred_quadratic, col="red", lwd = 2, lty = 4)
    lines(X, object$predictions$pred_cubic, col="purple", lwd = 2, lty = 5)
    legend("topright",                           
           legend = c("Order 2 (degree=1)", "Order 3 (degree=2)", "Order 4 (degree=3)"),    
           col = c("green4", "red", "purple"),
           lty = c(1, 4, 5),
           lwd = c(2, 2, 2),                                   
           cex = 1.5,
           bty="n")
    
  }
}


#' @export
visualize_boosting <- visualize_boosting.GeDSboost


################################################################################
############################ BASE LEARNER IMPORTANCE ###########################
################################################################################
#' @title Base Learner Importance for GeDSboost Objects
#' @name bl_imp.GeDSboost
#' @description
#' S3 method for \code{"GeDSboost"} class objects that calculates the
#' in-bag risk reduction ascribable to each base-learner of an FGB-GeDS model.
#' Essentially, it measures and aggregates the decrease in the empirical risk
#' attributable to each base-learner for every time it is selected across the
#' boosting iterations. This provides a measure on how much each base-learner
#' contributes to the overall improvement in the model's accuracy, as reflected
#' by the decrease in the empiral risk (loss function). This function is adapted
#' from \code{\link[mboost]{varimp}} and is compatible with the available
#' \code{\link[mboost]{mboost-package}} methods for \code{\link[mboost]{varimp}},
#' including \code{plot}, \code{print} and \code{as.data.frame}.
#' 
#' @param object An object of class \code{"GeDSboost"} class.
#' @param boosting_iter_only Logical value, if \code{TRUE} then base-learner
#' in-bag risk reduction is only computed across boosting iterations, i.e.,
#' without taking into account a potential initial GeDS learner.
#' @param ... Potentially further arguments.
#'
#' @details
#' See \code{\link[mboost]{varimp}} for details.
#' 
#' @return An object of class \code{varimp} with available \code{plot},
#' \code{print} and \code{as.data.frame} methods.
#' 
#' @references
#' Hothorn T., Buehlmann P., Kneib T., Schmid M. and Hofner B. (2022).
#' mboost: Model-Based Boosting. R package version 2.9-7, \url{https://CRAN.R-project.org/package=mboost}.
#' 
#' @examples
#' library(GeDS)
#' library(TH.data)
#' 
#' data("bodyfat", package = "TH.data")
#' N <- nrow(bodyfat); ratio <- 0.8
#' set.seed(123)
#' trainIndex <- sample(1:N, size = floor(ratio * N))
#' # Subset the data into training and test sets
#' train <- bodyfat[trainIndex, ]
#' test <- bodyfat[-trainIndex, ]
#' 
#' Gmodboost <- NGeDSboost(formula = DEXfat ~ f(hipcirc) + f(kneebreadth) + f(anthro3a),
#'                         data = train, phi = 0.7, initial_learner = FALSE)
#' MSE_Gmodboost_linear <- mean((test$DEXfat - predict(Gmodboost, newdata = test, n = 2))^2)
#' MSE_Gmodboost_quadratic <- mean((test$DEXfat - predict(Gmodboost, newdata = test, n = 3))^2)
#' MSE_Gmodboost_cubic <- mean((test$DEXfat - predict(Gmodboost, newdata = test, n = 4))^2)
#' 
#' # Print MSE
#' cat("\n", "TEST MEAN SQUARED ERROR", "\n",
#'     "Linear NGeDSboost:", MSE_Gmodboost_linear, "\n",
#'     "Quadratic NGeDSboost:", MSE_Gmodboost_quadratic, "\n",
#'     "Cubic NGeDSboost:", MSE_Gmodboost_cubic, "\n")
#' 
#' # Base Learner Importance
#' bl_imp <- bl_imp(Gmodboost)
#' print(bl_imp)
#' plot(bl_imp)
#' 
#' @rdname bl_imp
#' @aliases bl_imp bl_imp.GeDSboost
#' @method bl_imp GeDSboost
#' @export

bl_imp.GeDSboost <- function(object, boosting_iter_only = FALSE, ...)
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
  bl_selected <- lapply(object$models, function(model) model$best_bl$name)
  bl_indices <- match(bl_selected, bl_names)
  # Risk function
  risk <- object$args$family@risk
  
  # 2. In-bag risk
  # Calculate initial risk
  initial_risk <- risk(y = object$args$response[[response]],
                       f = 0, w = object$args$weights)
  # Get the riskdiff of model0
  model0_risk <- risk(y = object$args$response[[response]],
                      f = object$models[[paste0("model", 0)]]$F_hat, w = object$args$weights)
  first_riskdiff <- initial_risk - model0_risk
  
  # Get the number of models
  n_models <- length(object$models)
  if (!object$args$initial_learner) n_models <- n_models - 1
  # Initialize an empty vector to store the inbag risks
  inbag_risks <- vector(mode = "numeric", length = n_models)
  # Loop over each model from model1 onwards
  for (i in seq_len(n_models - 1)) {
    # Extract the selected predictor for each model
    F_hat <- object$models[[paste0("model", i)]]$F_hat
    # Calculate the inbag risk for this iteration
    inbag_risks[i] <- risk(y = object$args$response[[response]],
                           f = F_hat, w = object$args$weights)
  }
  inbag_risks <- c(model0_risk, inbag_risks)
  
  # 3. Calculate differences in inbag risk between consecutive boosting iterations
  riskdiff <- -1 * diff(inbag_risks)
  # Scale these differences by the number of observations
  riskdiff <- riskdiff / length(object$args$response[[response]])
  # Prepend the first risk difference to the riskdiff vector
  riskdiff <- c(first_riskdiff, riskdiff)
  
  # For offset initial-learner bl_imp is just measured within boosting iterations
  if(!object$args$initial_learner || boosting_iter_only){
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

