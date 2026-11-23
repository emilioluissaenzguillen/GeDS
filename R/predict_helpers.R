predict_all_orders <- function(object, newdata = NULL, type = c("response", "link")) {

  type <- match.arg(type)
  model <- object$final_model

  if (is.null(newdata)) {
    if (type == "response") {
      return(object$predictions)
    }

    return(list(
      pred_linear    = predict_training_data(object, model, 2L, "link"),
      pred_quadratic = predict_training_data(object, model, 3L, "link"),
      pred_cubic     = predict_training_data(object, model, 4L, "link")
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

  failed_order <- switch(as.character(n), "3" = "quadratic", "4" = "cubic")
  if (!is.null(failed_order) && !is.null(model$stageB_errors[[failed_order]]))
    return(rep(NA_real_, length(object$predictions$pred_linear)))

  if (type == "response") {
    key <- .prediction_key(n)
    pred <- object$predictions[[key]]

    if (is.null(pred)) {
      stop(sprintf("Prediction for '%s' not found.", key))
    }

    return(pred)
  }

  if (n == 2L && inherits(object, "GeDSboost") && !is.null(model$refinement)) {
    pred <- model$F_hat
    if (object$args$normalize_data && .prediction_family(object)$family != "binomial") {
      pred <- pred * object$args$Y_sd + object$args$Y_mean
    }
    return(as.numeric(pred))
  }

  # Linear GeDSgam fit: the local-scoring predictor itself (intercept and formula
  # offset included), matching the response-scale training predictions.
  if (n == 2L && inherits(object, "GeDSgam") && !is.null(model$Y_hat$eta)) {
    return(as.numeric(model$Y_hat$eta))
  }

  key <- .fit_key(n)
  if (!is.null(model[[key]]$linear.predictors))
    return(as.numeric(model[[key]]$linear.predictors))
  if (identical(model[[key]]$penalty$method, "REML"))
    return(as.numeric(.prediction_family(object)$linkfun(model[[key]]$predicted)))
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


.prediction_family <- function(object) {
  family <- object$args$family
  if (inherits(object, "GeDSboost")) {
    family <- get_mboost_family(family@name)
    if (!is.null(object$args$link)) family <- get(family$family)(link = object$args$link)
  }
  family
}

.prediction_clamp_linear <- function(family) {
  identical(family$family, "Gamma") && identical(family$link, "inverse")
}

# Evaluate formula offsets on the prediction rows, independently of the
# response and smooth terms. Boosting deliberately ignores formula offsets.
.gam_prediction_offset <- function(object, newdata) {
  mt <- stats::terms(object$formula)
  index <- attr(mt, "offset")
  result <- numeric(nrow(newdata))
  for (i in index) {
    expression <- attr(mt, "variables")[[i + 1L]][[2L]]
    value <- eval(expression, newdata, environment(object$formula))
    if (!is.numeric(value) || !length(value) %in% c(1L, nrow(newdata)) ||
        any(!is.finite(value))) {
      stop("Formula offsets must provide one finite numeric value per prediction row.",
           call. = FALSE)
    }
    result <- result + value
  }
  result
}

.validate_prediction_link <- function(eta, family) {
  if (.prediction_clamp_linear(family) &&
      any(!is.na(eta) & (!is.finite(eta) | eta <= 0))) {
    stop("Gamma predictions with the inverse link require a positive finite linear predictor; ",
         "check the parametric terms and offsets in 'newdata'.", call. = FALSE)
  }
  invisible(NULL)
}

.prediction_spline_basis <- function(knots, x, n, clamp_linear = FALSE) {
  x <- as.numeric(x)
  lower <- knots[n]
  upper <- knots[length(knots) - n + 1L]
  clipped <- pmax(lower, pmin(upper, x))
  basis <- splines::splineDesign(knots, clipped, ord = n, outer.ok = TRUE)
  outside <- x < lower | x > upper
  if (n == 2L && !clamp_linear && any(outside)) {
    # Use the one-sided segment slopes, including at the right endpoint.
    slope_x <- clipped
    slope_x[x < lower] <- (lower + min(knots[knots > lower])) / 2
    slope_x[x > upper] <- (upper + max(knots[knots < upper])) / 2
    slope <- splines::splineDesign(knots, slope_x, ord = n,
                                   derivs = 1L, outer.ok = TRUE)
    basis <- basis + slope * (x - clipped)
  }
  basis
}

predict_newdata_full_model <- function(object, model, newdata, n, type) {

  prediction_family <- .prediction_family(object)
  clamp_linear <- .prediction_clamp_linear(prediction_family)

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

  if (n == 2L && inherits(object, "GeDSboost") && !is.null(model$refinement) &&
      !isFALSE(model$refinement$prediction_active)) {
    anchor <- model
    anchor$base_learners <- model$refinement$anchor
    anchor$refinement <- NULL
    pred <- predict_newdata_full_model(object, anchor, newdata, n, "link")
    correction <- .boost_refinement_prediction(object, model, newdata)
    if (object$args$normalize_data && prediction_family$family != "binomial") {
      correction <- correction * object$args$Y_sd
    }
    pred <- pred + correction
    .validate_prediction_link(pred, prediction_family)
    return(as.numeric(if (type == "response") object$args$family@response(pred) else pred))
  }

  offset <- if (inherits(object, "GeDSgam")) .gam_prediction_offset(object, newdata) else rep(0, nobs)

  if (n == 2L) {

    # Extract predictor variables from newdata as a data.frame
    X_df <- newdata[, names(object$args$predictors), drop = FALSE]
    # 1.1 Extract linear base learners
    lin_bl <- base_learners[sapply(base_learners, function(x) x$type == "linear")]
    # 1.2 Extract univariate base learners
    univariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) == 1)]
    # 1.3 Extract bivariate base learners
    bivariate_bl <- base_learners[sapply(base_learners, function(x) x$type == "GeDS" & length(x$variables) >= 2)]

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
        pred_univ <- piecewise_multivar_linear_model(X = X_df, model, base_learners = univariate_bl,
                                                     clamp_linear = clamp_linear, extr = object$args$extr)
      }
      # 1.3 GeDS bivariate base learners
      if (length(bivariate_bl) != 0) {
        pred_biv <- bivariate_bl_linear_model(pred_vars = X_df, model,
                                              shrinkage, base_learners = bivariate_bl,
                                              extr = object$args$extr, clamp_linear = clamp_linear)
      }

      if (object$args$normalize_data && family_name != "binomial") {
        pred <- (pred0 + pred_lin + pred_univ + pred_biv)*Y_sd + Y_mean
      } else {
        pred <- pred0 + pred_lin + pred_univ + pred_biv
      }

      .validate_prediction_link(pred, prediction_family)
      pred <- if (type == "response") object$args$family@response(pred) else if (type == "link") pred

      return(as.numeric(pred))

      ##################
      ## 1.2. GeDSgam ##
      ##################
    } else if (inherits(object, "GeDSgam")) {

      # 1.0 Initial learner
      pred0 <- rep(.gam_linear_intercept(model$Y_hat), nrow(newdata)) + offset
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
        pred_univ <- piecewise_multivar_linear_model(X_df, model, base_learners = univariate_bl,
                                                     clamp_linear = clamp_linear, extr = object$args$extr)
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
                                              extr = object$args$extr, type = "gam", clamp_linear = clamp_linear)
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

      .validate_prediction_link(pred, prediction_family)
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
      int.knots <- "quadratic.int.knots"
      Fit <- "quadratic.fit"
    } else if (n == 4) {
      int.knots <- "cubic.int.knots"
      Fit <- "cubic.fit"
    }

    failed_order <- if (n == 3L) "quadratic" else "cubic"
    if (!is.null(model$stageB_errors[[failed_order]]))
      return(rep(NA_real_, nrow(newdata)))

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
      # extrListfit <- lapply(univariate_vars, function(var) range(object$args$predictors[[var]]))
      # extrListnew <- lapply(X_univ, range)
      # extrList <- mapply(function(var, range1, range2) {
      #   if(range2[1] < range1[1] || range2[2] > range1[2])
      #     warning(sprintf("Input values for variable '%s' exceed original boundary knots; extending boundary knots to cover new data range.", var))
      #   c(min(range1[1], range2[1]), max(range1[2], range2[2]))
      # }, var = univariate_vars, range1 = extrListfit, range2 = extrListnew, SIMPLIFY = FALSE)
      #
      extrList <- lapply(univariate_vars, function(var) range(object$args$predictors[[var]]))

      matrices_univ_list <- vector("list", length = ncol(X_univ))
      # Generate design matrices for each predictor
      for (j in 1:ncol(X_univ)) {
        matrices_univ_list[[j]] <- .prediction_spline_basis(
          sort(c(InterKnotsList_univ[[j]], rep(extrList[[j]], n))), X_univ[,j], n
        )
      }
      names(matrices_univ_list) <- names(univariate_learners)
    } else {
      matrices_univ_list <- NULL
    }
    # Joint tensor learners
    if (length(InterKnotsList_biv) != 0) {
      bivariate_learners <- base_learners[
        sapply(base_learners, function(bl) length(bl$variables)) >= 2
      ]
      matrices_biv_list <- list()

      for (learner_name in names(bivariate_learners)) {
        vars <- bivariate_learners[[learner_name]]$variables
        matrices_biv_list[[learner_name]] <- jointSplineBasisMultivar(
          X = X[, vars, drop = FALSE],
          base_learners = bivariate_learners[learner_name],
          InterKnotsList = InterKnotsList_biv[learner_name],
          n = n,
          extrList = lapply(object$args$predictors[vars], range), prediction = TRUE
        )$basis
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

    .validate_prediction_link(pred, family)
    pred <- if (type == "response") family$linkinv(pred) else if (type == "link") pred

    return(as.numeric(pred))

  }
}

predict_newdata_base_learner <- function(object, model, newdata, n,
                                         base_learner, type) {

  clamp_linear <- .prediction_clamp_linear(.prediction_family(object))

  # Single base-learner prediction
  bl_name <-  gsub(",\\s*", ", ", as.character(base_learner))
  bl <- object$args$base_learners[[bl_name]]
  if(is.null(bl)) stop(paste0(bl_name, " not found in the model."))

  # If GeDSboost and the requested base-learner was never selected by the
  # boosting algorithm, its contribution to the fitted predictor is zero.
  if (inherits(object, "GeDSboost") && !bl_name %in% names(model$base_learners)) {
    return(rep(0, nrow(newdata)))
  }

  if (n == 2L && inherits(object, "GeDSboost") && !is.null(model$refinement) &&
      !isFALSE(model$refinement$prediction_active)) {
    anchor <- model
    anchor$base_learners <- model$refinement$anchor
    anchor$linear.fit$theta <- model$refinement$anchor_theta
    anchor$refinement <- NULL
    anchor_object <- object
    anchor_object$final_model <- anchor
    pred <- predict_newdata_base_learner(anchor_object, anchor, newdata, n, bl_name, type)
    correction <- .boost_refinement_prediction(object, model, newdata, bl_name)
    if (bl$type == "GeDS" && length(bl$variables) == 1L &&
        !is.character(model$refinement$anchor_theta)) {
      reference <- .boost_refinement_prediction(object, model, object$args$predictors, bl_name)
      correction <- correction - mean(reference)
    }
    return(as.numeric(pred + correction))
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

  # For the linear GeDSgam fit, use the stored backfitting representation.
  # The combined B-spline coefficient vector is not unique for rank-deficient
  # models, so splitting that vector by base learner can be platform-dependent.
  # Centering each stored contribution separately is equivalent to the group
  # centering used by predict_newdata_full_model().
  if (n == 2L && inherits(object, "GeDSgam")) {
    selected_bl <- object$args$base_learners[bl_name]
    reference_X <- object$args$predictors[, bl$variables, drop = FALSE]

    if (object$args$normalize_data) {
      numeric_reference <- names(reference_X)[vapply(reference_X, is.numeric, logical(1))]
      if (length(numeric_reference) > 0) {
        reference_X[numeric_reference] <- scale(reference_X[numeric_reference])
      }
    }

    if (bl$type == "linear") {
      pred <- lin_model(X_df, model, selected_bl)
      alpha <- mean(lin_model(reference_X, model, selected_bl))
    } else if (length(bl$variables) == 1L) {
      pred <- piecewise_multivar_linear_model(X_df, model, selected_bl,
                                               clamp_linear = clamp_linear, extr = object$args$extr)
      alpha <- mean(piecewise_multivar_linear_model(reference_X, model, selected_bl))
    } else {
      pred <- bivariate_bl_linear_model(
        X_df, model, base_learners = selected_bl,
        extr = object$args$extr, type = "gam", clamp_linear = clamp_linear
      )
      alpha <- mean(bivariate_bl_linear_model(
        reference_X, model, base_learners = selected_bl,
        extr = object$args$extr, type = "gam"
      ))
    }

    return(as.numeric(pred - alpha))
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
      if (clamp_linear) {
        for (variable in bl$variables) {
          boundary <- range(object$args$predictors[[variable]])
          newdata[[variable]] <- pmax(boundary[1], pmin(boundary[2], newdata[[variable]]))
        }
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

      # # If new data exceeds boundary knots limits, redefine boundary knots
      # if(min(X_df) < extr[1] || max(X_df) > extr[2]) {
      #   extr <- range(c(extr, X_df))
      #   warning("Input values exceed original boundary knots; extending boundary knots to cover new data range.")
      # }

      # Create spline basis matrix using specified knots, evaluation points and order
      basisMatrix <- .prediction_spline_basis(sort(c(int.knt, rep(extr, n))), X_df[,1], n,
                                               clamp_linear = clamp_linear)
      # To recover backfitting predictions need de_mean
      pred <- basisMatrix %*% theta
      if (n == 2) {
        reference_x <- object$args$predictors[[bl$variables]]
        if (object$args$normalize_data) {
          reference_x <- (reference_x - object$args$X_mean[[bl$variables]]) /
            object$args$X_sd[[bl$variables]]
        }
        reference <- .prediction_spline_basis(sort(c(int.knt, rep(extr, n))), reference_x, n)
        pred <- pred - mean(reference %*% theta)
      }

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
    # # If new data exceeds boundary knots limits, redefine boundary knots
    # if(min(X_df[,1]) < Xextr[1] || max(X_df[,1]) > Xextr[2] || min(X_df[,2]) < Yextr[1] || max(X_df[,2]) > Yextr[2]) {
    #   Xextr <- range(c(Xextr, X_df[,1]))
    #   Yextr <- range(c(Yextr, X_df[,2]))
    #   warning("Input values exceed original boundary knots; extending boundary knots to cover new data range.")
    # }

    # Generate spline basis matrix for X and Y dimensions using object knots and given order
    basisMatrixX <- .prediction_spline_basis(sort(c(int.knt$ikX, rep(Xextr, n))), X_df[,1], n,
                                             clamp_linear = clamp_linear)
    basisMatrixY <- .prediction_spline_basis(sort(c(int.knt$ikY, rep(Yextr, n))), X_df[,2], n,
                                             clamp_linear = clamp_linear)
    # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
    basisMatrixbiv <- tensorProd(basisMatrixX, basisMatrixY)
    # Multiply the bivariate spline basis by model coefficients to get fitted values
    f_hat_XY_val <- basisMatrixbiv %*% theta[1:dim(basisMatrixbiv)[2]]

    return(as.numeric(f_hat_XY_val))

  } else if (NCOL(X_df) > 2 && bl$type == "GeDS") {
    boundary <- if (n == 2L) object$args$extr[bl$variables] else lapply(pred_vars, range)
    basis <- jointSplineBasisMultivar(
      X_df, object$args$base_learners[bl_name],
      setNames(list(int.knt), bl_name), n, boundary, prediction = TRUE,
      clamp_linear = clamp_linear
    )$basis
    return(as.numeric(basis %*% theta))
  }

}

.unfitted_order_reason <- function(object, n) {
  knot_key <- switch(as.character(n),
                     "3" = "quadratic.int.knots",
                     "4" = "cubic.int.knots",
                     NULL)
  if (is.null(knot_key)) return(NULL)
  knot_list <- object$internal_knots[[knot_key]]
  if (!length(knot_list)) return(NULL)

  joint_sizes <- vapply(names(knot_list), function(name) {
    learner <- object$args$base_learners[[name]]
    knots <- knot_list[[name]]
    if (is.null(learner) || length(learner$variables) < 2L ||
        !is.list(knots)) return(NA_real_)
    prod(as.double(lengths(knots) + n))
  }, numeric(1))
  n_informative <- sum(object$args$weights > 0)
  too_large <- joint_sizes[is.finite(joint_sizes) &
                             joint_sizes >= n_informative]
  if (!length(too_large)) return(NULL)

  paste0(
    " The joint tensor basis requires ",
    format(max(too_large), scientific = FALSE, big.mark = ","),
    " coefficients, but only ", n_informative,
    " observations have positive weight; the unpenalized Stage B fit was omitted."
  )
}
