#' @importFrom splines splineDesign
#' @importFrom stats lm coef .lm.fit lm.fit lm.wfit residuals
SplineReg_biv <- function(X ,Y , Z = NULL, W = NULL, offset = rep(0,length(X)), weights = rep(1,length(X)),
                          InterKnotsX, InterKnotsY, n, Xextr = range(X), Yextr = range(Y), flag = TRUE,
                          center = c(sum(Xextr)/2,sum(Yextr)/2), coefficients = NULL, fast = FALSE)
  {
  # Convert spline order to integer
  n <- as.integer(n)
  # Create spline basis matrix using specified knots, evaluation points and order
  basisMatrixX <- splineDesign(knots=sort(c(InterKnotsX,rep(Xextr,n))),derivs=rep(0,length(X)),
                           x=X,ord=n,outer.ok = T)
  basisMatrixY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))),derivs=rep(0,length(Y)),
                           x=Y,ord=n,outer.ok = T)

  basisMatrixbiv <- tensorProd(basisMatrixX,basisMatrixY)
  # basisMatrixbivbis <- tensorProd_R(basisMatrixX,basisMatrixY)
  # any(basisMatrixbiv-basisMatrixbivbis!=0)
  # basisMatrixXbis <- recoverXmat(basisMatrixbiv, basisMatrixY, basisMatrixX)
  # any(basisMatrixX-basisMatrixXbis!=0)

  # Combine spline basis with parametric design matrix (if provided)
  basisMatrixbiv2 <- cbind(basisMatrixbiv, W)

  if (fast && is.null(coefficients)) {
    if (is.null(Z)) {
      stop("'Z' is required when fitting a bivariate spline without supplied coefficients.",
           call. = FALSE)
    }
    Z0 <- Z - offset
    if (all(weights == 1)) {
      tmp <- .lm.fit(basisMatrixbiv2, Z0)
      if (tmp$rank < ncol(basisMatrixbiv2)) {
        tmp <- lm.fit(basisMatrixbiv2, Z0)
      }
      resid <- residuals(tmp)
    } else {
      tmp <- lm.wfit.light(basisMatrixbiv2, Z0, weights)
      resid <- tmp$residuals

      if (tmp$rank < ncol(basisMatrixbiv2)) {
        tmp <- lm.wfit(basisMatrixbiv2, Z0, weights)
        resid <- residuals(tmp)
      }
    }

    theta <- as.numeric(coef(tmp))
    if (anyNA(theta)) theta[is.na(theta)] <- 0
    predicted <- as.numeric(basisMatrixbiv2 %*% theta) + offset
    resid <- Z - predicted

    return(list("theta" = theta, "predicted" = predicted,
                "residuals" = resid, "rss" = .weighted_rss(resid, weights),
                "Xbasis" = basisMatrixX, "Ybasis" = basisMatrixY,
                "Xknots" = sort(c(InterKnotsX, rep(Xextr, n))),
                "Yknots" = sort(c(InterKnotsY, rep(Yextr, n))),
                "temporary" = tmp))
  }

  # If coefficients vector was provided, check whether this is conformable with the knots vectors, o.w. re-estimate the coefficients
  if (!is.null(coefficients)) {
    non_conformable <- NCOL(basisMatrixbiv2) != length(coefficients)
  }

  # 1) If coefficients are NOT provided or input vectors are NOT conformable, estimate the corresponding regression model
  if (is.null(coefficients) ||  non_conformable) {
    if (is.null(Z)) {
      stop("'Z' is required when fitting a bivariate spline without supplied coefficients, ",
           "or when supplied coefficients are not conformable with the spline basis.",
           call. = FALSE)
    }

    # Substract offset (if any) from Z
    Z0 <- Z - offset
    # Fit linear model without intercept, using weights
    tmp <- lm(Z0 ~ -1 + basisMatrixbiv2, weights = weights)
    # Extract fitted coefficients
    theta <- coef(tmp)
    # Avoid issues if there are NA values in the coefficients:
    if (any(is.na(theta))) theta[is.na(theta)] <- 0
    # Compute predicted values
    predicted <- basisMatrixbiv2 %*% theta + offset
    # Calculate residuals
    resid <- Z - predicted

  # 2) If coefficients are provided and conformable with InterKnotsX/InterKnotsY, compute the corresponding predicted values
  } else {
    tmp <- NULL
    theta <- coefficients
    # Compute predicted values
    predicted <- basisMatrixbiv2 %*% theta + offset
    # Calculate residuals
    resid <- NA
  }



  out <- list("theta"= theta,"predicted"= predicted,
              "residuals"= resid,"rss" = .weighted_rss(resid, weights),
              "Xbasis"= basisMatrixX, "Ybasis" = basisMatrixY,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,n))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,n))),
              "temporary"=tmp)
  return(out)
}

#' @importFrom splines splineDesign
#' @importFrom stats glm glm.fit coef
SplineReg_biv_GLM <- function(X, Y, Z, W = NULL, offset = rep(0,nobs), weights = rep(1,length(X)),
                              InterKnotsX, InterKnotsY, n, Xextr = range(X), Yextr = range(Y),
                              flag = TRUE, center = c(sum(Xextr)/2,sum(Yextr)/2),
                              family, mustart, inits = NULL, etastart = NULL,
                              coefficients = NULL, fast = FALSE)
  {
  # Check that 'n' (spline order) has length 1
  if(length(n) != 1) stop("'n' must have length 1")
  ord <- n # to avoid problem in use of family$initialize e.g. binomial()

  # Set required environment variables for family$initialize and IRLSfit
  y <- Z; nobs <- NROW(Z)

  # Create spline basis basisMatrixs using specified knots, order, and evaluation points
  basisMatrixX <- splineDesign(knots = sort(c(InterKnotsX,rep(Xextr,n))), derivs = rep(0,length(X)),
                           x = X, ord = n, outer.ok = T)
  basisMatrixY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))), derivs = rep(0,length(Y)),
                           x = Y, ord = n, outer.ok = T)

  basisMatrixbiv <- tensorProd(basisMatrixX, basisMatrixY)
  # Combine spline basis with parametric design matrix (if provided)
  basisMatrixbiv2 <- cbind(basisMatrixbiv,W)

  # If coefficients vector was provided, check whether this is conformable with the knots vectors, o.w. re-estimate the coefficients
  if(!is.null(coefficients)){
    non_conformable <- NCOL(basisMatrixbiv2) != length(coefficients)
  }

  # 1) If coefficients are NOT provided or input vectors are NOT conformable, estimate the corresponding regression model
  if (is.null(coefficients) ||  non_conformable) {

    # Initialize mustart based on input or defaults
    if(missing(mustart)||is.null(mustart)) {

      if (is.null(inits)) {
        if(is.null(etastart)) {
          # Set environment to parent frame
          env <- parent.frame()
          eval(family$initialize)
          mustart <- env$mustart
          } else {
            mustart <- family$linkinv(etastart)
            }

        } else {
          # Validate length of 'inits'
          if(length(inits)!= NCOL(basisMatrixbiv2)) stop("'inits' has wrong length")
          # Calculate initial mustart based on 'inits' (initial value for spline coefficients)
          mustart <- family$linkinv(basisMatrixbiv2 %*% inits)
        }
    }

    tmp <- tryCatch({
      if (fast) {
        glm.fit(x = basisMatrixbiv2, y = Z, family = family,
                weights = weights, mustart = mustart, offset = offset,
                intercept = FALSE)
      } else {
        glm(Z ~ -1 + basisMatrixbiv2, family = family, weights = weights,
            mustart = mustart, offset = offset)
      }
    }, error = function(e) {
      # If glm throws an error, fallback to IRLSfit
      IRLSfit(basisMatrixbiv2, Z, offset = offset,
              family = family, mustart = mustart, weights = weights)
    })

    # Extract fitted coefficients
    theta <- coef(tmp)
    # Avoid issues if there are NA values in the coefficients:
    if(any(is.na(theta))) theta[is.na(theta)] <- 0
    # Compute predicted values
    predicted <- family$linkinv(basisMatrixbiv2 %*% theta + offset)

  # 2) If coefficients are provided and conformable with InterKnotsX/InterKnotsY, compute the corresponding predicted values
  } else {
    tmp <- NULL
    theta <- coefficients
    predicted <- family$linkinv(basisMatrixbiv2 %*% theta + offset)
  }

  # Calculate residuals
  resid <- Z - predicted

  deviance <- if (is.null(tmp)) {
    sum(family$dev.resids(Z, predicted, weights))
  } else {
    tmp$deviance
  }

  out <- list("theta" = theta, "predicted" = predicted,
              "residuals" = resid, "rss" = deviance,
              "Xbasis" = basisMatrixX, "Ybasis" = basisMatrixY,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,ord))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,ord))),
              "temporary" = tmp )
  return(out)
}

