SplineReg_biv <- function(X ,Y , Z, W = NULL, offset = rep(0,length(X)), weights = rep(1,length(X)),
                          InterKnotsX, InterKnotsY, n, Xextr = range(X), Yextr = range(Y), flag = TRUE,
                          center = c(sum(Xextr)/2,sum(Yextr)/2), coefficients = NULL)
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
  
  # If coefficients vector was provided, check whether this is conformable with the knots vectors, o.w. re-estimate the coefficients
  if (!is.null(coefficients)) {
    non_conformable <- (length(InterKnotsX) + 2) *  (length(InterKnotsY) + 2) != length(coefficients)
  }
  
  # 1) If coefficients are NOT provided or input vectors are NOT conformable, estimate the corresponding regression model
  if (is.null(coefficients) ||  non_conformable) {
    
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
  
  
  
  out <- list("Theta"= theta,"Predicted"= predicted,
              "Residuals"= resid,"RSS" = as.numeric(crossprod(resid)),
              "XBasis"= basisMatrixX, "YBasis" = basisMatrixY,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,n))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,n))),
              "temporary"=tmp)
  return(out)
}


SplineReg_biv_GLM <- function(X, Y, Z, W = NULL, offset = rep(0,nobs), weights = rep(1,length(X)),
                              InterKnotsX, InterKnotsY, n, Xextr = range(X), Yextr = range(Y),
                              flag = TRUE, center = c(sum(Xextr)/2,sum(Yextr)/2), 
                              family, mustart, inits = NULL, etastart = NULL, coefficients = NULL)
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
    non_conformable <- (length(InterKnotsX) + 2) *  (length(InterKnotsY) + 2) != length(coefficients)
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
    
    tmp <- IRLSfit(basisMatrixbiv2, Z, offset = offset,
                   family = family, mustart = mustart, weights = weights)
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
  
  out <- list("Theta" = theta, "Predicted" = predicted,
              "Residuals" = resid,"RSS" = tmp$lastdeviance, "deviance" = tmp$deviance,
              "XBasis" = basisMatrixX, "YBasis" = basisMatrixY,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,ord))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,ord))),
              "temporary" = tmp )
  return(out)
}

