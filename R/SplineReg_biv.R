SplineReg_biv <- function(X,Y,Z,W=NULL,weights=rep(1,length(X)),InterKnotsX,InterKnotsY,n,
                          Xextr=range(X),Yextr=range(Y),flag=TRUE,center=c(sum(Xextr)/2,sum(Yextr)/2),
                          coefficients = NULL){
  matriceX <- splineDesign(knots=sort(c(InterKnotsX,rep(Xextr,n))),derivs=rep(0,length(X)),
                           x=X,ord=n,outer.ok = T)
  matriceY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))),derivs=rep(0,length(Y)),
                           x=Y,ord=n,outer.ok = T)
  matriceY_noint <- cut_int(matriceY)
  
  matricebiv <- tensorProd(matriceX,matriceY_noint)
  matricebiv2 <- cbind(matricebiv,W)
  
  # If coefficients vector was provided, check whether this is conformable with the knots vectors, o.w. re-estimate the coefficients
  if(!is.null(coefficients)){
    non_conformable <- (length(InterKnotsX) + 2) *  (length(InterKnotsY) + 2) != length(coefficients)
  }
  
  # 1) If coefficients are NOT provided or input vectors are NOT conformable, estimate the corresponding regression model
  if (is.null(coefficients) ||  non_conformable) {
    tmp <- lm.wfit(matricebiv2, Z, weights)
    if(any(is.na(coef(tmp)))) {
      warning("NAs in regression coefficients")
      tmp <- lm(Z~-1+matricebiv2,weights =  weights)
    }
    theta <- as.numeric(coef(tmp))
    # Avoid issues if there are NA values in the coefficients:
    if(any(is.na(theta))) {theta[is.na(theta)] <- 0}
    predicted <- matricebiv2%*%theta
    resid <- residuals(tmp)
  
  # 2) If coefficients are provided and conformable with InterKnotsX/InterKnotsY, compute the corresponding predicted values
  } else {
    tmp <- NULL
    theta <- coefficients
    predicted <- matricebiv2%*%theta
    resid <- NA
  }
  
  out <- list("Theta"= theta,"Predicted"= predicted,
              "Residuals"= resid,"RSS" = t(resid)%*%resid,
              "XBasis"= matriceX, "YBasis" = matriceY_noint,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,n))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,n))),
              "temporary"=tmp)
  return(out)
}


cut_int <- function(mat){
  d<-dim(mat)
  #we delete the intercept
  #otherwise tensor product basis is rank deficient
  mat_star <- t(rep(1,d[1])%*%mat)
  Q_star<-(qr.Q(qr(mat_star),complete=T))[,-1]
  #  return(mat%*%Q_star)
  return(mat)
}

SplineReg_biv_GLM <- function(X, Y, Z, W = NULL, offset = rep(0,nobs), weights = rep(1,length(Z)),
                              InterKnotsX, InterKnotsY, n, Xextr = range(X), Yextr = range(Y),
                              flag = TRUE, center = c(sum(Xextr)/2,sum(Yextr)/2), 
                              family, mustart, inits = NULL, etastart = NULL)
  {
  # Ensure X, Y, Z, InterKnots are numeric, and n is integer
  X           <- as.numeric(X)
  Y           <- as.numeric(Y)
  Z           <- as.numeric(Z)
  InterKnotsX <- as.numeric(InterKnotsX)
  InterKnotsY <- as.numeric(InterKnotsY)
  n           <- as.integer(n)
  
  # Check that 'n' (spline order) has length 1
  if(length(n) != 1) stop("'n' must have length 1")
  ord <- n # to avoid problem in use of family$initialize e.g. binomial()
  
  # Set required environment variables for family$initialize and IRLSfit
  y <- Z; nobs <- NROW(Z)
  
  # Create spline basis matrices using specified knots, order, and evaluation points
  matriceX <- splineDesign(knots = sort(c(InterKnotsX,rep(Xextr,n))), derivs = rep(0,length(X)),
                           x = X, ord = n, outer.ok = T)
  matriceY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))), derivs = rep(0,length(Y)),
                           x = Y, ord = n, outer.ok = T)
  matriceY_noint <- cut_int(matriceY)
  
  matricebiv <- tensorProd(matriceX,matriceY_noint)
  # Combine spline basis with parametric design matrix (if provided)
  matricebiv2 <- cbind(matricebiv,W)
  
  
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
      if(length(inits)!= NCOL(matricebiv2)) stop("'inits' has wrong length")
      # Calculate initial mustart based on 'inits' (initial value for spline coefficients)
      mustart <- family$linkinv(matricebiv2%*%inits)
    }
  }
  
  tmp <- IRLSfit(matricebiv2, Z, offset = offset,
                 family = family, mustart = mustart, weights = weights)
  theta <- coef(tmp)
  # Avoid issues if there are NA values in the coefficients:
  if(any(is.na(theta))) {theta[is.na(theta)] <- 0}
  predicted <- family$linkinv(matricebiv2%*%theta + offset)
  resid <- tmp$res2
  
  out <- list("Theta" = theta, "Predicted" = predicted,
              "Residuals" = resid,"RSS" = tmp$lastdeviance, "deviance" = tmp$deviance,
              "XBasis" = matriceX, "YBasis" = matriceY,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,ord))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,ord))),
              "temporary" = tmp )
  return(out)
}

