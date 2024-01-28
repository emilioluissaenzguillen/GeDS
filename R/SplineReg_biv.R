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
  if (is.null(coefficients) ||  non_conformable){
    tmp <- lm.wfit(matricebiv2, Z, weights)
    if(any(is.na(coef(tmp)))) {
      warning("NAs in regression coefficients")
      tmp <- lm(Z~-1+matricebiv2,weights =  weights)
    }
    resid <- residuals(tmp)
    theta <- as.numeric(coef(tmp))
    # Avoid issues if there are NA values in the coefficients:
    if(any(is.na(theta))) {theta[is.na(theta)] <- 0}
    predicted <- matricebiv2%*%theta
  
  # 2) If coefficients are provided and conformable with InterKnotsX/InterKnotsY, compute the corresponding predicted values
  } else {
    tmp <- NULL
    resid <- NA
    theta <- coefficients
    predicted <- matricebiv2%*%theta
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

