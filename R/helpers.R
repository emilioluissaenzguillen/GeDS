################################################################################

lm.wfit.light <- function (x, y, w, tol = 1e-07) {
  x.asgn <- attr(x, "assign")
  zero.weights <- any(w == 0)
  save.r <- y
  save.f <- x
  
  if (zero.weights) {
    save.w <- w
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x <- x[ok,  ,drop = FALSE]
    n <- nrow(x)
    y <-y[ok]
  }
  wts <- sqrt(w)
  out <- .lm.fit(x * as.numeric(wts), y * wts, tol)
  out$residuals <- save.r - save.f%*%out$coef
  return(out)
}

################################################################################

makeNewMatr <- function(matrice, tab, by.row=F){
  if(is.null(tab)){
    ret <- matrice
  } else {
    recurr <-if(by.row) c(t(tab)) else c(tab)
    recurr <- recurr[recurr!=0]
    ids <- cumsum(recurr)
    ids<- c(0,ids)
    newres <- numeric((length(ids)-1))
    newX <- numeric((length(ids)-1))
    newY <- numeric((length(ids)-1))
    for(i in 1:(length(ids)-1)){
      newres[i] <- sum(matrice[(ids[i]+1):ids[i+1],3])
      newX[i] <- matrice[i,1]
      newY[i] <- matrice[i,2]
    }
    ret <- cbind(newX,newY,newres)
  }
  return(ret)
}

################################################################################
# already fast - useless to do in C++
makeNewRes <- function(resold, recurr){
  ids <- cumsum(recurr)
  ids<- c(0,ids)
  newres <- numeric((length(ids)-1))
  for(i in 1:(length(ids)-1)){
    newres[i] <- sum(resold[(ids[i]+1):ids[i+1]])/recurr[i]
  }
  return(newres)
}

makeNewRes2 <- function(resold, recurr, weights){
  ids <- cumsum(recurr)
  ids<- c(0,ids)
  newres <- numeric((length(ids)-1))
  newweights <- numeric((length(ids)-1))
  
  for(i in 1:(length(ids)-1)){
    newres[i] <- sum(resold[(ids[i]+1):ids[i+1]])/recurr[i]
    newweights[i] <- sum(weights[(ids[i]+1):ids[i+1]])
  }
  newres <- newres*newweights
  return(newres)
}

################################################################################

SplineReg_fast_weighted_zed <- function(X, Y, Z, offset,
                                        weights = rep(1, length(X)), InterKnots,
                                        n, extr = range(X))
{
  matrice <- splineDesign(knots = sort(c(InterKnots, rep(extr,n))),
                          derivs = rep(0,length(X)), x = X, ord = n, outer.ok = TRUE)
  matrice2 <- cbind(matrice,Z)
  
  Y0 <- Y-offset
  tmp <-  if(all(weights==1)) .lm.fit(matrice2, Y0) else lm.wfit.light(matrice2, Y0, weights)
  theta <- coef(tmp)
  predicted <- matrice2 %*% theta + offset
  resid <- Y - predicted
  out <- list("Theta" = theta,"Predicted" = predicted,
              "Residuals" = resid, "RSS" = t(resid)%*%resid,
              "Basis" = matrice2,
              "temporary" = tmp)
  return(out)
}

################################################################################

SplineReg_fast_biv <- function(X, Y, Z, W=NULL, weights = rep(1, length(X)),
                               InterKnotsX, InterKnotsY, n, Xextr = range(X),
                               Yextr = range(Y), flag=TRUE, 
                               center = c(sum(Xextr)/2,sum(Yextr)/2))
  {
  matriceX <- splineDesign(knots = sort(c(InterKnotsX, rep(Xextr,n))), derivs = rep(0,length(X)),
                           x = X, ord = n, outer.ok = TRUE)
  matriceY <- splineDesign(knots = sort(c(InterKnotsY, rep(Yextr,n))), derivs = rep(0,length(Y)),
                           x = Y, ord = n, outer.ok = T)
  matriceY_noint <- cut_int(matriceY)
  matricebiv <- tensorProd(matriceX, matriceY_noint)
  matricebiv2 <- cbind(matricebiv, W)
  
  
  #fff <- !rankMatrix(matricebiv)==8
  #  Xknots<-makenewknots(sort(c(InterKnotsX,rep(Xextr,n-1))),deg=n)
  #  Yknots<-makenewknots(sort(c(InterKnotsY,rep(Yextr,n-1))),deg=n)
  if(all(weights==1)){
    tmp <- .lm.fit(matricebiv2, Z)
    if (tmp$rank<ncol(matricebiv2)){
      tmp <- lm.fit(matricebiv2, Z)
    }
    resid <- residuals(tmp)
  } else {
    tmp <- lm.wfit.light(matricebiv2, Z, weights) #ccc<-lm(Z ~ -1+matricebiv) #
    resid <- tmp$residuals
    
    
    if (tmp$rank<ncol(matricebiv2)){
      tmp <- lm.wfit(matricebiv2, Z, weights)
      resid <- residuals(tmp)
      
    }
  }
  theta <- as.numeric(coef(tmp))
  #  theta[is.na(theta)] <- 0
  out <- list("Theta" = theta, "Predicted" = matricebiv2%*%theta,
              "Residuals" = resid, "RSS" = t(resid)%*%resid,
              "XBasis" = matriceX, "YBasis" = matriceY_noint, #"Poly"=poly,
              "temporary" = tmp)
  return(out)
}

################################################################################

newknot.guess <- function(intknots, extr, guess, newknot) {
  # i. Determine the position of the new knot relative to existing internal knots
  newknot.position <- sum(intknots < as.numeric(newknot))
  
  # ii. Generate a spline design matrix for the new knot
  nk.design <- splineDesign(knots=sort(c(intknots, rep(extr,2))), derivs = 0,
                            x = newknot, ord = 2, outer.ok = TRUE)
  
  # iii. Calculate a new guess value (coefficient) for the spline at the new knot position
  pr.value <- sum(nk.design * guess)
  newguess <- pr.value
  
  # iv. Update the guess-coefficient vector based on the position of the new internal knot
  # (keep in mind the number of B-splines is p = k + 2, where k is the number of internal knots)
  if(newknot.position == 0) {
    # If the new knot is before all internal knots, insert the new guess just after the guess-coefficient for the lower boundary knot
    guess <- c(guess[1], newguess, guess[-1])
  } else if(newknot.position == length(intknots)) {
    # If the new knot is after all internal knots, insert the new guess just before the guess-coefficient for the upper boundary knot
    guess <- c(guess[1:(newknot.position+1)], newguess, guess[newknot.position+2])
  } else {
    # Otherwise, insert the new guess in its appropriate position
    guess <- c(guess[1:(newknot.position+1)], newguess, guess[-(1:(newknot.position+1))])
  }
  
  return(guess)
}

# newknot.guess_biv <- function(Dim, FixedDim, Dim.intknots, FixedDim.intknots, Dim.extr, FixedDim.extr, guess, Dim.newknot) {
#   # i. Determine the position of the new knot relative to existing internal knots
#   Dim.newknot_position <- sum(Dim.intknots < as.numeric(Dim.newknot))
#   
#   # ii. Generate a spline design matrix for the new knot
#   Dim.nk.design <- splineDesign(knots = sort(c(Dim.intknots,rep(Dim.extr, 2))), derivs = 0,
#                                 x = Dim.newknot, ord = 2, outer.ok = TRUE)
#   FixedDim.design <- splineDesign(knots = sort(c(FixedDim.intknots,rep(FixedDim.extr, 2))), derivs = rep(0,length(FixedDim)),
#                                   x = FixedDim, ord = 2, outer.ok = TRUE)
#   
#   nk.design <- tensorProd(Dim.nk.design, FixedDim.design)
#   
#   # iii. Calculate a new guess value (coefficient) for the spline at the new knot position
#   pr.value <- sum(nk.design * guess)
#   newguess <- pr.value
#   
#   # iv. Update the guess-coefficient vector based on the position of the new internal knot
#   # (keep in mind the number of B-splines is p = k + 2, where k is the number of internal knots)
#   if(newknot.position == 0) {
#     # If the new knot is before all internal knots, insert the new guess just after the guess-coefficient for the lower boundary knot
#     guess <- c(guess[1], newguess, guess[-1])
#   } else if(newknot.position == length(intknots)) {
#     # If the new knot is after all internal knots, insert the new guess just before the guess-coefficient for the upper boundary knot
#     guess <- c(guess[1:(newknot.position+1)], newguess, guess[newknot.position+2])
#   } else {
#     # Otherwise, insert the new guess in its appropriate position
#     guess <- c(guess[1:(newknot.position+1)], newguess, guess[-(1:(newknot.position+1))])
#   }
#   
#   return(guess)
# }

