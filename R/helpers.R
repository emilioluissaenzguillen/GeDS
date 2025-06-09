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

makeNewMatr <- function(basisMatrix, tab, by.row=F){
  if(is.null(tab)){
    ret <- basisMatrix
  } else {
    recurr <- if(by.row) c(t(tab)) else c(tab)
    recurr <- recurr[recurr!=0]
    ids <- cumsum(recurr)
    ids<- c(0,ids)
    newres <- numeric((length(ids)-1))
    newX <- numeric((length(ids)-1))
    newY <- numeric((length(ids)-1))
    for(i in 1:(length(ids)-1)){
      newres[i] <- sum(basisMatrix[(ids[i]+1):ids[i+1],3])
      newX[i] <- basisMatrix[i,1]
      newY[i] <- basisMatrix[i,2]
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
  basisMatrix <- splineDesign(knots = sort(c(InterKnots, rep(extr,n))),
                          derivs = rep(0,length(X)), x = X, ord = n, outer.ok = TRUE)
  basisMatrix2 <- cbind(basisMatrix,Z)
  
  Y0 <- Y-offset
  tmp <-  if(all(weights==1)) .lm.fit(basisMatrix2, Y0) else lm.wfit.light(basisMatrix2, Y0, weights)
  theta <- coef(tmp)
  predicted <- basisMatrix2 %*% theta + offset
  resid <- Y - predicted
  out <- list("Theta" = theta,"Predicted" = predicted,
              "Residuals" = resid, "RSS" = as.numeric(crossprod(resid)),
              "Basis" = basisMatrix2,
              "temporary" = tmp)
  return(out)
}

################################################################################

SplineReg_fast_biv <- function(X, Y, Z, W=NULL, weights = rep(1, length(X)),
                               InterKnotsX, InterKnotsY, n, Xextr = range(X),
                               Yextr = range(Y), flag=TRUE, 
                               center = c(sum(Xextr)/2,sum(Yextr)/2))
  {
  basisMatrixX <- splineDesign(knots = sort(c(InterKnotsX, rep(Xextr,n))), derivs = rep(0,length(X)),
                           x = X, ord = n, outer.ok = TRUE)
  basisMatrixY <- splineDesign(knots = sort(c(InterKnotsY, rep(Yextr,n))), derivs = rep(0,length(Y)),
                           x = Y, ord = n, outer.ok = T)
  basisMatrixbiv <- tensorProd(basisMatrixX, basisMatrixY)
  basisMatrixbiv2 <- cbind(basisMatrixbiv, W)
  
  
  #fff <- !rankMatrix(basisMatrixbiv)==8
  #  Xknots<-makenewknots(sort(c(InterKnotsX,rep(Xextr,n-1))),deg=n)
  #  Yknots<-makenewknots(sort(c(InterKnotsY,rep(Yextr,n-1))),deg=n)
  if(all(weights==1)){
    tmp <- .lm.fit(basisMatrixbiv2, Z)
    if (tmp$rank<ncol(basisMatrixbiv2)){
      tmp <- lm.fit(basisMatrixbiv2, Z)
    }
    resid <- residuals(tmp)
  } else {
    tmp <- lm.wfit.light(basisMatrixbiv2, Z, weights) #ccc<-lm(Z ~ -1+basisMatrixbiv) #
    resid <- tmp$residuals
    
    
    if (tmp$rank<ncol(basisMatrixbiv2)){
      tmp <- lm.wfit(basisMatrixbiv2, Z, weights)
      resid <- residuals(tmp)
      
    }
  }
  theta <- as.numeric(coef(tmp))
  #  theta[is.na(theta)] <- 0
  out <- list("Theta" = theta, "Predicted" = basisMatrixbiv2 %*% theta,
              "Residuals" = resid, "RSS" = as.numeric(crossprod(resid)),
              "XBasis" = basisMatrixX, "YBasis" = basisMatrixY, #"Poly"=poly,
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

################################################################################


CI <- function(tmp, resid, prob = 0.95, basisMatrix, basisMatrix2, predicted,
               n_obs = NROW(basisMatrix),
               type = "lm",
               huang = TRUE) {
  
  if (type == "lm") {
    # Residual standard error
    df <- if(!is.null(tmp)) tmp$df.residual else as.numeric(nrow(basisMatrix2) - rankMatrix(basisMatrix2)) # residual degrees of freedom
    sigma_hat <- sqrt(sum(resid^2)/df)
    # Adjust probability for two-tailed test
    prob <- 1-.5*(1-prob)
    # Diagonal of the hat matrix
    H_diag <- stats::hat(basisMatrix2, intercept = FALSE) # or influence(tmp)$hat
    # CI_j =\hat{y_j} ± t_{α/2,df}*\hat{σ}*\sqrt{H_{jj}}; H = X(X'X)^{−1}X'
    band <- qt(prob,df) * sigma_hat * H_diag^.5
    
    NCI = list("Upp" = predicted + band, "Low" = predicted - band)
    
    # Huang (2003) method for confidence band width (see Theorem 6.1)
    band_width_huang <- ACI <- NULL; dim_threshold = 1500
    if (huang && n_obs < dim_threshold && NCOL(basisMatrix) != 0) {
      # i. E_n[B(X)B^t(X)] = (1/n)*\sum_{i=1}^nB(X_i)B^t(X_i)
      matcb <- crossprod(basisMatrix) / n_obs
      matcbinv <- tryCatch({
        chol2inv(chol(matcb))  # Fastest if SPD
      }, error = function(e1) {
        message("SplineReg_LM, Huang CI: Matrix not SPD, using solve().")
        tryCatch({
          solve(matcb)
        }, error = function(e2) {
          message("SplineReg_LM, Huang CI: Matrix singular, using ginv().")
          MASS::ginv(matcb)
        })
      })
      # ii. Var(\hat{f} | X) = (1/n)*B^t(x) * E_n[B(X)B^t(X)]^-1 * B(x) * \hat{σ}^2
      S <- basisMatrix %*% matcbinv
      conditionalVariance <- (sigma_hat^2 / n_obs) * rowSums(S * basisMatrix)
      # iii. ± z_{1-α/2} * Var(\hat{f} | X)
      band_width_huang <- qnorm(prob) * sqrt(conditionalVariance)
      
      ACI = list("Upp" = predicted + band_width_huang,
                 "Low" = predicted - band_width_huang)
    }
    
    
    
  } else if (type == "glm") {
    
    if (is.numeric(tmp$coefficients)) {
      alpha <- 1 - prob
      z_val <- qnorm(1 - alpha / 2)  # For Wald-type CI
      
      # eta_hat <- predict(tmp, type = "link", se.fit = TRUE)
      eta_hat <- tryCatch(
        predict(tmp, type = "link", se.fit = TRUE),
        
        error = function(e) {
          
          eta <- predict(tmp, type = "link")
          
          matcb <- t(basisMatrix2) %*% diag(tmp$weights) %*% basisMatrix2
          Sigma <- summary(tmp)$dispersion * MASS::ginv(matcb)
          se_eta   <- sqrt(rowSums((basisMatrix2 %*% Sigma) * basisMatrix2))
          
          list(fit = eta, se.fit = se_eta)
        }
      )
      
      lower_eta <- eta_hat$fit - z_val * eta_hat$se.fit
      upper_eta <- eta_hat$fit + z_val * eta_hat$se.fit
      
      lower <- tmp$family$linkinv(lower_eta)
      upper <- tmp$family$linkinv(upper_eta)
      
      NCI = list("Upp" = upper, "Low" = lower)
      
      } else {
        # tmp$coefficients == "When using bivariate base-learners, the 'single spline representation' (in pp form or B-spline form) of the boosted fit is not available."
        NCI = NULL
      }
    
    ACI = NULL
    
  }
  
  return(list(NCI = NCI, ACI = ACI))
  
}

