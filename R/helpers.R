################################################################################
#' @importFrom stats .lm.fit
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

#' @importFrom splines splineDesign
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
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @importFrom stats hat qt qnorm
ci <- function(tmp, resid, prob = 0.95, basisMatrix, basisMatrix2, predicted,
               n_obs = NROW(basisMatrix),
               type = "lm",
               huang = TRUE) {

  if (type == "lm") {
    # Residual standard error
    df <- if(!is.null(tmp)) tmp$df.residual else as.numeric(nrow(basisMatrix2) - rankMatrix(basisMatrix2)) # residual degrees of freedom
    # Saturated fit (df <= 0): t-based bands are undefined; avoid spurious "NaNs produced".
    valid_df <- is.finite(df) && df > 0
    sigma_hat <- if (valid_df) sqrt(sum(resid^2)/df) else NA_real_
    # Adjust probability for two-tailed test
    prob <- 1-.5*(1-prob)
    # Diagonal of the hat matrix
    H_diag <- hat(basisMatrix2, intercept = FALSE) # or influence(tmp)$hat
    # CI_j =\hat{y_j} ± t_{α/2,df}*\hat{σ}*\sqrt{H_{jj}}; H = X(X'X)^{-1}X'
    band <- if (valid_df) qt(prob, df) * sigma_hat * H_diag^.5 else rep(NA_real_, length(predicted))

    nci = list("Upp" = predicted + band, "Low" = predicted - band)

    # Huang (2003) method for confidence band width (see Theorem 6.1)
    band_width_huang <- aci <- NULL; dim_threshold = 1500
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
          ginv(matcb)
        })
      })
      # ii. Var(\hat{f} | X) = (1/n)*B^t(x) * E_n[B(X)B^t(X)]^-1 * B(x) * \hat{σ}^2
      S <- basisMatrix %*% matcbinv
      conditionalVariance <- (sigma_hat^2 / n_obs) * rowSums(S * basisMatrix)
      # iii. ± z_{1-α/2} * Var(\hat{f} | X)
      # A variance is non-negative; clamp numerical noise (e.g. from ginv on an
      # ill-conditioned basis) to avoid spurious sqrt() NaNs.
      band_width_huang <- qnorm(prob) * sqrt(pmax(conditionalVariance, 0))

      aci = list("Upp" = predicted + band_width_huang,
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
          Sigma <- summary(tmp)$dispersion * ginv(matcb)
          se_eta   <- sqrt(rowSums((basisMatrix2 %*% Sigma) * basisMatrix2))

          list(fit = eta, se.fit = se_eta)
        }
      )

      lower_eta <- eta_hat$fit - z_val * eta_hat$se.fit
      upper_eta <- eta_hat$fit + z_val * eta_hat$se.fit

      lower <- tmp$family$linkinv(lower_eta)
      upper <- tmp$family$linkinv(upper_eta)

      nci = list("Upp" = upper, "Low" = lower)

      } else {
        # tmp$coefficients == "When using bivariate base-learners, the 'single spline representation' (in pp form or B-spline form) of the boosted fit is not available."
        nci = NULL
      }

    aci = NULL

  }

  return(list(nci = nci, aci = aci))

}

################################################################################
stopping_rule <- function(
    j, q, n_starting_intknots,
    rssnew,                           # numeric vector
    flag, phis,                       # logical, numeric vector
    stoptype = c("SR","RD","LR"),     # "Smoothed Ratio", "Ratio of Deviances", "Likelihood Ratio"
    intknots,                         # vector (for kappa = length(intknots))
    min.intknots,                     # scalar
    phi,                              # threshold: φ_exit (SR/RD) or χ^2 tail prob (LR)
    phis_star,                        # numeric vector to append (SR)
    oldintc, oldslp                   # numeric vectors to append (SR)
) {

  stoptype <- match.arg(stoptype)

  if (missing(n_starting_intknots)) n_starting_intknots <- 0
  if (missing(flag)) flag <- FALSE

  # Default return (no action)
  out <- list(
    should_break = FALSE,
    phis = phis,
    phis_star = phis_star,
    oldintc = oldintc,
    oldslp = oldslp,
    prnt = ""
  )

  # Guard: only active after enough steps
  if (j <= q + n_starting_intknots) return(out)

  # Early check: if ratio > 1 then stop (i.e. rss getting worse when adding more knots)
  ratio_j <- rssnew[j] / rssnew[j - q]
  if (ratio_j > 1) {
    out$should_break <- TRUE
    return(out)
  }

  # Adding the current ratio of deviances to the 'phis' vector
  if (flag) phis <- phis[1:(j - q - 1)]
  phis <- if (stoptype == "LR") {
    c(phis, rssnew[j - q] - rssnew[j])
  } else {
    c(phis, ratio_j)
  }

  # If not past min.intknots, continue with the iterations
  out$phis <- phis
  if ((j - q) <= min.intknots) return(out)

  # (I) Smoothed Ratio of deviances
  if(stoptype == "SR") {
    # \hat{φ}_κ = 1 - exp{\hat{γ}_0 + \hat{γ}_1*κ}
    # 1-\hat{φ}_κ = exp{\hat{γ}_0 + \hat{γ}_1*κ}
    # ln(1-\hat{φ}_κ) = \hat{γ}_0 + \hat{γ}_1*κ
    # Fit a linear model ln(1-φ) ~ \hat{γ}_0 + \hat{γ}_1*κ to the sample {φ_h, h}^κ_{h=q}
    phismod <- log(1-phis)
    kappa <- length(intknots)
    gamma <- .lm.fit(cbind(1, q:kappa), phismod)$coef

    # Calculate \hat{φ}_κ based on the estimated coefficients
    phi_kappa <- 1 - exp(gamma[1])*exp(gamma[2]*kappa)

    # Store \hat{φ}_κ and the estimated coefficients \hat{γ}_0 and \hat{γ}_1
    phis_star <- c(phis_star, phi_kappa)
    oldintc <- c(oldintc, gamma[1])
    oldslp <- c(oldslp, gamma[2])

    # Creating a print statement that shows the current adjusted phi value
    prnt <- paste0(", phi_hat = ", round(phi_kappa, 3))

    out$phis_star <- phis_star
    out$oldintc <- oldintc
    out$oldslp <- oldslp
    out$prnt <- prnt
    # Check if \hat{φ}_κ ≥ φ_{exit}
    if (phi_kappa >= phi) out$should_break <- TRUE
    return(out)

  # (II) Ratio of Deviances
  } else if (stoptype == "RD") {
    prnt <- paste0(", phi = ", round(ratio_j, 3))
    out$prnt <- prnt
    if (ratio_j >= phi) out$should_break <- TRUE
    return(out)

    # (III) Likelihood Ratio
  } else if (stoptype == "LR") {
    # stat = -(rss_j - rss_{j-q}); stop if stat < qchisq(phi, df=q)
    stat <- -(rssnew[j] - rssnew[j - q])
    prnt <- paste0(", p = ", round(pchisq(stat, df = q), 3))
    out$prnt <- prnt
    if (stat < qchisq(phi, df = q)) out$should_break <- TRUE
    return(out)
  }

  out
}

################################################################################
validate_GeDS_order <- function(n) {
  n <- suppressWarnings(as.integer(n))

  if (length(n) != 1L || is.na(n) || !(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.", call. = FALSE)
  }

  n
}

################################################################################
make_shape_constraint <- function(knots, n, p,
                                  shape_constraint = "none",
                                  eps = 0) {
  allowed <- c("none", "increasing", "decreasing", "convex", "concave")

  shape_constraint <- match.arg(
    shape_constraint,
    choices = allowed,
    several.ok = TRUE
  )

  if (length(shape_constraint) == 0 || identical(shape_constraint, "none")) {
    return(NULL)
  }

  if (any(!shape_constraint %in% allowed)) {
    stop("Unknown shape_constraint.", call. = FALSE)
  }

  if ("none" %in% shape_constraint && length(shape_constraint) > 1) {
    stop("'none' cannot be combined with other shape constraints.", call. = FALSE)
  }

  if (all(c("increasing", "decreasing") %in% shape_constraint)) {
    stop("Cannot impose both increasing and decreasing constraints.", call. = FALSE)
  }

  if (all(c("convex", "concave") %in% shape_constraint)) {
    stop("Cannot impose both convex and concave constraints.", call. = FALSE)
  }

  constraints <- list()

  # First derivative coefficient matrix
  # f'(x) = sum_i D1_i(theta) N_{i,n-1}(x)
  D1 <- matrix(0, nrow = p - 1, ncol = p)

  for (i in 2:p) {
    denom <- knots[i + n - 1] - knots[i]

    if (denom <= 0) {
      stop("Invalid knot sequence for derivative constraints.", call. = FALSE)
    }

    D1[i - 1, i]     <-  (n - 1) / denom
    D1[i - 1, i - 1] <- -(n - 1) / denom
  }

  if ("increasing" %in% shape_constraint) {
    constraints[["increasing"]] <- D1
  }

  if ("decreasing" %in% shape_constraint) {
    constraints[["decreasing"]] <- -D1
  }

  # Second derivative coefficient matrix
  # f''(x) = sum_i D2_i(theta) N_{i,n-2}(x)
  if (any(c("convex", "concave") %in% shape_constraint)) {
    if (n < 3) {
      stop("Convexity/concavity constraints require spline order n >= 3.",
           call. = FALSE)
    }

    D2 <- matrix(0, nrow = p - 2, ncol = p)

    for (i in 3:p) {
      denom <- knots[i + n - 2] - knots[i]

      if (denom <= 0) {
        stop("Invalid knot sequence for second-derivative constraints.",
             call. = FALSE)
      }

      D2[i - 2, ] <- (n - 2) * (D1[i - 1, ] - D1[i - 2, ]) / denom
    }

    if ("convex" %in% shape_constraint) {
      constraints[["convex"]] <- D2
    }

    if ("concave" %in% shape_constraint) {
      constraints[["concave"]] <- -D2
    }
  }

  C <- do.call(rbind, constraints)
  b <- rep(eps, nrow(C))

  list(C = C, b = b)
}

constrained_wls <- function(X, y, weights = rep(1, length(y)),
                            C, b, ridge = 1e-8) {

  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required for shape-constrained fitting.",
         call. = FALSE)
  }

  w_sqrt <- sqrt(as.numeric(weights))
  Xw <- X * w_sqrt
  yw <- y * w_sqrt

  Dmat <- crossprod(Xw)
  dvec <- as.numeric(crossprod(Xw, yw))

  # Ensure symmetry and numerical positive definiteness
  Dmat <- (Dmat + t(Dmat)) / 2
  Dmat <- Dmat + diag(ridge, ncol(Dmat))

  fit <- tryCatch(
    quadprog::solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = t(C),
      bvec = b,
      meq = 0
    ),
    error = function(e) {
      stop("Shape-constrained weighted least squares failed: ",
           conditionMessage(e),
           call. = FALSE)
    }
  )

  as.numeric(fit$solution)
}
