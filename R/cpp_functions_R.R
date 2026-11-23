################################################################################
# R reference implementations of the C++ functions in src/GeDS1.cpp.
#
# None of these is called by the fitters. They are kept for debugging and so
# that tests can compare them with the compiled versions. Each C++ function
# names its R reference in a comment.
#
#   C++                          R reference
#   Knotnew()                    Knotnew_R_clean() (current), Knotnew_R() (original)
#   findNewDimKnot()             findNewDimKnot_R()
#   tensorProd()                 tensorProd_R()
#   linear_spline_ls()           .linear_spline_ls_R()
#   geds_knot_step()             .geds_knot_step_R()
#   gd_neighbour_semivariance()  .gd_neighbour_semivariance_R()
#   gd_stage_a_path()            the Stage A loop of UnivariateFitter() (R/UnivariateFitter.R),
#                                which still runs whenever the C++ path does not apply
################################################################################

Knotnew_R_clean <- function(wht, restmp, x, dcm, oldknots, tol, support_order = 2) {

  # Number of residual clusters
  u <- length(dcm)
  if (u == 0 || length(wht) == 0) {
    return(list(newknot = NA_real_, cluster_index = NA_integer_))
  }

  candidate_wht <- as.numeric(wht)
  sw_sites <- sort(unique(x))
  # Number of old knots and internal knots
  noldint <- length(oldknots) - 2*support_order
  # Dimension of the x vector
  xdim <- length(x)
  # Initialize the index of the j_s-th cluster
  index <- NA
  # Initialize the new knot position
  newknot <- NA
  knot_found <- FALSE
  # Extract internal knots
  oldintknots <- if (noldint > 0) tail(oldknots, noldint) else numeric(0)

  for (kk in seq_len(u)) {

    # Index of the cluster with the maximum weight
    index <- which.max(candidate_wht)
    # If all weights are zero, assign triangular weights peaking at the centre
    if (candidate_wht[index] == 0) {
      n <- length(candidate_wht)
      center <- (n + 1) / 2
      candidate_wht <- 1 - abs(seq_len(n) - center) / center
      index <- which.max(candidate_wht)
    }

    # Cluster boundaries
    d_lower <- if (index == 1) 1 else dcm[index - 1] + 1
    d_upper <- dcm[index]
    inf <- x[d_lower]
    sup <- x[d_upper]

    # --- 1. Check for internal knot in the current cluster ---
    cluster_has_knot <- FALSE
    if (length(oldintknots) > 0) {
      if (inf != sup) {
        cluster_has_knot <- any(oldintknots >= inf & oldintknots <= sup)
      } else {
        cluster_has_knot <- any(abs(inf - oldintknots) < tol)
      }
    }

    # 2) If an internal knot already exists in this cluster, set its weight to zero
    if (cluster_has_knot) {
      candidate_wht[index] <- 0
      next
    }

    # --- 2. Propose new knot ---
    weights_cluster <- restmp[d_lower:d_upper]
    x_cluster <- x[d_lower:d_upper]
    newknot <- sum(weights_cluster * x_cluster) / sum(weights_cluster)

    sortedknots <- sort(c(oldknots, newknot))

    # --- 3. Minimum support followed by distinct-site SW matching ---
    support_valid <- all(vapply(
      seq_len(length(sortedknots) - support_order),
      function(i) {
        any(x > sortedknots[i] + tol &
              x < sortedknots[i + support_order] - tol)
      },
      logical(1)
    ))

    if (support_valid) {
      previous_site <- 0L
      n_basis <- length(sortedknots) - support_order
      for (i in seq_len(n_basis)) {
        interior <- sw_sites > sortedknots[i] + tol &
          sw_sites < sortedknots[i + support_order] - tol
        left_endpoint <- i == 1L && sortedknots[support_order] == sortedknots[1L]
        right_endpoint <- i == n_basis && sortedknots[n_basis + 1L] == tail(sortedknots, 1L)
        available <- interior |
          (left_endpoint & sw_sites == sortedknots[1L]) |
          (right_endpoint & sw_sites == tail(sortedknots, 1L))
        candidates <- which(seq_along(sw_sites) > previous_site & available)
        if (!length(candidates)) {
          support_valid <- FALSE
          break
        }
        previous_site <- candidates[1L]
      }
    }

    # --- 4. Check new knot is not a boundary knot ---
    newknot_is_internal <- isTRUE(
      newknot >= min(sortedknots, na.rm = TRUE) + tol &&
        newknot <= max(sortedknots, na.rm = TRUE) - tol
    )

    if (!support_valid || !newknot_is_internal) {
      candidate_wht[index] <- 0
    } else {
      knot_found <- TRUE
      break  # Valid knot found
    }

  }

  if (!knot_found) {
    return(list(newknot = NA_real_, cluster_index = NA_integer_))
  }

  return(list(newknot = newknot, cluster_index = index))
}

Knotnew_R <- function(wht, restmp, x, dcm, oldknots, tol) {

  # Number of residual clusters
  u <- length(dcm)
  # Number of old knots and internal knots
  noldint <- length(oldknots) - 6
  # Dimension of the x vector
  xdim <- length(x)
  # Initialize the index of the j_s-th cluster
  index <- NA
  # Initialize the new knot position
  newknot <- NA

  for (kk in seq_len(u)) {

    # Find the index of the cluster with the maximum weight
    index <- which.max(wht)

    # Get the boundaries of the current cluster
    if (index == 1) {
      d_lower <- 1
      } else {
        d_lower <- dcm[index - 1] + 1
        }
    d_upper <- dcm[index]
    inf <- x[d_lower]
    sup <- x[d_upper]

    # Initialize flag to check for internal knots in the cluster
    clusterFlag <- FALSE

    # 1) If there are internal knots, check whether they lie within the current cluster
    if (noldint > 0) {
      oldintknots <- get_internal_knots(oldknots, 3)  # Extract internal knots (ignoring the boundary knots)

      if (inf != sup) {
        # For each internal knot, check if it lies within the current cluster's boundaries
        for (jj in seq_len(noldint)) {
          clusterFlag <- clusterFlag || (inf <= oldintknots[jj] && oldintknots[jj] <= sup )
          }
        } else {
          # If inf == sup (single point), check if the singleton is close to previous knots
          for (jj in seq_len(noldint)) {
            clusterFlag <- clusterFlag || abs(inf - oldintknots[jj]) < tol
          }
        }
    }

    # 2) If an internal knot already exists in this cluster, set its weight to zero
    if (clusterFlag) {
      wht[index] <- 0
      } else {

        # 3) Compute the new knot as a weighted average of x-coordinates in the current cluster
        newknot <- sum(restmp[d_lower:d_upper] * x[d_lower:d_upper]) / sum(restmp[d_lower:d_upper])
        # Insert the new knot into the sorted list of old knots
        sortedknots <- c(oldknots, newknot)
        sortedknots <- sort(sortedknots)

        # For each consecutive set of 4 sortedknots, check whether there is at least one x that falls between the 1st and 4th knot in that set
        supportFlag <- TRUE  # initialize flag
        for (i in seq_len(length(sortedknots) - 3)) {
          temp <- FALSE
          for (jj in seq_len(xdim)) {
            temp <- (sortedknots[i] + tol < x[jj] && x[jj] < sortedknots[i + 3] - tol)
            if (temp) break
            }
          supportFlag <- temp
          if (!supportFlag) break
          }

        # If the new knot placement is invalid, set the weight of the cluster to zero
        if (!supportFlag) {
          wht[index] <- 0
          } else {
            # Break out of the loop if a valid knot is found
            break
          }
      }
  }

  # Return the new knot and the cluster index where it was inserted
  return(c(newknot, index))
}

#' @importFrom stats weighted.mean
findNewDimKnot_R <- function(dcumFixedDim.Dim, Dim.weights, Dim.oldknots,
                             Dim.values, residuals)
{

  # Loop through each cluster to find the optimal placement for a new Dim knot
  u <- length(dcumFixedDim.Dim) # total number of clusters
  flagDim <- FALSE              # flag to handle cases where all calculated weights are negative
  Dim.newknot <- weightDim <- NA_real_
  dcumInf <- dcumSup <- 0L

  if (length(Dim.weights) != u) {
    stop("Dim.weights and dcumFixedDim.Dim must have the same length.")
  }
  if (length(Dim.values) != length(residuals)) {
    stop("Dim.values and residuals must have the same length.")
  }

  for (i in seq_len(u)) {

    if (all(Dim.weights < 0)) {
      flagDim <- TRUE # No valid knot can be found when all remaining weights are negative
      break
    }

    # Find the index of the cluster with the highest weight
    indice <- which.max(Dim.weights)
    # Determine the (index) boundaries of the cluster with the highest weight
    if (indice == 1) {dcumInf = 1} else {dcumInf = dcumFixedDim.Dim[indice - 1] + 1}
    dcumSup <- dcumFixedDim.Dim[indice]
    # Calculate the superior and inferior Dim-bounds
    sup <- Dim.values[dcumSup]
    inf <- Dim.values[dcumInf]

    # (Step 7 - UnivariateFitter) Compute the new Dim knot as a weighted average of Dim values
    # within the selected cluster, weighted by their residuals
    # Dim.newknot <- residuals[dcumSup:dcumInf] %*%
    #   Dim.values[dcumSup:dcumInf] / sum(residuals[dcumSup:dcumInf])
    Dim.newknot <- weighted.mean(
      Dim.values[dcumSup:dcumInf],
      residuals[dcumSup:dcumInf]
      )

    # Check conditions to ensure the new knot is valid and does not conflict with existing knots
    # This involves ensuring there are no existing knots within the bounds of the selected cluster
    cond1 <- (dcumSup - dcumInf) != 0
    cond2 <- !any( Dim.oldknots >= inf & Dim.oldknots <= sup ) # no previous knot within the cluster
    cond3 <- !any( abs(inf - Dim.oldknots ) < as.double(1e-12) ) # no previous knot arbitrarily close to the singleton

    if ( cond1 && cond2 || !cond1 && cond3 ) {
      break # If conditions are met, exit the loop as a valid knot has been found
    } else {
      Dim.weights[indice] <- -Inf # Invalidate the current cluster by setting its weight to negative infinity and continue the search
    }

    # # Check conditions to ensure the new knot is valid and does not conflict with existing knots
    # # This involves
    # cond1 <- (dcumSup - dcumInf) != 0
    # cond2 <- !any((Dim.oldknots >= inf) & (Dim.oldknots <= sup)) # ensure there are no previous knots within the bounds of the selected cluster
    # cond3 <- dcumInf == 1 || dcumSup == length(FixedDim) # for the case in which the entire set is within one cluster
    # if ( (cond1 && cond2) || ( !cond1 && cond3 ) ) {
    #   break # If conditions are met, exit the loop as a valid knot has been found
    # } else {
    #   Dim.weights[indice] <- -Inf # Invalidate the current cluster by setting its weight to negative infinity and continue the search
    # }


  }

  # Check if all Dim.weights were turned to -Inf. In this case no admissible
  # cluster was found, matching the C++ return contract.
  if (all(Dim.weights < 0)) flagDim <- TRUE

  if (!flagDim) {
    weightDim <- Dim.weights[indice] # Store the selected cluster weight
  }


  return(list(Dim.newknot = Dim.newknot,
              weightDim = weightDim,
              flagDim = flagDim,
              dcumInf = dcumInf,
              dcumSup = dcumSup))

}


################################################################################
# Tensor product of two basis matrices, row by row.
tensorProd_R <- function(Xmat, Ymat) {
  # Get the number of rows and columns for both matrices
  ib <- nrow(Xmat)
  jx <- ncol(Xmat)
  jy <- ncol(Ymat)

  # Initialize the result matrix with the appropriate size
  ris <- matrix(0, nrow = ib, ncol = jx * jy)
  # Index kk to keep track of the column in the result matrix
  kk <- 1

  # Loop through columns of Xmat and Ymat
  for (jjx in 1:jx) {
    for (jjy in 1:jy) {
      for (ii in 1:ib) {
        # Fill in the result matrix with the product of elements from Xmat and Ymat
        ris[ii, kk] <- Xmat[ii, jjx] * Ymat[ii, jjy]
      }
      # Increment kk for the next column
      kk <- kk + 1
    }
  }

  return(ris)
}

################################################################################
# R reference for geds_knot_step() (Stage A steps 2 to 7 of UnivariateFitter).
# Inputs are prepared by UnivariateFitter() exactly as for the C++ call.
.geds_knot_step_R <- function(res.tmp, positive.weight, placement.weights,
                              Indicator, has_repeated_x, site_weight_scale,
                              geometry_tol, distinctX, beta, extr, intknots,
                              tol, support_order = 2) {
  weighted.contributions <- res.tmp[positive.weight] * placement.weights
  res.weighted <- if (has_repeated_x) {
    makeNewRes(resold = weighted.contributions, recurr = as.numeric(Indicator))
  } else weighted.contributions
  # Numerical zeros must not split geometric clusters after response centering.
  res.weighted[abs(res.weighted) <= geometry_tol * site_weight_scale] <- 0
  if (all(res.weighted == 0)) return(NA_real_)
  # Same-sign clusters, their mean residuals and their x-ranges
  d <- rle(sign(res.weighted))$lengths
  u <- length(d)
  dcum <- cumsum(d)
  cs <- cumsum(res.weighted)
  means <- abs((cs[dcum] - c(0, cs[dcum[-u]])) / d)
  left_idx <- c(1, dcum[-u] + 1)
  wc.range <- distinctX[dcum] - distinctX[left_idx]
  # Normalised cluster weights
  means <- means / max(means)
  if (max(wc.range) != 0) wc.range <- wc.range / max(wc.range)
  w <- beta * means + (1 - beta) * wc.range
  Knotnew(weights = w, residuals = res.weighted, x = distinctX, dcum = dcum,
          oldknots = c(rep(extr, support_order), intknots), tol = tol,
          support_order = support_order)[1]
}

################################################################################
# Weighted least-squares linear spline with breakpoints t (the boundary knots
# first and last, internal knots in between). The C++ version solves the
# tridiagonal normal equations and returns NULL for repeated breakpoints or an
# ill-conditioned system; this reference uses a QR fit, so results agree to
# rounding rather than bit for bit.
.linear_spline_ls_R <- function(x, y, w, t) {
  if (anyDuplicated(t)) return(NULL)
  B <- splines::splineDesign(c(t[1L], t, t[length(t)]), x, ord = 2L, outer.ok = TRUE)
  fit <- stats::lm.wfit(B, y, w)
  theta <- unname(fit$coefficients)
  list(theta = theta, fitted = as.vector(B %*% theta))
}

################################################################################
# Empirical semivariance of residuals over nearest neighbours, used by the GD
# noise estimator. nn is the p x n matrix of neighbour indices (1-based).
.gd_neighbour_semivariance_R <- function(r, nn) {
  colMeans((r - matrix(r[nn], nrow = length(r), byrow = TRUE))^2) / 2
}

