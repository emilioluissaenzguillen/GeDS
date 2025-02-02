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
          # If inf == sup (single point), check if the knot is close to the boundary
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
            temp <- (sortedknots[i] < x[jj] && x[jj] < sortedknots[i + 3])
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
