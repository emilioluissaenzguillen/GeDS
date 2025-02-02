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

recoverXmat <- function(ris, Ymat, Xmat) {
  # Get dimensions
  ib <- nrow(ris)
  jxjy <- ncol(ris)
  jy <- ncol(Ymat)
  jx <- jxjy / jy
  
  # Initialize the matrix to store the recovered Xmat
  Xmat <- matrix(0, nrow = ib, ncol = jx)
  
  # Recover Xmat by reversing the tensor product
  kk <- 1
  for (jjx in 1:jx) {
    for (jjy in 1:jy) {
      for (ii in 1:ib) {
        # Divide ris by Ymat to recover Xmat
        if (Ymat[ii, jjy] == 0) Ymat[ii, jjy] <- 1e-09
        Xmat[ii, jjx] <- ris[ii, kk] / Ymat[ii, jjy]
      }
      kk <- kk + 1
    }
  }
  
  return(Xmat)
}
