################################################################################
################################################################################
################################### PPolyRep ###################################
################################################################################
################################################################################
#' @title Piecewise Polynomial Spline Representation
#' @name PPolyRep
#' @description
#' The function converts a GeDS fit which has a  B-spline representation to a
#' piecewise polynomial form.
#' @param object  the \code{\link{GeDS-class}} where the GeDS fit to be
#' converted is found.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS fit which should be converted to a piecewise
#' polynomial form. By default equal to \code{3L}. Non-integer values will be
#' passed to the function \code{\link{as.integer}}.
#' 
#' @return An object that inherits from classes  \code{"spline"} and
#' \code{"polySpline"}. It is a list whose arguments are:
#' \item{knots}{ a vector of size  \eqn{k + 2} containing the complete set of 
#' knots (internal knots plus the limits of the interval) of the GeDS fit.}
#' \item{coefficients}{ a \eqn{(k + 2) \times n} matrix containing the
#' coefficients of the  polynomials in the required piecewise polynomial
#' representation. }
#' 
#' @details
#' This function converts a selected GeDS fit from a \code{\link{GeDS-class}}
#' object represented in terms of B-splines into an object where the fit is
#' represented in terms of piecewise polynomials.
#'
#' The function  wraps \code{\link[splines]{polySpline}} in order to let it 
#' accept \code{\link{GeDS-class}} objects as input. Hence the function provides
#' a useful link between the package \pkg{GeDS} and the package \pkg{splines},
#' allowing the user to take advantage of the functions provided in the
#' \pkg{splines} package.
#'
#' @examples
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression using NGeDS
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2))
#'
#' # construct the PP representation of the cubic GeDS fit
#' # and apply some functions of the package splines
#' Polymod <- PPolyRep(Gmod, 4)
#' require(splines)
#' class(Polymod)
#' splineKnots(Polymod)
#' knots(Gmod, n = 4)
#' plot(Polymod)
#'
#'
#' # Generate a plot showing the PP representation
#' # based on the same example
#' knt <- splineKnots(Polymod)
#' coeffs <- coef(Polymod)
#' plot(Gmod, n = 4, legend.pos = FALSE, main = "Cubic Curves")
#' cols <- sample(heat.colors(length(knt)), length(knt))
#' for(i in 1:(length(knt))){
#'   curve(coeffs[i,1] + coeffs[i,2]*(x - knt[i])+
#'           coeffs[i,3]*(x - knt[i])^2+
#'         coeffs[i,4]*(x - knt[i])^3,
#'         add = TRUE, col = cols[i])
#'   abline(v = knt[i])
#' }
#'
#' @export
#' 
#' @note Let us note that the first \eqn{k+1} rows of the matrix contain the
#' \code{n} coefficients of the \eqn{k+1} consecutive pieces of the piecewise
#' polynomial representation. The last \eqn{(k+2)}-th row is extraneous and it
#' appears as a result of the use of the function
#' \code{\link[splines]{polySpline}}.

PPolyRep <- function(object, n = 3)
  {
  if (!inherits(object, "GeDS")) stop("This function works only with GeDS class objects")
  n <- as.integer(n)
  kn <- knots(object, n = n, options="all")
  cf <- coef(object, n = n)
  newlist <- list(knots = kn, coefficients = cf, order = n)
  class(newlist) <- c("nbSpline", "bSpline", "spline")
  xname <- attr(object$terms,"specials")$f-1
  xname <- attr(object$terms,"term.labels")[xname]
  xname <- substr(xname,3,(nchar(xname)-1))
  yname <- rownames(attr(object$terms,"factors"))[1]
  fortmp <- paste0(yname," ~ ", xname)
  attr(newlist,"formula") <- as.formula(fortmp)
  out <- polySpline(newlist)
  return(out)
}


#####################
## Invert PPolyRep ##
#####################
PPolyInv <- function(ppoly, y_new)
{
  # Check if ppoly is of one of the specified classes
  if (!inherits(ppoly, "npolySpline") && !inherits(ppoly, "polySpline") && !inherits(ppoly, "spline")) {
    stop("ppoly must be of class 'npolySpline', 'polySpline', or 'spline'")
  }
  # Check if ppoly is invertible
  is_invppoly <- is_invertible(ppoly)
  if(!is_invppoly$is_inv) stop("ppoly needs to be stricly monotonic to be invertible!")
  
  # Check if y_new is a numeric vector
  if (!is.numeric(y_new)) {
    stop("y_new must be a numeric vector")
  }
  
  # Ensure y_new is a column matrix
  y_new <- matrix(y_new, ncol = 1)
  # Track original order with id variable
  y_new <- cbind(y_new, id = 1:length(y_new))
  # Ensure y_new remains a matrix after ordering
  y_new <- y_new[order(y_new[, 1]), , drop = FALSE]
  
  
  # Knots & Coefficients #
  # 1) ppoly$knots: a vector of size k + 2 containing the complete set of
  # knots (internal knots plus the limits of the interval) of the GeDS fit.
  knots <- ppoly$knots; k <- length(knots) - 2
  # 2) Let us note that the first k (= number of internal knots) + 1 rows of the
  # matrix contain the n coefficients of the k + 1 consecutive pieces of the
  # piecewise polynomial representation.
  coefficients <- ppoly$coefficients[1:(k+1), ]
  ylim <- ppoly$coefficients[,1]
  
  aux <- data.frame(
    start_x   = knots[1:(k+1)],
    end_x     = knots[2:(k+2)],
    start_y   = ylim[1:(k+1)],
    end_y     = ylim[2:(k+2)]
  )
  # Dynamically add polynomial terms based on the number of columns in `coefficients`
  poly_names <- c("constant", "linear", "quadratic", "cubic", "quartic")
  for (i in 1:ncol(coefficients)) {
    term_name <- poly_names[i]  
    aux[[term_name]] <- coefficients[, i]
  }
  
  
  x_new <- sapply(y_new[,1], solve_x, aux = aux, slope = is_invppoly$slp, knots = knots, k = k)
  
  # x_new <- vector("numeric", length(y_new[,1]))
  # 
  # for (i in 1:length(y_new[,1])) {
  #   
  #   # First and last value
  #   if (abs(y_new[,1][i] - aux[1, "start_y"]) < 1e-5) {
  #     x_new[i] <- knots[1]
  #     next
  #   }
  #   if (abs(y_new[,1][i] - aux[nrow(aux), "end_y"]) < 1e-5) {
  #     x_new[i] <- knots[k+2]
  #     next
  #   }
  #   
  #   interval_condition <- FALSE # Flag to track if the y interval is found
  #   
  #   for (j in 1:NROW(aux)) {
  #     
  #     if (is_invppoly$slp == "increasing") {
  #       interval_condition <- round(aux$start_y[j],5) <= round(y_new[,1][i],5) && y_new[,1][i] < aux$end_y[j]
  #     } else if (is_invppoly$slp == "decreasing") {
  #       interval_condition <- round(aux$start_y[j],5) >= round(y_new[,1][i],5) && y_new[,1][i] > aux$end_y[j]
  #     }
  #     
  #     if (interval_condition) {
  #       
  #       # p(x) = ax^2+bx+c = y ==> ax^2+bx+(c-y) = 0
  #       coef <- aux[,setdiff(names(aux), c("start_x", "end_x", "start_y", "end_y"))][j,]
  #       coef$constant <- coef$constant - y_new[,1][i]
  #       
  #       ########################
  #       ######## Linear ########
  #       ########################
  #       if (length(coef) == 2) {
  #         
  #         a <- coef$linear
  #         b <- coef$constant
  #         x_new[i]   <- -b/a + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #         
  #       ########################
  #       ###### Quadratic #######
  #       ########################  
  #       } else if (length(coef) == 3) {
  #         
  #         a <- coef$quadratic
  #         b <- coef$linear
  #         c <- coef$constant
  #         discriminant <- if (is_invppoly$slp == "increasing") {
  #           sqrt(b^2 - 4*a*c)
  #           } else {
  #             -sqrt(b^2 - 4*a*c)
  #             } 
  #         x_new[i]   <- (-b + discriminant) / (2*a) + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #       
  #       ########################
  #       #### Cubic/Quartic #####
  #       ########################   
  #       } else {
  #         roots <- polyroot(as.numeric(coef))
  #         
  #         # Filter only real roots
  #         roots <- roots[abs(Im(roots)) < 1e-5]
  #         suppressWarnings({
  #           real_roots <- as.numeric(roots)
  #         })
  #         # If there is more than one root, choose the min non-neg
  #         if (length(real_roots) > 1) {
  #           if (all(real_roots < 0)) {
  #             real_roots <- real_roots[which.min(abs(real_roots))]
  #           } else {
  #             real_roots <- real_roots[real_roots >= 0]
  #             real_roots <- real_roots[which.min(real_roots)]
  #           }
  #         }
  #         
  #         x_new[i]   <- real_roots + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #         
  #       }
  #     }
  #   }
  # }
  
  # Recover original ordering based on "id" column
  y_new <- cbind(y_new, x_new)
  y_new <- y_new[order(y_new[, "id"]), , drop = FALSE]
  # Extract the "x_new" column as a column vector
  x_new <- y_new[, "x_new", drop = FALSE]
  
  return(x_new)
}

# Helper function to check if the PPoly is invertible; for continuous functions
# strict monotonicity is a necessary and sufficient condition.
is_invertible <- function(ps, grid_points = 1000) {
  # Define a grid over the domain
  xs <- seq(min(ps$knots), max(ps$knots), length.out = grid_points)
  # Predict the spline values on the grid
  ys <- predict(ps, xs)$y
  # Compute differences between consecutive values
  diffs <- diff(ys)
  # Check if the spline is strictly increasing or strictly decreasing
  is_increasing <- all(diffs > -1e-6)
  is_decreasing <- all(diffs < 1e-6)
  is_invertible <- is_increasing || is_decreasing
  
  # Return a list containing the invertibility and monotonicity status
  return(list(is_inv = is_invertible,
              slp = if (is_increasing) "increasing" else if (is_decreasing) "decreasing" else "not monotonic"))
}

# Helper to calculate polynomial roots
solve_x <- function(y_value, aux, slope, knots, k) {
  
  # First and last value
  if (abs(y_value - aux$start_y[1]) < 1e-5) return(knots[1])
  if (abs(y_value - aux$end_y[nrow(aux)]) < 1e-5) return(knots[k+2])
  
  # Find the corresponding
  j <- which(if (slope == "increasing") {
    round(aux$start_y, 5) <= round(y_value, 5) & y_value < aux$end_y
  } else {
    round(aux$start_y, 5) >= round(y_value, 5) & y_value > aux$end_y
  })
  
  if (length(j) == 0) return(NA)  # If no interval is found, return NA
  if (length(j) > 1) stop("Overlapping intervals in ppoly!")
  
  
  # p(x) = ax^2+bx+c = y ==> ax^2+bx+(c-y) = 0
  coef <- aux[,setdiff(names(aux), c("start_x", "end_x", "start_y", "end_y"))][j,]
  coef$constant <- coef$constant - y_value
  # Order of polynomial
  n <- length(coef)
  ########################
  ######## Linear ########
  ########################
  if (n == 2) {
    a <- coef$linear
    b <- coef$constant
    return(-b/a + aux$start_x[j])
    
    ########################
    ###### Quadratic #######
    ########################  
    } else if (length(coef) == 3) {
      a <- coef$quadratic
      b <- coef$linear
      c <- coef$constant
      disc <- b^2 - 4 * a * c
      if(disc < 0) stop("Negative discriminant encountered.")
      sqrt_disc <- if (slope == "increasing") sqrt(disc) else -sqrt(disc)
      return( (-b + sqrt_disc) / (2*a) + aux$start_x[j] )
      
      ########################
      #### Cubic/Quartic #####
      ########################   
      } else {
        roots <- polyroot(as.numeric(coef))
        
        # Filter only real roots
        real_roots <- Re(roots[abs(Im(roots)) < 1e-5])
        
        # If there is more than one root, choose the min non-neg
        if (length(real_roots) > 1) {
          if (all(real_roots < 0)) {
            real_roots <- real_roots[which.min(abs(real_roots))]
            } else {
              real_roots <- real_roots[real_roots >= 0]
              real_roots <- real_roots[which.min(real_roots)]
            }
          }
        
        return(real_roots + aux$start_x[j])
      }
  
}


