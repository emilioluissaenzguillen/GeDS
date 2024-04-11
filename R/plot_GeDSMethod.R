################################################################################
################################################################################
################################ Plot_GeDSMethod ###############################
################################################################################
################################################################################
#' @title Plot method for GeDS objects.
#' @name plot,GeDS-method
#' @description
#' Plot method for GeDS objects. Plots GeDS fits.
#' @param x a \code{\link{GeDS-Class}} object from which the GeDS fit(s) should
#' be extracted.
#' @param f (optional) specifies the underlying function or generating process
#' to which the model was fit. This parameter is useful if the user wishes to
#' plot the specified function/process alongside the model fit and the data 
#' @param which a numeric vector specifying the iterations of stage A for which
#' the corresponding GeDS fits should be plotted.
#' It has to be a subset of  \code{1:nrow(x$stored)}. See details.
#' @param DEV logical variable specifying whether a plot representing the
#' deviance at each iteration of stage A should be produced or not.
#' @param ask logical variable specifying whether the user should be prompted
#' before changing the plot page.
#' @param main optional character string to be used as a title of the plot.
#' @param legend.pos the position of the legend within the panel. See
#' \link[graphics]{legend} for details.
#' @param new.window logical variable specifying whether the plot should be
#' shown in a new window or in the active one.
#' @param wait time, in seconds, the system should wait before plotting a new
#' page. Ignored if \code{ask = TRUE}.
#' @param ... further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS fit that should be plotted. By default equal to
#' \code{3L}. Non-integer values will be passed to the function
#' \code{\link{as.integer}}.
#' @param type character string specifying the type of plot required. Should be
#' set either to "\code{Polygon}" if the user wants to get also the control
#' polygon of the GeDS fit,  \code{"NCI"} or  \code{"ACI"} if 95\% confidence
#' bands for the predictions should be plotted (see details) or \code{"none"} if
#' only the fitted GeDS curve should be plotted. Applies only when plotting a
#' univariate spline regression.
#'
#' @details
#' This method is provided in order to allow the user to plot the GeDS  fits
#' contained in the \code{\link{GeDS-Class}} objects.
#' 
#' Since in Stage A of the GeDS algorithm the knots of a linear spline fit are
#' sequentially located, one at a time, the user may wish to visually inspect
#' this process using the argument \code{which}. The latter specifies a
#' particular iteration number (or a vector of such numbers) for which the
#' corresponding linear fit(s) should be plotted. The \code{ask} and \code{wait}
#' arguments can help the user to manage these pages.
#'
#' By means of \code{ask} the user can determine for how long each page should
#' appear on the screen. Pages are sequentially replaced by pressing the enter
#' button.
#'
#' Note that, in order to ensure stability, if the object was produced by the
#' function \code{\link{GGeDS}}, plotting intermediate fits of stage A is
#' allowed  only if \code{n = 2}, in contrast to objects produced by 
#' \code{\link{NGeDS}} for which plotting intermediate results is allowed also
#' for \code{n = }2 or 3 results.
#'
#' The confidence intervals obtained by setting \code{type = "NCI"} are
#' approximate local bands obtained considering the knots as fixed constants.
#' Hence the columns of the design matrix are seen as covariates and standard
#' methodology relying on the \code{se.fit} option of \code{predict.lm} or
#' \code{predict.glm} is applied.
#'
#' Setting \code{type = "ACI"}, asymptotic confidence intervals are plotted.
#' This option is applicable only if the canonical link function has been used
#' in the fitting procedure.
#'
#' @examples
#' ###################################################
#' # Generate a data sample for the response variable
#' # Y and the single covariate X, assuming Normal noise
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression using NGeDS
#' (Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Plot the final quadratic GeDS fit (red solid line)
#' # with its control polygon (blue dashed line)
#' plot(Gmod)
#'
#' # Plot the quadratic fit obtained from the linear fit at the 10th
#' # iteration of stage A i.e. after 9 internal knots have been inserted
#' # by the GeDS procedure
#' plot(Gmod, which=10)
#'
#' # Generate plots of all the intermediate fits obtained
#' # by running the GeDS procedure
#' \dontrun{
#' plot(Gmod, which=1:16)
#' }
#'
#' ###################################################
#' # Generate a data sample for the response variable Y and the covariate
#' # X assuming Poisson distributed error and a log link function
#'
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- exp(f_1(X))
#' # Generate Poisson distributed Y according to the mean model
#' Y <- rpois(N,means)
#'
#' # Fit a Poisson GeDS regression model using GGeDS
#' (Gmod2 <- GGeDS(Y ~ f(X), beta = 0.2, phi = 0.995, family = poisson(),
#'                 Xextr = c(-2,2)))
#'
#' # similar plots as before, but for the linear fit
#' plot(Gmod2, n = 2)
#' plot(Gmod2, which = 10, n = 2)
#' \dontrun{
#' plot(Gmod2, which = 1:16, n = 2)
#' plot(Gmod2, which = 1:16, n = 2, ask = T)
#' }
#'
#' @seealso \code{\link{NGeDS}} and \code{\link{GGeDS}};
#' \code{\link[graphics]{plot}}.
#' 
#' @export 
#' @importFrom plot3D persp3D points3D
#' 
#' @aliases plot.GeDS plot,GeDS-method plot,GeDS,ANY-method

setMethod("plot", signature(x = "GeDS"), function(x, f = NULL, which, DEV = FALSE, ask = FALSE,
                                                  main, legend.pos = "topright",
                                                  new.window = FALSE, wait = 0.5,
                                                  n=3L, type = c("none", "Polygon", "NCI", "ACI"), ...)
{
  results <- list()
  results$terms <- x$terms
  
  # Position of the legend within the panel
  draw.legend <- !(is.na(legend.pos))
  
  # Check length of DEV/ask/new.window/wait
  if(length(DEV)!= 1 || length(ask)!= 1 || length(new.window)!= 1 || length(wait)!= 1 )
    stop("Please, check the length of the parameters")
  # Logical variable specifying whether a plot representing the deviance at each iteration of stage A should be produced or not
  DEV <- as.logical(DEV)
  # Logical variable specifying whether the user should be prompted before changing the plot page
  ask <- as.logical(ask)
  # Logical variable specifying whether the plot should be shown in a new window or in the active one
  new.window <- as.logical(new.window)
  # Time, in seconds, the system should wait before plotting a new page
  wait <- as.numeric(wait)
  
  if (ask) {
    # Prompt before new page
    oask <- devAskNewPage(TRUE)
    # Restore the original setting
    on.exit(devAskNewPage(oask))
    # Since the plot waits for user input, set the sleep/wait time to 0
    slp <- 0
  } else {
    slp <- wait
  }
  
  # Extract arguments
  X <- x$Args$X; Y <- x$Args$Y; Z <- x$Args$Z
  weights <- x$Args$weights; q <- x$Args$q; offset <- x$Args$offset 
  
  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # Set toprint
  if (n == 2L) {
    toprint = "Linear"
  } else if (n==3L){
    toprint = "Quadratic"
  } else if (n==4L) {
    toprint = "Cubic"
  }
  
  maxim <- nrow(x$Stored) # number of iterations in stage A
  others <- list(...) # other arguments passed to the function
  
  # Data arguments
  col_data <- if (is.null(others$col)) "black" else others$col
  pch_data <- if (is.null(others$pch)) 1 else others$pch
  
  ########################
  ## 1. Univariate GeDS ##
  ########################
  if (x$Type == "LM - Univ" || x$Type == "GLM - Univ") {
    
    # Set plot color, default to "red" if not specified in additional arguments
    col_lines <- if ("col_lines" %in% names(others)) others$col_lines else "red"
    others$col_lines = NULL # remove 'col_lines' from additional arguments to prevent conflicts
    
    # Set default iteration(s) to plot to k + 1 if not specified
    if (missing(which)) which <- x$Nintknots + 1
    # If which == "all", plot all stage A iterations
    if (length(which) == 1 && which == "all") which <- 1:(maxim)
    # Independent variable extremes
    extr <- x$Args$extr
    
    # Validate ("Polygon", "NCI" or "ACI") and set plot type
    type <- match.arg(type)
    # Format independent variable name for plotting
    xname <- attr(x$terms,"specials")$f-1
    xname <- attr(x$terms,"term.labels")[xname]
    xname <- substr(xname,3,(nchar(xname)-1))
    # Extract dependent variable name for plotting
    yname <- rownames(attr(x$terms,"factors"))[1]
    # Set y-axis limits
    yylim <- range(c(Y, x$Linear$Predicted)) + 0.05 * c(-1, 1) *  range(c(Y, x$Linear$Predicted))
    
    # Determine the maximum number of iterations
    maxim <- nrow(x$Stored)
    
    # Validate 'which' parameter
    if ((!is.numeric(which) && !is.null(which)) || any(which < 1) || any(which > maxim)) {
      stop(sprintf("'which' must be between 1 and %d", maxim), domain = NA)
    }
    if (!is.null(which)) {
      which <- as.integer(which)
      last <- max(which)
    } else {
      last <-  0
    }
    # Warn if plotting iterations discarded by the algorithm
    if (any(which > maxim - q)) {
      warning("Plotting also iterations discarded by the algorithm")
    }
    
    # If generating function is also to be plotted
    if (!is.null(f)) {
      if (!is.function(f)) {
        stop("f must be a function")
      }
      # Check the number of arguments that f takes
      num_args <- length(formals(f))
      # If the function takes more than one argument, display a message
      if (num_args > 1) {
        stop("The function f should take only one argument. Instead of f(X, Y), you should provide f(X), with X having two columns.")
      }
      
      f_legend <- "f"; f_lty <- 1; f_col <- "black"; f_pch <- NA; f_lwd <- 2
      
    } else {
      f_legend <- f_lty <- f_col <- f_pch <- f_lwd <- NA
    }
    
    ###########################
    ## 1.1. Univariate NGeDS ##
    ###########################
    if (x$Type == "LM - Univ") {
      
      # Loop over specified iterations for plotting
      for (i in which) {
        
        # Logical variable specifying whether the plot should be shown in a new window or in the active one
        if(new.window) dev.new()
        if(slp > 0) {
          Sys.sleep(slp)
        }
        
        # Set main plot title
        main0 <- if(missing(main)) paste0( i-1, " internal knots") else main
        plot(X, Y, main = main0, xlab = xname, ylab = yname, ylim = yylim, ...)
        if (!is.null(f)) lines(X, f(X), col = "black", lwd = 2)
        
        # Obtain stage A knots and perform spline regression
        ik <- na.omit(x$Stored[i,-c(1,2,(i+2),(i+3))])
        # Stage B.1 (averaging knot location)
        knt <- if ( i> 1) makenewknots(ik, n) else NULL
        # Stage B.2
        temp <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                             InterKnots = knt, n = n)
        
        # Update results with predicted values
        results$X <- X
        results$predicted <- temp$Predicted
        
        # Plot the spline
        lines(X, temp$Predicted, col = col_lines)
        rug(c(knt, rep(extr,n))) # add a rug plot for knots
        
        ## Each branch now adds specific elements to the plot based on the selected type
        # 1) Polygon
        if (type == "Polygon") {
          results$Polykn <- temp$Poly$Kn
          results$Polyth <- temp$Poly$Thetas
          
          lines(temp$Poly$Kn,temp$Poly$Thetas, col = "blue",lty=2)
          points(temp$Poly$Kn,temp$Poly$Thetas, col = "blue")
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "Polygon", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, "blue", f_col), pch = c(pch_data, NA, 1, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 2) Normal Confidence Intervals
        } else if(type == "NCI") {
          lines(X,temp$NCI$Upp, col = "darkgrey", lty = 2)
          lines(X,temp$NCI$Low, col = "darkgrey", lty = 2)
          
          results$CIupp <- temp$NCI$Upp
          results$CIlow <- temp$NCI$Low
          
          if(draw.legend) legend(legend.pos,c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines,"darkgrey", f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 3) Asymptotic Confidence Intervals
        } else if (type == "ACI") {
          lines(X,temp$ACI$Upp, col = "darkgrey", lty = 2)
          lines(X,temp$ACI$Low, col = "darkgrey", lty = 2)
          
          results$CIupp <- temp$ACI$Upp
          results$CIlow <- temp$ACI$Low
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, "darkgrey", f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
        } else {
          if (draw.legend) legend(legend.pos, c("Data", toprint, f_legend), lty = c(NA, 1, f_lty),
                                  col = c(col_data, col_lines, f_col), pch = c(pch_data, NA, f_pch),
                                  lwd = c(NA, 1, f_lwd), bty = "n")
        }
        
      }
      
      ###########################
      ## 1.2. Univariate GGeDS ##
      ###########################
    } else if (x$Type == "GLM - Univ") {
      
      # Extract GLM family from model arguments
      family <- x$Args$family
      
      # Restrict plotting to linear spline in case which does not correspond to the final model's iteration
      if (n != 2L && which != maxim - q) {
        which <- maxim - q
        warning("Stage A iterations can be plotted only for the linear spline")
      }
      
      # Loop over specified iterations for plotting
      for(i in which) {
        
        # Logical variable specifying whether the plot should be shown in a new window or in the active one
        if(new.window) dev.new()
        if(slp > 0) {
          Sys.sleep(slp)
        }
        
        # Set main plot title
        main0 <- if(missing(main)) paste0(i-1, " internal knots") else main
        plot(X, Y, main = main0, xlab = xname, ylab = yname, ylim = yylim, ...)
        if (!is.null(f)) lines(X, f(X), col = "black", lwd = 2)
        
        # Obtain stage A knots and perform spline regression
        ik <- na.omit(x$Stored[i,-c(1,2,(i+2),(i+3))])
        # Stage B.1 (averaging knot location)
        knt <- makenewknots(ik, n)
        # Stage B.2
        temp <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                              InterKnots = knt, n = n, family = family, mustart = x$Linear$Predicted)
        
        # Update results with predicted values
        results$X <- X
        results$pred <- temp$Predicted
        
        # Plot the spline
        lines(X, temp$Predicted, col = col_lines)
        rug(c(knt, rep(extr,n))) # add a rug plot for knots
        
        ## Each branch now adds specific elements to the plot based on the selected type
        # 1) Polygon
        if (type == "Polygon") {
          results$Polykn <- temp$Poly$Kn
          results$Polyth <- family$linkinv(temp$Poly$Thetas)
          lines(results$Polykn, results$Polyth, col = "blue", lty = 2)
          points(results$Polykn, results$Polyth, col = "blue")
          if(draw.legend) legend(legend.pos, c("Data", toprint, "Polygon", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, "blue", f_col), pch = c(pch_data, NA, 1, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 2) Normal Confidence Intervals
        } else if(type=="NCI") {
          yy <- Y
          xx <- if(n!=2L) {
            temp$Basis
          } else {
            x$Linear$Basis
          }
          
          matrice <- splineDesign(knots = sort(c(knt,rep(extr,2))), derivs = rep(0,length(X)),
                                  x = X, ord = n, outer.ok = TRUE)
          
          temp_nci <- glm(yy ~ -1+xx, offset = offset, weights = weights, family = family, start = temp$Thetas)
          pred <- predict.glm(temp_nci, newdata = data.frame(xx = matrice, offset = 0), se.fit = T, type = "response")
          CIupp <- pred$fit + qnorm(.975)*pred$se.fit
          CIlow <- pred$fit - qnorm(.975)*pred$se.fit
          
          lines(X, CIupp, col = "darkgrey", lty = 2)
          lines(X, CIlow, col = "darkgrey", lty = 2)
          
          results$CIupp <- CIupp
          results$CIlow <- CIlow
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, "darkgrey", f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 3) Asymptotic Confidence Intervals
        } else if (type=="ACI") {
          CI <- confint.GeDS(object = x, n = n)
          CIupp <- CI[,2]
          CIlow <- CI[,1]
          results$CIupp <- CIupp
          results$CIlow <- CIlow
          
          lines(X, CIupp, col="darkgrey", lty = 2)
          lines(X, CIlow, col="darkgrey", lty = 2)
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col, "darkgrey", f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
        } else {
          if (draw.legend) legend(legend.pos, c("Data", toprint, f_legend), lty = c(NA, 1, f_lty),
                                  col = c(col_data, col_lines, f_col), pch = c(pch_data, NA, f_pch),
                                  lwd = c(NA, 1, f_lwd), bty = "n")
          
        }
      }
    }
    
    ## If DEV = TRUE, produce plot representing the deviance at each iteration of stage A
    if (x$Type == "LM - Univ") {
      print <- "RSS"
    } else if (x$Type == "GLM - Univ") {
      print <- "DEV"
    }
    
    if (DEV) {
      # i. RSS/DEV
      main0 <- if(missing(main)) print else main
      plot(seq(0, x$iters-1), x$RSS/length(X), main = main0, xlab = "Internal knots", ylab = "", ...)
      rug(seq(1, x$iters-1 - q))
      if (slp > 0) {
        Sys.sleep(slp)
      }
      # ii. \sqrt{RSS} / \sqrt{DEV}
      main0 <- if(missing(main)) bquote(sqrt(.(print))) else main
      plot(seq(0, x$iters-1), (x$RSS/length(X))^.5, main = main0, xlab = "Knots", ylab = "", ...)
      rug(seq(1, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
      # iii. \phi
      main0 <- if(missing(main)) expression(phi) else main
      plot(seq(q, x$iters-1), (x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)]), main = main0, xlab="Knots", ylab = "", ...)
      rug(seq(2, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
      # iv. \sqrt{\phi}
      main0 <- if(missing(main)) expression(sqrt(phi)) else main
      plot(seq(q, x$iters-1), (x$RSS[(1+q):x$iters]/x$RSS[(1):(x$iters-q)])^.5, main = main0, xlab="Knots", ylab = "", ...)
      rug(seq(2, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
    }
    
    ##############################
    ## 2. Bivariate NGeDS/GGeDS ##
    ##############################
  } else if (x$Type == "LM - Biv" || x$Type == "GLM - Biv") {
    
    if(n == 2L) {
      obj <- x$Linear
    } else if (n == 3L) {
      obj <- x$Quadratic
    } else if (n == 4L) {
      obj <- x$Cubic
    }
    
    # Extract parametric component of the predictor model + X, Y extremes
    W <- x$Args$W; Xextr <- x$Args$Yextr; Yextr <- x$Args$Yextr
    
    if (!is.null(f)) {
      if (!is.function(f)) {
        stop("f must be a function")
      }
      # Check the number of arguments that f takes
      num_args <- length(formals(f))
      # If the function takes more than one argument, display a message
      if (num_args > 1) {
        stop("The function f should take only one argument. Instead of f(X, Y), you should provide f(X), with X having two columns.")
      }
      
      # Create a grid of x, y values
      seq_valX <- seq(Xextr[1], Xextr[2], length.out = sqrt(length(X)))
      seq_valY <- seq(Xextr[1], Xextr[2], length.out = sqrt(length(Y)))
      grid_val <- expand.grid(x = seq_valX, y = seq_valY)
      # Compute z for each x, y in the grid
      f_XY_val <- matrix(f(grid_val), nrow = sqrt(length(X)))
      # Title based on the number of internal knots in X and Y
      title <- paste(x$Nintknots$X, "linear internal knots on X,",
                     x$Nintknots$Y, "linear internal knots on Y")
      
      # Plot
      persp3D(x = seq_valX, y = seq_valY, z = f_XY_val, main = title, phi = 25, theta = 50,
              xlab = 'X', ylab = 'Y', zlab = "f(X,Y)", zlim = range(f_XY_val),
              ticktype = "detailed", expand = 0.5, colkey = FALSE, border = "black", ...)
      # Add the data and predicted points to the plot
      points3D(x = X, y = Y, z = Z, col = "black", pch = 19, add = TRUE)
      points3D(x = X, y = Y, z = obj$Predicted, col = "red", pch = 19, add = TRUE)
      legend("topright",
             legend = c("Data", "Quadratic Fit"),
             col = c("black", "red"),
             pch = 19,
             bg = 'white',
             bty = "n")
      
    } else {
      
      # Create a sequence of sqrt(obs) evenly spaced numbers over the range of X/Y
      newX <- seq(from = Xextr[1], to = Xextr[2], length.out = sqrt(length(X)))
      newY <- seq(from = Yextr[1], to = Yextr[2], length.out = sqrt(length(Y)))
      # Create a grid data frame from all combinations of newX and newY
      grid.data <- expand.grid(newX, newY)
      
      # Generate spline basis matrix for X and Y dimensions using object knots and given order
      matriceX <- splineDesign(knots = obj$Xknots, derivs = rep(0,length(grid.data[,1])),
                               x = grid.data[,1], ord = n, outer.ok = T)
      matriceY <- splineDesign(knots = obj$Yknots, derivs = rep(0,length(grid.data[,2])),
                               x = grid.data[,2], ord = n, outer.ok = T)
      # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
      matricebiv <- tensorProd(matriceX, matriceY)
      # Multiply the bivariate spline basis by model coefficients to get fitted values
      f_hat_XY_val <- matricebiv %*% obj$Theta[1:dim(matricebiv)[2]]
      # Reshape the fitted values to a square matrix for plotting
      f_hat_XY_val <- matrix(f_hat_XY_val, nrow = sqrt(length(X)))
      if (x$Type == "GLM - Biv") f_hat_XY_val <- x$Args$family$linkinv(f_hat_XY_val)
      # Title based on the number of internal knots in X and Y
      title <- paste0(x$Nintknots$X, " X internal knots and ", x$Nintknots$Y ," Y internal knots" )
      
      # Plot the perspective 3D surface defined by newX, newY, and the fitted values
      persp3D(x = newX, y = newY, z = f_hat_XY_val, main = title, phi = 25, theta = 50,
              xlab = 'X', ylab = 'Y', zlab = "f_hat(X,Y)", zlim = range(f_hat_XY_val),
              ticktype = "detailed", expand = 0.5, colkey = FALSE, border = "black", ...)
      
      # Adjust Z by subtracting the parametric component of the predictor model if it exists
      if (!is.null(W)) {
        Z <- Z - W %*% x$Linear$Theta[-(1:dim(matricebiv)[2])]
      }
      # Identify points where Z is greater/smaller than the model predicted value
      tmp <- (obj$Predicted - Z > 0)
      points3D(x = X[!tmp], y = Y[!tmp], z = Z[!tmp], col = "red", pch = 19, add = TRUE)
      points3D(x = X[tmp], y = Y[tmp], z = Z[tmp], col = "blue", pch = 19, add = TRUE)
      legend("topright",
             legend = c("Z >= f_hat(X,Y)", "Z  <  f_hat(X,Y)"),
             col = c("red", "blue"),
             pch = 19,
             bg = 'white',
             bty = "n")
    }
    
  } else {
    stop("Type not recognized")
  }
  
  x$Plotinfo <- results
  invisible(x)
}
)




