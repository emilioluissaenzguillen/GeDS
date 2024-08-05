################################################################################
################################################################################
################################ Plot_GeDSMethod ###############################
################################################################################
################################################################################
#' @title Plot method for GeDS objects.
#' @name plot,GeDS-method
#' @description
#' Plot method for GeDS objects. Plots GeDS fits.
#' @param x a \code{\link{GeDS-class}} object from which the GeDS fit(s) should
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
#' @importFrom graphics abline
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
        int.knt <- if ( i> 1) makenewknots(ik, n) else NULL
        # Stage B.2
        temp <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                             InterKnots = int.knt, n = n)
        
        # Update results with predicted values
        results$X <- X
        results$predicted <- temp$Predicted
        
        # Plot the spline
        lines(X, temp$Predicted, col = col_lines, lwd = 2)
        # add vertical lines for knots
        if (length(int.knt) < 20) {
          for(ik in c(int.knt, extr)) {
            abline(v = ik, col = "gray", lty = 2)
          }
        } else {
          rug(c(int.knt, extr))
        }
        
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
          lines(X,temp$NCI$Upp, col = col_lines, lwd = 2, lty = 2)
          lines(X,temp$NCI$Low, col = col_lines, lwd = 2, lty = 2)
          
          results$CIupp <- temp$NCI$Upp
          results$CIlow <- temp$NCI$Low
          
          if(draw.legend) legend(legend.pos,c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, col_lines, f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 3) Asymptotic Confidence Intervals
        } else if (type == "ACI") {
          lines(X,temp$ACI$Upp, col = col_lines, lwd = 2, lty = 2)
          lines(X,temp$ACI$Low, col = col_lines, lwd = 2, lty = 2)
          
          results$CIupp <- temp$ACI$Upp
          results$CIlow <- temp$ACI$Low
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, col_lines, f_col), pch = c(pch_data, NA, NA, f_pch),
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
        int.knt <- makenewknots(ik, n)
        # Stage B.2
        temp <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                              InterKnots = int.knt, n = n, family = family, mustart = x$Linear$Predicted)
        
        # Update results with predicted values
        results$X <- X
        results$pred <- temp$Predicted
        
        # Plot the spline
        lines(X, temp$Predicted, col = col_lines, lwd = 2)
        # add vertical lines for knots
        if (length(int.knt) < 20) {
          for(ik in c(int.knt, extr)) {
            abline(v = ik, col = "gray", lty = 2)
          }
        } else {
          rug(c(int.knt, extr))
        }
        
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
          
          matrice <- splineDesign(knots = sort(c(int.knt,rep(extr,2))), derivs = rep(0,length(X)),
                                  x = X, ord = n, outer.ok = TRUE)
          
          temp_nci <- glm(yy ~ -1+xx, offset = offset, weights = weights, family = family, start = temp$Thetas)
          pred <- predict.glm(temp_nci, newdata = data.frame(xx = matrice, offset = 0), se.fit = T, type = "response")
          CIupp <- pred$fit + qnorm(.975)*pred$se.fit
          CIlow <- pred$fit - qnorm(.975)*pred$se.fit
          
          lines(X, CIupp, col = "black", lwd = 2, lty = 2)
          lines(X, CIlow, col = "black", lwd = 2, lty = 2)
          
          results$CIupp <- CIupp
          results$CIlow <- CIlow
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col_lines, "black", f_col), pch = c(pch_data, NA, NA, f_pch),
                                 lwd = c(NA, 1, 1, f_lwd), bty = "n")
          
          # 3) Asymptotic Confidence Intervals
        } else if (type == "ACI") {
          CI <- confint.GeDS(object = x, n = n)
          CIupp <- CI[,2]
          CIlow <- CI[,1]
          results$CIupp <- CIupp
          results$CIlow <- CIlow
          
          lines(X, CIupp, col = col_lines, lwd = 2, lty = 2)
          lines(X, CIlow, col = col_lines, lwd = 2, lty = 2)
          
          if(draw.legend) legend(legend.pos, c("Data", toprint, "CI", f_legend), lty = c(NA, 1, 2, f_lty),
                                 col = c(col_data, col, "black", f_col), pch = c(pch_data, NA, NA, f_pch),
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
      seq_valX <- seq(Xextr[1], Xextr[2], length.out = round(sqrt(length(X))))
      seq_valY <- seq(Xextr[1], Xextr[2], length.out = round(sqrt(length(Y))))
      grid_val <- expand.grid(x = seq_valX, y = seq_valY)
      # Compute z for each x, y in the grid
      f_XY_val <- matrix(f(grid_val), nrow = round(sqrt(length(X))))
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
      newX <- seq(from = Xextr[1], to = Xextr[2], length.out = round(sqrt(length(X))))
      newY <- seq(from = Yextr[1], to = Yextr[2], length.out = round(sqrt(length(Y))))
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
      f_hat_XY_val <- matrix(f_hat_XY_val, nrow = round(sqrt(length(X))))
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


################################################################################
################################################################################
############################# Plot_GeDSboostMethod #############################
################################################################################
################################################################################
#' @title Plot method for GeDSboost objects.
#' @name plot,GeDSboost-method
#' @description
#' Plot method for GeDSboost object with a single base-learner (i.e. non-additive model). 
#' @param x a \code{\link{GeDSboost-class}} object from which the GeDSboost fit(s) should
#' be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the FGB-GeDS fit.
#' @param ... further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' 
#' @export 
#' 
#' @aliases plot.GeDSboost

setMethod("plot", signature(x = "GeDSboost"), function(x, n = 3L,...)
{
  
  # Check if x is of class "GeDSboost"
  if(!inherits(x, "GeDSboost")) {
    stop("The input 'x' must be of class 'GeDSboost'")
  }
  # Check if x has only one predictor
  if(length(x$args$predictors) > 1) {
    stop("Only available for models with a single predictor")
  }
  
  Y <- x$args$response[[1]]; X_mat <- x$args$predictors[[1]] 
  int.knots <- x$internal_knots$linear.int.knots[[1]]
  
  if (n == 2) {
    fit <- x$predictions$pred_linear
    legend <- c("Data", "Linear")
  } else if (n == 3) {
    fit <- x$predictions$pred_quadratic
    legend <- c("Data", "Quadratic")
  } else if (n == 4) {
    fit <- x$predictions$pred_cubic
    legend <- c("Data", "Cubic")
  }
  
  y_range <- range(Y, x$predictions$pred_linear, x$predictions$pred_quadratic, x$predictions$pred_cubic)
  y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
  
  # Capture additional parameters
  additional_params <- list(...)
  # Set the main title
  main0 <- if("main" %in% names(additional_params)) {
    additional_params$main
  } else {
    paste0(length(int.knots), " internal knots")
  }
  
  plot(X_mat, Y, main = main0, ylim = y_range, ...)
  
  pred <- data.frame(X_mat, fit)
  pred <- pred[order(pred$X_mat),]
  lines(pred, col = "red", lwd = 2, lty = 1)
  
  legend("topright",                           
         legend = legend,    
         col = c("black", "red"),
         lty = c(NA, 1),
         lwd = c(NA, 2),                                   
         pch = c(1, NA),
         bty = "n")
  if (length(int.knots) < 20) {
    for(int.knot in c(int.knots, rep(range(X_mat),n))) {
      abline(v = int.knot, col = "gray", lty = 2)
    }
  } else {
    rug(c(int.knots, rep(range(X_mat),n)))
  }
  
}
)


################################################################################
################################################################################
############################## Plot_GeDSgamMethod ##############################
################################################################################
################################################################################
#' @title Plot method for GeDSgam objects.
#' @name plot,GeDSgam-method
#' @description
#' Plots the component functions of a GeDSgam object fitted using \code{\link{NGeDSgam}}.

#' @param x a \code{\link{GeDSgam-class}} object from which the GeDSgam fit(s) should
#' be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GAM-GeDS fit.
#' @param ... further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' 
#' @export 
#' @importFrom plot3D persp3D segments3D
#' 
#' @aliases plot.GeDSgam

setMethod("plot", signature(x = "GeDSgam"), function(x, n = 3L, ...)
{
  
  # Check if x is of class "GeDSgam"
  if(!inherits(x, "GeDSgam")) {
    stop("The input 'x' must be of class 'GeDSgam'")
  }
  
  base_learners <- x$args$base_learners
  
  Y <- x$args$response[[1]]; pred_vars <- x$args$predictors
  
  if (n == 2) {
    Theta <- x$final_model$Linear.Fit$Theta
    int.knots <- x$internal_knots$linear.int.knots
  } else if (n == 3) {
    Theta <- x$final_model$Quadratic.Fit$Theta
    int.knots <- x$internal_knots$quadratic.int.knots
  } else if (n == 4) {
    Theta <- x$final_model$Cubic.Fit$Theta
    int.knots <- x$internal_knots$cubic.int.knots
  }
  
  for (bl_name in names(base_learners)) {
    
    bl <- base_learners[[bl_name]]
    X_mat <- pred_vars[, intersect(bl$variables, colnames(pred_vars))]
    
    int.knt <- int.knots[[bl_name]]
    
    pattern <- paste0("^", gsub("([()])", "\\\\\\1", bl_name))
    theta <- Theta[grep(pattern, names(Theta))]
    # Replace NA values with 0
    theta[is.na(theta)] <- 0
    
    # 1. Univariate learners
    if (NCOL(X_mat) == 1) {
      
      if (bl$type == "GeDS") {
        # Including int.knt and giving more values enhances visualization
        X_mat <- seq(from = min(X_mat), to = max(X_mat), length.out = 1000)
        X_mat <- sort(c(X_mat, int.knt))
        # Create spline basis matrix using specified knots, evaluation points and order
        basisMatrix <- splineDesign(knots = sort(c(int.knt,rep(range(X_mat),n))),
                                    x = X_mat, ord = n, derivs = rep(0,length(X_mat)),
                                    outer.ok = T)
        # To recover backfitting predictions need de_mean
        if (n == 2) de_mean <- TRUE else de_mean <- FALSE
        # To recover backfitting predictions need de_mean
        Predicted <- if (de_mean) basisMatrix %*% theta - mean(basisMatrix %*% theta) else basisMatrix %*% theta 
        ylab <-  bl_name
      } else if (bl$type == "linear") {
        # Linear
        if (!is.factor(X_mat)) {
          Predicted <- theta * X_mat
          # Factor
        } else {
          names(theta) <- levels(X_mat)[-1]
          theta[levels(X_mat)[1]] <- 0 # set baseline coef to 0
        }
        ylab <- bquote(beta[1] %*% .(bl_name))
      }
      
      if  (!is.factor(X_mat)) {
        
        plot_data <- data.frame(X_mat, Predicted)
        plot_data <- plot_data[order(plot_data$X_mat),]
        
        plot(plot_data, type = "l", col = "steelblue",
             xlab = bl$variables, ylab = ylab,
             xlim = range(X_mat), ylim = range(Predicted), ...)
        
        if (length(int.knt) < 20) {
          knt <- c(int.knt, range(X_mat))
          for(k in knt) {
            abline(v = k, col = "gray", lty = 2)
          }
        } else {
          rug(int.knt)
        }
        
      } else {
        par(mar = c(7.1, 4.1, 4.1, 2.1))
        barplot(theta,
                las = 2,
                col = "steelblue",
                main = bl_name,
                ylab = bquote(beta))
        par(mar = c(5.1, 4.1, 4.1,2.1))
      }
      
      # 2. Bivariate learners
    } else if (NCOL(X_mat) == 2) {
      
      Xextr <- range(X_mat[,1])
      Yextr <- range(X_mat[,2])
      
      # Create a sequence of sqrt(obs) evenly spaced numbers over the range of X/Y
      newX <- seq(from = Xextr[1], to = Xextr[2], length.out = round(sqrt(length(X_mat[,1]))))
      newY <- seq(from = Yextr[1], to = Yextr[2], length.out = round(sqrt(length(X_mat[,2]))))
      
      # Create a grid data frame from all combinations of newX and newY
      grid.data <- expand.grid(newX, newY)
      
      # Generate spline basis matrix for X and Y dimensions using object knots and given order
      matriceX <- splineDesign(knots = sort(c(int.knt$ikX,rep(Xextr,n))), derivs = rep(0,length(grid.data[,1])),
                               x = grid.data[,1], ord = n, outer.ok = T)
      matriceY <- splineDesign(knots = sort(c(int.knt$ikY,rep(Yextr,n))), derivs = rep(0,length(grid.data[,2])),
                               x = grid.data[,2], ord = n, outer.ok = T)
      # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
      matricebiv <- tensorProd(matriceX, matriceY)
      # Multiply the bivariate spline basis by model coefficients to get fitted values
      f_hat_XY_val <- matricebiv %*% theta[1:dim(matricebiv)[2]]
      # Reshape the fitted values to a square matrix for plotting
      f_hat_XY_val <- matrix(f_hat_XY_val, nrow = round(sqrt(length(X_mat[,1]))))
      f_hat_XY_val <- x$args$family$linkinv(f_hat_XY_val)
      # Title based on the number of internal knots in X and Y
      title <- paste0(length(int.knt$ikX)," ", bl$variables[1]," internal knots and ",
                      length(int.knt$ikY), " ", bl$variables[2]," internal knots" )
      
      # Plot the perspective 3D surface defined by newX, newY, and the fitted values
      persp3D(x = newX, y = newY, z = f_hat_XY_val, main = title, phi = 25, theta = 50,
              xlab = paste0("\n", bl$variables[1]), ylab = paste0("\n",bl$variables[2]),
              zlab = paste0("\n", bl_name), zlim = range(f_hat_XY_val),
              ticktype = "detailed", expand = 0.5, colkey = FALSE, border = "black", ...)
      
      # Add rug plots to the x axis
      kntX <- c(int.knt$ikX, Xextr)
      segments3D(x0 = kntX, y0 = rep(min(newY), length(kntX)), z0 = rep(min(f_hat_XY_val), length(kntX)),
                 x1 = kntX, y1 = rep(min(newY), length(kntX)), z1 = rep(max(f_hat_XY_val), length(kntX)),
                 add = TRUE, col = "black", lty = 2)
      
      # Add rug plots to the y axis
      kntY <- c(int.knt$ikY, Yextr)
      segments3D(x0 = rep(max(newX), length(kntY)), y0 = kntY, z0 = rep(min(f_hat_XY_val), length(kntY)),
                 x1 = rep(max(newX), length(kntY)), y1 = kntY, z1 = rep(max(f_hat_XY_val), length(kntY)),
                 add = TRUE, col = "black", lty = 2)
    }
    
    
  }
  
}
)

