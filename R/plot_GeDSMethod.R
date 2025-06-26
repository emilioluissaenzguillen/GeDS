################################################################################
################################################################################
############################### plot.GeDS Method ###############################
################################################################################
################################################################################
#' @title Plot Method for GeDS Objects
#' @name plot.GeDS
#' @description
#' Plot method for GeDS objects. Plots GeDS fits.
#' @param x An object of class \code{"GeDS"}, as returned by \code{NGeDS()} or
#' \code{GGeDS()}.
#' @param f (Optional) specifies the underlying function or generating process
#' to which the model was fit. This parameter is useful if the user wishes to
#' plot the specified function/process alongside the model fit and the data 
#' @param which A numeric vector specifying the iterations of stage A for which
#' the corresponding GeDS fits should be plotted.
#' It has to be a subset of  \code{1:nrow(x$stored)}. See Details.
#' @param dev Logical variable specifying whether a plot representing the
#' deviance at each iteration of stage A should be produced or not.
#' @param ask Logical variable specifying whether the user should be prompted
#' before changing the plot page.
#' @param main An optional character string used as the plot title. If set to
#' \code{"detail"}, the knots vector will be displayed on the plot.
#' @param legend.pos The position of the legend within the panel. See
#' \link[graphics]{legend} for details.
#' @param legend.text A character vector specifying the legend text.
#' @param new.window Logical variable specifying whether the plot should be
#' shown in a new window or in the active one.
#' @param wait Time, in seconds, the system should wait before plotting a new
#' page. Ignored if \code{ask = TRUE}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GeDS fit that should be plotted. By default equal to
#' \code{3L}. Non-integer values will be passed to the function
#' \code{\link{as.integer}}.
#' @param type Character string specifying the type of plot required. Should be
#' set either to \code{"polygon"} if the user wants to get also the control
#' polygon of the GeDS fit,  \code{"nci"} or  \code{"aci"} if 95\% confidence
#' bands for the predictions should be plotted (see Details) or \code{"none"} if
#' only the fitted GeDS curve should be plotted. Applies only when plotting a
#' univariate spline regression.
#' @param ... Further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#'
#' @details
#' This method is provided in order to allow the user to plot the GeDS  fits
#' contained in the \code{"GeDS"} class objects.
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
#' The confidence intervals obtained by setting \code{type = "nci"} are
#' approximate local bands obtained considering the knots as fixed constants.
#' Hence the columns of the design matrix are seen as covariates and standard
#' methodology relying on the \code{se.fit} option of \code{predict.lm} or
#' \code{predict.glm} is applied.
#'
#' Setting \code{type = "aci"}, asymptotic confidence intervals are plotted.
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
#' plot(Gmod, main = "detail")
#' plot(Gmod, type = "nci")
#' plot(Gmod, type = "aci")
#'
#' # Plot the quadratic fit obtained from the linear fit at the 10th
#' # iteration of stage A i.e. after 9 internal knots have been inserted
#' # by the GeDS procedure
#' plot(Gmod, which=10)
#'
#' # Generate plots of all the intermediate fits obtained
#' # by running the GeDS procedure
#' \dontrun{
#' plot(Gmod, which=1:(Gmod$Nintknots + Gmod$args$q + 1))
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
#' iters <- (Gmod2$Nintknots + Gmod2$args$q + 1)
#' plot(Gmod2, which = 1:iters, n = 2)
#' plot(Gmod2, which = 1:iters, n = 2, ask = T)
#' }
#'
#' @seealso \code{\link{NGeDS}} and \code{\link{GGeDS}};
#' \code{\link[graphics]{plot}}.
#' 
#' @method plot GeDS
#' @importFrom plot3D persp3D points3D
#' @importFrom grDevices dev.new devAskNewPage
#' @importFrom graphics plot lines legend rug points abline mtext
#' @importFrom stats predict.glm
#' @export 

plot.GeDS <- function(x, f = NULL, which, dev = FALSE, ask = FALSE,
                      main, legend.pos = "topright", legend.text = NULL,
                      new.window = FALSE, wait = 0.5,
                      n = 3L, type = c("none", "polygon", "nci", "aci"), ...)
{
  results <- list()
  results$terms <- x$terms
  
  # Position of the legend within the panel
  draw.legend <- !(is.na(legend.pos))
  
  # Check length of dev/ask/new.window/wait
  if(length(dev)!= 1 || length(ask)!= 1 || length(new.window)!= 1 || length(wait)!= 1 )
    stop("Please, check the length of the parameters")
  # Logical variable specifying whether a plot representing the deviance at each iteration of stage A should be produced or not
  dev <- as.logical(dev)
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
  X <- x$args$X; Y <- x$args$Y; Z <- x$args$Z
  weights <- x$args$weights; q <- x$args$q; offset <- x$args$offset 
  
  # Check if order is correctly set
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  
  # Set toprint
  if (n == 2L) {
    toprint = "Linear"
  } else if (n == 3L) {
    toprint = "Quadratic"
  } else if (n == 4L) {
    toprint = "Cubic"
  }
  
  maxim <- nrow(x$stored) # number of iterations in stage A
  others <- list(...) # other arguments passed to the function
  
  # Data arguments
  col_data <- if (is.null(others$col)) "black" else others$col
  pch_data <- if (is.null(others$pch)) 1 else others$pch
  
  ########################
  ## 1. Univariate GeDS ##
  ########################
  if (x$type == "LM - Univ" || x$type == "GLM - Univ") {
    
    # Set plot color, default to "red" if not specified in additional arguments
    col_lines <- if ("col_lines" %in% names(others)) others$col_lines else "red"
    others$col_lines = NULL # remove 'col_lines' from additional arguments to prevent conflicts
    
    # Set default iteration(s) to plot to k + 1 if not specified
    if (missing(which)) which <- x$Nintknots + 1
    # If which == "all", plot all stage A iterations
    if (length(which) == 1 && which == "all") which <- 1:(maxim)
    # Independent variable extremes
    extr <- x$args$extr
    
    # Validate ("polygon", "nci" or "aci") and set plot type
    type <- match.arg(type)
    
    # Format independent variable name for plotting
    var_labels <- attr(x$terms, "variables")  # Full list of variables
    xname <- as.character(var_labels[[3]])[2] # Exclude response variable
    # Extract dependent variable name for plotting
    yname <- rownames(attr(x$terms,"factors"))[1]
    # Use provided labels or defaults
    xlab <- if (!is.null(others$xlab)) others$xlab else if ( !missing(main) && main=="detail" ) "" else xname
    ylab <- if (!is.null(others$ylab)) others$ylab else yname
    
    # Set y-axis limits
    yylim <- range(c(Y, x$linear.fit$predicted)) + 0.05 * c(-1, 1) *  range(c(Y, x$linear.fit$predicted))

    # Determine the maximum number of iterations
    maxim <- nrow(x$stored)
    
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
    
    ###############################
    ## 1. Univariate NGeDS/GGeDS ##
    ###############################
    if (x$type == "GLM - Univ") {
      # Extract GLM family from model arguments
      family <- x$args$family
      # Restrict plotting to linear spline in case which does not correspond to the final model's iteration
      if (n != 2L && which != maxim - q) {
        which <- maxim - q
        warning("Stage A iterations can be plotted only for the linear spline")
      }
    }
    
    # Loop over specified iterations for plotting
    for (i in which) {
      
      # Logical variable specifying whether the plot should be shown in a new window or in the active one
      if(new.window) dev.new()
      if(slp > 0) {
        Sys.sleep(slp)
      }
      # Set main plot title
      main0 <- if (missing(main)) {
        paste0(i-1, " internal knots")
      } else if (main == "detail") {
        ""
      } else {
        main
      }
      # Plot
      # if (!is.null(others$ylim)) {
      #   plot(X, Y, main = main0, ...)  
      # } else {
      #   plot(X, Y, main = main0, ylim = yylim, ...)
      # }
      plot_args <- list(x = X, y = Y, main = main0)
      if (!is.null(xlab)) plot_args$xlab <- xlab
      if (is.null(others$ylim)) plot_args$ylim <- yylim
      # Protect against weird y-axis titles
      if (is.null(others$ylab) && is.null(list(...)$ylab)) plot_args$ylab <- "Y"
      plot_args <- c(plot_args, list(...))
      do.call(plot, plot_args)
      
      if (!is.null(f)) lines(X, f(X), col = "black", lwd = 2)
      
      # Obtain stage A knots and perform spline regression
      ik <- na.omit(x$stored[i,-c(1,2,(i+2),(i+3))])
      # Stage B.1 (averaging knot location)
      int.knt <- if (i > 1) makenewknots(ik, n) else NULL
      
      ###########################
      ## 1.1. Univariate NGeDS ##
      ###########################
      if (x$type == "LM - Univ") {
        # Stage B.2
        temp <- SplineReg_LM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                             InterKnots = int.knt, n = n)
        
        ###########################
        ## 1.2. Univariate GGeDS ##
        ###########################
      } else if (x$type == "GLM - Univ") {
        # Stage B.2
        temp <- SplineReg_GLM(X = X, Y = Y, Z = Z, offset = offset, weights = weights, extr = extr,
                              InterKnots = int.knt, n = n, family = family, mustart = x$linear.fit$predicted)
      }
      
      # Update results with predicted values
      results$X <- X
      results$predicted <- temp$predicted
      
      # Plot the spline
      lines(X, temp$predicted, col = col_lines, lwd = 2)
      
      # Add vertical lines for knots
      if (length(int.knt) < 20) {
        for(ik in c(int.knt, extr)) {
          abline(v = ik, col = "gray", lty = 2)
        }
      } else {
        rug(c(int.knt, extr))
      }
      # Detailed knots vector
      if (!missing(main) && main == "detail") {
        # knot
        knot_values <- round(knots(x, n=2), 2)
        intknt_text <- if (!is.null(get_internal_knots(knot_values, 2))) {
          paste0(", ", get_internal_knots(knot_values, 2))
        } else {
          ""
        }
        knt_text <- paste0(
          "a = ",
          knot_values[1],
          paste(intknt_text, collapse = ""),
          ", b = ",
          knot_values[length(knot_values)]
        )
        if (nchar(knt_text) > 50) {
          
          knt_text_split <- strsplit(knt_text, ",\\s*")[[1]]
          mid <- ceiling(length(knt_text_split) / 2)
          knt_text1 <- paste0(paste(knt_text_split[1:mid], collapse = ", "), ",")
          knt_text2 <- paste(knt_text_split[(mid + 1):length(knt_text_split)], collapse = ", ")
          
          knots_text1 <- bquote(italic(bold(t)[1 * ";" * k[1] * "," * 2]) == italic(.(paste0("{", knt_text1))))
          knots_text2 <- bquote(italic(.(paste0(knt_text2, "}"))))
        } else {
          knots_text1 <- bquote(italic(bold(t)[1 * ";" * k[1] * "," * 2]) == italic(.(paste0("{", knt_text,"}"))))
          knots_text2 <- ""
        }
        
        mtext(knots_text1, side = 1, line = 3, cex = 0.9)
        mtext(knots_text2, side = 1, line = 4, cex = 0.9)
      }
      
      ## Each branch now adds specific elements to the plot based on the selected type
      # 1) polygon
      if (type == "polygon") {
        results$Polykn <- temp$Poly$Kn
        results$Polyth <- if (x$type == "LM - Univ") temp$Poly$thetas else if (x$type == "GLM - Univ") family$linkinv(temp$Poly$thetas)
        
        lines(results$Polykn, results$Polyth, col = "blue", lty = 2)
        points(results$Polykn, results$Polyth, col = "blue")
        
        if (missing(legend.text)) legend.text <- c("Data", toprint, "polygon", f_legend)
        
        if(draw.legend) legend(legend.pos, legend.text, lty = c(NA, 1, 2, f_lty),
                               col = c(col_data, col_lines, "blue", f_col),
                               pch = c(pch_data, NA, 1, f_pch),
                               lwd = c(NA, 1, 1, f_lwd), bty = "n")
        
        ## nci/aci
      } else if(type == "nci" || type == "aci") {
        if (missing(legend.text)) legend.text <- c("Data", toprint, "CI", f_legend)
        
        if (x$type == "LM - Univ") {
          
          # 2) Normal Confidence Intervals
          if (type == "nci") {
            lines(X,temp$nci$Upp, col = col_lines, lwd = 2, lty = 2)
            lines(X,temp$nci$Low, col = col_lines, lwd = 2, lty = 2)
            
            results$CIupp <- temp$nci$Upp
            results$CIlow <- temp$nci$Low
            
            if(draw.legend) legend(legend.pos, legend.text, lty = c(NA, 1, 2, f_lty),
                                   col = c(col_data, col_lines, col_lines, f_col),
                                   pch = c(pch_data, NA, NA, f_pch),
                                   lwd = c(NA, 1, 1, f_lwd), bty = "n")
            
            # 3) Asymptotic Confidence Intervals
          } else if (type == "aci") {
            lines(X,temp$aci$Upp, col = col_lines, lwd = 2, lty = 2)
            lines(X,temp$aci$Low, col = col_lines, lwd = 2, lty = 2)
            results$CIupp <- temp$aci$Upp
            results$CIlow <- temp$aci$Low
            
            if(draw.legend) legend(legend.pos, legend.text, lty = c(NA, 1, 2, f_lty),
                                   col = c(col_data, col_lines, col_lines, f_col),
                                   pch = c(pch_data, NA, NA, f_pch),
                                   lwd = c(NA, 1, 1, f_lwd), bty = "n")
          }
          
        } else if (x$type == "GLM - Univ") {
          
          # 2) Normal Confidence Intervals
          if (type == "nci") {
            yy <- Y
            xx <- if (n != 2L) {
              temp$basis
            } else {
              x$linear.fit$basis
            }
            
            matrice <- splineDesign(knots = sort(c(int.knt,rep(extr,2))), derivs = rep(0,length(X)),
                                    x = X, ord = n, outer.ok = TRUE)
            
            temp_nci <- glm(yy ~ -1+xx, offset = offset, weights = weights, family = family, start = temp$thetas)
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
          } else if (type == "aci") {
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
          }
        }
        
      } else {
        if (missing(legend.text)) legend.text <- c("Data", toprint, f_legend)
        if (draw.legend) legend(legend.pos, legend.text, lty = c(NA, 1, f_lty),
                                col = c(col_data, col_lines, f_col), pch = c(pch_data, NA, f_pch),
                                lwd = c(NA, 1, f_lwd), bty = "n")
      }
    }
    
    ## If dev = TRUE, produce plot representing the deviance at each iteration of stage A
    if (x$type == "LM - Univ") {
      print <- "rss"
    } else if (x$type == "GLM - Univ") {
      print <- "dev"
    }
    
    if (dev) {
      # i. rss/dev
      main0 <- if(missing(main)) print else main
      plot(seq(0, x$iters-1), x$rss/length(X), main = main0, xlab = "Internal knots", ylab = "", ...)
      rug(seq(1, x$iters-1 - q))
      if (slp > 0) {
        Sys.sleep(slp)
      }
      # ii. \sqrt{rss} / \sqrt{dev}
      main0 <- if(missing(main)) bquote(sqrt(.(print))) else main
      plot(seq(0, x$iters-1), (x$rss/length(X))^.5, main = main0, xlab = "Knots", ylab = "", ...)
      rug(seq(1, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
      # iii. \phi
      main0 <- if(missing(main)) expression(phi) else main
      plot(seq(q, x$iters-1), (x$rss[(1+q):x$iters]/x$rss[(1):(x$iters-q)]), main = main0,
           xlab = "Knots", ylab = "", ...)
      rug(seq(2, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
      # iv. \sqrt{\phi}
      main0 <- if(missing(main)) expression(sqrt(phi)) else main
      plot(seq(q, x$iters-1), (x$rss[(1+q):x$iters]/x$rss[(1):(x$iters-q)])^.5, main = main0, xlab="Knots", ylab = "", ...)
      rug(seq(2, x$iters-1 - q))
      if(slp > 0) {
        Sys.sleep(slp)
      }
    }
    
    ##############################
    ## 2. Bivariate NGeDS/GGeDS ##
    ##############################
  } else if (x$type == "LM - Biv" || x$type == "GLM - Biv") {
    
    if(n == 2L) {
      obj <- x$linear.fit
    } else if (n == 3L) {
      obj <- x$quadratic.fit
    } else if (n == 4L) {
      obj <- x$cubic.fit
    }
    
    # Extract parametric component of the predictor model + X, Y extremes
    W <- x$args$W; Xextr <- x$args$Xextr; Yextr <- x$args$Yextr
    
    # Default labels
    var_labels <- as.character(attr(x$terms, "variables"))[3]
    var_labels <- regmatches(var_labels, regexec("f\\(([^,]+), ([^)]+)\\)", var_labels))
    xname <- var_labels[[1]][2]; yname <- var_labels[[1]][3]
    xlab <- if (!is.null(others$xlab)) others$xlab else xname
    ylab <- if (!is.null(others$ylab)) others$ylab else yname
    zlab <- if (!is.null(others$zlab)) {
      others$zlab
    } else if (!is.null(f)) {
      var_labels[[1]][1]
    } else {
      gsub("f", "f_hat", var_labels[[1]][1], fixed = TRUE)
    }
    
    # Ensure they are removed from the 'others' list before passing '...'
    others <- others[!names(others) %in% c("xlab", "ylab", "zlab")]
    
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
      if(missing(main)) {
        main <- paste0(x$Nintknots$X, " internal knots in X and ",
                       x$Nintknots$Y ," internal knots in Y" )
      }
      
      # Plot
      persp3D(x = seq_valX, y = seq_valY, z = f_XY_val, phi = 25, theta = 50,
              zlim = range(f_XY_val), ticktype = "detailed", expand = 0.5,
              colkey = FALSE, border = "black",
              xlab = xlab, ylab = ylab, zlab = zlab)
      # Add the data and predicted points to the plot
      points3D(x = X, y = Y, z = Z, col = "black", pch = 19, add = TRUE)
      points3D(x = X, y = Y, z = obj$predicted, col = "red", pch = 19, add = TRUE)
      legend("topright",
             legend = c("Data", "Quadratic Fit"),
             col = c("black", "red"),
             pch = 19,
             bg = 'white',
             bty = "n")
      
    } else {
      
      df <- data.frame(X,Y)
      
      # Create a sequence of sqrt(obs) evenly spaced numbers over the range of X/Y
      newX <- seq(from = Xextr[1], to = Xextr[2], length.out = round(sqrt(length(X))))
      newY <- seq(from = Yextr[1], to = Yextr[2], length.out = round(sqrt(length(Y))))
      
      # Create a grid data frame from all combinations of newX and newY
      grid.data <- expand.grid(newX, newY)
      
      # Generate spline basis matrix for X and Y dimensions using object knots and given order
      basisMatrixX <- splineDesign(knots = obj$Xknots, derivs = rep(0,length(grid.data[,1])),
                                   x = grid.data[,1], ord = n, outer.ok = T)
      basisMatrixY <- splineDesign(knots = obj$Yknots, derivs = rep(0,length(grid.data[,2])),
                                   x = grid.data[,2], ord = n, outer.ok = T)
      # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
      basisMatrixBiv <- tensorProd(basisMatrixX, basisMatrixY)
      # Multiply the bivariate spline basis by model coefficients to get fitted values
      f_XY_hat_val <- basisMatrixBiv %*% obj$theta[1:dim(basisMatrixBiv)[2]]
      # Reshape the fitted values to a square matrix for plotting
      f_XY_hat_val <- matrix(f_XY_hat_val, nrow = length(newX) )
      if (x$type == "GLM - Biv") f_XY_hat_val <- x$args$family$linkinv(f_XY_hat_val)
      # Title based on the number of internal knots in X and Y
      if(missing(main)) {
        main <- paste0(x$Nintknots$X, " internal knots in X and ",
                       x$Nintknots$Y ," internal knots in Y" )
      }
      
      # Plot the perspective 3D surface defined by newX, newY, and the fitted values
      persp3D(x = newX, y = newY, z = f_XY_hat_val, phi = 25, theta = 50,
              zlim = range(f_XY_hat_val), ticktype = "detailed", expand = 0.7,
              colkey = FALSE, border = "black",
              xlab = xlab, ylab = ylab, zlab = zlab)
      
      if (main == "detail") {
        # KnotsX
        knotsX_values <- round(knots(x, n=2)$Xk,2)
        intkntX_text <- if (!is.null(get_internal_knots(knotsX_values,2))) {
          paste0(", ", get_internal_knots(knotsX_values,2))
        } else {
          ""
        }
        kntX_text <- paste0(
          "a = ",
          knotsX_values[1],
          paste(intkntX_text, collapse = ""),
          ", b = ",
          knotsX_values[length(knotsX_values)]
        )
        if (length(kntX_text) > 18) {
          kntX_text1 <- paste(kntX_text[1:round(length(kntX_text) / 2)], collapse = ", ")
          kntX_text2 <- paste(kntX_text[(round(length(kntX_text) / 2) + 1):length(kntX_text)], collapse = ", ")
          
          knotsX_text1 <- bquote(italic(bold(t)[1 * ";" * k[1] * "," * 2]) == italic(.(paste0("{", kntX_text1))))
          knotsX_text2 <- bquote(italic(.(paste0(knt_text2, "}"))))
        } else {
          knotsX_text1 <- bquote(italic(bold(t)[1 * ";" * k[1] * "," * 2]) == italic(.(paste0("{", kntX_text,"}"))))
          knotsX_text2 <- ""
        }
        
        mtext(knotsX_text1, side = 3, line = 2, adj = 0.5, cex = 0.7)
        mtext(knotsX_text2, side = 3, line = 1, adj = 0.5, cex = 0.7)
        
        # Knots Y
        knotsY_values <- round(knots(x, n=2)$Yk,2)
        intkntY_text <- if (!is.null(get_internal_knots(knotsY_values,2))) {
          paste0(", ", get_internal_knots(knotsY_values,2))
        } else {
          ""
        }
        kntY_text <- paste0(
          "a = ",
          knotsY_values[1],
          paste(intkntY_text, collapse = ""),
          ", b = ",
          knotsY_values[length(knotsY_values)]
        )
        
        if (length(kntY_text) > 18) {
          kntY_text1 <- paste(kntY_text[1:round(length(kntY_text) / 2)], collapse = ", ")
          kntY_text2 <- paste(kntY_text[(round(length(kntY_text) / 2) + 1):length(kntY_text)], collapse = ", ")
          
          knotsY_text1 <- bquote(italic(bold(t)[1 * ";" * k[2] * "," * 2]) == italic(.(paste0("{", kntY_text1))))
          knotsY_text2 <- bquote(italic(.(paste0(knt_text2, "}"))))
        } else {
          knotsY_text1 <- bquote(italic(bold(t)[2 * ";" * k[2] * "," * 2]) == italic(.(paste0("{", kntY_text,"}"))))
          knotsY_text2 <- ""
        }
        mtext(knotsY_text1, side = 3, line = 0, adj = 0.5, cex = 0.7)
        mtext(knotsY_text2, side = 3, line = -1, adj = 0.5, cex = 0.7)
      }
      
      
      # Adjust Z by subtracting the parametric component of the predictor model if it exists
      if (!is.null(W)) {
        Z <- Z - W %*% x$linear.fit$theta[-(1:dim(basisMatrixBiv)[2])]
      }
      # Identify points where Z is greater/smaller than the model predicted value
      tmp <- (obj$predicted - Z > 0)
      points3D(x = X[!tmp], y = Y[!tmp], z = Z[!tmp], col = "red", pch = 19, add = TRUE)
      points3D(x = X[tmp], y = Y[tmp], z = Z[tmp], col = "blue", pch = 19, add = TRUE)
      
      legend.text <- if (!is.null(legend.text)) legend.text else c(expression(z >= hat(f)(x, y)), expression(z < hat(f)(x, y)))
      legend("topright",
             legend = legend.text,
             col = c("red", "blue"),
             pch = 19,
             bg = 'white',
             bty = "n")
      
    }
    
  } else {
    stop("type not recognized")
  }
  
  x$Plotinfo <- results
  invisible(x)
}


################################################################################
################################################################################
############################ plot.GeDSboost method #############################
################################################################################
################################################################################
#' @title Plot Method for GeDSboost Objects
#' @name plot.GeDSboost
#' @description
#' Plots the component functions of a GeDSboost object fitted using
#' \code{\link{NGeDSboost}}. If the model has a single base-learner, the plot
#' will be returned on the response scale. Otherwise, plots are produced on the
#' linear predictor scale. Note that only univariate base-learner plots are
#' returned, as the representation of the boosted model as a single spline model
#' is available only for univariate base-learners (see Dimitrova et al. (2025)).
#' Additionally, since component-wise gradient boosting inherently performs
#' base-learner selection, plots will only be generated for the base-learners
#' selected during the boosting iterations.
#' @param x A GeDSboost object, as returned by \code{NGeDSboost()}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the FGB-GeDS fit to be extracted.
#' @param ... Further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' 
#' @references
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. L. (2025).
#' \pkg{GeDS}: An \proglang{R} Package for Regression, Generalized Additive
#' Models and Functional Gradient Boosting, based on Geometrically Designed
#' (GeD) Splines. \emph{Manuscript submitted for publication.}
#' @examples
#' data(mtcars)
#' # Convert specified variables to factors
#' categorical_vars <- c("cyl", "vs", "am", "gear", "carb")
#' mtcars[categorical_vars] <- lapply(mtcars[categorical_vars], factor)
#' N <- nrow(mtcars); ratio <- 0.8
#' set.seed(123)
#' trainIndex <- sample(1:N, size = floor(ratio * N))
#' # Subset the data into training and test sets
#' train <- mtcars[trainIndex, ]
#' test <- mtcars[-trainIndex, ]
#' Gmodboost <- NGeDSboost(mpg ~ cyl + f(drat) + f(wt) + f(hp) + vs + am,
#'                         data = train, phi = 0.7, shrinkage = 0.9, initial_learner = FALSE)
#' 
#' par(mfrow = c(2,3))
#' plot(Gmodboost, n = 2)
#' 
#' @seealso \code{\link{NGeDSboost}}
#' @method plot GeDSboost
#' @importFrom graphics plot lines legend abline rug barplot par
#' @importFrom splines splineDesign 
#' @export

plot.GeDSboost <- function(x, n = 3L,...)
{
  
  # Check if x is of class "GeDSboost"
  if(!inherits(x, "GeDSboost")) {
    stop("The input 'x' must be of class 'GeDSboost'")
  }
  
  base_learners <- x$args$base_learners
  Y <- x$args$response[[1]]
  if (x$args$family@name == "Negative Binomial Likelihood (logit link)") {
    Y <- (Y + 1)/2
  }
  
  # 1) FGB (one single base-learner)
  if(length(base_learners) == 1) {
    
    if (length(base_learners[[1]]$variables) == 2) {
      stop("Single spline model representation of the boosted model is only available for univariate base-learners.")
    }
    
    X_mat <- x$args$predictors[[1]]
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
    
    if (!is.null(additional_params$ylim)) {
      plot(X_mat, Y, main = main0, ...)  
    } else {
      plot(X_mat, Y, main = main0, ylim = y_range, ...)
    }
    
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
    
    # 2) Componentwise FGB (more than one base-learner)
  } else {
    
    # Check if the length of the variables element is equal to 2
    lapply(base_learners, function(learner) {
      if (length(learner$variables) == 2) {
        message("Single spline model representation of the boosted model is only available for univariate base-learners.")
      }
    })
    
    base_learners <- Filter(function(learner) length(learner$variables) == 1, base_learners)
    
    # Plot only base-learners that were selected
    bl_selected <- unique(unlist(lapply(x$models, function(model) model$best_bl$name)))
    # Keep the order of base_learners, but only include those that are in bl_selected
    base_learners <- base_learners[names(sapply(base_learners, function(learner) learner$name %in% bl_selected))]
    
    pred_vars <- x$args$predictors
    
    if (n == 2) {
      theta <- x$final_model$linear.fit$theta
      int.knots <- x$internal_knots$linear.int.knots
    } else if (n == 3) {
      theta <- x$final_model$quadratic.fit$theta
      int.knots <- x$internal_knots$quadratic.int.knots
    } else if (n == 4) {
      theta <- x$final_model$cubic.fit$theta
      int.knots <- x$internal_knots$cubic.int.knots
    }
    
    for (bl_name in names(base_learners)) {
      
      bl <- base_learners[[bl_name]]
      X_mat <- pred_vars[, intersect(bl$variables, colnames(pred_vars))]
      
      int.knt <- int.knots[[bl_name]]
      
      pattern <- paste0("^", gsub("([()])", "\\\\\\1", bl_name))
      theta_bl <- theta[grep(pattern, names(theta))]
      # Replace NA values with 0
      theta_bl[is.na(theta_bl)] <- 0
      
      if (bl$type == "GeDS") {
        # Including int.knt and giving more values enhances visualization
        X_mat <- seq(from = min(X_mat), to = max(X_mat), length.out = 1000)
        X_mat <- sort(c(X_mat, int.knt))
        # Create spline basis matrix using specified knots, evaluation points and order
        basisMatrix <- splineDesign(knots = sort(c(int.knt,rep(range(X_mat),n))),
                                    x = X_mat, ord = n, derivs = rep(0,length(X_mat)),
                                    outer.ok = T)
        predicted <- basisMatrix %*% theta_bl
        ylab <-  bl_name
      } else if (bl$type == "linear") {
        # Linear
        if (!is.factor(X_mat)) {
          predicted <- theta_bl * X_mat
          # Factor
        } else {
          names(theta_bl) <- levels(X_mat)[-1]
          theta_bl[levels(X_mat)[1]] <- 0 # set baseline coef to 0
        }
        ylab <- bquote(beta[1] %*% .(bl_name))
      }
      
      if (!is.factor(X_mat)) {
        plot_data <- data.frame(X_mat, predicted)
        plot_data <- plot_data[order(plot_data$X_mat),]
        
        plot(plot_data, type = "l", col = "steelblue",
             xlab = bl$variables, ylab = ylab,
             xlim = range(X_mat), ylim = range(predicted), ...)
        
        if (length(int.knt) < 20) {
          knt <- c(int.knt, range(X_mat))
          for (k in knt) {
            abline(v = k, col = "gray", lty = 2)
          }
        } else {
          rug(int.knt)
        }
        
      } else {
        par(mar = c(7.1, 4.1, 4.1, 2.1))
        barplot(theta_bl,
                las = 2,
                col = "steelblue",
                main = bl_name,
                ylab = bquote(beta))
        par(mar = c(5.1, 4.1, 4.1,2.1))
      }
    }
    
  }
}

################################################################################
################################################################################
############################# plot.GeDSgam method ##############################
################################################################################
################################################################################
#' @title Plot Method for GeDSgam Objects
#' @name plot.GeDSgam
#' @description
#' Plots on the linear predictor scale the component functions of a GeDSgam
#' object fitted using \code{\link{NGeDSgam}}.
#' 
#' @param x A GeDSgam object, as returned by \code{NGeDSgam()}.
#' @param n Integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the GAM-GeDS fit.
#' @param base_learners Either NULL or a vector of character string specifying
#' the base-learners of the model for which predictions should be plotted. Note
#' that single base-learner predictions are provided on the linear predictor scale.
#' @param f (Optional) specifies the underlying component function or generating
#' process to which the model was fit. This parameter is useful if the user wishes
#' to plot the specified function/process alongside the model fit and the data.
#' @param ... Further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' @examples
#' ## Gu and Wahba 4 univariate term example ##
#' # Generate a data sample for the response variable
#' # y and the covariates x0, x1 and x2; include a noise predictor x3
#' set.seed(123)
#' N <- 400
#' f_x0x1x2 <- function(x0,x1,x2) {
#'   f0 <- function(x0) 2 * sin(pi * x0)
#'   f1 <- function(x1) exp(2 * x1)
#'   f2 <- function(x2) 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
#'   f <- f0(x0) + f1(x1) + f2(x2)
#'   return(f)
#'}
#' x0 <- runif(N, 0, 1)
#' x1 <- runif(N, 0, 1)
#' x2 <- runif(N, 0, 1)
#' x3 <- runif(N, 0, 1)
#' # Specify a model for the mean of y
#' f <- f_x0x1x2(x0 = x0, x1 = x1, x2 = x2)
#' # Add (Normal) noise to the mean of y
#' y <- rnorm(N, mean = f, sd = 0.2)
#' data <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2, x3 = x3)
#' 
#' # Fit GAM-GeDS model
#' Gmodgam <- NGeDSgam(y ~ f(x0) + f(x1) + f(x2) + f(x3), data = data)
#' 
#' f0 <- function(x0) 2 * sin(pi * x0)
#' f1 <- function(x1) exp(2 * x1)
#' f2 <- function(x2) 0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * x2)^3 * (1 - x2)^10
#' fs <- list(f0, f1, f2)
#' main_f0 <- expression(
#'   f[0](x[0]) == 2 * sin(pi * x[0])
#'   )
#' main_f1 <- expression(
#'   f[1](x[1]) == e^(2 * x[1])
#'   )
#' main_f2 <- expression(
#'   f[2](x[2]) == 0.2 * x[2]^11 * (10 * (1 - x[2]))^6 +
#'     10 * (10 * x[2])^3 * (1 - x[2])^10
#'   )
#' mains <- list(main_f0, main_f1, main_f2)
#' 
#' # Create and display the plot
#' par(mfrow = c(1,3), mar = c(5.1, 5.1, 4.1, 2.1))
#' for (i in 1:3) {
#'   # Plot the base learner
#'   plot(Gmodgam, n = 3, base_learners = paste0("f(x", i-1, ")"), f = fs[[i]],
#'   main = mains[[i]], col = "seagreen",
#'   cex.lab = 2, cex.axis = 2, cex.main = 1.5)
#' # Add legend
#' if (i == 2) {
#'   position <- "topleft"
#' } else if (i == 3) {
#'   position <- "topright" 
#' } else {
#'   position <- "bottom"
#' }  
#'   legend(position, legend = c("NGeDSgam quad.", paste0("f(x", i-1, ")")),
#'          col = c("seagreen", "darkgray"),
#'          lwd = c(2, 2),
#'          bty = "n",
#'          cex = 1.5)
#' }
#' 
#' @seealso \code{\link{NGeDSgam}}
#' @method plot GeDSgam
#' @importFrom plot3D persp3D segments3D
#' @importFrom graphics plot lines abline rug barplot par
#' @importFrom splines splineDesign
#' @export

plot.GeDSgam <- function(x, base_learners = NULL,
                         f = NULL, n = 3L,...)
{
  
  # Check if x is of class "GeDSgam"
  if(!inherits(x, "GeDSgam")) {
    stop("The input 'x' must be of class 'GeDSgam'")
  }
  
  if (is.null(base_learners)) {
    base_learners <- x$args$base_learners
  } else {
    base_learners <- x$args$base_learners[names(x$args$base_learners) %in% base_learners]
  }
  
  # Other arguments
  others <- list(...)
  
  
  Y <- x$args$response[[1]]; pred_vars <- x$args$predictors
  
  if (n == 2) {
    theta <- x$final_model$linear.fit$theta
    int.knots <- x$internal_knots$linear.int.knots
  } else if (n == 3) {
    theta <- x$final_model$quadratic.fit$theta
    int.knots <- x$internal_knots$quadratic.int.knots
  } else if (n == 4) {
    theta <- x$final_model$cubic.fit$theta
    int.knots <- x$internal_knots$cubic.int.knots
  }
  
  
  # If you provide more/less f than base-learners stop
  if (!is.null(f) && length(f) != length(base_learners)) {
    stop(sprintf("Base learners to plot are %d, but %d functions were provided", 
                 length(base_learners), length(f)))
  }
  
  
  # Plot each base-learner
  for (bl_name in names(base_learners)) {
    
    bl <- base_learners[[bl_name]]
    X_mat <- pred_vars[, intersect(bl$variables, colnames(pred_vars))]
    
    int.knt <- int.knots[[bl_name]]
    
    pattern <- paste0("^", gsub("([()])", "\\\\\\1", bl_name))
    theta <- theta[grep(pattern, names(theta))]
    # Replace NA values with 0
    theta[is.na(theta)] <- 0
    
    ##########################
    # 1. Univariate learners #
    ##########################
    if (NCOL(X_mat) == 1) {
      
      # 1.1. Univariate GeDS
      if (bl$type == "GeDS") {
        # Including int.knt and giving more values enhances visualization
        X_mat <- seq(from = min(X_mat), to = max(X_mat), length.out = 1000)
        X_mat <- sort(c(X_mat, int.knt))
        # Create spline basis matrix using specified knots, evaluation points and order
        basisMatrix <- splineDesign(knots = sort(c(int.knt,rep(range(X_mat),n))),
                                    x = X_mat, ord = n, derivs = rep(0,length(X_mat)),
                                    outer.ok = T)
        # To recover backfitting predictions need de_mean
        predicted <- if (n == 2) basisMatrix %*% theta - mean(basisMatrix %*% theta) else basisMatrix %*% theta 
        ylab <-  bl_name
        
        # 1.2. Univariate Linear
      } else if (bl$type == "linear") {
        # Linear
        if (!is.factor(X_mat)) {
          predicted <- theta * X_mat
          # Factor
        } else {
          names(theta) <- levels(X_mat)[-1]
          theta[levels(X_mat)[1]] <- 0 # set baseline coef to 0
        }
        ylab <- bquote(beta[1] %*% .(bl_name))
      }
      
      # Non-factor
      if  (!is.factor(X_mat)) {
        
        if (!is.null(f)) predicted <- predicted + mean(f(X_mat)-predicted)
        
        plot_data <- data.frame(X_mat, predicted)
        plot_data <- plot_data[order(plot_data$X_mat),]
        
        ylim <- if (is.null(others$ylim)) range(predicted) else others$ylim
        col_fit <- if (is.null(others$col)) "red" else others$col
        
        plot(plot_data, type = "l",
             xlab = bl$variables, ylab = ylab,
             xlim = range(X_mat), lwd = 2,...)
        if (!is.null(f)){
          f_data <- data.frame(X_mat, f(X_mat))
          f_data <- f_data[order(f_data$X_mat),]
          lines(f_data, lwd = 2, col = "darkgray")
        }
        
        if (length(int.knt) < 20) {
          knt <- c(int.knt, range(X_mat))
          for(k in knt) {
            abline(v = k, col = "gray", lty = 2)
          }
        } else {
          rug(int.knt)
        }
        
        # Factor
      } else {
        par(mar = c(7.1, 4.1, 4.1, 2.1))
        barplot(theta,
                las = 2,
                col = col_fit,
                main = bl_name,
                ylab = bquote(beta))
        par(mar = c(5.1, 4.1, 4.1,2.1))
      }
      
      #########################
      # 2. Bivariate learners #
      #########################
    } else if (NCOL(X_mat) == 2) {
      
      if (!is.null(f)) warning("f not available for bivariate learners")
      
      Xextr <- range(X_mat[,1])
      Yextr <- range(X_mat[,2])
      
      # Create a sequence of sqrt(obs) evenly spaced numbers over the range of X/Y
      adj_rangeX <- diff(Xextr)/10; adj_rangeY <- diff(Yextr)/10 # constraining the range prevents extremes in the edges
      newX <- seq(from = Xextr[1]+adj_rangeX, to = Xextr[2]-adj_rangeX, length.out = round(sqrt(length(X_mat[,1]))))
      newY <- seq(from = Yextr[1]+adj_rangeY, to = Yextr[2]-adj_rangeY, length.out = round(sqrt(length(X_mat[,2]))))
      
      # Create a grid data frame from all combinations of newX and newY
      grid.data <- expand.grid(newX, newY)
      
      # Generate spline basis matrix for X and Y dimensions using object knots and given order
      basisMatrixX <- splineDesign(knots = sort(c(int.knt$ikX,rep(Xextr,n))), derivs = rep(0,length(grid.data[,1])),
                                   x = grid.data[,1], ord = n, outer.ok = T)
      basisMatrixY <- splineDesign(knots = sort(c(int.knt$ikY,rep(Yextr,n))), derivs = rep(0,length(grid.data[,2])),
                                   x = grid.data[,2], ord = n, outer.ok = T)
      # Calculate the tensor product of X and Y spline matrices to create a bivariate spline basis
      basisMatrixBiv <- tensorProd(basisMatrixX, basisMatrixY)
      # Multiply the bivariate spline basis by model coefficients to get fitted values
      f_XY_hat_val <- basisMatrixBiv %*% theta[1:dim(basisMatrixBiv)[2]]
      # Reshape the fitted values to a square matrix for plotting
      f_XY_hat_val <- matrix(f_XY_hat_val, nrow = round(sqrt(length(X_mat[,1]))))
      f_XY_hat_val <- x$args$family$linkinv(f_XY_hat_val)
      # Title based on the number of internal knots in X and Y
      title <- paste0(length(int.knt$ikX)," internal knots in ", bl$variables[1]," and ",
                      length(int.knt$ikY), " internal knots in ", bl$variables[2] )
      
      # Plot the perspective 3D surface defined by newX, newY, and the fitted values
      persp3D(x = newX, y = newY, z = f_XY_hat_val, main = title, phi = 25, theta = 50,
              xlab = paste0("\n", bl$variables[1]), ylab = paste0("\n",bl$variables[2]),
              zlab = paste0("\n", bl_name), zlim = range(f_XY_hat_val),
              ticktype = "detailed", expand = 0.5, colkey = FALSE, border = "black", ...)
      
      # Add rug plots to the x axis
      kntX <- c(int.knt$ikX, Xextr)
      segments3D(x0 = kntX, y0 = rep(min(newY), length(kntX)), z0 = rep(min(f_XY_hat_val), length(kntX)),
                 x1 = kntX, y1 = rep(min(newY), length(kntX)), z1 = rep(max(f_XY_hat_val), length(kntX)),
                 add = TRUE, col = "black", lty = 2)
      
      # Add rug plots to the y axis
      kntY <- c(int.knt$ikY, Yextr)
      segments3D(x0 = rep(max(newX), length(kntY)), y0 = kntY, z0 = rep(min(f_XY_hat_val), length(kntY)),
                 x1 = rep(max(newX), length(kntY)), y1 = kntY, z1 = rep(max(f_XY_hat_val), length(kntY)),
                 add = TRUE, col = "black", lty = 2)
    }
    
    
  }
  
}



