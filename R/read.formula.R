################################################################################
################################################################################
################################### formula ####################################
################################################################################
################################################################################

#######################################################
# Function for reading the formula of NGeDS and GGeDS #
#######################################################
#' @importFrom stats terms model.matrix model.frame model.response na.omit
read.formula <- function(formula, data, weights, offset)
  {
  # Determine model terms based on 'formula', using 'data' if provided
  mt <- if (missing(data)) {
    terms(formula, "f")
    } else terms(formula, "f", data = data)
  # Remove intercept from model terms
  if (attr(mt,"intercept") == 1) {
    attr(mt,"intercept") <- 0
    } else {
      # In case the intercept was explicitly removed in the formula
      warning("An intercept will be included in the basis functions")
      }
  data <- as.list(data)
  data$f <- function(x, xx = NULL, ...) {
    coordinates <- c(list(x), if (!is.null(xx)) list(xx), list(...))
    do.call(cbind, coordinates)
  }

  # Locate f(X)/f(X,Y) in model terms
  spec <- attr(mt,"specials")$f
  if(length(spec)!= 1) stop("Formula incorrectly specified. Read documentation for further information.")

  # Generate a model matrix based on the model terms and data
  mm <- model.matrix(mt, data)
  # Create a model frame based on the model terms and data, omitting rows with NAs
  mf <- model.frame(mt, data, na.action = na.omit)

  # Extract response variable
  if (is.factor(model.response(mf))) {
    # encode factor response to 0/1
    Y <- as.numeric(model.response(mf)) - 1
    if (any(Y != 0 & Y != 1)) stop("The factor variable is not binary.")
    } else {
      Y <- model.response(mf, type = "numeric")
    }
  attr(Y,"names")<- NULL
  # Extract GeDS covariates
  X <- mf[,spec]
  spline_vars <- all.vars(attr(mt, "variables")[[spec + 1L]])
  X <- as.matrix(X)
  colnames(X) <- spline_vars
  # Extract linear covariates
  if(ncol(mm)>ncol(X)) {
    Z <- mf[, -c(spec, attr(mt, "response"), attr(mt, "offset")), drop = TRUE]
  } else {
    Z <- NULL
  }

  # Guard against a parametric term collinear with the spline component:
  # a GeD spline basis already reproduces constant + linear functions, so a
  # bare linear term in the spline variable is not separately identifiable.
  if (!is.null(Z)) {
    f_idx        <- grep("^f\\(", attr(mt, "term.labels"))
    param_labels <- attr(mt, "term.labels")[-f_idx]
    param_vars   <- unique(unlist(lapply(param_labels,
                                         function(l) all.vars(parse(text = l)[[1]]))))
    clash <- intersect(spline_vars, param_vars)
    if (length(clash)) {
      stop(sprintf(
        "Variable '%s' enters both the spline component f(%s) and the parametric part.\n  The linear term is collinear with the GeD spline basis, so the two are not separately identifiable.\n  Fit 'Y ~ f(%s)' instead.",
        clash[1], spline_vars[1], spline_vars[1]), call. = FALSE)
    }
  }

  # Check if formula variables make sense
  response_is_coordinate <- any(vapply(
    seq_len(ncol(X)),
    function(j) all(Y == X[, j]),
    logical(1)
  ))
  if (response_is_coordinate || (!is.null(Z) && all(Y == Z))) {
    stop("Response variable cannot be equal to a covariate.")
  }
  if (ncol(X) > 1L) {
    equal_coordinates <- utils::combn(
      seq_len(ncol(X)),
      2L,
      FUN = function(index) all(X[, index[1L]] == X[, index[2L]])
    )
    if (any(equal_coordinates)) {
      stop("Covariates in the GeDS function must be distinct.")
    }
  }

  # Initialize an offset vector with zeros
  offset <- rep(0, nrow(X))
  # Sum up offset variables included in model terms
  if (!is.null(off.num <- attr(mt, "offset"))) {
    for (i in off.num) offset <- offset + eval(attr(mt,"variables")[[i + 1]], data)
  }

  out <- list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "terms" = mt, "model.matrix" = mm)
  return(out)
}

###############################################################
# Function for reading the formula of NGeDSboost and NGeDSgam #
###############################################################
# Recursive function to extract all calls from RHS
get_terms <- function(expr) {
  if (is.call(expr) && expr[[1]] == as.name("+")) {
    # If it's a sum, recurse on both sides
    c(get_terms(expr[[2]]), get_terms(expr[[3]]))
  } else {
    # Otherwise, return the term as a string
    deparse(expr)
  }
}
#' @importFrom stats reformulate setNames
read.formula.gam <- read.formula.boost <- function(formula, data,
                                                   type = c("gam", "boost")){


  # Collect raw terms
  terms_raw <- get_terms(formula[[3]])
  # Detect duplicates
  dups <- terms_raw[duplicated(terms_raw)]
  if (length(dups) > 0) {
    stop(paste("Duplicate learner(s) detected:", paste(unique(dups), collapse = ", ")))
  }

  formula <- as.formula(formula)

  # Check for offset
  offset_index <- attr(terms(formula, data = data), "offset")

  if (!is.null(offset_index)) {
    # Extract offset variable name
    offset_var <- as.character(attr(terms(formula), "variables")[[offset_index + 1]])[2]
    offset <- data[[offset_var]]

    # Inform the user
    if (type == "boost") {
      message("Note: An offset term ('", offset_var, "') was detected in the formula and will be ignored. ",
              "If needed, consider adjusting the response variable manually.")
    }

    # Remove offset from formula
    formula <- reformulate(attr(terms(formula), "term.labels"),
                           response = as.character(formula[[2]]))
  } else {
    offset <- NULL
  }

  # Parse formula
  terms <- all.vars(formula)
  response <- terms[1]
  predictors <- terms[-1]

  # Split the formula into LHS and RHS
  formula_string <- paste(deparse(formula), collapse = " ")
  parts <- strsplit(formula_string, " ~ ")[[1]]
  lhs <- parts[1]
  rhs <- parts[2]
  # Split the RHS into individual terms
  rhs_parts <- trimws(unlist(strsplit(rhs, " \\+ ")))
  # Reformat parts for validation: keep the first part as is, append others to LHS
  parts_to_validate <- c(lhs, paste(lhs, "~", rhs_parts))[-1]
  # Adjusted pattern for validation
  pattern <- "^[a-zA-Z0-9_.]+ ~ (\\.|(f\\([^)]+\\)|[a-zA-Z0-9_.]+))$"
  # Validate each part
  valid_parts <- sapply(parts_to_validate, function(part) {
    grepl(pattern, part)
  })
  # Check if all parts are valid
  all_parts_valid <- all(valid_parts)
  if (!all_parts_valid) {
    stop("Formula incorrectly specified. Read documentation for further information.")
  }

  # Predictors terms
  bl_part <- parts[2]
  bl_elements <- unlist(strsplit(bl_part, "\\+"))
  bl_elements <- trimws(bl_elements) # remove any leading or trailing white spaces from each element
  if(any(predictors == ".")) {
    # Handling "Y ~ ." case
    predictors <- predictors[predictors != "."]
    rest_predictors <- setdiff(names(data), c(response, predictors))
    predictors <- c(predictors, rest_predictors)
    bl_elements <- bl_elements[bl_elements != "."]
    rest_bl_elements <- gsub("(x\\.\\d+)", "f(\\1)", rest_predictors)
    bl_elements <- c(bl_elements, rest_bl_elements)
  }

  # Extract predictor variable names from each base-learner element
  base_learners <- setNames(lapply(bl_elements, function(bl) {
    # Determine the type of the base learner
    type <- ifelse(grepl("f\\(.*\\)", bl), "GeDS", "linear")
    # Extract the predictor variables
    variables <- trimws(unlist(strsplit(gsub("f\\((.*?)\\)", "\\1", bl), ",")))
    # Return a list containing variable names and the type
    return(list(variables = variables, type = type))
  }), bl_elements)

  # Check if response coincides with any of the covariates
  for (learner in base_learners) {
    if (response %in% learner$variables) stop("The response variable cannot be used as part of the predictors in the base learners.")
    }

  return(list(terms = terms, response = response, predictors = predictors,
              base_learners = base_learners, offset = offset))
}

#' @title Defining the Covariates for the Spline Component in a GeDS Formula
#' @name f
#' @description
#' In general the GeDS predictor model may include a GeD spline regression
#' component with respect to one or more independent variables and a parametric
#' component in which the remaining covariates may enter as additive terms.
#' GAM-GeDS and FGB-GeDS models may include more than one GeD spline regression
#' component.
#'
#' The function \code{f} is to be used in the
#' \code{\link[=formula.GeDS]{formula}} argument of \code{\link{NGeDS}},
#' \code{\link{GGeDS}}, \code{\link{NGeDSgam}} or \code{\link{NGeDSboost}} in
#' order to specify which independent variables (covariates) should be included
#' in the GeD spline regression component of the predictor model.
#' @param x Numeric vector containing \eqn{N} sample values of the covariate
#' chosen to enter the spline
#' regression component of the predictor model.
#' @param xx Numeric vector containing \eqn{N} sample values for the second
#' covariate (in case \code{\link{NGeDS}}/\code{\link{GGeDS}} is run for two
#' dimensions). It has to be either \code{NULL} (the default) or a vector of size
#' \eqn{N}, same as \code{x}.
#' @param ... Further covariates for an experimental multivariate Normal GeDS
#' spline. Generalized and additive-model smoothers currently remain limited to
#' one or two covariates per spline component.
#'
#' @examples
#' # Generate a data sample for the response variable Y and
#' # the covariates X, reg1, reg2 and off
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' reg1 <- runif(500, min = -0.1, max = 0.1)
#' reg2 <- runif(500, min = -0.2, max = 0.2)
#' off <- runif(500, min = -1, max = 1)
#' # Specify a model for the mean of Y to include a component non linear
#' # in X defined by the function f_1 and a linear one in the other covariates
#' means <- f_1(X) + 2*reg1 + 0.5*reg2 + off
#' # Add Normal noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Specify a formula that will be used to model Y as a
#' # function of X, reg1, reg2 and off.
#' # The covariate X is for the spline component modeled as GeDS,
#' # reg1 and reg2 enter linearly, off is an offset, i.e. no coefficient
#' # will be estimated for it
#' formula <- Y ~ f(X) + reg1 + reg2 + offset(off)
#'
#' # Fit a GeDS model specified in formula using NGeDS
#' (Gmod <- NGeDS(formula, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#'
#' @seealso \code{\link[=formula.GeDS]{formula}}; \link{NGeDS}; \link{GGeDS};
#' \link{NGeDSgam}; \link{NGeDSboost}
#'
#' @note This function is intended to be used only as part of the
#' \code{\link[=formula.GeDS]{formula}} in a GeDS model via
#' \code{\link{NGeDS}}, \code{\link{GGeDS}}, \code{\link{NGeDSgam}} or
#' \code{\link{NGeDSboost}} and not to be called in other cases by the user.

f <- function(x, xx = NULL, ...) {
  coordinates <- c(list(x), if (!is.null(xx)) list(xx), list(...))
  do.call(cbind, coordinates)
}


# this is to get the names of the Z variable(s)
getZnames <- function(out){
  names <- colnames(out$model.matrix)
  id <- attr(out$model.matrix,"assign")
  spec <- attr(out$terms,"specials")$f-1
  znames <- names[id!=spec]
  znames
}


