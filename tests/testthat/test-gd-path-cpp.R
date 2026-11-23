gd_path_outputs <- function(X, Y, w, off, cap, init, fit_init, staged) {
  base <- list(X = X, Y = Y, Z = NULL, offset = off, weights = w, beta = .5, phi = .5,
               min.intknots = 0, max.intknots = cap, q = 2, extr = range(X[w > 0]),
               show.iters = FALSE, higher_order = FALSE, intknots_init = init,
               fit_init = fit_init, only_pred = FALSE)
  if (!staged) {
    base$stoptype <- ".gamma_dynamic_path"
    return(list(suppressWarnings(do.call(UnivariateFitter, base))))
  }
  state <- structure(new.env(parent = emptyenv()), class = "geds_gd_path")
  lapply(c(2, 5, 11), function(k) {
    base$max.intknots <- k
    base$stoptype <- state
    list(suppressWarnings(do.call(UnivariateFitter, base)), state$restart)
  })
}

test_that("C++ Stage A path loop matches the R loop", {
  set.seed(77)
  cases <- lapply(seq_len(120), function(r) {
    n <- sample(c(8, 30, 250), 1)
    X <- sort(if (r %% 3 == 0) round(runif(n), 2) else runif(n))
    w <- if (r %% 4 == 0) replace(rexp(n), sample(n, n %/% 4), 0) else rep(1, n)
    if (!any(w > 0)) w[1] <- 1
    off <- if (r %% 5 == 0) rnorm(n) else rep(0, n)
    Y <- sin(7 * X) + rnorm(n, sd = .2)
    init <- if (r %% 6 == 0) sort(runif(2, min(X), max(X))) else NULL
    fit_init <- NULL
    if (r %% 7 == 0) {
      b <- splines::splineDesign(sort(c(init, rep(range(X[w > 0]), 2))), X, ord = 2L,
                                 outer.ok = TRUE)
      fit_init <- list(intknots = init, pred = as.numeric(b %*% qr.coef(qr(b), Y - off)) + off)
    }
    list(X = X, Y = Y, w = w, off = off, cap = sample(c(1, 4, 15), 1), init = init,
         fit_init = fit_init, staged = r %% 2 == 0)
  })
  run <- function(a) tryCatch(do.call(gd_path_outputs, a), error = conditionMessage)
  fast <- lapply(cases, run)
  local_mocked_bindings(gd_stage_a_path = function(...) NULL, .package = "GeDS")
  reference <- lapply(cases, run)
  for (r in seq_along(cases)) expect_identical(fast[[r]], reference[[r]], info = paste("case", r))
})
