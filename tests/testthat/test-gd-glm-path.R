gd_glm_case <- function(kind, seed, n = 300) {
  set.seed(seed)
  x <- sort(runif(n, -2, 2))
  f <- (10 * x / (1 + 100 * x^2)) * 4 + 4
  w <- rep(1, n)
  offset <- rep(0, n)
  switch(kind,
    poisson = list(y = rpois(n, exp(f / 2)), family = poisson()),
    poisson_weighted = {
      offset <- rnorm(n, 0, 0.2)
      w <- sample(c(0, 0.5, 1, 2), n, TRUE, c(0.1, 0.3, 0.3, 0.3))
      w[c(1, n)] <- 1
      list(y = rpois(n, exp(f / 2 + offset)), family = poisson())
    },
    poisson_repeated = {
      x <- sort(round(runif(n, -2, 2), 1))
      list(y = rpois(n, exp(((10 * x / (1 + 100 * x^2)) * 4 + 4) / 2)), family = poisson())
    },
    gamma = list(y = rgamma(n, shape = 10, rate = 10 / exp(f / 2)), family = Gamma("log")),
    binomial = {
      w <- sample(1:5, n, TRUE)
      list(y = rbinom(n, w, plogis(f - 4)) / w, family = binomial())
    },
    inverse_gaussian = list(y = GeDS:::.gd_draw_inverse_gaussian(exp(f / 4), 0.05),
                            family = inverse.gaussian())
  ) -> out
  c(out, list(x = x, w = w, offset = offset))
}

gd_glm_loss <- function(case, cap, compiled, state = ".gamma_dynamic_path") {
  fit <- function() suppressWarnings(GenUnivariateFitter(case$x, case$y,
    offset = case$offset, weights = case$w, family = case$family, max.intknots = cap,
    extr = range(case$x[case$w > 0]), stoptype = state))
  if (compiled) fit() else with_mocked_bindings(fit(), .gd_glm_fast_codes = function(family) NULL)
}

test_that("the compiled GLM GD path reproduces the R loop", {
  for (kind in c("poisson", "poisson_weighted", "poisson_repeated", "gamma", "binomial",
                 "inverse_gaussian")) {
    case <- gd_glm_case(kind, 7)
    compiled <- (gd_glm_loss(case, 12L, TRUE))
    reference <- (gd_glm_loss(case, 12L, FALSE))
    expect_identical(length(compiled), length(reference), info = kind)
    expect_equal(compiled, reference, tolerance = 1e-10, info = kind)
  }
})

test_that("the compiled GLM GD path continues staged paths", {
  case <- gd_glm_case("poisson", 3)
  state <- structure(new.env(parent = emptyenv()), class = "geds_gd_path")
  first <- gd_glm_loss(case, 6L, TRUE, state)
  expect_length(first, 7L)
  expect_identical(state$restart$j, 7L)
  staged <- gd_glm_loss(case, 12L, TRUE, state)
  expect_equal(staged, (gd_glm_loss(case, 12L, TRUE)), tolerance = 1e-12)
})

test_that("unsupported GLM families stay with the R loop", {
  expect_null(.gd_glm_fast_codes(quasipoisson()))
  expect_null(.gd_glm_fast_codes(poisson(link = "sqrt")))
  expect_identical(.gd_glm_fast_codes(binomial()), c(2L, 2L))
  calls <- 0L
  native <- gd_glm_path
  local_mocked_bindings(gd_glm_path = function(...) {
    calls <<- calls + 1L
    native(...)
  })
  case <- gd_glm_case("poisson", 5)
  case$family <- quasipoisson()
  suppressWarnings(GenUnivariateFitter(case$x, case$y, family = case$family,
    max.intknots = 8L, extr = range(case$x), stoptype = ".gamma_dynamic_path"))
  expect_identical(calls, 0L)
})

test_that("steps handed back by the compiled GLM path are refitted in R", {
  case <- gd_glm_case("gamma", 11)
  native <- gd_glm_path
  handbacks <- 0L
  # Place each knot in compiled code, then hand the fit back to the R loop.
  local_mocked_bindings(gd_glm_path = function(..., j, max_intknots, last_iter) {
    if (j >= last_iter) return(native(..., j = j, max_intknots = max_intknots, last_iter = last_iter))
    out <- native(..., j = j, max_intknots = max_intknots, last_iter = j)
    if (out$status == 1L) {
      handbacks <<- handbacks + 1L
      out$status <- 3L
      out$j <- j + 1L
    }
    out
  })
  handed <- gd_glm_loss(case, 10L, TRUE)
  expect_gt(handbacks, 5L)
  expect_equal(handed, gd_glm_loss(case, 10L, FALSE), tolerance = 1e-12)
})

test_that("custom family callbacks retain their R path semantics", {
  case <- gd_glm_case("poisson", 5, n = 80)
  deviance <- case$family$dev.resids
  case$family$dev.resids <- function(y, mu, wt) 3 * deviance(y, mu, wt)
  expect_null(.gd_glm_fast_codes(case$family))
  local_mocked_bindings(gd_glm_path = function(...) {
    stop("A custom family must not enter the compiled path")
  })
  actual <- gd_glm_loss(case, 6L, TRUE)
  reference <- gd_glm_loss(case, 6L, FALSE)
  expect_equal(actual, reference, tolerance = 1e-12)
  expect_gt(length(actual), 2L)
})

test_that("compiled dispatch checks callbacks and their captured bindings", {
  for (make in list(gaussian, poisson, binomial, Gamma, inverse.gaussian)) {
    family <- make()
    expect_length(.gd_glm_fast_codes(family), 2L)
    for (key in c("variance", "dev.resids", "linkfun", "linkinv", "mu.eta",
                  "validmu", "valideta")) {
      custom <- family
      custom[[key]] <- function(x) x
      expect_null(.gd_glm_fast_codes(custom), info = paste(family$family, key))
    }
  }
  custom <- poisson()
  # Identical function body, but a different log binding in its environment.
  environment(custom$dev.resids) <- list2env(list(log = function(x) 2 * base::log(x)),
                                           parent = environment(custom$dev.resids))
  expect_null(.gd_glm_fast_codes(custom))
})
