test_that("reaching one accumulated component cap stops the whole boost", {
  d <- droplevels(subset(iris, Species != "virginica"))
  fit <- suppressMessages(suppressWarnings(NGeDSboost(
    Species ~ f(Sepal.Length) + f(Sepal.Width) + f(Petal.Length) + f(Petal.Width),
    data = d, family = mboost::Binomial(), min_iterations = 10L,
    initial_learner = TRUE, boosting_with_memory = FALSE, int.knots_boost = 98L,
    penalty = "none"
  )))
  stop <- fit$args$knot_cap_stop
  expect_identical(stop$learner, "f(Petal.Length)")
  expect_equal(stop$cap, 26L)
  expect_equal(stop$iteration, length(fit$models) - 1L)
  expect_lt(stop$iteration, 10L)
  final_knots <- GeDS:::get_internal_knots(tail(fit$models, 1L)[[1L]]$base_learners[[stop$learner]]$knots)
  expect_equal(length(final_knots), stop$cap)
  for (model in fit$models) {
    for (bl in names(fit$args$component_knot_caps)) {
      count <- length(GeDS:::get_internal_knots(model$base_learners[[bl]]$knots))
      expect_lte(count, fit$args$component_knot_caps[[bl]])
    }
  }
  for (n in 2:4) {
    expect_true(all(is.finite(predict(fit, n = n))))
    expect_equal(predict(fit, n = n), predict(fit, newdata = d, n = n))
  }
})

test_that("gradient limits respect remaining room with and without memory", {
  d <- data.frame(x = seq(0, 1, length.out = 40))
  d$y <- sin(20 * d$x)
  limits <- integer()
  original <- GeDS:::.gamma_dynamic_fit
  local_mocked_bindings(.gamma_dynamic_fit = function(fitter, args, ...) {
    limits <<- c(limits, args$max.intknots)
    original(fitter, args, ...)
  }, .package = "GeDS")
  learners <- list("f(x)" = list(type = "GeDS", variables = "x"))
  cache <- list("f(x)" = list(xs = d$x, ord = seq_len(40),
                            w = rep(1, 40), extr = range(d$x)))
  for (memory in c(FALSE, TRUE)) {
    result <- GeDS:::componentwise_fit(
      bl_name = "f(x)", response = "y", data = d, model_formula_template = "y ~ ",
      weights = rep(1, 40), base_learners = learners, m = 2L,
      internal_knots = 20L, beta = .5, phi = .99, q = 2L,
      starting_intknots = if (memory) c(.2, .4, .6, .8) else NULL,
      univ_geds_cache = cache, default_spline_cap = TRUE,
      component_knot_cap = 5L, accumulated_knots = 4L
    )
    expect_false(result$error)
  }
  expect_equal(limits, c(1L, 5L))
})

test_that("zero-knot caps stop after an accepted linear update", {
  d <- data.frame(x = 1:4, y = c(1, 3, 2, 5))
  for (initial in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(NGeDSboost(
      y ~ f(x), data = d, initial_learner = initial, higher_order = FALSE,
      shrinkage = .1, phi_boost_exit = NULL, max_iterations = 5L
    )))
    expect_equal(fit$args$component_knot_caps[["f(x)"]], 0L)
    expect_equal(fit$args$knot_cap_stop$iteration, if (initial) 0L else 1L)
    expect_equal((length(fit$models) - 1L), if (initial) 0L else 1L)
    expect_true(all(is.finite(predict(fit, n = 2))))
  }
})

test_that("caps use positive-weight support and leave joint learners uncapped", {
  set.seed(735)
  d <- data.frame(x = rep(1:6, each = 4), z = runif(24), y = rnorm(24))
  w <- as.numeric(d$x <= 4)
  fit <- suppressMessages(suppressWarnings(NGeDSboost(
    y ~ f(x) + f(x, z), data = d, weights = w, initial_learner = FALSE,
    max_iterations = 1L, higher_order = FALSE
  )))
  expect_equal(fit$args$component_knot_caps, c("f(x)" = 2L))
})

test_that("a nonbinding accumulated cap leaves all spline orders unchanged", {
  set.seed(736)
  d <- data.frame(x = seq(-2, 2, length.out = 100))
  d$y <- sin(d$x) + rnorm(100, sd = .1)
  set.seed(737)
  fit <- suppressWarnings(NGeDSboost(y ~ f(x), data = d, int.knots_boost = 8L,
                                     max_iterations = 2L))
  original <- GeDS:::.geds_spline_cap
  local_mocked_bindings(.geds_spline_cap = function(...) {
    cap <- original(...)
    cap$max.intknots <- 1000L
    cap
  }, .package = "GeDS")
  set.seed(737)
  reference <- suppressWarnings(NGeDSboost(y ~ f(x), data = d, int.knots_boost = 8L,
                                           max_iterations = 2L))
  expect_null(fit$args$knot_cap_stop)
  expect_null(reference$args$knot_cap_stop)
  for (n in 2:4) expect_equal(predict(fit, n = n), predict(reference, n = n))
})

test_that("parallel candidates respect the same accumulated cap", {
  skip_on_cran()
  d <- droplevels(subset(iris, Species != "virginica"))
  formula <- Species ~ f(Sepal.Length) + f(Sepal.Width) + f(Petal.Length) + f(Petal.Width)
  serial <- suppressWarnings(NGeDSboost(formula, data = d,
                                        family = mboost::Binomial(), higher_order = FALSE))
  local_mocked_bindings(.ngedsboost_worker_count = function(...) 2L, .package = "GeDS")
  parallel_fit <- suppressWarnings(NGeDSboost(formula, data = d,
                                              family = mboost::Binomial(), higher_order = FALSE))
  expect_equal(parallel_fit$args$knot_cap_stop, serial$args$knot_cap_stop)
  expect_equal(N.boost.iter(parallel_fit), N.boost.iter(serial))
  expect_equal(predict(parallel_fit, n = 2), predict(serial, n = 2))
})
