test_that("iteration counts and plots include coefficient updates", {
  set.seed(123)
  d <- data.frame(x = sort(runif(80, -2, 2)))
  d$y <- 4 + 40 * d$x / (1 + 100 * d$x^2) + rnorm(80, sd = .2)
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 16, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  for (normalize in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(NGeDSboost(
      y ~ f(x), data = d, normalize_data = normalize, shrinkage = .25,
      max_iterations = 3L, phi_boost_exit = NULL
    )))
    ref <- fit$final_model$refinement
    expect_gt(ref$iterations, 0)
    expect_equal(fit$iters, length(fit$models) - 1L + ref$iterations)
    expect_equal(N.boost.iter(fit), fit$iters)
    expect_length(ref$coefficient_path, ref$iterations + 1L)
    expect_equal(ref$coefficient_path[[ref$selected$step + 1L]], ref$selected$delta)
    data <- d
    if (normalize) data$x <- (data$x - fit$args$X_mean[["x"]]) / fit$args$X_sd[["x"]]
    B <- GeDS:::.boost_fixed_basis(data, fit$args$base_learners[names(ref$anchor)],
                                  ref$knots, fit$args$extr, ref$levels, ref$contrasts)$B
    y <- if (normalize) (d$y - fit$args$Y_mean) / fit$args$Y_sd else d$y
    deviances <- vapply(ref$coefficient_path, function(delta) {
      sum((y - ref$initial_F_hat - as.numeric(B %*% delta))^2)
    }, numeric(1))
    expect_equal(deviances, ref$trace$deviance, tolerance = 1e-8)
    expect_equal(GeDS:::n_intknots(fit)$n_iterations, length(fit$models))
    expect_error(visualize_boosting(fit, 0:fit$iters, final_fits = TRUE), NA)
    expect_error(visualize_boosting(fit, fit$iters + 1L), "whole numbers")
    expect_error(visualize_boosting(fit, c(-1, .5, NA)), "whole numbers")
    legacy <- fit
    legacy$iters <- length(legacy$models) - 1L
    expect_equal(N.boost.iter(legacy), fit$iters)
    legacy$final_model$refinement$coefficient_path <- NULL
    expect_error(visualize_boosting(legacy, fit$iters), "refit")
  }
})

test_that("zero-update fits have no extra counted iterations", {
  d <- data.frame(x = seq(-1, 1, length.out = 40))
  d$y <- 2 + d$x
  fit <- suppressMessages(NGeDSboost(y ~ x, data = d, initial_learner = FALSE,
                                     shrinkage = 1, max_iterations = 3L))
  expect_equal(fit$final_model$refinement$iterations, 0L)
  expect_equal(fit$iters, 1L)
  expect_equal(N.boost.iter(fit), 1L)
})

test_that("coefficient-update panels retain knot titles and sequential numbers", {
  set.seed(123)
  d <- data.frame(x = sort(runif(80, -2, 2)))
  d$y <- 4 + 40 * d$x / (1 + 100 * d$x^2) + rnorm(80, sd = .2)
  fit <- suppressMessages(suppressWarnings(NGeDSboost(
    y ~ f(x), data = d, shrinkage = .25,
    max_iterations = 3L, phi_boost_exit = NULL
  )))
  expect_gt(fit$final_model$refinement$iterations, 0L)
  labels <- numeric()
  titles <- list()
  local_mocked_bindings(
    legend = function(x, legend, ...) {
      if (is.numeric(legend)) labels <<- c(labels, legend)
    },
    mtext = function(text, ...) { titles[[length(titles) + 1L]] <<- text },
    .package = "GeDS"
  )
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 16, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(N.boost.iter(fit) + 1L, 2L))
  visualize_boosting(fit, 0:N.boost.iter(fit), final_fits = TRUE)
  expect_equal(labels, seq.int(0L, 2L * N.boost.iter(fit)))
  knot_titles <- Filter(is.language, titles)
  expect_length(knot_titles, length(labels))
  expect_true(all(vapply(knot_titles, function(x) {
    grepl("Delta|delta", paste(deparse(x), collapse = ""))
  }, logical(1))))
  expect_equal(unname(graphics::par("mfg")[1:2]), c(N.boost.iter(fit) + 1L, 2L))
})
