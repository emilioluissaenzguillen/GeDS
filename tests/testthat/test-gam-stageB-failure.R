test_that("a failed GAM Stage B order preserves the linear and other-order fits", {
  skip_if_not_installed("mgcv")
  d <- data.frame(x = seq(0, 1, length.out = 60))
  d$y <- sin(5 * d$x) + .1 * cos(29 * d$x)
  args <- list(formula = y ~ f(x), data = d, stoptype = "RD",
               internal_knots = 2L, max_iterations = 3L)
  reference <- suppressMessages(do.call(NGeDSgam, args))
  native <- .stageB_penalised_fit
  for (failed in 3:4) {
    local_mocked_bindings(.stageB_penalised_fit = function(..., n) {
      # n is normally passed positionally by stageB_fit.
      values <- list(...)
      order <- if (missing(n)) values[[8L]] else n
      if (order == failed) stop("Injected REML failure")
      if (missing(n)) do.call(native, values) else native(..., n = n)
    }, .package = "GeDS")
    expect_warning(model <- suppressMessages(do.call(NGeDSgam, args)),
                   class = "geds_stageB_fit_warning")
    expect_identical(predict(model, n = 2), predict(reference, n = 2))
    expect_identical(predict(model, newdata = d, n = 2),
                     predict(reference, newdata = d, n = 2))
    expect_identical(knots(model, n = 2), knots(reference, n = 2))
    other <- if (failed == 3L) 4L else 3L
    expect_identical(predict(model, n = other), predict(reference, n = other))
    for (type in c("link", "response")) {
      expect_true(all(is.na(predict(model, n = failed, type = type))))
      expect_true(all(is.na(predict(model, newdata = d, n = failed, type = type))))
    }
    key <- if (failed == 3L) "quadratic" else "cubic"
    expect_identical(model$final_model$stageB_errors[[key]], "Injected REML failure")
    expect_null(model$final_model[[paste0(key, ".fit")]])
  }
})

test_that("two failed GAM Stage B orders still return a usable linear model", {
  d <- data.frame(x = seq(0, 1, length.out = 40))
  d$y <- 1 + d$x + .1 * sin(20 * d$x)
  args <- list(formula = y ~ f(x), data = d, stoptype = "RD",
               internal_knots = 0L, max_iterations = 3L)
  ref <- suppressMessages(do.call(NGeDSgam, c(args, list(higher_order = FALSE))))
  local_mocked_bindings(.stageB_penalised_fit = function(...) stop("Injected failure"),
                       .package = "GeDS")
  notes <- list()
  model <- withCallingHandlers(suppressMessages(do.call(NGeDSgam, args)),
    warning = function(w) {
      notes[[length(notes) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(notes, 2L)
  expect_true(all(vapply(notes, inherits, logical(1), "geds_stageB_fit_warning")))
  expect_identical(predict(model, n = 2), predict(ref, n = 2))
  expect_true(all(is.finite(predict(model, newdata = d, n = 2))))
  expect_true(all(is.na(predict(model, n = 3))))
  expect_true(all(is.na(predict(model, n = 4))))
  for (type in c("link", "response")) {
    all <- predict(model, n = "all", type = type)
    expect_identical(all$pred_linear, predict(model, n = 2, type = type))
    expect_true(all(is.na(all$pred_quadratic)))
    expect_true(all(is.na(all$pred_cubic)))
    all_new <- predict(model, newdata = d, n = "all", type = type)
    expect_true(all(is.finite(all_new$pred_linear)))
    expect_true(all(is.na(all_new$pred_quadratic)))
    expect_true(all(is.na(all_new$pred_cubic)))
  }
})
