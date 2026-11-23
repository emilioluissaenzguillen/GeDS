test_that("boosting retains successful orders after a Stage B failure", {
  d <- data.frame(x = seq(0, 1, length.out = 60))
  d$z <- cos(7 * d$x)
  d$y <- sin(5 * d$x) + .2 * d$z
  native <- GeDS:::stageB_fit
  for (formula in list(y ~ f(x), y ~ f(x) + f(z))) {
    args <- list(formula = formula, data = d, stoptype = "RD",
                 int.knots_boost = 2L, max_iterations = 3L,
                 initial_learner = FALSE)
    local_mocked_bindings(stageB_fit = native, .package = "GeDS")
    reference <- suppressMessages(do.call(NGeDSboost, args))
    for (failed in list(3L, 4L, c(3L, 4L))) {
      local_mocked_bindings(stageB_fit = function(n, ...) {
        if (n %in% failed) stop("Injected Stage B failure")
        native(n, ...)
      }, .package = "GeDS")
      notes <- list()
      model <- withCallingHandlers(suppressMessages(do.call(NGeDSboost, args)),
        warning = function(w) {
          if (inherits(w, "geds_stageB_fit_warning")) {
            notes[[length(notes) + 1L]] <<- w
            invokeRestart("muffleWarning")
          }
        })
      expect_length(notes, length(failed))
      expect_equal(vapply(notes, `[[`, numeric(1), "order"), failed)
      expect_identical(predict(model, n = 2), predict(reference, n = 2))
      expect_identical(predict(model, newdata = d, n = 2),
                       predict(reference, newdata = d, n = 2))
      expect_identical(knots(model, n = 2), knots(reference, n = 2))
      expect_identical(model$iters, reference$iters)
      for (n in 3:4) {
        key <- if (n == 3L) "quadratic" else "cubic"
        if (n %in% failed) {
          expect_identical(model$final_model$stageB_errors[[key]],
                           "Injected Stage B failure")
          expect_null(model$final_model[[paste0(key, ".fit")]])
          for (type in c("link", "response")) {
            expect_true(all(is.na(predict(model, n = n, type = type))))
            expect_true(all(is.na(predict(model, newdata = d, n = n, type = type))))
            pred <- predict(model, newdata = d, n = "all", type = type)
            expect_identical(pred$pred_linear, predict(reference, newdata = d,
                                                       n = 2, type = type))
            expect_true(all(is.na(pred[[paste0("pred_", key)]])))
          }
        } else {
          expect_identical(predict(model, n = n), predict(reference, n = n))
          expect_identical(predict(model, newdata = d, n = n),
                           predict(reference, newdata = d, n = n))
          expect_identical(knots(model, n = n), knots(reference, n = n))
        }
      }
    }
  }
})
