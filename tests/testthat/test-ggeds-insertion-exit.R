test_that("invalid generalized knot candidates never reach coefficient insertion", {
  d <- data.frame(x = seq(0, 1, length.out = 50))
  d$z <- sin(seq_len(nrow(d)))
  d$off <- .1 + .1 * d$x
  d$y <- 1 / (2 + .4 * d$x + .1 * d$z + d$off + .2 * sin(5 * d$x))

  for (candidate in c(NA_real_, NaN, Inf, -Inf, 0, 1, -.1, 1.1)) {
    local({
      calls <- 0L
      insertions <- 0L
      original <- GeDS:::newknot.guess
      local_mocked_bindings(
        Knotnew = function(...) {
          calls <<- calls + 1L
          c(if (calls == 1L) .5 else candidate, 1)
        },
        newknot.guess = function(...) {
          insertions <<- insertions + 1L
          original(...)
        }, .package = "GeDS"
      )
      fit <- suppressWarnings(GGeDS(
        y ~ f(x) + z + offset(off), data = d, family = Gamma(),
        min.intknots = 4L, max.intknots = 6L, q = 2L, stoptype = "RD"
      ))
      expect_equal(calls, 2L)
      expect_equal(insertions, 1L)
      expect_equal(fit$iters, 2L)
      expect_equal(fit$Nintknots, 1L)
      expect_equal(fit$linear.intknots, .5)
      reference <- glm(y ~ splines::bs(x, knots = .5, degree = 1) + z + offset(off),
                       data = d, family = Gamma())
      expect_equal(predict(fit, n = 2), unname(predict(reference, type = "response")),
                   tolerance = 1e-7)
      for (n in 2:4) {
        expect_true(all(is.finite(predict(fit, n = n))))
        expect_equal(predict(fit, n = n), predict(fit, newdata = d, n = n), tolerance = 1e-8)
      }
    })
  }
})

test_that("generalized insertion exhaustion can retain the initial fit", {
  d <- data.frame(x = seq(0, 1, length.out = 40))
  d$y <- exp(1 + d$x + .2 * sin(4 * d$x))
  local_mocked_bindings(
    Knotnew = function(...) c(NA_real_, NA_real_),
    newknot.guess = function(...) stop("Invalid candidate reached coefficient insertion"),
    .package = "GeDS"
  )
  fit <- suppressWarnings(GGeDS(y ~ f(x), data = d, family = Gamma(link = "log"),
                                min.intknots = 3L, max.intknots = 5L, q = 2L))
  reference <- glm(y ~ x, data = d, family = Gamma(link = "log"))
  expect_equal(fit$iters, 1L)
  expect_equal(fit$Nintknots, 0L)
  expect_length(fit$linear.intknots, 0L)
  expect_equal(predict(fit, n = 2), unname(predict(reference, type = "response")), tolerance = 1e-7)
})

test_that("ordinary generalized stopping retains its q-step model selection", {
  d <- data.frame(x = seq(0, 1, length.out = 50))
  d$y <- exp(1 + d$x + .2 * sin(4 * d$x))
  original <- GeDS:::stopping_rule
  local_mocked_bindings(
    Knotnew = function(...) c(.5, 1),
    stopping_rule = function(...) {
      input <- list(...)
      result <- original(...)
      result$should_break <- input$j == 2L
      result
    }, .package = "GeDS"
  )
  fit <- suppressWarnings(GGeDS(y ~ f(x), data = d, family = Gamma(link = "log"),
                                max.intknots = 5L, q = 1L, stoptype = "SR"))
  expect_equal(fit$iters, 2L)
  expect_equal(fit$Nintknots, 0L)
  expect_length(fit$linear.intknots, 0L)
})
