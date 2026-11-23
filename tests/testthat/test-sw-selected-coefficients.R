test_that("SW rejects shared supporting sites and tries the next cluster", {
  x <- c(.5, 2, 3.4, 3.6, 3.8)
  oldknots <- c(0, 4, 0, 4, 1, 3)
  proposed <- sort(c(oldknots, 2))
  B <- splines::splineDesign(proposed, x, ord = 2)
  expect_true(all(colSums(B > 0) > 0))
  expect_lt(stats::.lm.fit(B, rep(1, length(x)))$rank, ncol(B))

  result <- GeDS:::Knotnew(c(.1, 1, .8, .2, .1), rep(1, length(x)), x,
                          seq_along(x), oldknots, 1e-12)
  reference <- GeDS:::Knotnew_R_clean(c(.1, 1, .8, .2, .1), rep(1, length(x)), x,
                                      seq_along(x), oldknots, 1e-12)
  expect_equal(result, c(3.6, 4))
  expect_equal(c(reference$newknot, reference$cluster_index), result)
  accepted <- splines::splineDesign(sort(c(oldknots, result[1])), x, ord = 2)
  expect_equal(stats::.lm.fit(accepted, rep(1, length(x)))$rank, ncol(accepted))
})

test_that("SW applies tolerance to each matched site", {
  tol <- 1e-12
  near <- 2 + tol / 4
  x <- c(.5, 2, near, 3.4, 3.6, 3.8)
  w <- c(.1, .2, 1, .8, .2, .1)
  oldknots <- c(0, 4, 0, 4, 1, 3)
  for (margin in c(0, tol)) {
    result <- GeDS:::Knotnew(w, rep(1, length(x)), x, seq_along(x), oldknots, margin)
    reference <- GeDS:::Knotnew_R_clean(w, rep(1, length(x)), x,
                                        seq_along(x), oldknots, margin)
    expect_equal(result, if (margin == 0) c(near, 3) else c(3.4, 4))
    expect_equal(c(reference$newknot, reference$cluster_index), result)
  }
})

test_that("SW handles clamped endpoints without counting repeated sites twice", {
  oldknots <- c(0, 4, 0, 4, 1, 3)
  cases <- list(
    endpoints = list(x = c(0, .5, 2, 3.5, 4), w = c(.1, .2, 1, .2, .1), expected = c(2, 3)),
    repeated = list(x = c(.5, .5, 2, 3.4, 3.6, 3.8),
                    w = c(.1, .1, 1, .8, .2, .1), expected = c(3.6, 5))
  )
  for (case in cases) {
    result <- GeDS:::Knotnew(case$w, rep(1, length(case$x)), case$x,
                            seq_along(case$x), oldknots, 1e-12)
    reference <- GeDS:::Knotnew_R_clean(case$w, rep(1, length(case$x)), case$x,
                                        seq_along(case$x), oldknots, 1e-12)
    expect_equal(result, case$expected)
    expect_equal(c(reference$newknot, reference$cluster_index), result)
  }
})

test_that("C++ and R SW checks agree for higher spline orders", {
  x <- seq(0, 4, by = .25)
  w <- rep(.1, length(x))
  w[x == 2] <- 1
  for (n in 2:4) {
    oldknots <- c(rep(c(0, 4), n), 1, 3)
    result <- GeDS:::Knotnew(w, rep(1, length(x)), x, seq_along(x), oldknots,
                            1e-12, support_order = n)
    reference <- GeDS:::Knotnew_R_clean(w, rep(1, length(x)), x, seq_along(x), oldknots,
                                        1e-12, support_order = n)
    expect_equal(result, c(2, which(x == 2)))
    expect_equal(c(reference$newknot, reference$cluster_index), result)
  }
})

test_that("linear learner conversion uses the selected refit, not its history", {
  dat <- data.frame(x = seq(0, 1, length.out = 60))
  dat$y <- 2 + sin(7 * dat$x)
  fit <- suppressWarnings(NGeDS(y ~ f(x), data = dat, phi = .9))
  expected <- as.numeric(predict(fit, newdata = dat, n = 2))
  theta <- as.numeric(fit$linear.fit$theta)
  fit$stored[] <- NA_real_
  fit$coefficients[] <- -999
  converted <- GeDS:::predict_GeDS_linear(fit, X = dat$x)
  expect_equal(as.numeric(converted$theta), theta)
  expect_equal(as.numeric(converted$int.knt), as.numeric(fit$linear.intknots))
  expect_lt(max(abs(as.numeric(converted$Y_hat) - expected)), 1e-8)
})

test_that("positive-weight sites determine spline identifiability", {
  dat <- data.frame(x = rep(seq(0, 1, length.out = 40), each = 2))
  dat$y <- sin(9 * dat$x) + .03 * cos(seq_len(nrow(dat)))
  w <- rep(c(0, 1, .5, 2), length.out = nrow(dat))
  fit <- suppressWarnings(NGeDS(y ~ f(x), data = dat, weights = w, phi = .95))
  ik <- as.numeric(unlist(knots(fit, n = 2, options = "internal")))
  x <- sort(unique(dat$x[w > 0]))
  B <- splines::splineDesign(sort(c(ik, rep(fit$args$extr, 2))), x, ord = 2)
  expect_equal(stats::.lm.fit(B, rep(1, length(x)))$rank, ncol(B))
})

test_that("diffraction fits agree before the boosted Stage B safeguard", {
  data("BaFe2As2", package = "GeDS")
  dat <- data.frame(X = BaFe2As2$angle, Y = BaFe2As2$intensity)
  cases <- data.frame(seed = c(124L, 124L, 127L), fold = c(4L, 9L, 1L))
  removed <- matrix(0L, nrow(cases), 2L)
  for (i in seq_len(nrow(cases))) {
    set.seed(cases$seed[i])
    id <- integer(nrow(dat))
    blocks <- split(order(dat$X), ceiling(seq_len(nrow(dat)) / 10))
    for (block in blocks) id[block] <- sample(1:10, length(block))
    train <- dat[id != cases$fold[i], ]
    test <- dat[id == cases$fold[i], ]
    invisible(capture.output({
      plain <- NGeDS(Y ~ f(X), data = train, beta = .6, phi = .995, q = 3,
                     stoptype = "RD")
      boost <- NGeDSboost(Y ~ f(X), data = train, beta = .6, phi = .995, q = 3,
                          initial_learner = FALSE, shrinkage = 1, penalty = "none",
                          stoptype = "RD")
      # Isolate the selected-coefficient check from intentional Stage B deletion.
      unguarded <- local({
        local_mocked_bindings(.stageB_stabilise = function(X, Z, weights, learners,
          internal, n, ...) list(internal = internal, stability = NULL), .package = "GeDS")
        NGeDSboost(Y ~ f(X), data = train, beta = .6, phi = .995, q = 3,
                    initial_learner = FALSE, shrinkage = 1, penalty = "none",
                    stoptype = "RD")
      })
    }))
    for (n in 2:4) {
      kp <- as.numeric(unlist(knots(plain, n = n, options = "internal")))
      kb <- as.numeric(unlist(knots(boost, n = n, options = "internal")))
      ku <- as.numeric(unlist(knots(unguarded, n = n, options = "internal")))
      expect_length(ku, length(kp))
      expect_equal(ku, kp, tolerance = 1e-9)
      expect_lt(max(abs(as.numeric(predict(plain, newdata = test, n = n)) -
                         as.numeric(predict(unguarded, newdata = test, n = n)))), 1e-6)
      if (n == 2L) {
        expect_equal(kb, ku, tolerance = 1e-9)
        expect_equal(predict(boost, newdata = test, n = n),
                     predict(unguarded, newdata = test, n = n), tolerance = 1e-6)
      } else {
        key <- c("quadratic.fit", "cubic.fit")[n - 2L]
        diagnostic <- boost$final_model[[key]]$stability
        expect_null(plain[[key]]$stability)
        expect_equal(as.numeric(unlist(diagnostic$original_knots)), ku, tolerance = 1e-9)
        expect_equal(as.numeric(unlist(diagnostic$knots)), kb, tolerance = 1e-9)
        expect_equal(length(ku) - length(kb), diagnostic$removed)
        expect_true(all(kb %in% ku))
        expect_true(diagnostic$target_met)
        expect_lte(diagnostic$max_ratio, 1 + 1e-10)
        expect_true(all(is.finite(predict(boost, newdata = test, n = n))))
        removed[i, n - 2L] <- diagnostic$removed
        if (diagnostic$removed == 0L) {
          expect_equal(predict(boost, newdata = test, n = n),
                       predict(unguarded, newdata = test, n = n), tolerance = 1e-6)
        }
      }
      for (ik in list(kp, ku, kb)) {
        B <- splines::splineDesign(sort(c(ik, rep(range(train$X), n))), train$X, ord = n)
        expect_equal(stats::.lm.fit(B, train$Y)$rank, ncol(B))
      }
    }
  }
  expect_equal(removed[, 1L], c(0L, 0L, 1L))
  expect_equal(removed[, 2L], c(0L, 0L, 0L))
})
