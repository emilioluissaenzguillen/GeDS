test_that("univariate knot-cap warnings are explicit and selectively suppressible", {
  set.seed(617)
  d <- data.frame(x = seq(0, 1, length.out = 40), y = rnorm(40))
  for (fun in list(NGeDS, GGeDS)) {
    expect_warning(
      fit <- fun(y ~ f(x), data = d, min.intknots = 2, max.intknots = 2,
                 phi = 1, q = 1, higher_order = FALSE),
      "Internal-knot cap reached.*max.intknots = 2",
      class = "geds_knot_cap_warning"
    )
    expect_no_warning(
      quiet <- suppressWarnings(
        fun(y ~ f(x), data = d, min.intknots = 2, max.intknots = 2,
            phi = 1, q = 1, higher_order = FALSE),
        classes = "geds_knot_cap_warning"
      )
    )
    expect_equal(predict(fit, n = 2), predict(quiet, n = 2))
  }
  expect_warning(
    suppressWarnings(warning("Unrelated warning"), classes = "geds_knot_cap_warning"),
    "Unrelated warning"
  )
})
