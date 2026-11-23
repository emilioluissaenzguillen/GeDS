test_that("Gamma GGeDS defaults to beta 0.3 for every supported link", {
  set.seed(273)
  dat <- data.frame(x = seq(0, 1, length.out = 60))
  for (link in c("log", "inverse", "identity")) {
    family <- Gamma(link = link)
    mu <- family$linkinv(1.5 + .1 * dat$x)
    dat$y <- rgamma(nrow(dat), shape = 30, rate = 30 / mu)
    args <- list(formula = y ~ f(x), data = dat, family = family,
                 phi = .8, q = 1, max.intknots = 5, higher_order = FALSE)
    default <- suppressWarnings(do.call(GGeDS, args))
    explicit <- suppressWarnings(do.call(GGeDS, c(args, list(beta = .3))))
    overridden <- suppressWarnings(do.call(GGeDS, c(args, list(beta = .1))))

    expect_equal(default$args$beta, .3)
    expect_equal(overridden$args$beta, .1)
    expect_equal(knots(default, n = 2), knots(explicit, n = 2))
    expect_equal(predict(default, n = 2), predict(explicit, n = 2),
                 tolerance = 1e-10)
  }
})

test_that("other GGeDS beta defaults and the boosting default are unchanged", {
  set.seed(274)
  dat <- data.frame(x = seq(0, 1, length.out = 60))
  families <- list(gaussian(), poisson(), quasipoisson(), binomial(), quasibinomial())
  expected <- c(.5, .2, .2, .1, .1)
  for (i in seq_along(families)) {
    family <- families[[i]]
    mu <- family$linkinv(.5 + .1 * dat$x)
    dat$y <- switch(family$family,
      gaussian = rnorm(nrow(dat), mu, .1),
      poisson = rpois(nrow(dat), mu), quasipoisson = rpois(nrow(dat), mu),
      binomial = rbinom(nrow(dat), 1, mu), quasibinomial = rbinom(nrow(dat), 1, mu))
    fit <- suppressWarnings(GGeDS(y ~ f(x), data = dat, family = family,
      phi = .8, q = 1, max.intknots = 5, higher_order = FALSE))
    expect_equal(fit$args$beta, expected[i])
  }
  expect_equal(formals(NGeDSboost)$beta, .5)
})
