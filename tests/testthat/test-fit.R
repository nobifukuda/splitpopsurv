make_test_data <- function(n = 150, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
  t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
  censor_time <- rexp(n, rate = 0.2)
  time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
  event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
  data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
}

test_that("all five models converge on well-behaved simulated data", {
  dat <- make_test_data()

  fit_ll <- fit_splitpop_loglogistic(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  fit_w  <- fit_splitpop_weibull(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  fit_ln <- fit_splitpop_lognormal(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  fit_g  <- fit_splitpop_gamma(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  fit_gg <- fit_splitpop_ggamma(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")

  for (fit in list(fit_ll, fit_w, fit_ln, fit_g, fit_gg)) {
    expect_s3_class(fit, "maxLik")
    expect_equal(fit$code, 0)
    expect_true(is.finite(fit$maximum))
    expect_false(any(is.na(coef(fit))))
  }
})

test_that("parameter names follow the H_regression/P_regression convention", {
  dat <- make_test_data()
  fit <- fit_splitpop_weibull(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  nm <- names(coef(fit))
  expect_true(any(grepl("^H_regression:", nm)))
  expect_true(any(grepl("^P_regression:", nm)))
  expect_true("ln_sigma" %in% nm)
})

test_that("the corrected likelihood recovers known Weibull parameters (regression test)", {
  skip_on_cran()  # larger n than needed for CRAN's routine test run; keeps CI fast elsewhere
  set.seed(42)
  n <- 5000
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  beta0 <- 0.5; beta1 <- -0.4; shape <- 1.3
  alpha0 <- 0.2; alpha1 <- 0.8
  theta1 <- beta0 + beta1 * x1
  cured <- rbinom(n, 1, plogis(alpha0 + alpha1 * x2))
  t_latent <- (-log(runif(n)) / exp(theta1))^(1 / shape)
  censor_time <- rexp(n, rate = 0.05)
  time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-4)
  event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
  dat <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)

  fit <- fit_splitpop_weibull(~x1, ~x2, dat, "time", "event", "group", method = "BFGS")
  est <- coef(fit)

  # Absolute-difference checks: expect_equal(..., tolerance=) is a RELATIVE
  # tolerance under testthat 3e, which is too tight for small-magnitude true
  # values like alpha0 = 0.2 (15% relative = only 0.03 absolute). We want an
  # absolute bound here, since simulation noise doesn't scale with |true|.
  expect_lt(abs(unname(est["H_regression:(Intercept)"]) - beta0), 0.1)
  expect_lt(abs(unname(est["H_regression:x1"]) - beta1), 0.1)
  expect_lt(abs(unname(est["P_regression:(Intercept)"]) - alpha0), 0.1)
  expect_lt(abs(unname(est["P_regression:x2"]) - alpha1), 0.15)
})
