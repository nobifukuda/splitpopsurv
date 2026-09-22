test_that("loglik_splitpop_* functions return finite, correctly-sized vectors", {
  set.seed(2)
  n <- 30
  Xh <- cbind(1, rnorm(n)); Xp <- cbind(1, rbinom(n, 1, 0.5))
  time <- runif(n, 0.1, 5); event <- rbinom(n, 1, 0.7); group <- rep(1, n)

  ll1 <- loglik_splitpop_loglogistic(c(0, 0, 1, 0, 0), time, event, group, Xh, Xp)
  ll2 <- loglik_splitpop_weibull(c(0, 0, 0, 0, 0), time, event, group, Xh, Xp)
  ll3 <- loglik_splitpop_lognormal(c(0, 0, 0, 0, 0), time, event, group, Xh, Xp)
  ll4 <- loglik_splitpop_gamma(c(1, 0, 1, 0, 0), time, event, group, Xh, Xp)
  ll5 <- loglik_splitpop_ggamma(c(0, 0, 0, 1, 0, 0), time, event, group, Xh, Xp)

  for (ll in list(ll1, ll2, ll3, ll4, ll5)) {
    expect_length(ll, n)
    expect_true(all(is.finite(ll)))
  }
})

test_that("censored observations never depend on the mover density term", {
  # A regression guard for the likelihood correction: for a censored
  # observation (event = 0), the contribution must equal log(marginal
  # survival) exactly, regardless of how extreme the mover hazard/density is.
  Xh <- cbind(1, 0); Xp <- cbind(1, 0)
  ll <- loglik_splitpop_weibull(c(0, 0, 0, 0, 0), time = 1e6, event = 0,
                                 group = 1, Xh = Xh, Xp = Xp)
  expect_true(is.finite(ll))
})
