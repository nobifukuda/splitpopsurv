# ============================================================================
# splitpopsurv: Split-Population (Cure / Mover-Stayer) Survival Models
#
# Author:  Nobutaka Fukuda <nobutaka.fukuda@tohoku.ac.jp>
# License: MIT (see LICENSE)
#
# An R translation of five Stata `ml` programs (SphLog, SphWieb, SphNom,
# SphGam, SphGGam) for Log-logistic, Weibull, Log-normal, Gamma and
# Generalized Gamma split-population survival models.
#
# Statistical background: Schmidt & Witte (1989) "split population" survival
# models and Yamaguchi (1992, 1998) mover-stayer / accelerated-failure-time
# models with a logit regression on the surviving ("stayer"/"cured") fraction.
# See Yamaguchi & Ferguson (1995, American Sociological Review 60(2)) eqs.
# (3)-(6) and Yamaguchi (1998, Sociological Methodology 28) eq. (1)-(2) for
# the underlying model. Full references and derivations: docs/manual.html.
#
# Stata -> R correspondence used throughout:
#   $ML_y1            -> time    (duration/survival time)
#   $ML_y2            -> event   (1 = failure observed, 0 = censored)
#   $ML_y3            -> group   (0/1 indicator that flips the cure-prob eq,
#                                  exactly as in the original Stata code)
#   invlogit(x)        -> plogis(x)
#   norm(x)            -> pnorm(x)
#   gammap(a,x)         -> pgamma(x, shape = a)      (regularized lower incomplete gamma)
#   lngamma(x)          -> lgamma(x)
#   _pi                 -> pi
#   ml model lf/d0 ...   -> maxLik::maxLik(..., method = "NR"/"BHHH")
#
# Each fit_splitpop_*() function mirrors one Stata `ml model` block: it takes
# a formula for the hazard/location regression (H_regression), a formula for
# the cure-probability regression (P_regression), plus the time/event/group
# variables, builds the design matrices, and maximizes the log-likelihood.
#
# Ancillary shape parameters (Stata's unnamed "()" or "ln_sigma:"/"kappa:"
# equations, i.e. intercept-only equations) are represented as single free
# scalars in the parameter vector, exactly as in the original.
#
# ----------------------------------------------------------------------------
# CORRECTION relative to the original Stata code (see docs/manual.html for
# the full derivation and a simulation-based verification):
#
# For a split-population model with mover survival S_m(t)/density f_m(t) and
# cure probability p, the marginal survival and density are
#     S(t) = (1-p) S_m(t) + p            f(t) = (1-p) f_m(t)
# and the standard hazard-based log-likelihood contribution for person i is
#     delta_i * log( f(t_i) )  +  (1 - delta_i) * log( S(t_i) )
# equivalently  delta_i * log( h(t_i) ) + log( S(t_i) )  where h(t) = f(t)/S(t)
# is the TRUE marginal hazard (Yamaguchi 1998, eq. 2).
#
# The original Stata code computes hz2/mhz = (1-p)*hz_m(t)*S_m(t), which
# equals the marginal DENSITY f(t), not the marginal hazard h(t) = f(t)/S(t)
# (it never divides by S(t)), then evaluates lnf = $ML_y2*ln(hz2) + ln(sv2)
# for every observation -- double-counting log(S(t_i)) for every observed
# event. Fitting simulated data with known true parameters confirmed this:
# the as-translated formula gives visibly biased estimates (especially for
# the cure-probability/P_regression coefficients), while the corrected
# formula below recovers the true parameters accurately.
#
# This package implements the CORRECTED formula:
#     loglik_i = delta_i * log(f(t_i))  +  (1 - delta_i) * log(S(t_i))
# written as ifelse(event==1, log(marginal density), log(marginal survival)),
# which is both statistically correct and numerically robust (it avoids ever
# forming hazard = density / survival, so there is no 0/0 or Inf*0 risk when
# survival underflows to exactly 0 in floating point at extreme time values).
# ============================================================================

## ---------------------------------------------------------------------
## Helpers (internal, not exported)
## ---------------------------------------------------------------------

# Build a design matrix from a one-sided formula (~ x1 + x2), i.e. the
# right-hand side used in Stata's (eqname: depvar = indepvars) equations.
.design_matrix <- function(formula, data) {
  formula <- update(formula, NULL ~ .)
  mf <- model.frame(formula, data, na.action = na.pass)
  model.matrix(formula, mf)
}

# p = invlogit(theta_p) if group==1, 1-invlogit(theta_p) if group==0
# (kept exactly as in the Stata code: $ML_y3 flips which side of the
# logistic curve is used for the cure probability)
.cure_prob <- function(theta_p, group) {
  p <- plogis(theta_p)
  p[group == 0] <- 1 - p[group == 0]
  p
}

.name_params <- function(k_h, extra_names, k_p, Xh, Xp) {
  c(paste0("H_regression:", colnames(Xh)),
    extra_names,
    paste0("P_regression:", colnames(Xp)))
}

# Correct split-population log-likelihood contribution:
#   delta * log(marginal density) + (1 - delta) * log(marginal survival)
# mdens = (1-p)*f_m(t)  [marginal/mixture density]
# msurv = (1-p)*S_m(t)+p [marginal/mixture survival]
.split_loglik <- function(event, mdens, msurv) {
  ifelse(event == 1, log(mdens), log(msurv))
}

# Safe starting values for the H_regression part of the Gamma / generalized-
# gamma models, where theta1 = Xh %*% b_h is used directly (no exp()
# transform) and must stay positive: start with intercept = 1, all slopes
# = 0, so theta1 is a positive constant regardless of covariate values.
.safe_start_h <- function(Xh, intercept_value = 1) {
  b_h0 <- rep(0, ncol(Xh))
  icol <- which(colnames(Xh) == "(Intercept)")
  if (length(icol) == 1) b_h0[icol] <- intercept_value else b_h0[] <- intercept_value
  b_h0
}

## ---------------------------------------------------------------------
## 1. Log-logistic distribution  (Stata: SphLog)
## ---------------------------------------------------------------------

#' Log-likelihood: split-population log-logistic model
#'
#' Internal log-likelihood used by [fit_splitpop_loglogistic()]. Exposed so
#' it can be evaluated by hand at candidate starting values (the R analogue
#' of Stata's `ml check`).
#'
#' @param par Numeric parameter vector: `H_regression` coefficients, the
#'   ancillary scalar `theta2`, then `P_regression` coefficients.
#' @param time,event,group Numeric vectors: duration, 0/1 event indicator,
#'   0/1 group indicator (see [fit_splitpop_loglogistic()]).
#' @param Xh,Xp Design matrices for the H_regression and P_regression parts.
#' @return A numeric vector of per-observation log-likelihood contributions.
#' @examples
#' Xh <- cbind(1, rnorm(20)); Xp <- cbind(1, rbinom(20, 1, 0.5))
#' par0 <- c(0, 0, 1, 0, 0)
#' loglik_splitpop_loglogistic(par0, time = runif(20, 0.1, 5),
#'                              event = rbinom(20, 1, 0.7),
#'                              group = rep(1, 20), Xh = Xh, Xp = Xp)
#' @export
loglik_splitpop_loglogistic <- function(par, time, event, group, Xh, Xp) {
  k_h <- ncol(Xh); k_p <- ncol(Xp)
  b_h    <- par[1:k_h]
  theta2 <- par[k_h + 1]                       # ancillary scalar ("()")
  b_p    <- par[(k_h + 2):(k_h + 1 + k_p)]

  theta1 <- as.vector(Xh %*% b_h)
  theta3 <- as.vector(Xp %*% b_p)

  p     <- .cure_prob(theta3, group)
  lam   <- exp(-theta1)
  sigma <- 1 / theta2
  sv    <- 1 / (1 + (lam * time)^sigma)                 # S_m(t)
  hz    <- (lam * sigma * (lam * time)^(sigma - 1)) / (1 + (lam * time)^sigma)  # h_m(t)
  msv   <- (1 - p) * sv + p                              # S(t)
  mdens <- (1 - p) * hz * sv                             # f(t) = (1-p) h_m(t) S_m(t)

  .split_loglik(event, mdens, msv)
}

#' Fit a split-population log-logistic survival model
#'
#' Maximum-likelihood estimation of a split-population (cure) survival model
#' with a log-logistic baseline for the "mover" (susceptible) population and
#' a logistic regression on the probability of being a "stayer" (immune to
#' the event). R translation of the Stata program `SphLog`; see
#' `docs/manual.html` for the model, formulas, and a note on a likelihood
#' correction relative to the original Stata code.
#'
#' @param hform One-sided formula (e.g. `~ x1 + x2`) for the H_regression
#'   equation: covariates for the log-logistic timing distribution.
#' @param pform One-sided formula for the P_regression equation: covariates
#'   for the logit on the stayer ("cure") probability.
#' @param data A `data.frame` containing `time`, `event`, `group`, and every
#'   variable referenced in `hform`/`pform`.
#' @param time Character: name of the duration column ($ML_y1 in Stata).
#' @param event Character: name of the 0/1 event column, 1 = failure
#'   observed, 0 = censored ($ML_y2 in Stata).
#' @param group Character: name of a 0/1 column that flips which side of the
#'   cure-probability logit is used ($ML_y3 in Stata) -- see the manual.
#'   If there is no such distinction in your data, pass a column of all 1s.
#' @param start Optional numeric starting vector; if `NULL`, a default is
#'   constructed automatically.
#' @param method Optimizer passed to [maxLik::maxLik()]; default `"NR"`
#'   (Newton-Raphson). Try `"BFGS"` if that fails to converge.
#' @return A `maxLik` object; use `summary()`, `coef()`, `logLik()`, `vcov()`.
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
#' t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
#' censor_time <- rexp(n, rate = 0.2)
#' time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
#' event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
#' mydata <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
#'
#' fit <- fit_splitpop_loglogistic(~x1, ~x2, mydata,
#'                                  "time", "event", "group", method = "BFGS")
#' summary(fit)
#' @export
fit_splitpop_loglogistic <- function(hform, pform, data, time, event, group,
                                      start = NULL, method = "NR") {
  Xh <- .design_matrix(hform, data)
  Xp <- .design_matrix(pform, data)
  y1 <- data[[time]]; y2 <- data[[event]]; y3 <- data[[group]]

  if (is.null(start)) start <- c(rep(0, ncol(Xh)), 1, rep(0, ncol(Xp)))
  names(start) <- .name_params(ncol(Xh), "theta2", ncol(Xp), Xh, Xp)

  maxLik::maxLik(loglik_splitpop_loglogistic, start = start, method = method,
                  time = y1, event = y2, group = y3, Xh = Xh, Xp = Xp)
}

## ---------------------------------------------------------------------
## 2. Weibull distribution  (Stata: SphWieb)
## ---------------------------------------------------------------------

#' Log-likelihood: split-population Weibull model
#' @inheritParams loglik_splitpop_loglogistic
#' @return A numeric vector of per-observation log-likelihood contributions.
#' @examples
#' Xh <- cbind(1, rnorm(20)); Xp <- cbind(1, rbinom(20, 1, 0.5))
#' par0 <- c(0, 0, 0, 0, 0)
#' loglik_splitpop_weibull(par0, time = runif(20, 0.1, 5),
#'                          event = rbinom(20, 1, 0.7),
#'                          group = rep(1, 20), Xh = Xh, Xp = Xp)
#' @export
loglik_splitpop_weibull <- function(par, time, event, group, Xh, Xp) {
  k_h <- ncol(Xh); k_p <- ncol(Xp)
  b_h      <- par[1:k_h]
  theta2   <- par[k_h + 1]                     # "ln_sigma:" scalar
  b_p      <- par[(k_h + 2):(k_h + 1 + k_p)]

  theta1 <- as.vector(Xh %*% b_h)
  theta3 <- as.vector(Xp %*% b_p)

  p     <- .cure_prob(theta3, group)
  sigma <- exp(-theta1)
  lam   <- exp(theta2)
  sv    <- exp(-exp(theta1) * (time^lam))                # S_m(t)
  hz    <- lam * (time^(lam - 1)) * exp(theta1)           # h_m(t)
  msv   <- (1 - p) * sv + p                               # S(t)
  mdens <- (1 - p) * hz * sv                              # f(t)

  .split_loglik(event, mdens, msv)
}

#' Fit a split-population Weibull survival model
#'
#' Maximum-likelihood estimation of a split-population (cure) survival model
#' with a Weibull baseline for the "mover" population. R translation of the
#' Stata program `SphWieb`; see `docs/manual.html` for the model and a note
#' on a likelihood correction relative to the original Stata code.
#'
#' @inheritParams fit_splitpop_loglogistic
#' @return A `maxLik` object; use `summary()`, `coef()`, `logLik()`, `vcov()`.
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
#' t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
#' censor_time <- rexp(n, rate = 0.2)
#' time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
#' event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
#' mydata <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
#'
#' fit <- fit_splitpop_weibull(~x1, ~x2, mydata,
#'                              "time", "event", "group", method = "BFGS")
#' summary(fit)
#' @export
fit_splitpop_weibull <- function(hform, pform, data, time, event, group,
                                  start = NULL, method = "NR") {
  Xh <- .design_matrix(hform, data)
  Xp <- .design_matrix(pform, data)
  y1 <- data[[time]]; y2 <- data[[event]]; y3 <- data[[group]]

  if (is.null(start)) start <- c(rep(0, ncol(Xh)), 0, rep(0, ncol(Xp)))
  names(start) <- .name_params(ncol(Xh), "ln_sigma", ncol(Xp), Xh, Xp)

  maxLik::maxLik(loglik_splitpop_weibull, start = start, method = method,
                  time = y1, event = y2, group = y3, Xh = Xh, Xp = Xp)
}

## ---------------------------------------------------------------------
## 3. Log-normal distribution  (Stata: SphNom)
## ---------------------------------------------------------------------

#' Log-likelihood: split-population log-normal model
#' @inheritParams loglik_splitpop_loglogistic
#' @return A numeric vector of per-observation log-likelihood contributions.
#' @examples
#' Xh <- cbind(1, rnorm(20)); Xp <- cbind(1, rbinom(20, 1, 0.5))
#' par0 <- c(0, 0, 0, 0, 0)
#' loglik_splitpop_lognormal(par0, time = runif(20, 0.1, 5),
#'                            event = rbinom(20, 1, 0.7),
#'                            group = rep(1, 20), Xh = Xh, Xp = Xp)
#' @export
loglik_splitpop_lognormal <- function(par, time, event, group, Xh, Xp) {
  k_h <- ncol(Xh); k_p <- ncol(Xp)
  b_h    <- par[1:k_h]
  theta2 <- par[k_h + 1]                       # "ln_sigma:" scalar
  b_p    <- par[(k_h + 2):(k_h + 1 + k_p)]

  theta1 <- as.vector(Xh %*% b_h)
  theta3 <- as.vector(Xp %*% b_p)

  p     <- .cure_prob(theta3, group)
  sigma <- exp(theta2)
  lam   <- log(time) - theta1
  sv    <- 1 - pnorm(lam / sigma)                                  # S_m(t)
  pdf_m <- exp(-(lam^2) / (2 * sigma^2)) / (sqrt(2 * pi) * sigma * time)  # f_m(t), direct (no /sv)
  msv   <- (1 - p) * sv + p                                        # S(t)
  mdens <- (1 - p) * pdf_m                                         # f(t)

  .split_loglik(event, mdens, msv)
}

#' Fit a split-population log-normal survival model
#'
#' Maximum-likelihood estimation of a split-population (cure) survival model
#' with a log-normal baseline for the "mover" population. R translation of
#' the Stata program `SphNom`; see `docs/manual.html` for the model and a
#' note on a likelihood correction relative to the original Stata code.
#'
#' @inheritParams fit_splitpop_loglogistic
#' @return A `maxLik` object; use `summary()`, `coef()`, `logLik()`, `vcov()`.
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
#' t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
#' censor_time <- rexp(n, rate = 0.2)
#' time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
#' event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
#' mydata <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
#'
#' fit <- fit_splitpop_lognormal(~x1, ~x2, mydata,
#'                                "time", "event", "group", method = "BFGS")
#' summary(fit)
#' @export
fit_splitpop_lognormal <- function(hform, pform, data, time, event, group,
                                    start = NULL, method = "NR") {
  Xh <- .design_matrix(hform, data)
  Xp <- .design_matrix(pform, data)
  y1 <- data[[time]]; y2 <- data[[event]]; y3 <- data[[group]]

  if (is.null(start)) start <- c(rep(0, ncol(Xh)), 0, rep(0, ncol(Xp)))
  names(start) <- .name_params(ncol(Xh), "ln_sigma", ncol(Xp), Xh, Xp)

  maxLik::maxLik(loglik_splitpop_lognormal, start = start, method = method,
                  time = y1, event = y2, group = y3, Xh = Xh, Xp = Xp)
}

## ---------------------------------------------------------------------
## 4. Gamma distribution  (Stata: SphGam)
## ---------------------------------------------------------------------

#' Log-likelihood: split-population gamma model
#' @inheritParams loglik_splitpop_loglogistic
#' @return A numeric vector of per-observation log-likelihood contributions.
#' @examples
#' Xh <- cbind(1, rnorm(20)); Xp <- cbind(1, rbinom(20, 1, 0.5))
#' par0 <- c(1, 0, 1, 0, 0)  # positive H_regression intercept: theta1 must stay > 0
#' loglik_splitpop_gamma(par0, time = runif(20, 0.1, 5),
#'                        event = rbinom(20, 1, 0.7),
#'                        group = rep(1, 20), Xh = Xh, Xp = Xp)
#' @export
loglik_splitpop_gamma <- function(par, time, event, group, Xh, Xp) {
  k_h <- ncol(Xh); k_p <- ncol(Xp)
  b_h    <- par[1:k_h]
  theta2 <- par[k_h + 1]                       # "kappa:" scalar
  b_p    <- par[(k_h + 2):(k_h + 1 + k_p)]

  theta1 <- as.vector(Xh %*% b_h)
  theta3 <- as.vector(Xp %*% b_p)

  p     <- .cure_prob(theta3, group)
  l     <- theta1 * time
  k     <- theta2
  cdf   <- 1 - pgamma(l, shape = k)            # gammap(k, l)
  gam   <- exp(lgamma(k))
  pdf_m <- (theta1 * (l^(k - 1)) * exp(-l)) / gam    # f_m(t), direct (no /sv)
  sv    <- 1 - cdf                                    # S_m(t)
  msv   <- (1 - p) * sv + p                           # S(t)
  mdens <- (1 - p) * pdf_m                            # f(t)

  .split_loglik(event, mdens, msv)
}

#' Fit a split-population gamma survival model
#'
#' Maximum-likelihood estimation of a split-population (cure) survival model
#' with a gamma baseline for the "mover" population. R translation of the
#' Stata program `SphGam`; see `docs/manual.html` for the model and a note
#' on a likelihood correction relative to the original Stata code.
#'
#' Note: unlike the other four models, the H_regression linear index here is
#' used directly as a rate (no `exp()` transform), so it must stay positive;
#' the default `start` sets a positive intercept and zero slopes to keep it
#' safe regardless of covariate values.
#'
#' @inheritParams fit_splitpop_loglogistic
#' @return A `maxLik` object; use `summary()`, `coef()`, `logLik()`, `vcov()`.
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
#' t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
#' censor_time <- rexp(n, rate = 0.2)
#' time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
#' event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
#' mydata <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
#'
#' fit <- fit_splitpop_gamma(~x1, ~x2, mydata,
#'                            "time", "event", "group", method = "BFGS")
#' summary(fit)
#' @export
fit_splitpop_gamma <- function(hform, pform, data, time, event, group,
                                start = NULL, method = "NR") {
  Xh <- .design_matrix(hform, data)
  Xp <- .design_matrix(pform, data)
  y1 <- data[[time]]; y2 <- data[[event]]; y3 <- data[[group]]

  if (is.null(start)) start <- c(.safe_start_h(Xh), 1, rep(0, ncol(Xp)))
  names(start) <- .name_params(ncol(Xh), "kappa", ncol(Xp), Xh, Xp)

  maxLik::maxLik(loglik_splitpop_gamma, start = start, method = method,
                  time = y1, event = y2, group = y3, Xh = Xh, Xp = Xp)
}

## ---------------------------------------------------------------------
## 5. Generalized Gamma distribution  (Stata: SphGGam / KYmodel5, d0 method)
## ---------------------------------------------------------------------

#' Log-likelihood: split-population generalized gamma model
#'
#' @param par Numeric parameter vector: `H_regression` coefficients,
#'   `ln_sigma`, `kappa`, then `P_regression` coefficients.
#' @inheritParams loglik_splitpop_loglogistic
#' @return A numeric vector of per-observation log-likelihood contributions.
#' @examples
#' Xh <- cbind(1, rnorm(20)); Xp <- cbind(1, rbinom(20, 1, 0.5))
#' par0 <- c(0, 0, 0, 1, 0, 0)
#' loglik_splitpop_ggamma(par0, time = runif(20, 0.1, 5),
#'                         event = rbinom(20, 1, 0.7),
#'                         group = rep(1, 20), Xh = Xh, Xp = Xp)
#' @export
loglik_splitpop_ggamma <- function(par, time, event, group, Xh, Xp) {
  k_h <- ncol(Xh); k_p <- ncol(Xp)
  b_h      <- par[1:k_h]
  theta2   <- par[k_h + 1]                     # "ln_sigma:" scalar
  theta3s  <- par[k_h + 2]                     # "kappa:" scalar
  b_p      <- par[(k_h + 3):(k_h + 2 + k_p)]

  theta1 <- as.vector(Xh %*% b_h)
  theta4 <- as.vector(Xp %*% b_p)

  p <- .cure_prob(theta4, group)

  k <- theta3s
  if (abs(k) < 0.01) k <- sign(k) * 0.01       # replace k=sign(k)*.01 if abs(k)<0.01

  s <- exp(theta2)
  l <- (abs(k))^(-2)
  z <- sign(k) * (log(time) - theta1) / s
  u <- l * exp(abs(k) * z)

  if (abs(k) < 0.01) {
    cdf <- pnorm(z)
  } else if (k >= 0.01) {
    cdf <- pgamma(u, shape = l)                # gammap(l, u)
  } else {
    cdf <- 1 - pgamma(u, shape = l)
  }

  gam <- exp(lgamma(l))
  if (abs(k) < 0.01) {
    pdf_m <- exp(-(z^2) / 2) / (s * time * sqrt(2 * pi))
  } else {
    pdf_m <- (l^l * exp(z * sqrt(l) - u)) / (s * time * sqrt(l) * gam)   # f_m(t), direct (no /sv)
  }

  sv    <- 1 - cdf                              # S_m(t)
  msv   <- (1 - p) * sv + p                     # S(t)
  mdens <- (1 - p) * pdf_m                      # f(t)

  .split_loglik(event, mdens, msv)
}

#' Fit a split-population generalized gamma survival model
#'
#' Maximum-likelihood estimation of a split-population (cure) survival model
#' with a generalized gamma baseline (Prentice 1974 / Yamaguchi & Ferguson
#' 1995, note 10) for the "mover" population -- the family that nests the
#' Weibull (kappa=1), log-normal (kappa=0), and gamma (sigma=1) models above.
#' R translation of the Stata program `SphGGam` (`d0` method); see
#' `docs/manual.html` for the model and a note on a likelihood correction
#' relative to the original Stata code.
#'
#' Note: `|kappa| < 0.01` switches to a log-normal-equivalent branch, a
#' genuine kink in the likelihood surface (preserved from the Stata guard).
#' If `summary()` reports a non-finite standard error, the optimizer likely
#' landed near that threshold -- try a different starting `kappa`.
#'
#' @inheritParams fit_splitpop_loglogistic
#' @return A `maxLik` object; use `summary()`, `coef()`, `logLik()`, `vcov()`.
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' cured <- rbinom(n, 1, plogis(0.2 + 0.5 * x2))
#' t_latent <- (-log(runif(n)) / exp(0.5 - 0.4 * x1))^(1 / 1.3)
#' censor_time <- rexp(n, rate = 0.2)
#' time  <- pmax(ifelse(cured == 1, censor_time, pmin(t_latent, censor_time)), 1e-3)
#' event <- ifelse(cured == 1, 0, as.numeric(t_latent <= censor_time))
#' mydata <- data.frame(time = time, event = event, x1 = x1, x2 = x2, group = 1)
#'
#' fit <- fit_splitpop_ggamma(~x1, ~x2, mydata,
#'                             "time", "event", "group", method = "BFGS")
#' summary(fit)
#' @export
fit_splitpop_ggamma <- function(hform, pform, data, time, event, group,
                                 start = NULL, method = "NR") {
  Xh <- .design_matrix(hform, data)
  Xp <- .design_matrix(pform, data)
  y1 <- data[[time]]; y2 <- data[[event]]; y3 <- data[[group]]

  if (is.null(start)) start <- c(rep(0, ncol(Xh)), 0, 1, rep(0, ncol(Xp)))
  names(start) <- c(paste0("H_regression:", colnames(Xh)),
                     "ln_sigma", "kappa",
                     paste0("P_regression:", colnames(Xp)))

  maxLik::maxLik(loglik_splitpop_ggamma, start = start, method = method,
                  time = y1, event = y2, group = y3, Xh = Xh, Xp = Xp)
}

## ---------------------------------------------------------------------
## Notes on the pieces that don't carry over from Stata:
##
##  - `ml check` / `ml search`: maxLik has no direct equivalent. Check your
##    log-likelihood at the starting values manually (call the loglik_*
##    function once) and/or try several `start` vectors / optimizers
##    (method = "NR", "BFGS", "BHHH") if convergence fails.
##  - `ml graph`: plot fitted survival/hazard curves yourself, e.g. with
##    ggplot2, after extracting coef(fit) and plugging back into the
##    sv/hz formulas above.
##  - Use summary(fit) after fitting for coefficients, standard errors,
##    z-values and a log-likelihood value, analogous to `ml maximize` output.
## ---------------------------------------------------------------------
