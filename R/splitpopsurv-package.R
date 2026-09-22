#' splitpopsurv: Split-Population (Cure / Mover-Stayer) Survival Models
#'
#' Maximum-likelihood estimation of split-population survival models, which
#' combine an accelerated failure-time regression for event timing among
#' "movers" with a logistic regression on the probability of belonging to
#' the immune "stayer" population. See Schmidt & Witte (1989) and Yamaguchi
#' (1992, 1998) for the underlying theory, and the package README for a full
#' manual including a likelihood correction relative to the Stata `ml`
#' programs this package translates.
#'
#' @section Main functions:
#' [fit_splitpop_loglogistic()], [fit_splitpop_weibull()],
#' [fit_splitpop_lognormal()], [fit_splitpop_gamma()],
#' [fit_splitpop_ggamma()].
#'
#' @keywords internal
#' @importFrom stats model.frame model.matrix na.pass pgamma plogis pnorm update
"_PACKAGE"
