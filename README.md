# splitpopsurv

Split-population (cure / mover–stayer) survival models in R: an accelerated
failure-time regression for event timing among "movers", combined with a
logistic regression on the probability of belonging to the immune "stayer"
population.

Five baseline timing distributions are provided — **log-logistic**,
**Weibull**, **log-normal**, **gamma**, and the **generalized gamma** that
nests the other four — following Schmidt & Witte (1989) and Yamaguchi
(1992, 1998). This package is an R translation of a set of Stata `ml`
programs, with the log-likelihood corrected to match the published model
and verified by simulation against known parameters.

📖 **Full manual** (theory, formulas, the likelihood correction, and a
complete function reference): [`docs/manual.html`](https://github.com/nobifukuda/splitpopsurv/blob/main/docs/manual.html)

## Author

**Nobutaka Fukuda**, Tohoku University — <nobutaka.fukuda@tohoku.ac.jp>

## Installation

```r
# install.packages("remotes")
remotes::install_github("nobifukuda/splitpopsurv")
```

## Usage

```r
library(splitpopsurv)

# mydata needs: time, event (0/1), group (0/1), and your covariates
fit <- fit_splitpop_weibull(
  hform = ~ x1,           # H_regression: covariates for timing
  pform = ~ x2,           # P_regression: covariates for cure probability
  data  = mydata,
  time  = "time",
  event = "event",
  group = "group",
  method = "BFGS"
)

summary(fit)   # coefficients, SEs, z-values, log-likelihood
coef(fit)      # named parameter vector
```

The other four distributions use the same signature:
`fit_splitpop_loglogistic()`, `fit_splitpop_lognormal()`,
`fit_splitpop_gamma()`, `fit_splitpop_ggamma()`.

## A note on the original Stata code

All five Stata programs this package translates compute a log-likelihood
term that turns out to be the marginal *density* where the formula requires
the marginal *hazard* (density divided by survival) — a discrepancy from
Yamaguchi's own published model. This was confirmed by fitting simulated
data with known parameters: the as-translated formula gives visibly biased
estimates (especially for the cure-probability coefficients), while the
corrected formula implemented here recovers the true parameters accurately.
See [`docs/manual.html`](https://github.com/nobifukuda/splitpopsurv/blob/main/docs/manual.html#correction) for the full
derivation and the simulation results.

A companion Stata command implementing the same corrected models is also
available; see the author's contact details below.

## License

MIT — see [LICENSE](LICENSE).

## References

Yamaguchi, K., & Ferguson, L. R. (1995). The stopping and spacing of
childbirths and their birth-history predictors: Rational-choice theory and
event-history analysis. *American Sociological Review*, 60(2), 272–298.

Yamaguchi, K. (1998). Mover-stayer models for analyzing event nonoccurrence
and event timing with time-dependent covariates: An application to an
analysis of remarriage. *Sociological Methodology*, 28(1), 327–361.

Schmidt, P., & Witte, A. D. (1989). Predicting criminal recidivism using
"split population" survival time models. *Journal of Econometrics*, 40(1),
141–159.
