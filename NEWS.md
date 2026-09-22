# splitpopsurv 0.1.0

* Initial release.
* R translation of five Stata `ml` split-population survival programs
  (log-logistic, Weibull, log-normal, gamma, generalized gamma).
* Corrected the log-likelihood relative to the original Stata source: the
  Stata code used the marginal density where the standard hazard-based
  formula requires the marginal hazard (density divided by survival),
  which double-counted the survival term for every observed event. The
  correction was verified by simulation against known parameters (see
  `docs/manual.html`).
