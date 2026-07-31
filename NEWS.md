# blueprint 0.0.0.9000

## Design grammar

* Added a modular design grammar: all new design classes share the
  `trial_design` contract, and data-generating assumptions live in `scenario`
  objects created with `assume_response()` and `assume_survival()`, or derived
  from a design with `assume_null()` and `assume_alternative()`.
* Added `design_simon()` to construct Simon (1989) two-stage designs from
  known parameters, and `search_simon_designs()` plus `pick_design()` to find
  optimal and minimax designs meeting target error rates. Exact operating
  characteristics are computed in closed form with no dependency on clinfun.
* Added generics `evaluate()` (operating characteristics, exact or simulated),
  `verify()` (exact versus simulated operating characteristics side by side),
  `draft_protocol_text()` (protocol-ready sample size paragraph), and
  broom-style `tidy()`/`glance()` methods.
* `simulate_trial()` is now a true S3 generic. Legacy design classes
  (`two_stage_design`, `single_stage_design`, `basket_trial_design`) keep
  their existing behavior through delegating methods.

## Initial Development

### New Features

* Added `two_stage_single_arm_tte()` function for calculating optimal two-stage single-arm phase II trial designs with time-to-event endpoints
  - Based on Wu et al. (2020) methodology
  - Supports Weibull, log-normal, gamma, log-logistic, and logspline distributions
  - Uses one-sample log-rank test with exact variance estimates
  - Allows for restricted follow-up periods
  - Enables early stopping for futility at interim analysis

* Added S3 methods for `two_stage_design` objects:
  - `print.two_stage_design()`: Display design results
  - `summary.two_stage_design()`: Summarize design characteristics
  - `print.summary.two_stage_design()`: Print summary output

### Documentation

* Added comprehensive vignette "Two-Stage Single-Arm Designs for Time-to-Event Endpoints"
* Added example scripts in `inst/examples/`
* Created extensive function documentation with examples

### Testing

* Added testthat unit tests for `two_stage_single_arm_tte()`
* Tests cover multiple distributions and input validation

## References

Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Optimal two-stage phase II survival trial design. *Pharmaceutical Statistics*, 19(3), 214-229. https://doi.org/10.1002/pst.1983
