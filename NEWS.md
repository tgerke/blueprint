# blueprint 0.0.0.9000

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
