
<!-- README.md is generated from README.Rmd. Please edit that file -->

# blueprint <img src="man/figures/logo.jpg" align="right" height="114" alt="" />

<!-- badges: start -->

<!-- badges: end -->

> **⚠️ Warning: Active Development**  
> This package is currently under active development and has not yet
> reached a stable release (v1.0). The API may change without notice,
> and the package should **not be used for production use-cases** or
> critical decision-making until version 1.0 is released. Use at your
> own risk.

The **blueprint** package provides a suite of functions to calculate
sample size, power, and simulate data according to various assumptions
for clinical trial designs.

## Installation

You can install the development version of blueprint from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("tgerke/blueprint")
```

## Features

### Two-Stage Single-Arm Designs with Time-to-Event Endpoints

The package implements optimal two-stage designs for single-arm phase II
clinical trials with time-to-event endpoints, based on the methodology
described in:

> Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Optimal
> two-stage phase II survival trial design. *Pharmaceutical Statistics*,
> 19(3), 214-229. <https://doi.org/10.1002/pst.1983>

Key features:

- **Flexible survival distributions**: Supports Weibull, log-normal,
  gamma, log-logistic, and non-parametric logspline distributions
- **One-sample log-rank test**: Uses exact variance estimates for
  adequate power
- **Early stopping**: Allows for futility stopping at the first stage
- **Restricted follow-up**: Accommodates realistic trial constraints

## Usage

### Example: Small-Cell Lung Cancer Trial

This example reproduces the design from Wu et al. (2020) for a
single-arm phase II trial evaluating immunotherapy as second-line
treatment for small-cell lung cancer patients.

``` r
library(blueprint)

# Design parameters:
# - Null hypothesis: median PFS = 3.5 months (Weibull distribution)
# - Alternative: median PFS = 5 months (hazard ratio = 0.5913)
# - Restricted follow-up: 5 months
# - Accrual rate: 2 patients per month
# - Type I error: 5%, Power: 80%

design <- two_stage_single_arm_tte(
  shape = 1.47327,      # Weibull shape parameter
  S0 = 0.5,             # Survival probability at x0 under H0
  x0 = 3.5,             # Time point for S0 (median)
  hr = 0.5913,          # Hazard ratio under H1
  tf = 5,               # Follow-up time
  rate = 2,             # Accrual rate (patients/month)
  alpha = 0.05,         # Type I error (one-sided)
  beta = 0.2,           # Type II error (80% power)
  dist = "WB",          # Weibull distribution
  restricted = FALSE    # Unrestricted follow-up (Wu et al. default)
)

# View the design
print(design)
```

### Interpreting Results

The function returns a list with three components:

1.  **`param`**: Input parameters
2.  **`Single_stage`**: Single-stage design for comparison
    - `nsingle`: Total sample size
    - `tasingle`: Total accrual time
    - `csingle`: Critical value
3.  **`Two_stage`**: Optimal two-stage design
    - `n1`: Stage 1 sample size
    - `c1`: Stage 1 critical value (stop for futility if Z₁ ≤ c₁)
    - `n`: Total sample size
    - `c`: Final critical value (reject H₀ if Z \> c)
    - `t1`: Calendar time for interim analysis
    - `MTSL`: Maximum total study length
    - `ES`: Expected sample size under H₀
    - `PS`: Probability of early stopping under H₀

### Conducting the Trial

Based on the example design:

1.  **Stage 1**: Enroll `n1` patients. After all patients are enrolled,
    conduct an interim analysis when each patient has been followed for
    up to `tf` months or until an event occurs.
    - If Z₁ ≤ c₁, stop the trial for futility
    - Otherwise, continue to Stage 2
2.  **Stage 2**: Enroll additional patients to reach total sample size
    `n`. Conduct final analysis after the last patient has been followed
    for `tf` months.
    - If Z ≤ c, conclude no efficacy
    - If Z \> c, conclude the treatment is promising

## Supported Distributions

The package supports five survival distributions under the null
hypothesis:

| Code | Distribution | Survival Function                             |
|------|--------------|-----------------------------------------------|
| `WB` | Weibull      | S(t) = exp(-(t/b)^a)                          |
| `LN` | Log-normal   | S(t) = 1 - Φ((log(t) - μ)/σ)                  |
| `GM` | Gamma        | S(t) = 1 - I_a(t/b)                           |
| `LG` | Log-logistic | S(t) = 1/(1 + (t/b)^a)                        |
| `SP` | Logspline    | Non-parametric (requires `logspline` package) |

For parametric distributions, you need to specify only the shape
parameter; the scale parameter is derived from the specified survival
probability (`S0`) at time point `x0`.
