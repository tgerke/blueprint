
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

### 1. Two-Stage Single-Arm Designs with Time-to-Event Endpoints

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

### 2. Two-Stage Basket Trial Designs with Binary Endpoints

The package also implements two-stage basket trial designs for
evaluating a single treatment across multiple cancer types or genetic
subtypes, based on:

> Jing, Y., Qin, R., & Liu, S. (2022). Two-stage basket trial design
> with time-to-event endpoint allowing for early termination of
> individual baskets and efficient information borrowing across baskets.
> *Statistics in Medicine*, 41(27), 5427-5443.
> <https://doi.org/10.1002/sim.9576>

Key features:

- **Multiple indications**: Simultaneously evaluate treatment across K
  different cancer types/subtypes
- **Information borrowing**: Bayesian approach borrows strength across
  baskets when appropriate
- **Individual basket early termination**: Each indication can stop for
  futility at interim
- **Flexible settings**: Supports both homogeneous (all baskets
  identical) and heterogeneous (baskets differ) scenarios

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
2.  **`single_stage`**: Single-stage design for comparison
    - `nsingle`: Total sample size
    - `tasingle`: Total accrual time
    - `csingle`: Critical value
3.  **`two_stage`**: Optimal two-stage design
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

### Example: Basket Trial

This example designs a basket trial to evaluate a targeted therapy
across three different cancer types with a shared biomarker.

``` r
# Design a basket trial:
# - 3 cancer types (baskets/indications)
# - Null response rate: 10% per basket
# - Target response rate: 30% per basket
# - Type I error: 5% per basket, Power: 80% per basket

basket_design <- two_stage_basket_trial(
  K = 3,                    # Number of indications
  p0 = 0.10,                # Null response rate
  pa = 0.30,                # Alternative response rate
  constraint_alpha = 0.05,  # Type I error per basket
  constraint_beta = 0.20,   # Type II error per basket (80% power)
  min_S = 45,               # Minimum sample size per basket to search
  max_S = 70,               # Maximum sample size per basket to search
  alpha1_grid = seq(0.05, 0.95, by = 0.05),  # Grid for interim futility threshold
  alpha2_grid = seq(0.05, 0.95, by = 0.05),  # Grid for final efficacy threshold
  nsim_oc = 1000,           # Simulations per design evaluation
  parallel = TRUE,          # Use parallel processing
  ncores = 4                # Number of cores
)

# View the design
print(basket_design$design)

# View operating characteristics by basket
print(basket_design$operating_characteristics)
```

The basket trial design returns:

- **`S`**: Sample size per basket (stage 2 total)
- **`alpha1`**: Interim futility threshold (Bayesian posterior
  probability)
- **`alpha2`**: Final efficacy threshold (Bayesian posterior
  probability)
- **Operating characteristics** for each basket:
  - Type I error rate
  - Power
  - Expected sample size under H₀
  - Probability of early stopping

**Key advantage**: Information borrowing across baskets reduces total
sample size compared to running separate trials for each indication.

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

## Simulation and Validation

The package includes comprehensive simulation capabilities for
validating designs and understanding trial behavior:

``` r
# Simulate a single trial realization
trial_result <- simulate_trial(
  design = design,  # Any design object (TTE or basket)
  shape = 1.47327,  # For TTE designs
  scale = scale_h1  # For TTE designs
  # p_true = c(0.30, 0.30, 0.30)  # For basket trials
)

# Validate operating characteristics with Monte Carlo simulation
oc_results <- simulate_operating_characteristics(
  design = design,
  shape = 1.47327,
  scale_h0 = scale_h0,  # Under null hypothesis
  scale_h1 = scale_h1,  # Under alternative hypothesis
  n_sims = 1000,
  seed = 123
)

# Check Type I error and power
cat("Type I Error:", round(oc_results$type1_error, 3), "\n")
cat("Power:       ", round(oc_results$power, 3), "\n")
```

See the package vignettes for detailed examples and validation against
published results.

## Documentation

Comprehensive documentation is available:

- **Vignettes**:
  - `vignette("two-stage-tte")` - Introduction to two-stage TTE designs
  - `vignette("single-stage-simulations")` - Single-stage design
    simulations
  - `vignette("wu-validation")` - Validation against Wu et al. (2020)
  - `vignette("jing-validation")` - Validation against Jing et
    al. (2022)
- **Function help**: `?two_stage_single_arm_tte`,
  `?two_stage_basket_trial`, `?simulate_trial`

## References

Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Optimal
two-stage phase II survival trial design. *Pharmaceutical Statistics*,
19(3), 214-229. <https://doi.org/10.1002/pst.1983>

Jing, Y., Qin, R., & Liu, S. (2022). Two-stage basket trial design with
time-to-event endpoint allowing for early termination of individual
baskets and efficient information borrowing across baskets. *Statistics
in Medicine*, 41(27), 5427-5443. <https://doi.org/10.1002/sim.9576>
