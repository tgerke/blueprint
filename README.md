
<!-- README.md is generated from README.Rmd. Please edit that file -->

# blueprint <img src="man/figures/logo.jpg" align="right" height="114" alt="" />

<!-- badges: start -->

<!-- badges: end -->

> **⚠️ Warning: Active Development** This package is currently under
> active development and has not yet reached a stable release (v1.0).
> The API may change without notice, and the package should **not be
> used for production use-cases** or critical decision-making until
> version 1.0 is released. Use at your own risk.

blueprint calculates sample size and power for clinical trial designs,
then verifies the operating characteristics by simulation. The two
halves are the point: you design to hit target operating
characteristics, and you confirm by simulating that the design actually
hits them.

The package is built around a small grammar rather than one monolithic
function per design:

- `search_simon_designs()` and friends find designs meeting your
  constraints; `design_simon()` constructs a design directly from known
  parameters, for example when reproducing a design written into a
  protocol.
- `assume_null()`, `assume_alternative()`, and `assume_response()`
  describe data-generating truths, kept separate from the design itself.
- `simulate_trial()` and `evaluate()` produce tidy tibbles of trial
  realizations and operating characteristics.
- `verify()` compares exact and simulated operating characteristics side
  by side.
- `tidy()`, `glance()`, and `draft_protocol_text()` turn a design into
  report tables and protocol language.

## Installation

You can install the development version of blueprint from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("tgerke/blueprint")
```

## A Simon two-stage design in four verbs

``` r
library(blueprint)

# 1. Search for designs: uninteresting response rate 20%, target 40%
candidates <- search_simon_designs(p0 = 0.20, pa = 0.40, alpha = 0.05,
                                   power = 0.80, n_max = 60)

# 2. Pick the optimal design (minimum expected N under the null)
design <- pick_design(candidates, "optimal")
design
#> <Simon two-stage design: optimal> 
#>   Hypotheses: H0 p = 0.2 vs. H1 p = 0.4 (one-sided)
#>   Stage 1: enroll 13; stop for futility if <= 3 responses
#>   Stage 2: enroll 30 more (43 total); reject H0 if >= 13 responses
#>   Exact operating characteristics:
#>     type I error 0.0496 | power 0.8002
#>     PET(H0) 0.747 | E[N | H0] 20.6

# 3. Verify the operating characteristics by simulation
verify(design, n_sims = 10000, seed = 2026)
#> <verification: exact vs. simulated (10,000 trials/scenario)>
#> # A tibble: 6 × 6
#>   scenario    metric            exact simulated difference    mc_se
#> * <chr>       <chr>             <dbl>     <dbl>      <dbl>    <dbl>
#> 1 null        prob_success     0.0496    0.0493  -0.000281  0.00217
#> 2 null        prob_early_stop  0.747     0.751    0.00388   0.00435
#> 3 null        expected_n      20.6      20.5     -0.116    NA      
#> 4 alternative prob_success     0.800     0.803    0.00289   0.00400
#> 5 alternative prob_early_stop  0.169     0.167   -0.00128   0.00374
#> 6 alternative expected_n      37.9      38.0      0.0384   NA      
#> 
#> All probability metrics are within 3.5 Monte Carlo standard errors of their exact values.

# 4. Ask what happens at a truth the design was not built for
evaluate(design, under = assume_response(0.30))
#> # A tibble: 1 × 6
#>   scenario p_true prob_success prob_early_stop expected_n method
#>   <chr>     <dbl>        <dbl>           <dbl>      <dbl> <chr> 
#> 1 custom      0.3        0.408           0.421       30.4 exact
```

And when the protocol is due:

``` r
cat(draft_protocol_text(design))
```

A Simon optimal two-stage design (Simon, 1989) will be used. The null
hypothesis that the true response rate is 20% will be tested against a
one-sided alternative. In the first stage, 13 patients will be accrued.
If there are 3 or fewer responses in these 13 patients, the study will
be stopped. Otherwise, 30 additional patients will be accrued for a
total of 43. The null hypothesis will be rejected if 13 or more
responses are observed in 43 patients. This design yields a type I error
rate of 0.050 and power of 0.800 when the true response rate is 40%.
Under the null hypothesis, the probability of early termination is 0.75
and the expected sample size is 20.6 patients.

See `vignette("design-grammar")` for the full workflow, including gt
tables of operating characteristics across scenarios.

## Other design families

These predate the design grammar and keep their original interfaces for
now; they will migrate to the grammar in a future release.

- **Two-stage single-arm designs with time-to-event endpoints**
  (`two_stage_single_arm_tte()`, `single_stage_single_arm_tte()`),
  following Wu et al. (2020). Supports Weibull, log-normal, gamma,
  log-logistic, and logspline survival distributions, with futility
  stopping at stage 1 and restricted follow-up. See
  `vignette("two-stage-tte")`.
- **Simple basket trial design** (`simple_basket_trial()`), the modified
  Simon approach of Cunanan et al. (2017): a Fisher’s exact
  heterogeneity test at interim decides whether baskets are pooled or
  analyzed separately. See `vignette("simple-basket")`.
- **Bayesian basket trial design** (`two_stage_basket_trial()`),
  following Jing et al. (2022): information borrowing across baskets
  with per-basket early termination. See `vignette("basket-trial")`.

## References

Simon, R. (1989). Optimal two-stage designs for phase II clinical
trials. *Controlled Clinical Trials*, 10(1), 1-10.
<https://doi.org/10.1016/0197-2456(89)90015-9>

Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Optimal
two-stage phase II survival trial design. *Pharmaceutical Statistics*,
19(3), 214-229. <https://doi.org/10.1002/pst.1983>

Cunanan, K. M., Iasonos, A., Shen, R., Begg, C. B., & Gönen, M. (2017).
An efficient basket trial design. *Statistics in Medicine*, 36(10),
1568-1579. <https://doi.org/10.1002/sim.7227>

Jing, Y., Qin, R., & Liu, S. (2022). Two-stage basket trial design with
time-to-event endpoint allowing for early termination of individual
baskets and efficient information borrowing across baskets. *Statistics
in Medicine*, 41(27), 5427-5443. <https://doi.org/10.1002/sim.9576>
