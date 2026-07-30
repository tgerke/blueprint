# blueprint development notes

## Testing

- `tests/testthat/test-basket-simple.R` holds extended validation tests that
  replicate Cunanan et al. (2017) Tables 2-4 by full grid search. Each takes
  hours on the 2 cores R CMD check allows, so they skip unless
  `BLUEPRINT_EXTENDED_TESTS=true` is set (see
  `tests/testthat/helper-extended.R`). Run them on purpose with:
  `BLUEPRINT_EXTENDED_TESTS=true Rscript -e 'devtools::test(filter = "basket-simple")'`
- Routine coverage for the same functions lives in
  `tests/testthat/test-basket-simple-quick.R`, which uses narrow grids with
  `parallel = FALSE` and runs under R CMD check.

## Parallel processing

- R CMD check sets `_R_CHECK_LIMIT_CORES_` and errors if
  `parallel::mclapply()` spawns more than 2 processes. Any function that
  picks a default core count must cap at 2 when that variable is set;
  `simple_basket_trial()` has the reference implementation.
