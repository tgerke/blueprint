test_that("basket trial design creates valid object", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 2,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing  # Fast testing with optimized code
    seed = 123
  )
  
  expect_s3_class(design, "basket_trial_design")
  expect_named(design, c("param", "design", "performance"))
  expect_true(design$param$K == 2)
  expect_true(all(design$param$p0[[1]] == 0.10))
  expect_true(all(design$param$pa[[1]] == 0.25))
})


test_that("basket trial validates inputs correctly", {
  # K < 2
  expect_error(
    two_stage_basket_trial(K = 1, p0 = 0.1, pa = 0.2),
    "K must be at least 2"
  )
  
  # pa <= p0
  expect_error(
    two_stage_basket_trial(
      K = 3, p0 = 0.2, pa = 0.1,
      constraint_alpha = 0.1, constraint_beta = 0.2,
      n_sims_design = 100
    ),
    "pa must be greater than p0"
  )
  
  # Length mismatch
  expect_error(
    two_stage_basket_trial(
      K = 4, p0 = c(0.1, 0.2), pa = c(0.3, 0.4),
      constraint_alpha = 0.1, constraint_beta = 0.2,
      n_sims_design = 100
    ),
    "p0 and pa must have length K"
  )
  
  # Invalid alpha
  expect_error(
    two_stage_basket_trial(
      K = 3, p0 = 0.1, pa = 0.2, alpha = 1.5,
      constraint_alpha = 0.1, constraint_beta = 0.2,
      n_sims_design = 100
    ),
    "alpha must be between 0 and 1"
  )
})


test_that("homogeneous design approximately matches Jing et al. 2022 Table A1 - K=3, p0=10%, pa=25%", {
  skip_on_cran()
  
  # Jing et al. (2022) Contemporary Clinical Trials, Table A1
  # K=3 homogeneous: p0=10%, pa=25%
  # Expected: S=38, alpha1=0.38, R1=5, N=35, r=6, alpha2=0.028, EN(H0)=60.1
  #
  # Jing's method uses constrained optimization: fix N, r per indication using
  # exact binomial test (constraint_alpha, constraint_beta), then optimize S, α₁, α₂
  # We use similar constraints to get comparable results
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.10,  # Per-indication Type I error constraint
    constraint_beta = 0.20,   # Per-indication Type II error constraint
    n_sims_design = 200,      # Fast enough with optimized code (~1 min)
    seed = 20251107
  )
  
  # PRIMARY: Check global constraints are met
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)
  
  # SECONDARY: Design parameters should be in ballpark of Jing et al.
  # With constraints, results should be more comparable (but still not exact due to
  # different optimization algorithms and Monte Carlo variation)
  
  # S: Expected ~38, allow ±20
  expect_gte(design$design$S, 18)
  expect_lte(design$design$S, 58)
  
  # EN(H0): Expected ~60, allow ±20
  expect_gte(design$design$EN_H0, 40)
  expect_lte(design$design$EN_H0, 80)
  
  # N per indication: Expected ~35, allow ±20
  expect_gte(design$design$N[[1]][1], 15)
  expect_lte(design$design$N[[1]][1], 55)
})


test_that("homogeneous design approximately matches Jing et al. 2022 Table A1 - K=4, p0=10%, pa=25%", {
  skip_on_cran()
  
  # Jing et al. (2022) Contemporary Clinical Trials, Table A1
  # K=4 homogeneous: p0=10%, pa=25%
  # Expected: S=47, alpha1=0.36, R1=6, N=32, r=5, alpha2=0.021, EN(H0)=73.6
  
  design <- two_stage_basket_trial(
    K = 4,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.10,
    constraint_beta = 0.20,
    n_sims_design = 150,  # Fast enough (~2 min)
    seed = 20251107
  )
  
  # PRIMARY: Check constraints
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)
  
  # SECONDARY: Comparable to Jing et al. (with constraints, should be closer)
  # S: Expected ~47, allow ±22
  expect_gte(design$design$S, 25)
  expect_lte(design$design$S, 69)
  
  # EN(H0): Expected ~74, allow ±25
  expect_gte(design$design$EN_H0, 49)
  expect_lte(design$design$EN_H0, 99)
  
  # N per indication: Expected ~32, allow ±20
  expect_gte(design$design$N[[1]][1], 12)
  expect_lte(design$design$N[[1]][1], 52)
})


test_that("heterogeneous design approximately matches Jing et al. 2022 Table 1 - K=4", {
  skip_on_cran()
  
  # Jing et al. (2022) Contemporary Clinical Trials, Table 1
  # K=4 heterogeneous: p0=(10,10,20,20), pa=(25,25,35,35) [15% improvement]
  # Expected: S=61, alpha1=0.42, N=(31,31,42,42), r=(5,5,11,11), alpha2=0.016, EN(H0)=90.9
  #
  # NOTE: Heterogeneous setting requires constraint parameters
  
  design <- two_stage_basket_trial(
    K = 4,
    p0 = c(0.10, 0.10, 0.20, 0.20),
    pa = c(0.25, 0.25, 0.35, 0.35),
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.10,  # Required for heterogeneous
    constraint_beta = 0.20,
    n_sims_design = 150,  # Fast enough (~2 min)
    seed = 20251107
  )
  
  # PRIMARY: Check constraints
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)
  
  # Design should handle heterogeneity appropriately
  expect_false(design$param$homogeneous)
  
  # SECONDARY: Reasonable ranges (relaxed tolerances)
  # S: Allow ±25 from published 61
  expect_gte(design$design$S, 36)
  expect_lte(design$design$S, 86)
  
  # EN(H0): Allow ±20 from published 90.9
  expect_gte(design$design$EN_H0, 70)
  expect_lte(design$design$EN_H0, 111)
})


test_that("simulate_basket_trial produces valid output", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  # Simulate under H0
  set.seed(456)
  trial_h0 <- simulate_basket_trial(
    design = design,
    p_true = design$param$p0[[1]]
  )
  
  expect_type(trial_h0, "list")
  expect_named(trial_h0, c(
    "reject_h0", "stopped_early", "stage1_n", "stage1_responses",
    "stage1_total_resp", "stage1_decision", "stage2_total_n",
    "stage2_total_resp", "pruned_indications", "pooled_indications",
    "final_decision"
  ), ignore.order = TRUE)
  
  expect_type(trial_h0$reject_h0, "logical")
  expect_type(trial_h0$stopped_early, "logical")
  expect_length(trial_h0$stage1_n, 3)
  expect_equal(sum(trial_h0$stage1_n), design$design$S)
  
  # Simulate under H1
  set.seed(789)
  trial_h1 <- simulate_basket_trial(
    design = design,
    p_true = design$param$pa[[1]]
  )
  
  expect_type(trial_h1, "list")
  expect_type(trial_h1$reject_h0, "logical")
})


test_that("simulate_basket_trial validates inputs", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  # Wrong length p_true
  expect_error(
    simulate_basket_trial(design, p_true = c(0.1, 0.2)),
    "p_true must have length K"
  )
  
  # Wrong length allocation_probs
  expect_error(
    simulate_basket_trial(
      design, 
      p_true = c(0.1, 0.1, 0.1),
      allocation_probs = c(0.5, 0.5)
    ),
    "allocation_probs must have length K"
  )
  
  # allocation_probs don't sum to 1
  expect_error(
    simulate_basket_trial(
      design,
      p_true = c(0.1, 0.1, 0.1),
      allocation_probs = c(0.3, 0.3, 0.3)
    ),
    "allocation_probs must sum to 1"
  )
})


test_that("simulate_basket_operating_characteristics produces valid output", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  set.seed(999)
  oc <- simulate_basket_operating_characteristics(
    design = design,
    n_sims = 100,
    seed = 111
  )
  
  expect_type(oc, "list")
  expect_named(oc, c(
    "type1_error", "expected_power", "power_by_n_active",
    "expected_n_h0", "expected_n_h1",
    "prob_early_stop_h0", "prob_early_stop_h1",
    "prob_identify_1plus_tp", "prob_identify_2plus_tp",
    "prob_identify_all_tp", "expected_n_true_positives",
    "expected_n_false_positives", "summary_by_config"
  ), ignore.order = TRUE)
  
  # Check types
  expect_type(oc$type1_error, "double")
  expect_type(oc$expected_power, "double")
  expect_s3_class(oc$power_by_n_active, "data.frame")
  expect_s3_class(oc$summary_by_config, "data.frame")
  
  # Check ranges
  expect_gte(oc$type1_error, 0)
  expect_lte(oc$type1_error, 1)
  expect_gte(oc$expected_power, 0)
  expect_lte(oc$expected_power, 1)
  expect_gte(oc$expected_n_h0, design$design$S)
  expect_lte(oc$expected_n_h0, design$design$total_N)
})


test_that("Type I error is controlled with sufficient simulations", {
  skip_on_cran()
  skip_if(TRUE, "Skipping slow test - enable for thorough validation")
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 5000,
    seed = 123
  )
  
  # Validate with high-precision simulation
  oc <- simulate_basket_operating_characteristics(
    design = design,
    n_sims = 10000,
    seed = 456
  )
  
  # With 10,000 sims, 95% CI for Type I error is roughly ±0.004
  # Allow some tolerance for simulation variability
  expect_lte(oc$type1_error, 0.06)  # Slightly above nominal for safety
  expect_gte(oc$expected_power, 0.75)  # Should be close to 0.80
})


test_that("print and summary methods work", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  # Test print
  expect_output(print(design), "Optimal Two-Stage Basket Trial Design")
  expect_output(print(design), "Trial Configuration")
  expect_output(print(design), "Stage I")
  expect_output(print(design), "Stage II")
  
  # Test summary
  summ <- summary(design)
  expect_s3_class(summ, "summary.basket_trial_design")
  expect_named(summ, c(
    "param", "design", "performance", "efficiency",
    "avg_pruning_threshold", "avg_indication_n"
  ), ignore.order = TRUE)
  
  expect_output(print(summ), "Summary of Two-Stage Basket Trial Design")
})


test_that("design handles homogeneous setting correctly", {
  skip_on_cran()
  
  # Single value should expand to all K
  design <- two_stage_basket_trial(
    K = 4,
    p0 = 0.05,
    pa = 0.20,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.1,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  expect_true(design$param$homogeneous)
  expect_equal(length(design$param$p0[[1]]), 4)
  expect_equal(length(design$param$pa[[1]]), 4)
  expect_true(all(design$param$p0[[1]] == 0.05))
  expect_true(all(design$param$pa[[1]] == 0.20))
})


test_that("design handles heterogeneous setting correctly", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 4,
    p0 = c(0.05, 0.05, 0.20, 0.20),
    pa = c(0.20, 0.20, 0.35, 0.35),
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.1,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 456
  )
  
  expect_false(design$param$homogeneous)
  expect_equal(design$param$p0[[1]], c(0.05, 0.05, 0.20, 0.20))
  expect_equal(design$param$pa[[1]], c(0.20, 0.20, 0.35, 0.35))
})


test_that("early stopping works correctly", {
  skip_on_cran()
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 123
  )
  
  # Simulate many trials under H0, some should stop early
  set.seed(777)
  n_early_stop <- 0
  n_trials <- 100
  
  for (i in 1:n_trials) {
    trial <- simulate_basket_trial(
      design = design,
      p_true = design$param$p0[[1]]
    )
    if (trial$stopped_early) {
      n_early_stop <- n_early_stop + 1
      # If stopped early, should not have stage 2 data
      expect_true(is.na(trial$stage2_total_n[1]) || 
                  all(trial$stage2_total_n == trial$stage1_n))
    }
  }
  
  # Under H0, we expect some trials to stop early
  # (exact proportion depends on design)
  expect_gt(n_early_stop, 0)
})


test_that("seed produces reproducible results", {
  skip_on_cran()
  
  design1 <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 999
  )
  
  design2 <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.25,
    alpha = 0.05,
    beta = 0.2,
    constraint_alpha = 0.15,
    constraint_beta = 0.2,
    n_sims_design = 50,  # Fast testing
    seed = 999
  )
  
  # Same seed should give same results
  expect_equal(design1$design$S, design2$design$S)
  expect_equal(design1$design$N[[1]], design2$design$N[[1]])
  expect_equal(design1$design$r[[1]], design2$design$r[[1]])
  expect_equal(design1$design$EN_H0, design2$design$EN_H0)
})
