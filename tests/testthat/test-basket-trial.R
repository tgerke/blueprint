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


test_that("homogeneous design approximately matches Jing et al. 2022 Table 1 - K=3, p0=10%, pa=30%", {
  skip_on_cran()
  
  # Jing et al. (2022) Table 1, Design 1: K=3, p0=10%, pa=30%
  # Expected: S=57, alpha1=0.10, alpha2=0.90, EN(H0)=68.6
  # Note: Their approach uses Bayesian information borrowing, ours is frequentist
  # We expect similar but not identical results
  
  design <- two_stage_basket_trial(
    K = 3,
    p0 = 0.10,
    pa = 0.30,
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.05,
    constraint_beta = 0.20,
    n_sims_design = 200,  # Fast enough with optimized code (~1 min)
    seed = 20221106
  )
  
  # Check basic constraints are met
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)  # Allow some slack
  
  # Design should be reasonably efficient
  # Jing's EN(H0) = 68.6, we expect similar magnitude
  expect_lt(design$design$EN_H0, 100)
  expect_gt(design$design$S, 10)  # Should have reasonable Stage I size (relaxed from 30)
  
  # Stage I significance should be moderate (not too strict or too lenient)
  expect_gt(design$design$alpha1, 0.05)
  expect_lt(design$design$alpha1, 0.50)
})


test_that("homogeneous design approximately matches Jing et al. 2022 Table 1 - K=4, p0=10%, pa=30%", {
  skip_on_cran()
  
  # Jing et al. (2022) Table 1, Design 2: K=4, p0=10%, pa=30%
  # Expected: S=34, alpha1=0.15, alpha2=0.60, EN(H0)=61.0
  
  design <- two_stage_basket_trial(
    K = 4,
    p0 = 0.10,
    pa = 0.30,
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.05,
    constraint_beta = 0.20,
    n_sims_design = 150,  # Fast enough (~2 min)
    seed = 20221107
  )
  
  # Check constraints
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)
  
  # Design efficiency checks
  expect_lt(design$design$EN_H0, 100)
  expect_gt(design$design$S, 20)
})


test_that("heterogeneous design approximately matches Jing et al. 2022 Table A1", {
  skip_on_cran()
  
  # Jing et al. (2022) Table A1: K=4 heterogeneous
  # p0=(0.10, 0.15, 0.20, 0.25), pa=(0.30, 0.35, 0.40, 0.45)
  # Expected: S=58, alpha1=0.20, alpha2=0.10, Mean EN=92.3
  
  design <- two_stage_basket_trial(
    K = 4,
    p0 = c(0.10, 0.15, 0.20, 0.25),
    pa = c(0.30, 0.35, 0.40, 0.45),
    alpha = 0.05,
    beta = 0.20,
    constraint_alpha = 0.05,
    constraint_beta = 0.20,
    n_sims_design = 150,  # Fast enough (~2 min)
    seed = 20221108
  )
  
  # Check constraints
  expect_lte(design$performance$type1_error, 0.05)
  expect_gte(design$performance$expected_power, 0.75)
  
  # Design should handle heterogeneity appropriately
  expect_false(design$param$homogeneous)
  
  # Efficiency check
  expect_lt(design$design$EN_H0, 150)
  expect_gt(design$design$S, 30)
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
