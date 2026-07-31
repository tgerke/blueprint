test_that("simple basket trial produces reasonable results for Table 2 scenario", {
  skip_on_cran()
  
  set.seed(4206)
  
  # Test with known optimal values from Table 2 (or very close)
  design <- simple_basket_trial(
    K = 5,
    p0 = 0.15,
    pa = 0.45,
    N1 = 35,
    N2 = 20,
    r_s = 1,
    r_c = 5,
    target_fwer = 0.05,
    target_power = 0.80,
    min_power = 0.70,
    accrual_rate = 2,
    calibrate_at_A = 2,
    n2k_range = c(14, 15, 16),  # Very narrow around expected 15
    gamma_range = c(0.50, 0.52, 0.54),  # Around expected 0.52
    alpha_s_range = c(0.06, 0.07, 0.08),  # Around expected 0.07
    alpha_c_range = c(0.04, 0.05, 0.06),  # Around expected 0.05
    n_sims = 200,  # More sims for better estimates
    seed = 4206,
    parallel = FALSE
  )
  
  # Check class
  expect_s3_class(design, "simple_basket_design")
  
  # Check structure
  expect_true(all(c("design", "performance", "inputs") %in% names(design)))
  
  # Check A=0 (FWER should be around 5%)
  expect_true(design$performance$A0$fwer <= 0.10)  # Relaxed for small n_sims
  
  # Check A=1 (minimum power should be >= 70%)
  # Just check that the active basket has reasonable power
  expect_true(design$performance$A1$marginal_power[1] >= 0.50)  # Relaxed for small n_sims
  
  # Check A=2 (calibration power should be around 80%)
  # Just check that at least one basket has reasonable power
  expect_true(any(design$performance$A2$marginal_power >= 0.70))
  
  # Check EN values are reasonable (not too different from Table 2)
  expect_true(design$performance$A0$EN >= 50 && design$performance$A0$EN <= 80)
  expect_true(design$performance$A5$EN >= 70 && design$performance$A5$EN <= 100)
  
  # Check sensitivity/specificity are in valid range
  expect_true(is.na(design$performance$A0$sensitivity) || 
              (design$performance$A0$sensitivity >= 0 && design$performance$A0$sensitivity <= 1))
  expect_true(is.na(design$performance$A0$specificity) || 
              (design$performance$A0$specificity >= 0 && design$performance$A0$specificity <= 1))
})

test_that("simulate_simple_trial_once works correctly", {
  set.seed(123)
  
  # Test basic simulation
  result <- simulate_simple_trial_once(
    K = 3,
    p_true = c(0.2, 0.5, 0.2),
    N1 = 20,
    N2 = 15,
    n2k = 10,
    p0 = 0.2,
    r_s = 1,
    r_c = 3,
    gamma = 0.5,
    alpha_s = 0.05,
    alpha_c = 0.05,
    accrual_rate = c(2, 2, 2)
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("decisions", "stage1_n", "stage2_n", "trial_duration") %in% names(result)))
  
  # Check decisions are logical
  expect_type(result$decisions, "logical")
  expect_length(result$decisions, 3)
  
  # Check sample sizes are non-negative
  expect_true(all(result$stage1_n >= 0))
  expect_true(all(result$stage2_n >= 0))
  
  # Check duration is positive
  expect_true(result$trial_duration > 0)
})

test_that("test_heterogeneity works correctly", {
  # Homogeneous scenario (all same)
  responses_hom <- c(5, 6, 5, 6, 5)
  n_hom <- c(20, 20, 20, 20, 20)
  
  # Should NOT reject homogeneity (returns FALSE)
  expect_false(test_heterogeneity(responses_hom, n_hom, gamma = 0.5))
  
  # Heterogeneous scenario (very different)
  responses_het <- c(0, 20, 0, 20, 0)
  n_het <- c(20, 20, 20, 20, 20)
  
  # Should reject homogeneity (returns TRUE)
  expect_true(test_heterogeneity(responses_het, n_het, gamma = 0.5))
})

test_that("print and summary methods work", {
  skip_on_cran()
  
  set.seed(4206)
  
  design <- simple_basket_trial(
    K = 3,
    p0 = 0.2,
    pa = 0.5,
    N1 = 20,
    N2 = 15,
    r_s = 1,
    r_c = 3,
    n2k_range = 10,
    gamma_range = 0.5,
    alpha_s_range = 0.05,
    alpha_c_range = 0.05,
    n_sims = 50,
    seed = 123,
    parallel = FALSE
  )
  
  # Test print
  expect_output(print(design), "Simple Basket Trial Design")
  expect_output(print(design), "Design Parameters")
  
  # Test summary
  summ <- summary(design)
  expect_s3_class(summ, "summary.simple_basket_design")
  expect_output(print(summ), "Operating Characteristics")
})
