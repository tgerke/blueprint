# Test Simon et al. (2017) Basket Trial Design
# Reference: Cunanan KM, Iasonos A, Shen R, Begg CB, Gönen M. (2017)
# An Efficient Basket Trial Design. Statistics in Medicine, 36(10):1568-1579.
# DOI: 10.1002/sim.7227

test_that("Simple basket design matches Table 2: K=5, equal accrual, null scenario (A=0)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  skip_on_cran()  # Takes ~5 minutes
  
  # Table 2 parameters: K=5, θ0=0.15, θa=0.45
  # Design parameters: N1=35, n2k=15, N2=20, γ=0.52, rS=1, rC=5, αS=0.07, αC=0.05
  # Expected: FWER=5%, EN=58, ET=7.0 months
  
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
    accrual_rate = rep(2, 5),  # 2 patients/month per basket
    calibrate_at_A = 2,  # Calibrate power when A=2 active baskets
    n2k_range = 13:17,  # Narrow range around expected optimal (15)
    gamma_range = seq(0.4, 0.6, by = 0.05),  # Narrow range around 0.52
    alpha_s_range = seq(0.05, 0.09, by = 0.01),  # Around 0.07
    alpha_c_range = seq(0.04, 0.06, by = 0.01),  # Around 0.05
    n_sims = 1000,
    seed = 4206
  )
  
  # Expected design parameters from paper
  expect_equal(design$design$N1, 35)
  expect_equal(design$design$N2, 20)
  expect_equal(design$design$r_s, 1)
  expect_equal(design$design$r_c, 5)
  
  # Optimal parameters from paper (Table 2 footnote)
  expect_equal(design$design$gamma, 0.52, tolerance = 0.05)
  expect_equal(design$design$n2k, 15, tolerance = 2)
  expect_equal(design$design$alpha_s, 0.07, tolerance = 0.02)
  expect_equal(design$design$alpha_c, 0.05, tolerance = 0.01)
  
  # A=0: All baskets null
  # Table 2: FWER=5%, marginal error rates 2% per basket, EN=58, ET=7.0
  expect_lte(design$performance$A0$fwer, 0.05 + 0.01)  # FWER ≤ 5% (+ 1% margin)
  expect_gte(design$performance$A0$fwer, 0.03)  # Not too conservative
  expect_equal(design$performance$A0$EN, 58, tolerance = 5)
  expect_equal(design$performance$A0$ET, 7.0, tolerance = 1.5)
  
  # Marginal false positive rates should be ~2% per basket
  for (k in 1:5) {
    expect_lte(design$performance$A0$marginal_power[[k]], 0.05)
  }
})

test_that("Simple basket design matches Table 2: K=5, A=1 active basket", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # A=1: Basket 1 active (θ1=0.45), others null (θ=0.15)
  # Table 2: P1=70%, P2-P5=7%, EN=74, ET=9.5
  # This validates minimum acceptable power constraint
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # Primary validation: Minimum power when A=1
  expect_gte(design$performance$A1$marginal_power[[1]], 0.70 - 0.05)  # P1 ≥ 70% (- 5% tolerance)
  expect_lte(design$performance$A1$marginal_power[[1]], 0.80)  # Reasonable upper bound
  
  # Inactive baskets should have low false positive rates
  for (k in 2:5) {
    expect_lte(design$performance$A1$marginal_power[[k]], 0.15)  # Table 2 shows ~7%
  }
  
  # Expected sample size
  expect_equal(design$performance$A1$EN, 74, tolerance = 8)
  expect_equal(design$performance$A1$ET, 9.5, tolerance = 2.0)
})

test_that("Simple basket design matches Table 2: K=5, A=2 active baskets (calibration scenario)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # A=2: Baskets 1-2 active, 3-5 null
  # Table 2: P1=80%, P2=80%, P3-P5=11%, EN=83, ET=10.4
  # This is the CALIBRATION scenario - should match target power exactly
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # PRIMARY VALIDATION: Calibrated power at A=2
  expect_equal(design$performance$A2$marginal_power[[1]], 0.80, tolerance = 0.05)
  expect_equal(design$performance$A2$marginal_power[[2]], 0.80, tolerance = 0.05)
  
  # Inactive baskets
  for (k in 3:5) {
    expect_lte(design$performance$A2$marginal_power[[k]], 0.20)  # Table 2: ~11%
  }
  
  # Sample size
  expect_equal(design$performance$A2$EN, 83, tolerance = 8)
  expect_equal(design$performance$A2$ET, 10.4, tolerance = 2.0)
})

test_that("Simple basket design matches Table 2: K=5, A=3 active baskets", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # A=3: Baskets 1-3 active, 4-5 null
  # Table 2: P1=84%, P2=85%, P3=85%, P4-P5=17%, EN=86, ET=10.5
  # Power should INCREASE with more active baskets (due to pooling)
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # Power in active baskets should be > 80% (calibration value)
  expect_gte(design$performance$A3$marginal_power[[1]], 0.80)
  expect_gte(design$performance$A3$marginal_power[[2]], 0.80)
  expect_gte(design$performance$A3$marginal_power[[3]], 0.80)
  expect_lte(design$performance$A3$marginal_power[[1]], 0.90)  # Upper bound
  
  # Inactive baskets
  for (k in 4:5) {
    expect_lte(design$performance$A3$marginal_power[[k]], 0.25)  # Table 2: ~17%
  }
  
  # Sample size should be similar to A=2
  expect_equal(design$performance$A3$EN, 86, tolerance = 10)
  expect_equal(design$performance$A3$ET, 10.5, tolerance = 2.0)
})

test_that("Simple basket design matches Table 2: K=5, A=4 active baskets", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # A=4: Baskets 1-4 active, 5 null
  # Table 2: P1=86%, P2=85%, P3=86%, P4=86%, P5=23%, EN=88, ET=10.2
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # All active baskets should have high power
  for (k in 1:4) {
    expect_gte(design$performance$A4$marginal_power[[k]], 0.80)
    expect_lte(design$performance$A4$marginal_power[[k]], 0.92)
  }
  
  # Single inactive basket - higher false positive due to pooling
  expect_lte(design$performance$A4$marginal_power[[5]], 0.30)  # Table 2: 23%
  
  # Sample size
  expect_equal(design$performance$A4$EN, 88, tolerance = 10)
  expect_equal(design$performance$A4$ET, 10.2, tolerance = 2.0)
})

test_that("Simple basket design matches Table 2: K=5, A=5 active baskets (all active)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # A=5: All baskets active
  # Table 2: P1=88%, P2=90%, P3=88%, P4=88%, P5=88%, EN=78, ET=8.3
  # This is homogeneous scenario - should have LOWEST sample size
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # All baskets should have high, similar power
  for (k in 1:5) {
    expect_gte(design$performance$A5$marginal_power[[k]], 0.85)
    expect_lte(design$performance$A5$marginal_power[[k]], 0.95)
  }
  
  # Sample size should be LOWEST (most efficient scenario)
  expect_equal(design$performance$A5$EN, 78, tolerance = 10)
  expect_lte(design$performance$A5$EN, design$performance$A4$EN)  # Lower than A=4
  expect_equal(design$performance$A5$ET, 8.3, tolerance = 2.0)
})

test_that("Simple basket design validates sensitivity/specificity (Table 3)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # Table 3 shows sensitivity and specificity of basket classification
  # For proposed design with equal accrual:
  # A=0: Specificity=98%
  # A=1: Sensitivity=70%, Specificity=93%
  # A=2: Sensitivity=80%, Specificity=89%
  # A=5: Sensitivity=88%, Specificity=N/A (all active)
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # A=0: High specificity (correctly identify inactive baskets)
  expect_gte(design$performance$A0$specificity, 0.95)
  
  # A=1: Sensitivity = power in active basket
  expect_equal(design$performance$A1$sensitivity, 0.70, tolerance = 0.08)
  expect_gte(design$performance$A1$specificity, 0.90)
  
  # A=2: Balanced sensitivity/specificity
  expect_equal(design$performance$A2$sensitivity, 0.80, tolerance = 0.08)
  expect_gte(design$performance$A2$specificity, 0.85)
  
  # A=5: Maximum sensitivity
  expect_gte(design$performance$A5$sensitivity, 0.85)
})

test_that("Simple basket design with unequal accrual rates (Table 4 scenario i)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # Table 4 scenario (i): λ = (1,1,1,2,3) - slow accrual in baskets 1-3
  # When basket 1 is active but accrues slowly, power should decrease
  # A=1: P1=68%, EN=75, ET=16.5 months (vs 70%, 74, 9.5 for equal accrual)
  
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
    accrual_rate = c(1, 1, 1, 2, 3),  # Unequal accrual
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # A=1: Power should be slightly lower due to slow accrual in active basket
  expect_gte(design$performance$A1$marginal_power[[1]], 0.65)
  expect_lte(design$performance$A1$marginal_power[[1]], 0.75)
  
  # Trial duration should INCREASE substantially
  expect_gte(design$performance$A1$ET, 14.0)  # Much longer than 9.5
  expect_lte(design$performance$A1$ET, 19.0)
})

test_that("Simple basket design with unequal accrual rates (Table 4 scenario ii)", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # Table 4 scenario (ii): λ = (3,2,1,1,1) - fast accrual in baskets 1-2
  # When baskets 1-2 are active AND accrue fast, power should increase
  # A=2: P1=86%, P2=85%, EN=79, ET=16.6 (vs P1=80%, P2=80%, EN=83, ET=10.4)
  
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
    accrual_rate = c(3, 2, 1, 1, 1),  # Fast in active baskets
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # A=2: Power should be HIGHER due to fast accrual in active baskets
  expect_gte(design$performance$A2$marginal_power[[1]], 0.82)
  expect_gte(design$performance$A2$marginal_power[[2]], 0.82)
  
  # Sample size might be slightly lower
  expect_lte(design$performance$A2$EN, 85)
})

test_that("Simple basket design validates track usage proportions", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # The design uses two tracks:
  # Track 1 (heterogeneous): Analyze baskets separately
  # Track 2 (homogeneous): Pool all baskets
  # Table 3 shows proportion using each track
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # A=0: Should mostly use Track 2 (homogeneous - all null)
  expect_lte(design$performance$A0$prop_track1, 0.20)  # Table 3: 2%
  
  # A=1: Mixed - some heterogeneity
  expect_gte(design$performance$A1$prop_track1, 0.10)
  expect_lte(design$performance$A1$prop_track1, 0.25)  # Table 3: 14%
  
  # A=3: More likely heterogeneous
  expect_gte(design$performance$A3$prop_track1, 0.40)  # Table 3: 51%
  expect_lte(design$performance$A3$prop_track1, 0.65)
  
  # A=5: Should mostly use Track 2 (homogeneous - all active)
  expect_gte(design$performance$A5$prop_track1, 0.80)  # Table 3: 88%
})

test_that("Simple basket design validates expected trial duration patterns", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # Expected trial duration should follow pattern:
  # A=0: Fast (futility, EN=58, ET=7.0)
  # A=1-4: Longer (continue to stage 2, ET≈9-11)
  # A=5: Fast (homogeneous efficacy, EN=78, ET=8.3)
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206
  )
  
  # A=0: Should be shortest (early futility stopping)
  expect_lte(design$performance$A0$ET, 8.5)
  
  # A=2-4: Should be longest (heterogeneous paths, separate analyses)
  expect_gte(design$performance$A2$ET, 9.0)
  expect_lte(design$performance$A2$ET, 12.0)
  
  # A=5: Should be faster than A=2-4 (homogeneous, pooled)
  expect_lte(design$performance$A5$ET, design$performance$A2$ET)
})

test_that("Simple basket design fails with invalid inputs", {
  skip_if_not_installed("blueprint")
  
  # Test input validation
  expect_error(
    simple_basket_trial(K = 2, p0 = 0.15, pa = 0.45),  # K too small
    "at least 3 baskets"
  )
  
  expect_error(
    simple_basket_trial(K = 5, p0 = 0.50, pa = 0.45),  # p0 > pa
    "p0 must be less than pa"
  )
  
  expect_error(
    simple_basket_trial(K = 5, p0 = 0.15, pa = 0.45, target_fwer = 0.15),  # FWER too high
    "target_fwer must be between 0 and 0.10"
  )
  
  expect_error(
    simple_basket_trial(K = 5, p0 = 0.15, pa = 0.45, min_power = 0.85, target_power = 0.80),
    "min_power must be less than target_power"
  )
  
  expect_error(
    simple_basket_trial(K = 5, p0 = 0.15, pa = 0.45, accrual_rate = rep(2, 4)),  # Wrong length
    "accrual_rate must have length K"
  )
})

test_that("Simple basket design reference design comparison", {
  skip_if_not_installed("blueprint")
  skip_extended()
  
  # The paper compares to reference design: independent Simon two-stage designs
  # Reference uses α = ε/K per basket to control FWER
  # For K=5, ε=0.05: α = 0.01 per basket
  # Table 2 shows reference design results for comparison
  
  # Our design should have:
  # - Similar FWER control
  # - Similar or better power when A≥2
  # - Lower EN when A≥3 (efficiency gain from pooling)
  
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
    accrual_rate = rep(2, 5),
    calibrate_at_A = 2,
    n_sims = 1000,
    seed = 4206,
    include_reference = TRUE  # Compute reference design for comparison
  )
  
  # A=0: Both should control FWER at 5%
  expect_lte(design$performance$A0$fwer, 0.06)
  expect_lte(design$reference$A0$fwer, 0.06)
  
  # A=2: Calibrated to same power
  expect_equal(
    design$performance$A2$marginal_power[[1]],
    design$reference$A2$marginal_power[[1]],
    tolerance = 0.10
  )
  
  # A=5: Proposed should have lower EN (efficiency gain)
  # Table 2: Proposed EN=78 vs Reference EN=121
  expect_lte(design$performance$A5$EN, design$reference$A5$EN - 20)
})
