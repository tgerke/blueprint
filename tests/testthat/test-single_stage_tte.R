test_that("single_stage_single_arm_tte validates against internal app: Case 1 (shape=1.0, tf=1.0)", {
  # From internal application output:
  # Power=0.9004, N=81, E=68, Ta=2.5, Tf=1.0
  # HR=0.667, M0=0.5, M1=0.75, k=1.00, Alpha=0.050, P1=0.8453
  
  design <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 81 / 2.5,  # N=81, Ta=2.5 -> rate=32.4
    alpha = 0.05,
    beta = 1 - 0.9004,  # Power = 0.9004
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  # Validate sample size
  expect_equal(design$Single_stage$n, 81, tolerance = 2)
  
  # Validate expected events (should be close to 68)
  expect_equal(design$Single_stage$E, 68, tolerance = 2)
  
  # Validate probability of event (should be close to 0.8453)
  expect_equal(design$Single_stage$P1, 0.8453, tolerance = 0.02)
  
  # Validate accrual time (should be close to 2.5)
  expect_equal(design$Single_stage$ta, 2.5, tolerance = 0.2)
})

test_that("single_stage_single_arm_tte validates against internal app: Case 2 (shape=1.25, tf=1.0)", {
  # Power=0.9026, N=52, E=47, Ta=2.5, Tf=1.0
  # HR=0.602, M0=0.5, M1=0.75, k=1.25, Alpha=0.050, P1=0.8957
  
  design <- single_stage_single_arm_tte(
    shape = 1.25,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.602,
    tf = 1.0,
    rate = 52 / 2.5,
    alpha = 0.05,
    beta = 1 - 0.9026,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_equal(design$Single_stage$n, 52, tolerance = 2)
  expect_equal(design$Single_stage$E, 47, tolerance = 2)
  expect_equal(design$Single_stage$P1, 0.8957, tolerance = 0.02)
})

test_that("single_stage_single_arm_tte validates against internal app: Case 3 (shape=1.50, tf=1.0)", {
  # Power=0.9043, N=37, E=34, Ta=2.5, Tf=1.0
  # HR=0.544, M0=0.5, M1=0.75, k=1.50, Alpha=0.050, P1=0.9284
  
  design <- single_stage_single_arm_tte(
    shape = 1.50,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.544,
    tf = 1.0,
    rate = 37 / 2.5,
    alpha = 0.05,
    beta = 1 - 0.9043,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_equal(design$Single_stage$n, 37, tolerance = 2)
  expect_equal(design$Single_stage$E, 34, tolerance = 2)
  expect_equal(design$Single_stage$P1, 0.9284, tolerance = 0.02)
})

test_that("single_stage_single_arm_tte validates against internal app: Case 4 (shape=2.00, tf=1.0)", {
  # Power=0.9035, N=22, E=21, Ta=2.5, Tf=1.0
  # HR=0.444, M0=0.5, M1=0.75, k=2.00, Alpha=0.050, P1=0.9628
  
  design <- single_stage_single_arm_tte(
    shape = 2.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.444,
    tf = 1.0,
    rate = 22 / 2.5,
    alpha = 0.05,
    beta = 1 - 0.9035,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_equal(design$Single_stage$n, 22, tolerance = 2)
  expect_equal(design$Single_stage$E, 21, tolerance = 1)
  expect_equal(design$Single_stage$P1, 0.9628, tolerance = 0.02)
})

test_that("single_stage_single_arm_tte validates against internal app: Extended Case (tf=1.5)", {
  # Power=0.9019, N=78, E=70, Ta=2.5, Tf=1.5
  # HR=0.667, M0=0.5, M1=0.75, k=1.00, Alpha=0.050, P1=0.9025
  
  design <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.5,
    rate = 78 / 2.5,
    alpha = 0.05,
    beta = 1 - 0.9019,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_equal(design$Single_stage$n, 78, tolerance = 2)
  expect_equal(design$Single_stage$E, 70, tolerance = 2)
  expect_equal(design$Single_stage$P1, 0.9025, tolerance = 0.02)
})

test_that("single_stage_single_arm_tte validates against internal app: Extended Case (tf=2.0)", {
  # Power=0.9006, N=76, E=71, Ta=2.5, Tf=2.0
  # HR=0.667, M0=0.5, M1=0.75, k=1.00, Alpha=0.050, P1=0.9386
  
  design <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 2.0,
    rate = 76 / 2.5,
    alpha = 0.05,
    beta = 1 - 0.9006,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_equal(design$Single_stage$n, 76, tolerance = 2)
  expect_equal(design$Single_stage$E, 71, tolerance = 2)
  expect_equal(design$Single_stage$P1, 0.9386, tolerance = 0.02)
})

test_that("single_stage_single_arm_tte works with one-sided test", {
  # Test that one-sided option works (should give different critical value)
  design_one_sided <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    two_sided = FALSE,
    dist = "WB",
    restricted = FALSE
  )
  
  # One-sided critical value should be qnorm(0.95) = 1.645
  expect_equal(design_one_sided$Single_stage$c, qnorm(0.95), tolerance = 0.01)
  
  # Should have smaller sample size than two-sided
  design_two_sided <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  expect_true(design_one_sided$Single_stage$n < design_two_sided$Single_stage$n)
})

test_that("single_stage_single_arm_tte returns correct structure", {
  design <- single_stage_single_arm_tte(
    shape = 1.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    two_sided = TRUE,
    dist = "WB",
    restricted = FALSE
  )
  
  # Check structure
  expect_type(design, "list")
  expect_s3_class(design, "single_stage_design")
  expect_named(design, c("param", "Single_stage"))
  
  # Check param data frame
  expect_s3_class(design$param, "data.frame")
  expect_true("shape" %in% names(design$param))
  expect_true("S0" %in% names(design$param))
  expect_true("hr" %in% names(design$param))
  expect_true("two_sided" %in% names(design$param))
  
  # Check Single_stage data frame
  expect_s3_class(design$Single_stage, "data.frame")
  expect_named(design$Single_stage, c("n", "ta", "c", "MTSL", "E", "P1"))
  
  # Check values are reasonable
  expect_true(design$Single_stage$n > 0)
  expect_true(design$Single_stage$ta > 0)
  expect_true(design$Single_stage$MTSL > design$Single_stage$ta)
  expect_true(design$Single_stage$E > 0)
  expect_true(design$Single_stage$P1 > 0 && design$Single_stage$P1 < 1)
})

test_that("single_stage_single_arm_tte validates input parameters", {
  expect_error(
    single_stage_single_arm_tte(
      shape = 1.0, S0 = 0.5, x0 = 0.5, hr = 0.667, tf = 1.0,
      rate = 32.4, alpha = 0.05, beta = 0.1, dist = "INVALID"
    ),
    "dist must be one of"
  )
  
  expect_error(
    single_stage_single_arm_tte(
      S0 = 0.5, x0 = 0.5, hr = 0.667, tf = 1.0,
      rate = 32.4, alpha = 0.05, beta = 0.1, dist = "WB"
    ),
    "shape, S0, and x0 are required"
  )
  
  expect_error(
    single_stage_single_arm_tte(
      shape = 1.0, S0 = 0.5, x0 = 0.5, hr = -0.5, tf = 1.0,
      rate = 32.4, alpha = 0.05, beta = 0.1, dist = "WB"
    ),
    "hr must be positive"
  )
})

test_that("single_stage_single_arm_tte works with different distributions", {
  # Log-normal
  design_ln <- single_stage_single_arm_tte(
    shape = 0.8,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    dist = "LN",
    restricted = FALSE
  )
  expect_s3_class(design_ln, "single_stage_design")
  
  # Log-logistic
  design_lg <- single_stage_single_arm_tte(
    shape = 1.5,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    dist = "LG",
    restricted = FALSE
  )
  expect_s3_class(design_lg, "single_stage_design")
  
  # Gamma
  design_gm <- single_stage_single_arm_tte(
    shape = 2.0,
    S0 = 0.5,
    x0 = 0.5,
    hr = 0.667,
    tf = 1.0,
    rate = 32.4,
    alpha = 0.05,
    beta = 0.1,
    dist = "GM",
    restricted = FALSE
  )
  expect_s3_class(design_gm, "single_stage_design")
})
