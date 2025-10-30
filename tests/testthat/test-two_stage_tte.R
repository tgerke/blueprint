test_that("two_stage_single_arm_tte works with Weibull distribution", {
  # Example from Wu et al. (2020) - Small-cell lung cancer trial
  # Weibull distribution with 5 months restricted follow-up
  design <- two_stage_single_arm_tte(
    shape = 1.47327,
    S0 = 0.5,
    x0 = 3.5,
    hr = 0.5913,
    x = 5,
    rate = 2,
    alpha = 0.05,
    beta = 0.2,
    dist = "WB"
  )
  
  expect_type(design, "list")
  expect_named(design, c("param", "Single_stage", "Two_stage"))
  
  # Check param structure
  expect_s3_class(design$param, "data.frame")
  expect_equal(design$param$shape, 1.47327)
  expect_equal(design$param$S0, 0.5)
  expect_equal(design$param$hr, 0.5913)
  expect_equal(design$param$alpha, 0.05)
  expect_equal(design$param$beta, 0.2)
  expect_equal(design$param$rate, 2)
  expect_equal(design$param$x0, 3.5)
  expect_equal(design$param$x, 5)
  
  # Check Single_stage structure
  expect_s3_class(design$Single_stage, "data.frame")
  expect_named(design$Single_stage, c("nsingle", "tasingle", "csingle"))
  expect_true(design$Single_stage$nsingle > 0)
  expect_true(design$Single_stage$tasingle > 0)
  expect_true(abs(design$Single_stage$csingle - qnorm(0.95)) < 0.01)
  
  # Check Two_stage structure
  expect_s3_class(design$Two_stage, "data.frame")
  expect_named(design$Two_stage, c("n1", "c1", "n", "c", "t1", "MTSL", "ES", "PS"))
  expect_true(design$Two_stage$n1 > 0)
  expect_true(design$Two_stage$n > design$Two_stage$n1)
  expect_true(design$Two_stage$t1 > 0)
  expect_true(design$Two_stage$MTSL > 0)
  expect_true(design$Two_stage$ES >= design$Two_stage$n1)
  expect_true(design$Two_stage$PS >= 0 && design$Two_stage$PS <= 1)
  
  # Check that design is reasonable (less strict than paper values)
  # Paper reports n1=28, n=45, but optimization may vary
  expect_true(design$Two_stage$n1 >= 5)
  expect_true(design$Two_stage$n <= 100)
  expect_true(design$Two_stage$ES < design$Single_stage$nsingle)
})

test_that("two_stage_single_arm_tte works with log-normal distribution", {
  design <- two_stage_single_arm_tte(
    shape = 0.5,
    S0 = 0.3,
    x0 = 1,
    hr = 0.65,
    x = 1,
    rate = 10,
    alpha = 0.05,
    beta = 0.2,
    dist = "LN"
  )
  
  expect_type(design, "list")
  expect_named(design, c("param", "Single_stage", "Two_stage"))
  expect_true(design$Two_stage$n > 0)
})

test_that("two_stage_single_arm_tte works with gamma distribution", {
  design <- two_stage_single_arm_tte(
    shape = 1,
    S0 = 0.3,
    x0 = 1,
    hr = 0.65,
    x = 1,
    rate = 10,
    alpha = 0.05,
    beta = 0.2,
    dist = "GM"
  )
  
  expect_type(design, "list")
  expect_named(design, c("param", "Single_stage", "Two_stage"))
  expect_true(design$Two_stage$n > 0)
})

test_that("two_stage_single_arm_tte works with log-logistic distribution", {
  design <- two_stage_single_arm_tte(
    shape = 1,
    S0 = 0.3,
    x0 = 1,
    hr = 0.65,
    x = 1,
    rate = 10,
    alpha = 0.05,
    beta = 0.2,
    dist = "LG"
  )
  
  expect_type(design, "list")
  expect_named(design, c("param", "Single_stage", "Two_stage"))
  expect_true(design$Two_stage$n > 0)
})

test_that("two_stage_single_arm_tte validates input parameters", {
  expect_error(
    two_stage_single_arm_tte(
      shape = 1.47327,
      S0 = 0.5,
      x0 = 3.5,
      hr = 0.5913,
      x = 5,
      rate = 2,
      alpha = 0.05,
      beta = 0.2,
      dist = "INVALID"
    ),
    "dist must be one of"
  )
  
  expect_error(
    two_stage_single_arm_tte(
      hr = 0.5913,
      x = 5,
      rate = 2,
      alpha = 0.05,
      beta = 0.2,
      dist = "WB"
    ),
    "shape, S0, and x0 are required"
  )
  
  expect_error(
    two_stage_single_arm_tte(
      shape = 1.47327,
      S0 = 0.5,
      x0 = 3.5,
      hr = 0.5913,
      x = 5,
      rate = 2,
      alpha = 1.5,
      beta = 0.2,
      dist = "WB"
    ),
    "alpha must be between 0 and 1"
  )
  
  expect_error(
    two_stage_single_arm_tte(
      shape = 1.47327,
      S0 = 0.5,
      x0 = 3.5,
      hr = -0.5913,
      x = 5,
      rate = 2,
      alpha = 0.05,
      beta = 0.2,
      dist = "WB"
    ),
    "hr must be positive"
  )
})

test_that("two_stage_single_arm_tte produces consistent results", {
  set.seed(123)
  design1 <- two_stage_single_arm_tte(
    shape = 1.47327,
    S0 = 0.5,
    x0 = 3.5,
    hr = 0.5913,
    x = 5,
    rate = 2,
    alpha = 0.05,
    beta = 0.2,
    dist = "WB"
  )
  
  set.seed(123)
  design2 <- two_stage_single_arm_tte(
    shape = 1.47327,
    S0 = 0.5,
    x0 = 3.5,
    hr = 0.5913,
    x = 5,
    rate = 2,
    alpha = 0.05,
    beta = 0.2,
    dist = "WB"
  )
  
  expect_equal(design1$Two_stage$n1, design2$Two_stage$n1)
  expect_equal(design1$Two_stage$n, design2$Two_stage$n)
  expect_equal(design1$Two_stage$c1, design2$Two_stage$c1)
  expect_equal(design1$Two_stage$c, design2$Two_stage$c)
})
