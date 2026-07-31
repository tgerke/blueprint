test_that("logrank_stat_wu calculates correct statistic for Weibull", {
  # Simple test case with known values
  set.seed(123)
  
  # Simulate some data under Weibull with shape=1.5, scale=2
  shape <- 1.5
  scale <- 2
  n <- 50
  
  # Generate event times from Weibull
  event_times <- rweibull(n, shape = shape, scale = scale)
  obs_times <- pmin(event_times, 3)  # Censor at time 3
  event_ind <- event_times <= 3
  
  # Calculate statistic
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "WB",
    shape_h0 = shape,
    scale_h0 = scale
  )
  
  # Should return a numeric value
  expect_type(Z, "double")
  expect_length(Z, 1)
  expect_false(is.na(Z))
  expect_false(is.infinite(Z))
  
  # When data matches H0 exactly, Z should be close to 0
  # (though with randomness it won't be exactly 0)
  expect_true(abs(Z) < 5)  # Reasonable range for Z-score
})

test_that("logrank_stat_wu works with log-normal distribution", {
  set.seed(456)
  
  # Log-normal parameters
  meanlog <- 0.5
  sdlog <- 0.8
  n <- 40
  
  # Generate data
  event_times <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  obs_times <- pmin(event_times, 5)
  event_ind <- event_times <= 5
  
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "LN",
    shape_h0 = sdlog,
    scale_h0 = meanlog
  )
  
  expect_type(Z, "double")
  expect_length(Z, 1)
  expect_false(is.na(Z))
})

test_that("logrank_stat_wu works with log-logistic distribution", {
  set.seed(789)
  
  shape <- 2
  scale <- 1.5
  n <- 30
  
  # For log-logistic: generate using flexsurv approach
  u <- runif(n)
  event_times <- scale * ((1 - u)^(-1/shape) - 1)^(1/shape)
  obs_times <- pmin(event_times, 4)
  event_ind <- event_times <= 4
  
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "LG",
    shape_h0 = shape,
    scale_h0 = scale
  )
  
  expect_type(Z, "double")
  expect_length(Z, 1)
  expect_false(is.na(Z))
})

test_that("logrank_stat_wu works with gamma distribution", {
  set.seed(101)
  
  shape <- 2.5
  scale <- 1.2
  n <- 35
  
  # Generate gamma data
  event_times <- rgamma(n, shape = shape, scale = scale)
  obs_times <- pmin(event_times, 6)
  event_ind <- event_times <= 6
  
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "GM",
    shape_h0 = shape,
    scale_h0 = scale
  )
  
  expect_type(Z, "double")
  expect_length(Z, 1)
  expect_false(is.na(Z))
})

test_that("logrank_stat_wu has correct sign convention", {
  # Fewer events than expected should give positive Z
  obs_times <- c(1, 2, 3, 4, 5)
  event_ind <- c(FALSE, FALSE, FALSE, FALSE, FALSE)  # All censored
  
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "WB",
    shape_h0 = 1,  # Exponential
    scale_h0 = 2
  )
  
  # No events, but some expected -> positive Z favoring H1
  expect_true(Z > 0)
  
  # More events than expected under very good survival should give negative Z
  obs_times2 <- rep(0.01, 10)  # Events very early
  event_ind2 <- rep(TRUE, 10)
  
  Z2 <- logrank_stat_wu(
    obs_times = obs_times2,
    event_ind = event_ind2,
    dist = "WB",
    shape_h0 = 1,
    scale_h0 = 100  # Excellent survival under H0
  )
  
  # Many events when few expected -> negative Z
  expect_true(Z2 < 0)
})

test_that("logrank_stat_wu handles edge cases", {
  # All events at same time
  Z1 <- logrank_stat_wu(
    obs_times = rep(2, 10),
    event_ind = rep(TRUE, 10),
    dist = "WB",
    shape_h0 = 1.5,
    scale_h0 = 2
  )
  expect_type(Z1, "double")
  expect_false(is.na(Z1))
  
  # Single observation
  Z2 <- logrank_stat_wu(
    obs_times = 1.5,
    event_ind = TRUE,
    dist = "WB",
    shape_h0 = 1,
    scale_h0 = 1
  )
  expect_type(Z2, "double")
  expect_false(is.na(Z2))
  
  # Very small observation times (near 0)
  Z3 <- logrank_stat_wu(
    obs_times = c(0.001, 0.002, 0.003),
    event_ind = c(TRUE, TRUE, FALSE),
    dist = "WB",
    shape_h0 = 1,
    scale_h0 = 2
  )
  expect_type(Z3, "double")
  expect_false(is.na(Z3))
})

test_that("logrank_stat_wu gives error for unsupported distribution", {
  expect_error(
    logrank_stat_wu(
      obs_times = c(1, 2, 3),
      event_ind = c(TRUE, FALSE, TRUE),
      dist = "UNSUPPORTED",
      shape_h0 = 1,
      scale_h0 = 1
    ),
    "Unsupported distribution"
  )
})

test_that("logrank_stat_wu matches expected calculation", {
  # Manual calculation for simple case
  obs_times <- c(1, 2, 3)
  event_ind <- c(TRUE, TRUE, FALSE)
  
  # Weibull with shape=1, scale=2 (exponential with rate=1/2)
  # H(t) = (t/2)^1 = t/2
  # Expected events: H(1) + H(2) + H(3) = 0.5 + 1.0 + 1.5 = 3
  # Observed events: 2
  # Z = (3 - 2) / sqrt(3) ≈ 0.577
  
  Z <- logrank_stat_wu(
    obs_times = obs_times,
    event_ind = event_ind,
    dist = "WB",
    shape_h0 = 1,
    scale_h0 = 2
  )
  
  expected_Z <- (3 - 2) / sqrt(3)
  expect_equal(Z, expected_Z, tolerance = 1e-10)
})
