test_that("assume_response() validates probabilities", {
  expect_s3_class(assume_response(0.3), c("response_scenario", "scenario"))
  expect_equal(assume_response(c(0.2, 0.4))$p, c(0.2, 0.4))
  expect_error(assume_response(0), "strictly between 0 and 1")
  expect_error(assume_response(1.2), "strictly between 0 and 1")
  expect_error(assume_response("a"), "strictly between 0 and 1")
})

test_that("assume_survival() validates parameters", {
  scenario <- assume_survival("weibull", shape = 1.5, scale = 4.2)
  expect_s3_class(scenario, c("survival_scenario", "scenario"))
  expect_equal(scenario$dist, "weibull")
  expect_error(assume_survival("weibull", shape = -1, scale = 2),
               "`shape` must be a single positive number")
  expect_error(assume_survival("weibull", shape = 1, scale = 0),
               "`scale` must be a single positive number")
  expect_error(assume_survival("cauchy", shape = 1, scale = 1))
})

test_that("scenario print methods are stable", {
  expect_snapshot(print(assume_response(0.3)))
  expect_snapshot(print(assume_survival("weibull", shape = 1.5, scale = 4.2)))
})
