# Reference values cross-validated against clinfun::ph2simon() and
# Simon (1989); see dev-docs/adr-001-design-grammar.md.

test_that("design_simon() validates inputs", {
  expect_error(design_simon(n1 = 13, r1 = 13, n = 43, r = 12, p0 = 0.2, pa = 0.4),
               "`r1` must be less than `n1`")
  expect_error(design_simon(n1 = 13, r1 = 3, n = 13, r = 12, p0 = 0.2, pa = 0.4),
               "`n` must be greater than `n1`")
  expect_error(design_simon(n1 = 13, r1 = 3, n = 43, r = 43, p0 = 0.2, pa = 0.4),
               "`r` must satisfy")
  expect_error(design_simon(n1 = 13, r1 = 3, n = 43, r = 2, p0 = 0.2, pa = 0.4),
               "`r` must satisfy")
  expect_error(design_simon(n1 = 13.5, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4),
               "whole number")
  expect_error(design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.4, pa = 0.2),
               "`pa` must be greater than `p0`")
  expect_error(design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 1.4),
               "strictly between 0 and 1")
})

test_that("exact OCs match published Simon (1989) values: p0=0.2, pa=0.4 optimal", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  ocs <- evaluate(design, method = "exact")

  expect_s3_class(ocs, "tbl_df")
  expect_equal(ocs$scenario, c("null", "alternative"))
  null_row <- ocs[ocs$scenario == "null", ]
  alt_row <- ocs[ocs$scenario == "alternative", ]

  expect_equal(null_row$prob_success, 0.0496, tolerance = 1e-3)
  expect_equal(alt_row$prob_success, 0.8002, tolerance = 1e-3)
  expect_equal(null_row$prob_early_stop, 0.7473, tolerance = 1e-3)
  expect_equal(null_row$expected_n, 20.58, tolerance = 1e-3)
})

test_that("exact OCs match published Simon (1989) values: p0=0.05, pa=0.25 optimal", {
  design <- design_simon(n1 = 9, r1 = 0, n = 17, r = 2, p0 = 0.05, pa = 0.25)
  ocs <- evaluate(design, method = "exact")
  null_row <- ocs[ocs$scenario == "null", ]

  expect_equal(null_row$prob_early_stop, 0.6302, tolerance = 1e-3)
  expect_equal(null_row$expected_n, 11.96, tolerance = 1e-2)
  expect_lte(null_row$prob_success, 0.05)
  expect_gte(ocs$prob_success[ocs$scenario == "alternative"], 0.80)
})

test_that("search_simon_designs() reproduces published optimal and minimax designs", {
  candidates <- search_simon_designs(p0 = 0.2, pa = 0.4, alpha = 0.05,
                                     power = 0.8, n_max = 45)

  optimal <- candidates[grepl("optimal", candidates$criterion), ]
  expect_equal(
    unlist(optimal[, c("r1", "n1", "r", "n")]),
    c(r1 = 3, n1 = 13, r = 12, n = 43)
  )

  minimax <- candidates[grepl("minimax", candidates$criterion), ]
  expect_equal(
    unlist(minimax[, c("r1", "n1", "r", "n")]),
    c(r1 = 4, n1 = 18, r = 10, n = 33)
  )

  # every candidate meets the constraints
  expect_true(all(candidates$type1_error <= 0.05))
  expect_true(all(candidates$power >= 0.8))
})

test_that("search_simon_designs() errors when no design fits", {
  expect_error(
    search_simon_designs(p0 = 0.2, pa = 0.25, alpha = 0.05, power = 0.8,
                         n_max = 30),
    "Increase `n_max`"
  )
})

test_that("pick_design() constructs the flagged design with metadata", {
  candidates <- search_simon_designs(p0 = 0.05, pa = 0.25, n_max = 20)
  design <- pick_design(candidates, "minimax")

  expect_s3_class(design, "simon_design")
  expect_true(is_trial_design(design))
  expect_equal(design$design, list(n1 = 12L, r1 = 0L, n = 16L, r = 2L))
  expect_equal(design$criterion, "minimax")
  expect_equal(design$hypotheses$alpha, 0.05)
  expect_equal(design$hypotheses$power, 0.80)
})

test_that("simulate_trial() returns a coherent tidy tibble", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  sims <- simulate_trial(design, under = assume_alternative(design),
                         n_sims = 500, seed = 42)

  expect_s3_class(sims, "tbl_df")
  expect_equal(nrow(sims), 500)
  expect_named(sims, c("sim", "responses_stage1", "stopped_early",
                       "n_enrolled", "responses_total", "success"))
  expect_setequal(unique(sims$n_enrolled), c(13L, 43L))
  expect_equal(sims$stopped_early, sims$responses_stage1 <= 3)
  expect_equal(sims$n_enrolled == 13L, sims$stopped_early)
  expect_false(any(sims$success & sims$stopped_early))
  expect_equal(sims$success, !sims$stopped_early & sims$responses_total > 12)

  # reproducible with a seed
  sims2 <- simulate_trial(design, under = assume_alternative(design),
                          n_sims = 500, seed = 42)
  expect_identical(sims, sims2)
})

test_that("simulate_trial() requires a valid scenario", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  expect_error(simulate_trial(design, n_sims = 5), "Supply a scenario")
  expect_error(simulate_trial(design, under = list(p = 0.3)),
               "must be a response scenario")
  expect_error(simulate_trial(design, under = assume_response(c(0.3, 0.4))),
               "single-arm")
})

test_that("evaluate() by simulation agrees with exact values", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  exact <- evaluate(design, method = "exact")
  simulated <- evaluate(design, method = "simulation", n_sims = 5000,
                        seed = 123)

  expect_equal(simulated$prob_success, exact$prob_success, tolerance = 0.05)
  expect_equal(simulated$prob_early_stop, exact$prob_early_stop,
               tolerance = 0.05)
  expect_equal(simulated$expected_n, exact$expected_n, tolerance = 0.05)
  expect_equal(unique(simulated$method), "simulation")
})

test_that("evaluate() accepts custom scenarios", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)

  one <- evaluate(design, under = assume_response(0.3))
  expect_equal(one$scenario, "custom")
  expect_equal(one$p_true, 0.3)

  several <- evaluate(design, under = list(
    pessimistic = assume_response(0.25),
    optimistic = assume_response(0.45)
  ))
  expect_equal(several$scenario, c("pessimistic", "optimistic"))
  # rejection probability increases with the true response rate
  expect_lt(several$prob_success[1], several$prob_success[2])
})

test_that("verify() compares exact and simulated OCs", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  verification <- verify(design, n_sims = 4000, seed = 99)

  expect_s3_class(verification, "trial_verification")
  expect_named(verification, c("scenario", "metric", "exact", "simulated",
                               "difference", "mc_se"))
  expect_setequal(unique(verification$scenario), c("null", "alternative"))
  expect_true(all(abs(verification$difference) < 0.5))

  prob_rows <- verification[startsWith(verification$metric, "prob_"), ]
  expect_true(all(abs(prob_rows$difference) < 0.03))
  expect_false(any(is.na(prob_rows$mc_se)))
})

test_that("tidy() and glance() return tidy summaries", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4,
                         criterion = "optimal")

  tidied <- tidy(design)
  expect_s3_class(tidied, "tbl_df")
  expect_equal(nrow(tidied), 2)

  glanced <- glance(design)
  expect_equal(nrow(glanced), 1)
  expect_equal(glanced$n, 43L)
  expect_equal(glanced$criterion, "optimal")
  expect_equal(glanced$type1_error, 0.0496, tolerance = 1e-3)
  expect_equal(glanced$power, 0.8002, tolerance = 1e-3)
})

test_that("print and protocol text output are stable", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4,
                         criterion = "optimal")
  expect_snapshot(print(design))
  expect_snapshot(cat(draft_protocol_text(design)))
})

test_that("assume_null() and assume_alternative() carry the design hypotheses", {
  design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  expect_equal(assume_null(design)$p, 0.2)
  expect_equal(assume_alternative(design)$p, 0.4)
  expect_s3_class(assume_null(design), "response_scenario")
})
