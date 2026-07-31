# Contract tests: every modern design class must satisfy the shared
# trial_design contract. Add new design families to the list below.

modern_designs <- function() {
  list(
    simon = design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
  )
}

test_that("designs satisfy the trial_design structural contract", {
  for (design in modern_designs()) {
    expect_true(is_trial_design(design))
    expect_named(design, c("design", "hypotheses", "endpoint", "citation"),
                 ignore.order = TRUE, ignore.case = FALSE)
    expect_type(design$design, "list")
    expect_type(design$hypotheses, "list")
    expect_true(design$endpoint %in% c("binary", "tte"))
  }
})

test_that("designs support the scenario helpers", {
  for (design in modern_designs()) {
    expect_s3_class(assume_null(design), "scenario")
    expect_s3_class(assume_alternative(design), "scenario")
  }
})

test_that("designs support simulate_trial() with required columns", {
  for (design in modern_designs()) {
    sims <- simulate_trial(design, under = assume_null(design),
                           n_sims = 20, seed = 1)
    expect_s3_class(sims, "tbl_df")
    expect_equal(nrow(sims), 20)
    expect_contains(names(sims), c("sim", "n_enrolled", "success"))
  }
})

test_that("designs support evaluate() for exact and simulation methods", {
  for (design in modern_designs()) {
    for (method in c("exact", "simulation")) {
      ocs <- evaluate(design, method = method, n_sims = 100, seed = 1)
      expect_s3_class(ocs, "tbl_df")
      expect_contains(names(ocs),
                      c("scenario", "prob_success", "expected_n", "method"))
      expect_contains(ocs$scenario, c("null", "alternative"))
    }
  }
})

test_that("designs support verify()", {
  for (design in modern_designs()) {
    verification <- verify(design, n_sims = 200, seed = 1)
    expect_s3_class(verification, "trial_verification")
  }
})

test_that("generics fail informatively for unknown classes", {
  expect_error(simulate_trial(structure(list(), class = "mystery")),
               "No `simulate_trial\\(\\)` method")
  expect_error(evaluate(structure(list(), class = "mystery")),
               "No `evaluate\\(\\)` method")
  expect_error(verify(structure(list(), class = "mystery")),
               "No `verify\\(\\)` method")
  expect_error(draft_protocol_text(structure(list(), class = "mystery")),
               "No `draft_protocol_text\\(\\)` method")
  expect_error(assume_null(structure(list(), class = "mystery")),
               "No `assume_null\\(\\)` method")
})

test_that("contract: modern_designs() covers exactly one entry per modern class", {
  # When a new design family lands, register a representative design above so
  # it runs through the whole contract suite.
  expect_named(modern_designs(), "simon")
})
