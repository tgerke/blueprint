# Extended validation tests replicate published results (e.g. Cunanan et al.
# 2017 Tables 2-4) via large grid searches and are far too slow for routine
# R CMD check. Opt in with Sys.setenv(BLUEPRINT_EXTENDED_TESTS = "true").
skip_extended <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("BLUEPRINT_EXTENDED_TESTS"), "true"),
    "Extended validation tests run only when BLUEPRINT_EXTENDED_TESTS=true"
  )
}
