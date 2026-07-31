#' @describeIn simulate_trial Legacy two-stage TTE designs; delegates to
#'   [simulate_two_stage_tte()] (pass `shape` and `scale`).
#' @export
simulate_trial.two_stage_design <- function(design, ...) {
  simulate_two_stage_tte(design, ...)
}

#' @describeIn simulate_trial Legacy single-stage TTE designs; delegates to
#'   [simulate_single_stage_tte()] (pass `shape` and `scale`).
#' @export
simulate_trial.single_stage_design <- function(design, ...) {
  simulate_single_stage_tte(design, ...)
}

#' @describeIn simulate_trial Legacy Bayesian basket designs; delegates to
#'   [simulate_basket_trial()] (pass `p_true` and optionally
#'   `allocation_probs`).
#' @export
simulate_trial.basket_trial_design <- function(design, ...) {
  simulate_basket_trial(design, ...)
}
