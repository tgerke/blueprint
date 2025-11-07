#' Simulate a Trial from a Design Object
#'
#' @description
#' Generic function that simulates a single realization of a trial based on a 
#' design object. This wrapper automatically dispatches to the appropriate 
#' simulation function based on the design class:
#' 
#' \itemize{
#'   \item \code{two_stage_design}: Dispatches to \code{\link{simulate_two_stage_tte}}
#'   \item \code{single_stage_design}: Dispatches to \code{\link{simulate_single_stage_tte}}
#'   \item \code{basket_trial_design}: Dispatches to \code{\link{simulate_basket_trial}}
#' }
#'
#' @param design A design object from any blueprint design function
#' @param ... Additional arguments passed to the specific simulation function.
#'   For TTE designs: \code{shape} and \code{scale}.
#'   For basket trials: \code{p_true} and optionally \code{allocation_probs}.
#'
#' @return A list with simulation results. The exact structure depends on the 
#'   design type. See the specific simulation functions for details.
#'
#' @examples
#' \dontrun{
#' # Two-stage TTE design
#' design_two <- two_stage_single_arm_tte(
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5, hr = 0.5913,
#'   tf = 5, rate = 2, alpha = 0.05, beta = 0.2,
#'   dist = "WB", restricted = FALSE
#' )
#' scale_h0 <- 3.5 / (-log(0.5))^(1/1.47327)
#' trial_result <- simulate_trial(design_two, shape = 1.47327, scale = scale_h0)
#' 
#' # Single-stage TTE design
#' design_single <- single_stage_single_arm_tte(
#'   shape = 1.0, S0 = 0.5, x0 = 0.5, hr = 0.667,
#'   tf = 1.0, rate = 32.4, alpha = 0.05, beta = 0.10,
#'   two_sided = TRUE, dist = "WB", restricted = FALSE
#' )
#' scale_h1 <- scale_h0 / (0.667^(1/1.0))
#' trial_result <- simulate_trial(design_single, shape = 1.0, scale = scale_h1)
#' 
#' # Basket trial design
#' design_basket <- two_stage_basket_trial(
#'   K = 3, p0 = 0.10, pa = 0.25,
#'   alpha = 0.05, beta = 0.2,
#'   n_sims_design = 1000, seed = 123
#' )
#' trial_result <- simulate_trial(design_basket, p_true = c(0.10, 0.10, 0.10))
#' }
#'
#' @seealso 
#' \code{\link{simulate_two_stage_tte}}, 
#' \code{\link{simulate_single_stage_tte}},
#' \code{\link{simulate_basket_trial}}
#'
#' @export
simulate_trial <- function(design, ...) {
  
  # Dispatch based on design class
  if (inherits(design, "two_stage_design")) {
    return(simulate_two_stage_tte(design, ...))
  } else if (inherits(design, "single_stage_design")) {
    return(simulate_single_stage_tte(design, ...))
  } else if (inherits(design, "basket_trial_design")) {
    return(simulate_basket_trial(design, ...))
  } else {
    stop(
      "Unknown design class: ", class(design)[1], "\n",
      "simulate_trial() supports: two_stage_design, single_stage_design, basket_trial_design"
    )
  }
}
