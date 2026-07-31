#' Simulate Operating Characteristics of Two-Stage Design
#'
#' @description
#' Runs multiple simulations to estimate operating characteristics (Type I error,
#' power, expected sample size) for a two-stage single-arm trial design.
#' 
#' This function simulates many trial realizations under both the null hypothesis 
#' (H₀) and alternative hypothesis (H₁) to empirically estimate the design's 
#' operating characteristics. Results can be compared against the theoretical 
#' values to validate the design implementation.
#'
#' @param design A design object from `two_stage_single_arm_tte()` or similar
#'   blueprint design functions.
#' @param shape Numeric. Shape parameter for the survival distribution under 
#'   simulation (applies to both H₀ and H₁).
#' @param scale_h0 Numeric. Scale parameter under null hypothesis.
#' @param scale_h1 Numeric. Scale parameter under alternative hypothesis.
#' @param n_sims Integer. Number of simulations to run (default: 1000). Larger 
#'   values provide more precise estimates but take longer to compute.
#' @param seed Integer. Random seed for reproducibility (optional). If NULL, 
#'   no seed is set.
#'
#' @return A list containing:
#'   \item{type1_error}{Estimated Type I error rate (proportion of H₀ trials 
#'     that reject H₀)}
#'   \item{power}{Estimated power (proportion of H₁ trials that reject H₀)}
#'   \item{expected_n_h0}{Expected sample size under H₀ (mean enrolled patients)}
#'   \item{expected_n_h1}{Expected sample size under H₁ (mean enrolled patients)}
#'   \item{prob_early_stop_h0}{Probability of early stopping under H₀}
#'   \item{prob_early_stop_h1}{Probability of early stopping under H₁}
#'   \item{avg_interim_time_h0}{Average calendar time of interim analysis under H₀}
#'   \item{avg_interim_time_h1}{Average calendar time of interim analysis under H₁}
#'   \item{avg_final_time_h0}{Average calendar time of final analysis under H₀ 
#'     (excluding early-stopped trials)}
#'   \item{avg_final_time_h1}{Average calendar time of final analysis under H₁ 
#'     (excluding early-stopped trials)}
#'   \item{avg_interim_events_h0}{Average number of events at interim under H₀}
#'   \item{avg_interim_events_h1}{Average number of events at interim under H₁}
#'   \item{avg_final_events_h0}{Average number of events at final under H₀ 
#'     (excluding early-stopped trials)}
#'   \item{avg_final_events_h1}{Average number of events at final under H₁ 
#'     (excluding early-stopped trials)}
#'   \item{results_h0}{List of individual simulation results under H₀. Each 
#'     element is the output from `simulate_trial()`.}
#'   \item{results_h1}{List of individual simulation results under H₁. Each 
#'     element is the output from `simulate_trial()`.}
#'
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Validating that a design achieves its target Type I error and power
#'   \item Comparing simulated operating characteristics against theoretical values
#'   \item Understanding the distribution of sample sizes and analysis times
#'   \item Accessing individual trial data for detailed analysis
#' }
#' 
#' The function runs `n_sims` independent trials under each hypothesis. For each 
#' trial, it calls `simulate_trial()` which simulates patient accrual, event 
#' generation, and applies the two-stage decision rules.
#' 
#' **Computing Time**: Runtime scales linearly with `n_sims`. Typical values:
#' \itemize{
#'   \item 1,000 simulations: ~1-2 minutes
#'   \item 10,000 simulations: ~10-20 minutes
#'   \item 100,000 simulations: ~1-2 hours
#' }
#'
#' @examples
#' \dontrun{
#' # Create a design
#' design <- two_stage_single_arm_tte(
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5, hr = 0.5913,
#'   tf = 5, rate = 2, alpha = 0.05, beta = 0.2,
#'   dist = "WB", restricted = FALSE
#' )
#' 
#' # Calculate scale parameters
#' # Under H0: median = 3.5 months
#' scale_h0 <- 3.5 / (log(2)^(1/1.47327))
#' # Under H1: median = 5 months (hazard ratio = 0.5913)
#' scale_h1 <- 5 / (log(2)^(1/1.47327))
#' 
#' # Run simulations
#' sim_results <- simulate_operating_characteristics(
#'   design = design,
#'   shape = 1.47327,
#'   scale_h0 = scale_h0,
#'   scale_h1 = scale_h1,
#'   n_sims = 1000,
#'   seed = 123
#' )
#' 
#' # Check operating characteristics
#' cat("Type I Error:", sim_results$type1_error, "\n")
#' cat("Power:", sim_results$power, "\n")
#' cat("Expected N | H0:", sim_results$expected_n_h0, "\n")
#' cat("P(Early Stop | H0):", sim_results$prob_early_stop_h0, "\n")
#' 
#' # Access individual trial data
#' first_trial_h0 <- sim_results$results_h0[[1]]
#' head(first_trial_h0$interim_data)
#' }
#'
#' @seealso 
#' \code{\link{simulate_trial}} for simulating a single trial realization.
#' \code{\link{two_stage_single_arm_tte}} for creating two-stage designs.
#'
#' @export
simulate_operating_characteristics <- function(design, shape,
                                                scale_h0, scale_h1,
                                                n_sims = 1000,
                                                seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Simulate under H0
  results_h0 <- replicate(n_sims, {
    simulate_trial(design = design, shape = shape, scale = scale_h0)
  }, simplify = FALSE)
  
  # Simulate under H1
  results_h1 <- replicate(n_sims, {
    simulate_trial(design = design, shape = shape, scale = scale_h1)
  }, simplify = FALSE)
  
  # Calculate operating characteristics
  type1_error <- mean(sapply(results_h0, function(x) x$reject_h0))
  power <- mean(sapply(results_h1, function(x) x$reject_h0))
  expected_n_h0 <- mean(sapply(results_h0, function(x) x$enrolled_n))
  expected_n_h1 <- mean(sapply(results_h1, function(x) x$enrolled_n))
  prob_early_stop_h0 <- mean(sapply(results_h0, function(x) x$stopped_early))
  prob_early_stop_h1 <- mean(sapply(results_h1, function(x) x$stopped_early))
  
  # Calculate average analysis times and event counts
  avg_interim_time_h0 <- mean(sapply(results_h0, function(x) x$interim_time))
  avg_interim_time_h1 <- mean(sapply(results_h1, function(x) x$interim_time))
  avg_final_time_h0 <- mean(sapply(results_h0, function(x) ifelse(is.na(x$final_time), NA, x$final_time)), na.rm = TRUE)
  avg_final_time_h1 <- mean(sapply(results_h1, function(x) ifelse(is.na(x$final_time), NA, x$final_time)), na.rm = TRUE)
  
  avg_interim_events_h0 <- mean(sapply(results_h0, function(x) x$interim_events))
  avg_interim_events_h1 <- mean(sapply(results_h1, function(x) x$interim_events))
  avg_final_events_h0 <- mean(sapply(results_h0, function(x) ifelse(is.na(x$final_events), NA, x$final_events)), na.rm = TRUE)
  avg_final_events_h1 <- mean(sapply(results_h1, function(x) ifelse(is.na(x$final_events), NA, x$final_events)), na.rm = TRUE)
  
  list(
    type1_error = type1_error,
    power = power,
    expected_n_h0 = expected_n_h0,
    expected_n_h1 = expected_n_h1,
    prob_early_stop_h0 = prob_early_stop_h0,
    prob_early_stop_h1 = prob_early_stop_h1,
    avg_interim_time_h0 = avg_interim_time_h0,
    avg_interim_time_h1 = avg_interim_time_h1,
    avg_final_time_h0 = avg_final_time_h0,
    avg_final_time_h1 = avg_final_time_h1,
    avg_interim_events_h0 = avg_interim_events_h0,
    avg_interim_events_h1 = avg_interim_events_h1,
    avg_final_events_h0 = avg_final_events_h0,
    avg_final_events_h1 = avg_final_events_h1,
    results_h0 = results_h0,
    results_h1 = results_h1
  )
}
