#' Simulate a Single-Stage Time-to-Event Trial
#'
#' @description
#' Simulates a single realization of a single-stage time-to-event trial based on a 
#' design object created by `single_stage_single_arm_tte()`.
#'
#' @param design A single-stage design object from `single_stage_single_arm_tte()`
#' @param shape Shape parameter for the survival distribution under simulation
#' @param scale Scale parameter for the survival distribution under simulation
#'
#' @return A list containing:
#'   \item{reject_h0}{Logical; TRUE if null hypothesis is rejected}
#'   \item{stopped_early}{Logical; FALSE (single-stage doesn't have interim)}
#'   \item{enrolled_n}{Number of patients enrolled}
#'   \item{interim_data}{NULL (no interim for single-stage)}
#'   \item{final_data}{Data frame of patient data at final analysis}
#'   \item{interim_stat}{NA (no interim for single-stage)}
#'   \item{final_stat}{Test statistic value at final}
#'   \item{interim_time}{NA (no interim for single-stage)}
#'   \item{final_time}{Calendar time of final analysis}
#'   \item{interim_events}{NA (no interim for single-stage)}
#'   \item{final_events}{Number of events observed at final}
#'
#' @examples
#' \dontrun{
#' design <- single_stage_single_arm_tte(
#'   shape = 1.0, S0 = 0.5, x0 = 0.5, hr = 0.667,
#'   tf = 1.0, rate = 32.4, alpha = 0.05, beta = 0.10,
#'   two_sided = TRUE, dist = "WB", restricted = FALSE
#' )
#' 
#' # Simulate one trial under H1
#' scale_h0 <- 0.5 / (-log(0.5))^(1/1.0)
#' scale_h1 <- scale_h0 / (0.667^(1/1.0))
#' trial_result <- simulate_single_stage_tte(design, shape = 1.0, scale = scale_h1)
#' }
#'
#' @export
simulate_single_stage_tte <- function(design, shape, scale) {
  
  # Validate design object
  if (!inherits(design, "single_stage_design")) {
    stop("design must be a single_stage_design object from single_stage_single_arm_tte()")
  }
  
  # Extract parameters from single-stage design
  n <- design$single_stage$n
  ta <- design$single_stage$ta
  c <- design$single_stage$c  # Critical value
  
  # Extract from param
  tf <- design$param$tf
  rate <- design$param$rate
  restricted <- design$param$restricted
  two_sided <- design$param$two_sided
  
  # Determine distribution
  dist <- if ("dist" %in% names(design$param)) design$param$dist else "WB"
  
  # Extract H0 parameters for log-rank test
  S0 <- design$param$S0
  x0 <- design$param$x0
  shape_h0 <- design$param$shape
  
  # Calculate H0 scale parameter
  if (!is.null(S0) && !is.null(x0) && !is.null(shape_h0)) {
    scale_h0 <- .calculate_h0_scale(dist, S0, x0, shape_h0)
  } else {
    scale_h0 <- NULL
  }
  
  # Generate accrual times for all n patients
  accrual_period <- n / rate
  accrual_times <- sort(stats::runif(n, 0, accrual_period))
  
  # Generate event times based on distribution
  event_times <- switch(dist,
    "WB" = stats::rweibull(n, shape = shape, scale = scale),
    "LN" = stats::rlnorm(n, meanlog = log(scale), sdlog = shape),
    "GM" = stats::rgamma(n, shape = shape, scale = scale),
    "LG" = {
      # Log-logistic: use inverse CDF method
      u <- stats::runif(n)
      scale * ((u / (1 - u))^(1/shape))
    },
    stop("Unsupported distribution. Use 'WB', 'LN', 'GM', or 'LG'.")
  )
  
  # === FINAL ANALYSIS (Single-stage) ===
  # Time of final analysis
  if (restricted) {
    # Restricted: everyone followed for exactly tf
    t_final <- max(accrual_times) + tf
  } else {
    # Unrestricted: last patient gets tf follow-up
    t_final <- max(accrual_times) + tf
  }
  
  # Calculate observed times and event indicators at final
  final_obs_time <- pmin(event_times, t_final - accrual_times)
  final_event_ind <- event_times <= (t_final - accrual_times)
  
  final_data <- data.frame(
    patient_id = 1:n,
    accrual_time = accrual_times,
    event_time = event_times,
    obs_time = final_obs_time,
    event = final_event_ind
  )
  
  # Calculate one-sample log-rank statistic at final
  if (!is.null(scale_h0)) {
    final_stat <- logrank_stat_wu(
      obs_times = final_obs_time,
      event_ind = final_event_ind,
      dist = dist,
      shape_h0 = shape_h0,
      scale_h0 = scale_h0
    )
  } else {
    # Fallback to simple standardized count
    O_final <- sum(final_event_ind)
    E_final <- n * mean(final_event_ind)
    V_final <- n * stats::var(final_event_ind)
    
    if (V_final == 0 || is.na(V_final)) {
      final_stat <- O_final - E_final
    } else {
      final_stat <- (O_final - E_final) / sqrt(V_final)
    }
  }
  
  # Decision at final
  # For two-sided test: reject if |Z| > c (i.e., Z > c OR Z < -c)
  # For one-sided test: reject if Z > c
  if (two_sided) {
    reject <- !is.na(final_stat) && abs(final_stat) > c
  } else {
    reject <- !is.na(final_stat) && final_stat > c
  }
  
  return(list(
    reject_h0 = reject,
    stopped_early = FALSE,  # Single-stage doesn't have interim
    enrolled_n = n,
    interim_data = NULL,  # No interim for single-stage
    final_data = final_data,
    interim_stat = NA,  # No interim for single-stage
    final_stat = final_stat,
    interim_time = NA,  # No interim for single-stage
    final_time = t_final,
    interim_events = NA,  # No interim for single-stage
    final_events = sum(final_event_ind)
  ))
}


#' Calculate H0 Scale Parameter from S0 and x0
#'
#' @description
#' Internal helper function to calculate the scale parameter for the null hypothesis
#' distribution given survival probability S0 at time x0.
#'
#' @keywords internal
#' @noRd
.calculate_h0_scale <- function(dist, S0, x0, shape) {
  if (dist == "WB") {
    # Weibull: S(t) = exp(-(t/scale)^shape)
    # S0 = exp(-(x0/scale)^shape)
    # scale = x0 / (-log(S0))^(1/shape)
    scale <- x0 / (-log(S0))^(1/shape)
  } else if (dist == "LN") {
    # Log-normal: S(t) = 1 - Phi((log(t) - meanlog) / sdlog)
    # meanlog = scale (in our parameterization)
    # sdlog = shape
    scale <- log(x0) - shape * stats::qnorm(1 - S0)
  } else if (dist == "LG") {
    # Log-logistic: S(t) = 1 / (1 + (t/scale)^shape)
    # S0 = 1 / (1 + (x0/scale)^shape)
    # scale = x0 / (1/S0 - 1)^(1/shape)
    scale <- x0 / ((1/S0 - 1)^(1/shape))
  } else if (dist == "GM") {
    # Gamma: need to solve numerically
    root_fn <- function(s) 1 - stats::pgamma(x0, shape, s) - S0
    scale <- stats::uniroot(root_fn, c(0, 10))$root
  } else {
    stop("Unsupported distribution for H0 calculation")
  }
  return(scale)
}
