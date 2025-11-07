#' Simulate a Trial from a Design Object
#'
#' @description
#' Simulates a single realization of a trial based on a design object created by
#' `two_stage_single_arm_tte()`, `single_stage_single_arm_tte()`, or other design 
#' functions in this package.
#' 
#' This is a generic function that will work with any trial design created by
#' blueprint package functions. It automatically detects whether the design is
#' single-stage or two-stage based on the class attribute.
#'
#' @param design A design object (e.g., from `two_stage_single_arm_tte()` or 
#'   `single_stage_single_arm_tte()`)
#' @param shape Shape parameter for the survival distribution under simulation
#' @param scale Scale parameter for the survival distribution under simulation
#'
#' @return A list containing:
#'   \item{reject_h0}{Logical; TRUE if null hypothesis is rejected}
#'   \item{stopped_early}{Logical; TRUE if trial stopped at interim (two-stage only; FALSE for single-stage)}
#'   \item{enrolled_n}{Number of patients enrolled}
#'   \item{interim_data}{Data frame of patient data at interim analysis (two-stage only; NULL for single-stage)}
#'   \item{final_data}{Data frame of patient data at final analysis}
#'   \item{interim_stat}{Test statistic value at interim (two-stage only; NA for single-stage)}
#'   \item{final_stat}{Test statistic value at final}
#'   \item{interim_time}{Calendar time of interim analysis (two-stage only; NA for single-stage)}
#'   \item{final_time}{Calendar time of final analysis}
#'   \item{interim_events}{Number of events observed at interim (two-stage only; NA for single-stage)}
#'   \item{final_events}{Number of events observed at final}
#'
#' @examples
#' \dontrun{
#' # Two-stage design
#' design_two <- two_stage_single_arm_tte(
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5, hr = 0.5913,
#'   tf = 5, rate = 2, alpha = 0.05, beta = 0.2,
#'   dist = "WB", restricted = FALSE
#' )
#' 
#' # Simulate one trial under H0
#' scale_h0 <- 3.5 / (-log(0.5))^(1/1.47327)
#' trial_result <- simulate_trial(design_two, shape = 1.47327, scale = scale_h0)
#' 
#' # Single-stage design
#' design_single <- single_stage_single_arm_tte(
#'   shape = 1.0, S0 = 0.5, x0 = 0.5, hr = 0.667,
#'   tf = 1.0, rate = 32.4, alpha = 0.05, beta = 0.10,
#'   two_sided = TRUE, dist = "WB", restricted = FALSE
#' )
#' 
#' # Simulate one trial under H1
#' scale_h0 <- 0.5 / (-log(0.5))^(1/1.0)
#' scale_h1 <- scale_h0 / (0.667^(1/1.0))
#' trial_result <- simulate_trial(design_single, shape = 1.0, scale = scale_h1)
#' }
#'
#' @export
simulate_trial <- function(design, shape, scale) {
  
  # Extract parameters from design object
  # Check for valid design object (either two_stage_design or single_stage_design)
  is_valid <- inherits(design, c("two_stage_design", "single_stage_design")) || 
              (inherits(design, "list") && !is.null(design$param))
  
  if (!is_valid) {
    stop("design must be a valid design object from blueprint package")
  }
  
  # Check if this is a single-stage design
  is_single_stage <- inherits(design, "single_stage_design")
  
  if (is_single_stage) {
    return(.simulate_single_stage_trial(design, shape, scale))
  }
  
  # Two-stage design
  if (is.null(design$two_stage)) {
    stop("Two-stage design object must have a 'two_stage' component")
  }
  
  n <- design$two_stage$n
  n1 <- design$two_stage$n1
  t1 <- design$two_stage$t1  # Interim analysis time (calendar time from study start)
  c1 <- design$two_stage$c1  # Interim critical value
  c <- design$two_stage$c     # Final critical value
  
  # Extract from param
  tf <- design$param$tf
  rate <- design$param$rate
  restricted <- design$param$restricted
  
  # Determine distribution (check for "dist" column first, otherwise use "WB")
  dist <- if ("dist" %in% names(design$param)) design$param$dist else "WB"
  
  # Extract H0 parameters for log-rank test
  # These are needed to calculate expected events under null hypothesis
  S0 <- design$param$S0
  x0 <- design$param$x0
  shape_h0 <- design$param$shape
  
  # Calculate H0 scale parameter from S0 and x0
  if (!is.null(S0) && !is.null(x0) && !is.null(shape_h0)) {
    scale_h0 <- .calculate_h0_scale(dist, S0, x0, shape_h0)
  } else {
    scale_h0 <- NULL
  }
  
  # Generate accrual times for all n patients
  accrual_period <- n / rate
  accrual_times <- sort(runif(n, 0, accrual_period))
  
  # Generate event times based on distribution
  event_times <- switch(dist,
    "WB" = rweibull(n, shape = shape, scale = scale),
    "LN" = rlnorm(n, meanlog = log(scale), sdlog = shape),
    "GM" = rgamma(n, shape = shape, scale = scale),
    "LG" = {
      # Log-logistic: use inverse CDF method
      u <- runif(n)
      scale * ((u / (1 - u))^(1/shape))
    },
    stop("Unsupported distribution. Use 'WB', 'LN', 'GM', or 'LG'.")
  )
  
  # === INTERIM ANALYSIS (Stage 1) ===
  # Only first n1 patients have been enrolled
  interim_accrual <- accrual_times[1:n1]
  interim_events <- event_times[1:n1]
  
  # Time of interim analysis (calendar time from study start)
  # t1 is provided by the design
  # For unrestricted follow-up: patients followed from enrollment to t1
  # For restricted follow-up: patients followed for min(tf, t1-enrollment)
  
  # Calculate observed times and event indicators at interim
  # obs_time = min(event_time, time_from_enrollment_to_interim)
  # Following OneArm2stage: xt1 = pmax(0, pmin(w1, t1-u1))
  if (restricted) {
    # Restricted: follow each patient for min(tf, t1 - accrual_time)
    time_on_study <- pmax(0, pmin(tf, t1 - interim_accrual))
    interim_obs_time <- pmin(interim_events, time_on_study)
    interim_event_ind <- interim_events < time_on_study
  } else {
    # Unrestricted: follow each patient from enrollment to t1
    time_on_study <- pmax(0, t1 - interim_accrual)
    interim_obs_time <- pmin(interim_events, time_on_study)
    interim_event_ind <- interim_events < time_on_study
  }
  
  interim_data <- data.frame(
    patient_id = 1:n1,
    accrual_time = interim_accrual,
    event_time = interim_events,
    obs_time = interim_obs_time,
    event = interim_event_ind
  )
  
  # Calculate one-sample log-rank statistic at interim
  # Use proper expected events under H0 if parameters available
  if (!is.null(scale_h0)) {
    interim_stat <- logrank_stat_wu(
      obs_times = interim_obs_time,
      event_ind = interim_event_ind,
      dist = dist,
      shape_h0 = shape_h0,
      scale_h0 = scale_h0
    )
  } else {
    # Fallback to simple standardized count (less accurate)
    O_interim <- sum(interim_event_ind)
    E_interim <- n1 * mean(interim_event_ind)
    V_interim <- n1 * var(interim_event_ind)
    
    if (V_interim == 0 || is.na(V_interim)) {
      interim_stat <- O_interim - E_interim
    } else {
      interim_stat <- (O_interim - E_interim) / sqrt(V_interim)
    }
  }
  
  # Decision at interim
  if (is.na(interim_stat) || interim_stat <= c1) {
    # Stop early for futility
    return(list(
      reject_h0 = FALSE,
      stopped_early = TRUE,
      enrolled_n = n1,
      interim_data = interim_data,
      final_data = NULL,
      interim_stat = interim_stat,
      final_stat = NA,
      interim_time = t1,
      final_time = NA,
      interim_events = sum(interim_event_ind),
      final_events = NA
    ))
  }
  
  # === FINAL ANALYSIS (Stage 2) ===
  # Continue to enroll remaining patients
  
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
    V_final <- n * var(final_event_ind)
    
    if (V_final == 0 || is.na(V_final)) {
      final_stat <- O_final - E_final
    } else {
      final_stat <- (O_final - E_final) / sqrt(V_final)
    }
  }
  
  # Decision at final
  reject <- !is.na(final_stat) && final_stat > c
  
  return(list(
    reject_h0 = reject,
    stopped_early = FALSE,
    enrolled_n = n,
    interim_data = interim_data,
    final_data = final_data,
    interim_stat = interim_stat,
    final_stat = final_stat,
    interim_time = t1,
    final_time = t_final,
    interim_events = sum(interim_event_ind),
    final_events = sum(final_event_ind)
  ))
}


#' Calculate H0 Scale Parameter from S0 and x0
#'
#' @description
#' Internal helper to calculate the scale parameter for the null hypothesis
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
    scale <- log(x0) - shape * qnorm(1 - S0)
  } else if (dist == "LG") {
    # Log-logistic: S(t) = 1 / (1 + (t/scale)^shape)
    # S0 = 1 / (1 + (x0/scale)^shape)
    # scale = x0 / (1/S0 - 1)^(1/shape)
    scale <- x0 / ((1/S0 - 1)^(1/shape))
  } else if (dist == "GM") {
    # Gamma: need to solve numerically
    root_fn <- function(s) 1 - pgamma(x0, shape, s) - S0
    scale <- uniroot(root_fn, c(0, 10))$root
  } else {
    stop("Unsupported distribution for H0 calculation")
  }
  return(scale)
}


#' Simulate Single-Stage Trial
#'
#' @description
#' Internal helper function to simulate a single-stage trial.
#'
#' @keywords internal
#' @noRd
.simulate_single_stage_trial <- function(design, shape, scale) {
  
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
  accrual_times <- sort(runif(n, 0, accrual_period))
  
  # Generate event times based on distribution
  event_times <- switch(dist,
    "WB" = rweibull(n, shape = shape, scale = scale),
    "LN" = rlnorm(n, meanlog = log(scale), sdlog = shape),
    "GM" = rgamma(n, shape = shape, scale = scale),
    "LG" = {
      # Log-logistic: use inverse CDF method
      u <- runif(n)
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
    V_final <- n * var(final_event_ind)
    
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
