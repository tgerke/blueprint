#' Simulate a Trial from a Design Object
#'
#' @description
#' Simulates a single realization of a trial based on a design object created by
#' `two_stage_single_arm_tte()` or other design functions in this package.
#' 
#' This is a generic function that will work with any trial design created by
#' blueprint package functions.
#'
#' @param design A design object (e.g., from `two_stage_single_arm_tte()`)
#' @param shape Shape parameter for the survival distribution under simulation
#' @param scale Scale parameter for the survival distribution under simulation
#'
#' @return A list containing:
#'   \item{reject_h0}{Logical; TRUE if null hypothesis is rejected}
#'   \item{stopped_early}{Logical; TRUE if trial stopped at interim}
#'   \item{enrolled_n}{Number of patients enrolled}
#'   \item{interim_data}{Data frame of patient data at interim analysis}
#'   \item{final_data}{Data frame of patient data at final analysis (if applicable)}
#'   \item{interim_stat}{Test statistic value at interim}
#'   \item{final_stat}{Test statistic value at final (if applicable)}
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
#' # Simulate one trial under H0
#' # (Calculate scale_h0 from design parameters)
#' scale_h0 <- 3.5 / (-log(0.5))^(1/1.47327)
#' trial_result <- simulate_trial(design, shape = 1.47327, scale = scale_h0)
#' 
#' # Simulate one trial under H1  
#' scale_h1 <- scale_h0 / (0.5913^(1/1.47327))
#' trial_result <- simulate_trial(design, shape = 1.47327, scale = scale_h1)
#' }
#'
#' @export
simulate_trial <- function(design, shape, scale) {
  
  # Extract parameters from design object
  if (!inherits(design, "list") || is.null(design$Two_stage) || is.null(design$param)) {
    stop("design must be a valid design object from blueprint package")
  }
  
  n <- design$Two_stage$n
  n1 <- design$Two_stage$n1
  t1 <- design$Two_stage$t1  # Interim analysis time (calendar time from study start)
  c1 <- design$Two_stage$c1  # Interim critical value
  c <- design$Two_stage$c     # Final critical value
  
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
  
  .simulate_two_stage_trial_internal(
    n = n, n1 = n1, t1 = t1, c1 = c1, c = c, tf = tf, rate = rate,
    dist = dist, shape = shape, scale = scale, restricted = restricted,
    S0 = S0, x0 = x0, shape_h0 = shape_h0
  )
}


#' Internal Function: Simulate a Two-Stage Single-Arm Trial
#'
#' @description
#' Internal worker function that performs the actual simulation.
#' Users should call `simulate_trial()` instead.
#'
#' @keywords internal
#' @noRd
.simulate_two_stage_trial_internal <- function(n, n1, t1, c1, c, tf, rate, dist, shape, scale,
                                      restricted = FALSE, S0 = NULL, x0 = NULL, shape_h0 = NULL) {
  
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
    interim_stat <- .calculate_logrank_stat(
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
      final_stat = NA
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
    final_stat <- .calculate_logrank_stat(
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
    final_stat = final_stat
  ))
}


#' Simulate Operating Characteristics of Two-Stage Design
#'
#' @description
#' Runs multiple simulations to estimate operating characteristics (Type I error,
#' power, expected sample size) for a two-stage single-arm trial design.
#'
#' @param design A design object from `two_stage_single_arm_tte()` or similar
#' @param shape Shape parameter for survival distribution 
#' @param scale_h0 Scale parameter under null hypothesis
#' @param scale_h1 Scale parameter under alternative hypothesis
#' @param n_sims Number of simulations to run (default: 1000)
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list containing:
#'   \item{type1_error}{Estimated Type I error rate}
#'   \item{power}{Estimated power}
#'   \item{expected_n_h0}{Expected sample size under H0}
#'   \item{expected_n_h1}{Expected sample size under H1}
#'   \item{prob_early_stop_h0}{Probability of early stopping under H0}
#'   \item{prob_early_stop_h1}{Probability of early stopping under H1}
#'   \item{results_h0}{Individual simulation results under H0}
#'   \item{results_h1}{Individual simulation results under H1}
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
#' scale_h0 <- 3.5 / (-log(0.5))^(1/1.47327)
#' scale_h1 <- scale_h0 / (0.5913^(1/1.47327))
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
#' }
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
  
  list(
    type1_error = type1_error,
    power = power,
    expected_n_h0 = expected_n_h0,
    expected_n_h1 = expected_n_h1,
    prob_early_stop_h0 = prob_early_stop_h0,
    prob_early_stop_h1 = prob_early_stop_h1,
    results_h0 = results_h0,
    results_h1 = results_h1
  )
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


#' Calculate One-Sample Log-Rank Test Statistic
#'
#' @description
#' Computes the one-sample log-rank test statistic for comparing observed
#' survival times against a null hypothesis survival function.
#'
#' Based on Wu et al. (2020) and OneArm2stage implementation.
#' Uses cumulative hazard under H0 for expected events.
#'
#' @keywords internal
#' @noRd
.calculate_logrank_stat <- function(obs_times, event_ind, dist, shape_h0, scale_h0) {
  
  # Create cumulative hazard function under H0
  # H(t) = -log(S(t))
  H0 <- switch(dist,
    "WB" = function(t) {
      # H(t) = (t/scale)^shape for Weibull
      (t / scale_h0)^shape_h0
    },
    "LN" = function(t) {
      # H(t) = -log(1 - plnorm(...))
      S <- 1 - plnorm(t, meanlog = scale_h0, sdlog = shape_h0)
      -log(pmax(S, 1e-10))  # Avoid log(0)
    },
    "LG" = function(t) {
      # H(t) = -log(S(t)) where S(t) = 1/(1+(t/scale)^shape)
      -log(1 / (1 + (t/scale_h0)^shape_h0))
    },
    "GM" = function(t) {
      S <- 1 - pgamma(t, shape = shape_h0, scale = scale_h0)
      -log(pmax(S, 1e-10))
    },
    stop("Unsupported distribution")
  )
  
  # Expected events = sum of cumulative hazards at observation times
  # This matches OneArm2stage: M=H(shape,scale0,xt); E=sum(M)
  expected_events <- sum(H0(obs_times))
  
  # Observed events
  observed_events <- sum(event_ind)
  
  # Test statistic: Z = (E - O) / sqrt(E)
  # OneArm2stage uses variance ≈ E (Poisson approximation)
  # This is simpler than the full variance formula
  # NOTE: Sign convention - positive Z favors H1 (fewer events than expected)
  if (expected_events > 0) {
    Z <- (expected_events - observed_events) / sqrt(expected_events)
  } else {
    # Degenerate case
    Z <- 0
  }
  
  return(Z)
}
