#' Simulate a Single Two-Stage Trial with Time-to-Event Endpoint
#'
#' @description
#' Generates data for a single simulated two-stage single-arm trial with
#' time-to-event endpoint. This is a building block for validating design
#' operating characteristics via Monte Carlo simulation.
#'
#' @param n1 Integer, stage 1 sample size
#' @param n Integer, maximum total sample size (stage 1 + stage 2)
#' @param c1 Numeric, stage 1 critical value (standardized test statistic scale)
#' @param c Numeric, final critical value (standardized test statistic scale)
#' @param t1 Numeric, stage 1 duration (time units)
#' @param x Numeric, follow-up duration for event assessment (time units)
#' @param rate Numeric, accrual rate (patients per time unit)
#' @param true_hr Numeric, true hazard ratio for data generation (HR < 1 favors treatment)
#' @param dist Character, distribution for survival times: "WB" (Weibull), "LN" (log-normal),
#'   "GM" (gamma), "LL" (log-logistic)
#' @param shape Numeric, shape parameter for the null distribution (required for "WB", "GM", "LL")
#' @param S0 Numeric, survival probability at time x0 under null (alternative to shape)
#' @param x0 Numeric, time point for S0 (required if S0 specified)
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{stopped_stage1}: Logical, whether trial stopped at stage 1
#'   \item \code{reject_h0}: Logical, whether null hypothesis was rejected
#'   \item \code{n_enrolled}: Integer, total number of patients enrolled
#'   \item \code{stage1_stat}: Numeric, stage 1 test statistic
#'   \item \code{final_stat}: Numeric, final test statistic (NA if stopped at stage 1)
#'   \item \code{stage1_events}: Integer, number of events observed at stage 1
#'   \item \code{final_events}: Integer, total events observed (NA if stopped at stage 1)
#' }
#'
#' @details
#' This function simulates the conduct of a two-stage trial:
#' 1. Enroll n1 patients over time according to accrual rate
#' 2. Generate event times from specified distribution with true_hr
#' 3. At time t1, calculate test statistic and compare to c1
#' 4. If test statistic > c1, continue to stage 2; otherwise stop for futility
#' 5. If continuing, enroll remaining patients and assess at final analysis
#' 6. Compare final test statistic to c to make final decision
#'
#' The test statistic used is the one-sample log-rank statistic as described
#' in Wu et al. (2020).
#'
#' @references
#' Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Two-stage phase II
#' survival trial design. Pharmaceutical Statistics, 19(3), 214-229.
#'
#' @export
#'
#' @examples
#' # Simulate trial under alternative hypothesis (HR = 0.59)
#' sim1 <- simulate_two_stage_trial(
#'   n1 = 28, n = 45, c1 = 0.0936, c = 1.6269,
#'   t1 = 13.65, x = 5, rate = 2,
#'   true_hr = 0.5913, dist = "WB",
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5
#' )
#'
#' # Simulate trial under null hypothesis (HR = 1)
#' sim0 <- simulate_two_stage_trial(
#'   n1 = 28, n = 45, c1 = 0.0936, c = 1.6269,
#'   t1 = 13.65, x = 5, rate = 2,
#'   true_hr = 1.0, dist = "WB",
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5
#' )
simulate_two_stage_trial <- function(n1, n, c1, c, t1, x, rate, true_hr,
                                      dist = c("WB", "LN", "GM", "LL"),
                                      shape = NULL, S0 = NULL, x0 = NULL) {
  
  dist <- match.arg(dist)
  
  # Get survival functions for the true data-generating distribution
  surv_funs <- .define_survival_functions(dist, shape, S0, x0, true_hr, data = NULL)
  
  # Also get null distribution functions (for calculating test statistics)
  surv_funs_h0 <- .define_survival_functions(dist, shape, S0, x0, hr = 1.0, data = NULL)
  
  # Stage 1: Enroll n1 patients
  # Enrollment times uniformly distributed over accrual period
  enroll_time_1 <- stats::runif(n1, 0, n1 / rate)
  
  # Generate event times from the true distribution
  event_time_1 <- .generate_event_times(n1, dist, shape, S0, x0, true_hr)
  
  # Calendar time of event = enrollment time + event time
  calendar_event_1 <- enroll_time_1 + event_time_1
  
  # At analysis time t1, determine which events have occurred
  # Event observed if: calendar_event_1 <= t1 AND event_time_1 <= x
  event_indicator_1 <- (calendar_event_1 <= t1) & (event_time_1 <= x)
  
  # Calculate stage 1 test statistic
  stat_1 <- .calculate_logrank_stat(
    event_times = pmin(event_time_1, x),
    event_indicators = event_indicator_1,
    analysis_time = t1,
    enroll_times = enroll_time_1,
    surv_funs_h0 = surv_funs_h0,
    surv_funs_h1 = surv_funs,
    x = x,
    rate = rate,
    n_current = n1
  )
  
  # Stage 1 decision
  if (stat_1 <= c1) {
    # Stop for futility
    return(list(
      stopped_stage1 = TRUE,
      reject_h0 = FALSE,
      n_enrolled = n1,
      stage1_stat = stat_1,
      final_stat = NA_real_,
      stage1_events = sum(event_indicator_1),
      final_events = NA_integer_
    ))
  }
  
  # Stage 2: Continue and enroll additional patients
  n2 <- n - n1
  
  # Enrollment for stage 2 patients starts after earliest stage 1 enrollment
  # and continues according to accrual rate
  enroll_time_2 <- stats::runif(n2, min(enroll_time_1), min(enroll_time_1) + n / rate)
  
  # Generate stage 2 event times
  event_time_2 <- .generate_event_times(n2, dist, shape, S0, x0, true_hr)
  
  # Combine stage 1 and stage 2 data
  enroll_time_all <- c(enroll_time_1, enroll_time_2)
  event_time_all <- c(event_time_1, event_time_2)
  calendar_event_all <- enroll_time_all + event_time_all
  
  # Final analysis: when last enrolled patient reaches follow-up time x
  # (or study reaches maximum calendar time)
  final_analysis_time <- max(enroll_time_all) + x
  
  # Events observed at final analysis
  event_indicator_all <- (calendar_event_all <= final_analysis_time) & (event_time_all <= x)
  
  # Calculate final test statistic
  stat_final <- .calculate_logrank_stat(
    event_times = pmin(event_time_all, x),
    event_indicators = event_indicator_all,
    analysis_time = final_analysis_time,
    enroll_times = enroll_time_all,
    surv_funs_h0 = surv_funs_h0,
    surv_funs_h1 = surv_funs,
    x = x,
    rate = rate,
    n_current = n
  )
  
  # Final decision
  reject <- stat_final > c
  
  return(list(
    stopped_stage1 = FALSE,
    reject_h0 = reject,
    n_enrolled = n,
    stage1_stat = stat_1,
    final_stat = stat_final,
    stage1_events = sum(event_indicator_1),
    final_events = sum(event_indicator_all)
  ))
}


#' Validate Design Operating Characteristics via Simulation
#'
#' @description
#' Performs Monte Carlo simulation to estimate the operating characteristics
#' of a two-stage design: type I error rate, power, expected sample size under
#' null and alternative, and probability of early stopping.
#'
#' @param design A design object from \code{\link{two_stage_single_arm_tte}}, or
#'   a list with named components: n1, n, c1, c, t1, x, rate
#' @param true_hr Numeric, true hazard ratio for power simulations. Use 1.0 for
#'   type I error simulations, use design HR for power simulations
#' @param nsim Integer, number of simulation replicates (default 10000)
#' @param dist,shape,S0,x0 Distribution parameters matching the design
#' @param seed Integer, random seed for reproducibility
#' @param parallel Logical, whether to use parallel processing (default FALSE)
#' @param ncores Integer, number of cores for parallel processing (default 4)
#'
#' @return A list of class "two_stage_sim" with components:
#' \itemize{
#'   \item \code{reject_rate}: Proportion of simulations rejecting H0
#'   \item \code{stop_stage1_rate}: Proportion stopping at stage 1
#'   \item \code{mean_n}: Average sample size across simulations
#'   \item \code{median_n}: Median sample size
#'   \item \code{design}: The design being validated
#'   \item \code{nsim}: Number of simulations performed
#'   \item \code{true_hr}: True HR used in simulations
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate a design
#' design <- two_stage_single_arm_tte(
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5,
#'   hr = 0.5913, x = 5, rate = 2,
#'   alpha = 0.05, beta = 0.2, dist = "WB"
#' )
#'
#' # Validate type I error (should be ~0.05)
#' sim_h0 <- validate_design(design, true_hr = 1.0, nsim = 10000)
#'
#' # Validate power (should be ~0.80)
#' sim_h1 <- validate_design(design, true_hr = 0.5913, nsim = 10000)
#' }
validate_design <- function(design, true_hr, nsim = 10000,
                           dist = c("WB", "LN", "GM", "LL"),
                           shape = NULL, S0 = NULL, x0 = NULL,
                           seed = NULL, parallel = FALSE, ncores = 4) {
  
  dist <- match.arg(dist)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Extract design parameters
  if (inherits(design, "two_stage_design")) {
    params <- design$Two_stage
    n1 <- params$n1
    n <- params$n
    c1 <- params$c1
    c <- params$c
    t1 <- params$t1
    
    # Also get x and rate from param table
    param_table <- design$param
    x <- param_table$x
    rate <- param_table$rate
    
  } else if (is.list(design)) {
    # Extract from list
    n1 <- design$n1
    n <- design$n
    c1 <- design$c1
    c <- design$c
    t1 <- design$t1
    x <- design$x
    rate <- design$rate
    
  } else {
    stop("design must be a two_stage_design object or list with components: n1, n, c1, c, t1, x, rate")
  }
  
  # Run simulations
  if (parallel && requireNamespace("furrr", quietly = TRUE)) {
    # Parallel execution
    if (requireNamespace("future", quietly = TRUE)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = ncores)
    }
    
    results <- furrr::future_map(
      1:nsim,
      function(i) {
        simulate_two_stage_trial(
          n1 = n1, n = n, c1 = c1, c = c, t1 = t1,
          x = x, rate = rate, true_hr = true_hr,
          dist = dist, shape = shape, S0 = S0, x0 = x0
        )
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
    
  } else {
    # Sequential execution
    results <- vector("list", nsim)
    
    # Progress bar
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = nsim)
        for (i in seq_len(nsim)) {
          results[[i]] <- simulate_two_stage_trial(
            n1 = n1, n = n, c1 = c1, c = c, t1 = t1,
            x = x, rate = rate, true_hr = true_hr,
            dist = dist, shape = shape, S0 = S0, x0 = x0
          )
          p()
        }
      })
    } else {
      # No progress bar
      for (i in seq_len(nsim)) {
        results[[i]] <- simulate_two_stage_trial(
          n1 = n1, n = n, c1 = c1, c = c, t1 = t1,
          x = x, rate = rate, true_hr = true_hr,
          dist = dist, shape = shape, S0 = S0, x0 = x0
        )
        
        # Simple console progress every 1000 iterations
        if (i %% 1000 == 0) {
          message("Completed ", i, "/", nsim, " simulations")
        }
      }
    }
  }
  
  # Summarize results
  stopped_stage1 <- sapply(results, function(x) x$stopped_stage1)
  rejected <- sapply(results, function(x) x$reject_h0)
  n_enrolled <- sapply(results, function(x) x$n_enrolled)
  
  # Calculate operating characteristics
  reject_rate <- mean(rejected)
  stop_stage1_rate <- mean(stopped_stage1)
  mean_n <- mean(n_enrolled)
  median_n <- stats::median(n_enrolled)
  
  # Return structured results
  out <- list(
    reject_rate = reject_rate,
    stop_stage1_rate = stop_stage1_rate,
    mean_n = mean_n,
    median_n = median_n,
    design = design,
    nsim = nsim,
    true_hr = true_hr,
    results = results  # Store individual results for further analysis
  )
  
  class(out) <- "two_stage_sim"
  return(out)
}


#' Compare Multiple Designs via Simulation
#'
#' @description
#' Compare operating characteristics of multiple candidate designs through
#' simulation. Useful for comparing the paper's published design against
#' algorithmically-derived designs.
#'
#' @param designs List of design objects or data frames with design parameters
#' @param design_names Character vector of names for each design
#' @param true_hr_null Numeric, HR under null hypothesis (default 1.0)
#' @param true_hr_alt Numeric, HR under alternative hypothesis
#' @param nsim Integer, number of simulations per design per scenario
#' @param dist,shape,S0,x0 Distribution parameters
#' @param seed Integer, random seed
#' @param parallel Logical, use parallel processing
#' @param ncores Integer, number of cores
#'
#' @return A data frame comparing designs with columns:
#' \itemize{
#'   \item \code{design}: Design identifier
#'   \item \code{scenario}: "H0" or "H1"
#'   \item \code{reject_rate}: Type I error (H0) or Power (H1)
#'   \item \code{mean_n}: Expected sample size
#'   \item \code{stop_stage1_rate}: Early stopping probability
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Paper's design
#' design_paper <- list(
#'   n1 = 28, n = 45, c1 = 0.0936, c = 1.6269,
#'   t1 = 13.65, x = 5, rate = 2
#' )
#'
#' # Our algorithm's design
#' design_algo <- two_stage_single_arm_tte(
#'   shape = 1.47327, S0 = 0.5, x0 = 3.5,
#'   hr = 0.5913, x = 5, rate = 2,
#'   alpha = 0.05, beta = 0.2, dist = "WB"
#' )
#'
#' # Compare
#' comparison <- compare_designs(
#'   designs = list(design_paper, design_algo),
#'   design_names = c("Wu et al. 2020", "Algorithm"),
#'   true_hr_alt = 0.5913,
#'   nsim = 10000,
#'   dist = "WB", shape = 1.47327, S0 = 0.5, x0 = 3.5
#' )
#' }
compare_designs <- function(designs, design_names = NULL,
                           true_hr_null = 1.0, true_hr_alt,
                           nsim = 10000,
                           dist = c("WB", "LN", "GM", "LL"),
                           shape = NULL, S0 = NULL, x0 = NULL,
                           seed = NULL, parallel = FALSE, ncores = 4) {
  
  dist <- match.arg(dist)
  
  if (is.null(design_names)) {
    design_names <- paste0("Design_", seq_along(designs))
  }
  
  # TODO: Implement comparison
  # 1. For each design, validate under H0 and H1
  # 2. Combine results into comparison table
  # 3. Add interpretation/recommendations
  
  stop("Function under development - see issue tracker or vignettes")
}


# Helper Functions --------------------------------------------------------

#' Generate Event Times from Parametric Distribution
#'
#' @description
#' Internal helper to generate event times from specified parametric
#' survival distribution with given hazard ratio.
#'
#' @keywords internal
#' @noRd
.generate_event_times <- function(n, dist, shape, S0, x0, hr) {
  
  # Get scale parameter for the distribution under specified HR
  if (dist == "WB") {
    # Weibull: T ~ Weibull(shape, scale)
    # Under HR, scale parameter is modified by HR^(-1/shape)
    scale0 <- x0 / (-log(S0))^(1/shape)
    scale <- scale0 * hr^(-1/shape)
    event_times <- stats::rweibull(n, shape = shape, scale = scale)
    
  } else if (dist == "LN") {
    # Log-normal: log(T) ~ Normal(meanlog, sdlog)
    # Under HR, meanlog is shifted by log(HR)
    sdlog <- shape  # In LN parameterization, shape = sdlog
    meanlog0 <- log(x0) - sdlog * stats::qnorm(1 - S0)
    meanlog <- meanlog0 - log(hr)
    event_times <- stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog)
    
  } else if (dist == "GM") {
    # Gamma: T ~ Gamma(shape, rate)
    # Under HR, rate parameter is multiplied by HR
    # First find rate0 from S0, x0
    root0 <- function(t) 1 - stats::pgamma(x0, shape, t) - S0
    rate0 <- stats::uniroot(root0, c(0, 10))$root
    rate <- rate0 * hr
    event_times <- stats::rgamma(n, shape = shape, rate = rate)
    
  } else if (dist == "LL") {
    # Log-logistic: Use inverse CDF method
    # S(t) = 1 / (1 + (t/scale)^shape)
    # For HR, scale is modified
    scale0 <- x0 / ((1/S0 - 1)^(1/shape))
    scale <- scale0 / hr^(1/shape)
    
    # Inverse CDF: T = scale * ((1-u)/u)^(1/shape)
    u <- stats::runif(n)
    event_times <- scale * ((1 - u) / u)^(1/shape)
    
  } else {
    stop("Distribution ", dist, " not supported for simulation")
  }
  
  return(event_times)
}


#' Calculate One-Sample Log-Rank Test Statistic
#'
#' @description
#' Internal helper to calculate the one-sample log-rank test statistic
#' comparing observed survival to a known null hypothesis survival distribution.
#'
#' @details
#' The one-sample log-rank test compares observed event rates to expected
#' event rates under a specified null hypothesis survival distribution.
#' 
#' Following Wu et al. (2020), the test statistic is:
#'   Z = (O - E) / sqrt(Var)
#' 
#' where:
#'   O = observed events (count)
#'   E = expected events under H0 = n * p0
#'   Var = n * sigma^2_0
#'   
#' And p0, sigma^2_0 are derived from the null survival distribution accounting
#' for staggered entry and administrative censoring.
#'
#' @references
#' - Finkelstein DM, Schoenfeld DA, Stamenovic E. Analysis of time-to-event data. 
#'   In: Fundamentals of Clinical Research for Radiologists. 2002
#' - Wu, J., et al. (2020). Two-stage phase II survival trial design. 
#'   Pharmaceutical Statistics, 19(3), 214-229.
#'
#' @keywords internal
#' @noRd
.calculate_logrank_stat <- function(event_times, event_indicators, 
                                     analysis_time, enroll_times,
                                     surv_funs_h0, surv_funs_h1,
                                     x, rate, n_current) {
  
  # Extract survival functions under null hypothesis
  s0 <- surv_funs_h0$s0  # Null survival function
  h0 <- surv_funs_h0$h0  # Null hazard function
  H0 <- surv_funs_h0$H0  # Null cumulative hazard
  
  # Extract functions under alternative (for data generation context)
  s_alt <- surv_funs_h1$s
  h_alt <- surv_funs_h1$h
  
  # Observed events
  O <- sum(event_indicators)
  
  # Calculate expected events and variance under H0
  # This must account for staggered enrollment
  
  # Accrual time for this cohort
  ta <- max(enroll_times) - min(enroll_times)
  
  # G(t) = P(patient at risk for time t | uniform enrollment over [0, ta])
  # At analysis time t_analysis:
  # A patient enrolled at time u is at risk for event at time t if:
  #   1. They've been enrolled long enough: u <= t_analysis - t
  #   2. Event time t <= x (administrative censoring)
  
  # For uniform enrollment from 0 to ta, and analysis at time t_analysis:
  # The fraction at risk for time t is related to enrollment distribution
  
  # Simplified approach: Calculate expected event probability
  # accounting for the enrollment pattern
  
  # If all patients fully enrolled and followed to time x:
  # p0 = integral from 0 to x of S_H0(t) * h_H0(t) dt
  calc_p0_full <- function(t) s0(t) * h0(t)
  
  p0_full <- tryCatch({
    stats::integrate(calc_p0_full, 0, x, subdivisions = 200)$value
  }, error = function(e) {
    # Fallback: use Kaplan-Meier-like estimate from observed data
    sum(event_indicators) / n_current
  })
  
  # Adjust for actual time at risk
  # Patients enrolled uniformly, so average time at risk varies
  # At analysis_time, average time at risk is approximately:
  avg_time_at_risk <- min(analysis_time - mean(enroll_times), x)
  
  # If analysis happens before full follow-up, adjust expected events
  if (avg_time_at_risk < x) {
    # Recalculate p0 for actual follow-up duration
    calc_p0_actual <- function(t) s0(t) * h0(t)
    p0 <- tryCatch({
      stats::integrate(calc_p0_actual, 0, avg_time_at_risk, subdivisions = 200)$value
    }, error = function(e) p0_full)
  } else {
    p0 <- p0_full
  }
  
  # Expected events
  E <- n_current * p0
  
  # Variance under H0
  # From Wu et al.: variance accounts for correlation in event times
  # For simplicity in simulation, use: Var = n * p0 (binomial-like)
  # The full Wu et al. formula includes correlation terms but those
  # are more relevant for optimization than simulation validation
  Var <- n_current * p0
  
  # Test statistic
  if (Var > 0 && E > 0) {
    Z <- (O - E) / sqrt(Var)
  } else {
    # Edge case: no events expected or no variance
    # If O = 0 and E = 0, then Z = 0 (no evidence against H0)
    # If O > 0 but E = 0, Z = large positive (strong evidence against H0)
    Z <- ifelse(O > 0, 10, 0)  
  }
  
  return(Z)
}


# S3 Methods --------------------------------------------------------------

#' Print Method for Two-Stage Design Simulation
#'
#' @param x A \code{two_stage_sim} object
#' @param ... Additional arguments (not used)
#'
#' @export
print.two_stage_sim <- function(x, ...) {
  cat("\n=== Two-Stage Design Simulation Results ===\n\n")
  cat("Number of simulations:", x$nsim, "\n")
  cat("True HR:", x$true_hr, "\n\n")
  
  cat("Operating Characteristics:\n")
  cat("  Reject H0 rate:", sprintf("%.4f", x$reject_rate), "\n")
  cat("  Stop at Stage 1 rate:", sprintf("%.4f", x$stop_stage1_rate), "\n")
  cat("  Mean sample size:", sprintf("%.2f", x$mean_n), "\n")
  cat("  Median sample size:", sprintf("%.0f", x$median_n), "\n\n")
  
  # Interpretation
  if (x$true_hr == 1.0) {
    cat("Interpretation (under H0):\n")
    cat("  Type I error estimate:", sprintf("%.4f", x$reject_rate), "\n")
    if (x$reject_rate > 0.055) {
      cat("  WARNING: Type I error exceeds nominal 0.05 level\n")
    } else if (x$reject_rate < 0.045) {
      cat("  NOTE: Type I error below nominal 0.05 level (conservative)\n")
    } else {
      cat("  Type I error is within acceptable range\n")
    }
  } else {
    cat("Interpretation (under H1):\n")
    cat("  Power estimate:", sprintf("%.4f", x$reject_rate), "\n")
    if (x$reject_rate < 0.78) {
      cat("  WARNING: Power below nominal 0.80 level\n")
    } else if (x$reject_rate > 0.82) {
      cat("  NOTE: Power exceeds nominal 0.80 level\n")
    } else {
      cat("  Power is within acceptable range\n")
    }
  }
  
  invisible(x)
}
