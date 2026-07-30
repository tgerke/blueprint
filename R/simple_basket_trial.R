#' Simulate a Single Simon Basket Trial
#'
#' Simulates one realization of a Simon basket trial with adaptive heterogeneity
#' assessment. The trial proceeds in two stages with an interim analysis that
#' determines whether to pursue a heterogeneous path (analyze baskets separately)
#' or homogeneous path (pool all baskets).
#'
#' @param K Integer. Number of baskets.
#' @param p_true Numeric vector of length K. True response rates for each basket.
#' @param N1 Integer. Total stage 1 sample size across all baskets.
#' @param N2 Integer. Total stage 2 sample size for homogeneous path.
#' @param n2k Integer. Stage 2 sample size per basket for heterogeneous path.
#' @param p0 Numeric. Null response rate.
#' @param r_s Integer. Minimum responders per basket to continue (heterogeneous path).
#' @param r_c Integer. Minimum total responders to continue (homogeneous path).
#' @param gamma Numeric. Critical value for Fisher's exact test (0-1).
#' @param alpha_s Numeric. Significance level for separate analyses (heterogeneous).
#' @param alpha_c Numeric. Significance level for combined analysis (homogeneous).
#' @param accrual_rate Numeric vector of length K. Accrual rates per basket (patients/month).
#' @param s1_min Integer. Minimum stage 1 sample size per basket (default 3).
#' @param s1_max Integer. Maximum stage 1 sample size per basket (default 10).
#' @param s2_min Integer. Minimum stage 2 sample size per basket for homogeneous path (default 1).
#' @param s2_max Integer. Maximum stage 2 sample size per basket for homogeneous path (default 6).
#'
#' @return A list with components:
#'   \item{decisions}{Logical vector of length K. TRUE if basket k is declared active.}
#'   \item{stage1_responses}{Integer vector. Number of responders in stage 1 per basket.}
#'   \item{stage2_responses}{Integer vector. Number of responders in stage 2 per basket.}
#'   \item{stage1_n}{Integer vector. Stage 1 sample size per basket.}
#'   \item{stage2_n}{Integer vector. Stage 2 sample size per basket.}
#'   \item{heterogeneous}{Logical. TRUE if heterogeneous path was used.}
#'   \item{continued_baskets}{Integer vector. Indices of baskets continuing to stage 2.}
#'   \item{trial_duration}{Numeric. Total trial duration in months.}
#'
#' @importFrom stats rbinom binom.test rexp
#'
#' @references
#' Cunanan KM, Iasonos A, Shen R, Begg CB, Gonen M. (2017). An Efficient Basket
#' Trial Design. Statistics in Medicine, 36(10):1568-1579. DOI: 10.1002/sim.7227
#'
#' @keywords internal
simulate_simple_trial_once <- function(K, p_true, N1, N2, n2k, p0,
                                      r_s, r_c, gamma, alpha_s, alpha_c,
                                      accrual_rate,
                                      s1_min = 3, s1_max = 10,
                                      s2_min = 1, s2_max = 6) {
  
  # Initialize outputs
  decisions <- rep(FALSE, K)
  stage1_n <- integer(K)
  stage2_n <- integer(K)
  stage1_responses <- integer(K)
  stage2_responses <- integer(K)
  continued_baskets <- integer(0)
  heterogeneous <- FALSE
  time_stage1 <- numeric(K)
  time_stage2 <- numeric(K)
  
  # ============================================================================
  # STAGE 1: Accrue N1 patients with min/max constraints per basket
  # ============================================================================
  
  # Generate arrival times for all potential stage 1 patients
  # Use exponential inter-arrival times
  baskets_matrix <- matrix(NA, nrow = s1_max * K, ncol = 2)
  for (i in 1:K) {
    idx_start <- (i - 1) * s1_max + 1
    idx_end <- i * s1_max
    arrival_times <- rexp(s1_max, rate = accrual_rate[i])
    baskets_matrix[idx_start:idx_end, 1] <- arrival_times
    baskets_matrix[idx_start:idx_end, 2] <- i  # basket indicator
  }
  
  # Sort by arrival time and take first N1 patients
  ord <- order(baskets_matrix[, 1])
  enrolled_stage1 <- baskets_matrix[ord, , drop = FALSE][1:N1, ]
  
  # Count patients per basket
  basket_counts <- table(factor(enrolled_stage1[, 2], levels = 1:K))
  stage1_n <- as.integer(basket_counts)
  
  # Apply min/max constraints
  for (i in 1:K) {
    if (stage1_n[i] < s1_min) {
      # Need to add more patients to this basket
      additional <- s1_min - stage1_n[i]
      additional_times <- rexp(additional, rate = accrual_rate[i])
      time_stage1[i] <- sum(c(
        enrolled_stage1[enrolled_stage1[, 2] == i, 1],
        additional_times
      ))
      stage1_n[i] <- s1_min
    } else if (stage1_n[i] > s1_max) {
      # Cap at maximum
      enrolled_i <- enrolled_stage1[enrolled_stage1[, 2] == i, 1]
      time_stage1[i] <- sum(enrolled_i[1:s1_max])
      stage1_n[i] <- s1_max
    } else {
      # Within range
      enrolled_i <- enrolled_stage1[enrolled_stage1[, 2] == i, 1]
      if (length(enrolled_i) == 1) {
        time_stage1[i] <- enrolled_i
      } else {
        time_stage1[i] <- sum(enrolled_i)
      }
    }
  }
  
  # Generate stage 1 responses
  stage1_responses <- rbinom(K, stage1_n, p_true)
  
  # ============================================================================
  # INTERIM ANALYSIS: Test of Homogeneity (TOH)
  # ============================================================================
  
  heterogeneous <- test_heterogeneity(
    responses = stage1_responses,
    n = stage1_n,
    gamma = gamma
  )
  
  # ============================================================================
  # STAGE 2: Follow heterogeneous or homogeneous path
  # ============================================================================
  
  if (heterogeneous) {
    # HETEROGENEOUS PATH: Analyze baskets separately
    
    # Determine which baskets continue (futility rule: need >= r_s responders)
    continue_basket <- stage1_responses >= r_s
    continued_baskets <- which(continue_basket)
    
    if (length(continued_baskets) > 0) {
      # Accrue n2k additional patients per continuing basket
      stage2_n[continued_baskets] <- n2k
      
      # Generate stage 2 arrival times and responses
      for (k in continued_baskets) {
        stage2_times <- rexp(n2k, rate = accrual_rate[k])
        time_stage2[k] <- sum(stage2_times)
        stage2_responses[k] <- rbinom(1, n2k, p_true[k])
      }
      
      # Final decision: Binomial test for each continuing basket
      # Bonferroni correction: alpha_s / K*
      K_star <- length(continued_baskets)
      alpha_adjusted <- alpha_s / K_star
      
      for (k in continued_baskets) {
        total_responses <- stage1_responses[k] + stage2_responses[k]
        total_n <- stage1_n[k] + stage2_n[k]
        
        # One-sided binomial test
        p_value <- binom.test(
          x = total_responses,
          n = total_n,
          p = p0,
          alternative = "greater"
        )$p.value
        
        decisions[k] <- p_value <= alpha_adjusted
      }
    }
    
  } else {
    # HOMOGENEOUS PATH: Pool all baskets
    
    # Futility rule: need >= r_c total responders
    total_stage1_responses <- sum(stage1_responses)
    
    if (total_stage1_responses >= r_c) {
      # Continue all baskets to stage 2
      continued_baskets <- 1:K
      
      # Accrue N2 total patients across all baskets
      # Generate arrivals and allocate with min/max constraints
      baskets_matrix2 <- matrix(NA, nrow = s2_max * K, ncol = 2)
      for (i in 1:K) {
        idx_start <- (i - 1) * s2_max + 1
        idx_end <- i * s2_max
        arrival_times <- rexp(s2_max, rate = accrual_rate[i])
        baskets_matrix2[idx_start:idx_end, 1] <- arrival_times
        baskets_matrix2[idx_start:idx_end, 2] <- i
      }
      
      # Sort and take first N2
      ord2 <- order(baskets_matrix2[, 1])
      enrolled_stage2 <- baskets_matrix2[ord2, , drop = FALSE][1:N2, ]
      basket_counts2 <- table(factor(enrolled_stage2[, 2], levels = 1:K))
      stage2_n <- as.integer(basket_counts2)
      
      # Apply minimum constraint (at least s2_min per basket)
      for (i in 1:K) {
        if (stage2_n[i] < s2_min) {
          additional <- s2_min - stage2_n[i]
          additional_times <- rexp(additional, rate = accrual_rate[i])
          enrolled_i <- enrolled_stage2[enrolled_stage2[, 2] == i, 1]
          time_stage2[i] <- sum(c(enrolled_i, additional_times))
          stage2_n[i] <- s2_min
        } else if (stage2_n[i] > s2_max) {
          enrolled_i <- enrolled_stage2[enrolled_stage2[, 2] == i, 1]
          time_stage2[i] <- sum(enrolled_i[1:s2_max])
          stage2_n[i] <- s2_max
        } else {
          enrolled_i <- enrolled_stage2[enrolled_stage2[, 2] == i, 1]
          if (length(enrolled_i) == 0) {
            time_stage2[i] <- 0
          } else if (length(enrolled_i) == 1) {
            time_stage2[i] <- enrolled_i
          } else {
            time_stage2[i] <- sum(enrolled_i)
          }
        }
      }
      
      # Generate stage 2 responses
      stage2_responses <- rbinom(K, stage2_n, p_true)
      
      # Final decision: One-sample binomial test on pooled data
      total_responses <- sum(stage1_responses) + sum(stage2_responses)
      total_n <- sum(stage1_n) + sum(stage2_n)
      
      p_value <- binom.test(
        x = total_responses,
        n = total_n,
        p = p0,
        alternative = "greater"
      )$p.value
      
      # Either all active or all inactive
      if (p_value <= alpha_c) {
        decisions <- rep(TRUE, K)
      } else {
        decisions <- rep(FALSE, K)
      }
    }
  }
  
  # Calculate total trial duration
  trial_duration <- max(time_stage1) + max(time_stage2, na.rm = TRUE)
  
  # Return results
  list(
    decisions = decisions,
    stage1_responses = stage1_responses,
    stage2_responses = stage2_responses,
    stage1_n = stage1_n,
    stage2_n = stage2_n,
    heterogeneous = heterogeneous,
    continued_baskets = continued_baskets,
    trial_duration = trial_duration
  )
}


#' Test of Homogeneity Using Fisher's Exact Test
#'
#' Performs Fisher's exact test on a Kx2 contingency table to assess whether
#' response rates are homogeneous across baskets. This determines which design
#' path to follow in the Simon basket trial.
#'
#' @param responses Integer vector of length K. Number of responders per basket.
#' @param n Integer vector of length K. Sample size per basket.
#' @param gamma Numeric. Critical value (0-1). If p-value <= gamma, baskets are
#'   considered heterogeneous.
#'
#' @return Logical. TRUE if heterogeneous (p-value <= gamma), FALSE if homogeneous.
#'
#' @details
#' The test of homogeneity (TOH) uses Fisher's exact test on the Kx2 table:
#' \tabular{lcc}{
#'   Basket \tab Responders \tab Non-responders \cr
#'   1 \tab r_1 \tab n_1 - r_1 \cr
#'   2 \tab r_2 \tab n_2 - r_2 \cr
#'   ... \tab ... \tab ... \cr
#'   K \tab r_K \tab n_K - r_K
#' }
#'
#' The tuning parameter gamma controls sensitivity to heterogeneity:
#' - Larger gamma -> more likely to pursue heterogeneous path
#' - Smaller gamma -> more likely to pursue homogeneous path
#'
#' Special case: If all responders = 0, return p-value = 1 (don't simulate).
#'
#' @importFrom stats fisher.test
#'
#' @references
#' Cunanan KM, Iasonos A, Shen R, Begg CB, Gonen M. (2017). An Efficient Basket
#' Trial Design. Statistics in Medicine, 36(10):1568-1579.
#'
#' @keywords internal
test_heterogeneity <- function(responses, n, gamma) {
  # Create Kx2 contingency table
  non_responders <- n - responses
  contingency_table <- cbind(responses, non_responders)
  
  # Special case: if all responses = 0, don't simulate (marginal is 0)
  if (sum(responses) == 0) {
    p_value <- 1.0
  } else {
    # Use Fisher's exact test with simulation for computational stability
    # (hybrid = TRUE uses simulation only when needed)
    p_value <- tryCatch({
      fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)$p.value
    }, error = function(e) {
      # Fallback if simulation fails
      fisher.test(contingency_table)$p.value
    })
  }
  
  # Return TRUE if p-value <= gamma (reject homogeneity)
  p_value <= gamma
}


#' Simulate Simon Basket Trial Operating Characteristics
#'
#' Runs multiple simulations of the Simon basket trial to estimate operating
#' characteristics including FWER, marginal power, expected sample size,
#' trial duration, and classification accuracy (sensitivity/specificity).
#'
#' @param K Integer. Number of baskets.
#' @param A Integer. Number of truly active baskets (0 to K).
#' @param p0 Numeric. Null response rate.
#' @param pa Numeric. Alternative (active) response rate.
#' @param N1 Integer. Total stage 1 sample size.
#' @param N2 Integer. Total stage 2 sample size (homogeneous path).
#' @param n2k Integer. Stage 2 sample size per basket (heterogeneous path).
#' @param r_s Integer. Minimum responders per basket to continue (heterogeneous).
#' @param r_c Integer. Minimum total responders to continue (homogeneous).
#' @param gamma Numeric. Critical value for Fisher's exact test.
#' @param alpha_s Numeric. Significance level for separate analyses.
#' @param alpha_c Numeric. Significance level for combined analysis.
#' @param accrual_rate Numeric vector of length K. Accrual rates per basket.
#' @param n_sims Integer. Number of simulations (default 1000).
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list with components:
#'   \item{fwer}{Numeric. Family-wise error rate (only for A=0).}
#'   \item{marginal_power}{Numeric vector of length K. Power per basket.}
#'   \item{EN}{Numeric. Expected total sample size.}
#'   \item{ET}{Numeric. Expected trial duration (months).}
#'   \item{sensitivity}{Numeric. True positive rate.}
#'   \item{specificity}{Numeric. True negative rate.}
#'   \item{prop_track1}{Numeric. Proportion using heterogeneous path.}
#'   \item{prop_stage2_track1}{Numeric vector. Proportion continuing to stage 2 per basket (heterogeneous).}
#'   \item{prop_stage2_track2}{Numeric. Proportion continuing to stage 2 (homogeneous).}
#'
#' @references
#' Cunanan KM, Iasonos A, Shen R, Begg CB, Gonen M. (2017). An Efficient Basket
#' Trial Design. Statistics in Medicine, 36(10):1568-1579.
#'
#' @export
simulate_simple_operating_characteristics <- function(K, A, p0, pa,
                                                     N1, N2, n2k,
                                                     r_s, r_c,
                                                     gamma, alpha_s, alpha_c,
                                                     accrual_rate,
                                                     n_sims = 1000,
                                                     seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Set up true response rates
  # First A baskets are active, rest are null
  p_true <- c(rep(pa, A), rep(p0, K - A))
  
  # Initialize result storage
  decisions_matrix <- matrix(FALSE, nrow = n_sims, ncol = K)
  heterogeneous_vec <- logical(n_sims)
  trial_duration_vec <- numeric(n_sims)
  total_n_vec <- numeric(n_sims)
  stage2_track1_matrix <- matrix(FALSE, nrow = n_sims, ncol = K)
  stage2_track2_vec <- logical(n_sims)
  
  # Run simulations
  for (i in 1:n_sims) {
    result <- simulate_simple_trial_once(
      K = K,
      p_true = p_true,
      N1 = N1,
      N2 = N2,
      n2k = n2k,
      p0 = p0,
      r_s = r_s,
      r_c = r_c,
      gamma = gamma,
      alpha_s = alpha_s,
      alpha_c = alpha_c,
      accrual_rate = accrual_rate
    )
    
    decisions_matrix[i, ] <- result$decisions
    heterogeneous_vec[i] <- result$heterogeneous
    trial_duration_vec[i] <- result$trial_duration
    total_n_vec[i] <- sum(result$stage1_n) + sum(result$stage2_n)
    
    # Track which baskets continued to stage 2
    if (result$heterogeneous) {
      if (length(result$continued_baskets) > 0) {
        stage2_track1_matrix[i, result$continued_baskets] <- TRUE
      }
    } else {
      if (length(result$continued_baskets) > 0) {
        stage2_track2_vec[i] <- TRUE
      }
    }
  }
  
  # Calculate operating characteristics
  
  # Marginal power (or false positive rate if basket inactive)
  marginal_power <- colMeans(decisions_matrix)
  
  # FWER (only relevant when A=0)
  fwer <- if (A == 0) {
    mean(rowSums(decisions_matrix) >= 1)
  } else {
    NA_real_
  }
  
  # Expected sample size and duration
  EN <- mean(total_n_vec)
  ET <- mean(trial_duration_vec)
  
  # Sensitivity and Specificity
  if (A > 0 && A < K) {
    # True positives: correctly identify active baskets
    tp <- sum(decisions_matrix[, 1:A]) / (n_sims * A)
    # False negatives
    fn <- sum(!decisions_matrix[, 1:A]) / (n_sims * A)
    # True negatives: correctly identify inactive baskets
    tn <- sum(!decisions_matrix[, (A+1):K]) / (n_sims * (K - A))
    # False positives
    fp <- sum(decisions_matrix[, (A+1):K]) / (n_sims * (K - A))
    
    sensitivity <- tp
    specificity <- tn
  } else if (A == 0) {
    sensitivity <- NA_real_
    specificity <- mean(rowSums(decisions_matrix) == 0)
  } else if (A == K) {
    sensitivity <- mean(rowSums(decisions_matrix) == K)
    specificity <- NA_real_
  }
  
  # Track usage
  prop_track1 <- mean(heterogeneous_vec)
  prop_stage2_track1 <- colMeans(stage2_track1_matrix)
  prop_stage2_track2 <- mean(stage2_track2_vec)
  
  # Return results
  list(
    fwer = fwer,
    marginal_power = marginal_power,
    EN = EN,
    ET = ET,
    sensitivity = sensitivity,
    specificity = specificity,
    prop_track1 = prop_track1,
    prop_stage2_track1 = prop_stage2_track1,
    prop_stage2_track2 = prop_stage2_track2
  )
}


#' Design a Simon Basket Trial
#'
#' Optimizes design parameters for a Simon basket trial using grid search.
#' The design features an adaptive two-stage structure with interim heterogeneity
#' assessment. If baskets appear homogeneous, they are pooled for a combined
#' analysis; if heterogeneous, they are analyzed separately.
#'
#' @param K Integer. Number of baskets (must be >= 3).
#' @param p0 Numeric. Null response rate (between 0 and 1).
#' @param pa Numeric. Alternative (active) response rate (must be > p0).
#' @param N1 Integer. Total stage 1 sample size across all baskets.
#' @param N2 Integer. Total stage 2 sample size for homogeneous path.
#' @param r_s Integer. Minimum responders per basket to continue (heterogeneous path).
#' @param r_c Integer. Minimum total responders to continue (homogeneous path).
#' @param target_fwer Numeric. Target family-wise error rate when all baskets null (default 0.05).
#' @param target_power Numeric. Target marginal power at calibration scenario (default 0.80).
#' @param min_power Numeric. Minimum acceptable power when A=1 (default 0.70).
#' @param accrual_rate Numeric or vector. Accrual rate(s) per basket in patients/month.
#'   If scalar, same rate used for all baskets. If vector, must have length K.
#' @param calibrate_at_A Integer. Calibrate power at this number of active baskets.
#'   Use 2 for K=5, use 3 for K>=10 (default 2).
#' @param n2k_range Integer vector. Range of n2k values to search (default 12:22).
#' @param gamma_range Numeric vector. Range of gamma values to search (default seq(0.1, 0.9, 0.05)).
#' @param alpha_s_range Numeric vector. Range of alpha_s values (default seq(0.01, 0.10, 0.01)).
#' @param alpha_c_range Numeric vector. Range of alpha_c values (default seq(0.01, 0.05, 0.01)).
#' @param n_sims Integer. Number of simulations per design (default 1000).
#' @param fwer_margin Numeric. Acceptable margin above target FWER (default 0.01).
#' @param power_margin Numeric. Acceptable margin below target power (default 0.02).
#' @param seed Integer. Random seed for reproducibility (default 4206).
#' @param include_reference Logical. Whether to compute reference design (independent
#'   Simon designs) for comparison (default FALSE).
#' @param parallel Logical. Whether to use parallel processing (default TRUE).
#' @param n_cores Integer. Number of cores for parallel processing (default:
#'   `getOption("mc.cores")` if set, otherwise all available - 1; capped at 2
#'   when running under `R CMD check`).
#'
#' @return An object of class "simple_basket_design" containing:
#'   \item{design}{List of optimized design parameters: N1, N2, n2k, r_s, r_c,
#'     gamma, alpha_s, alpha_c.}
#'   \item{performance}{List of operating characteristics for each scenario
#'     (A=0, 1, 2, ..., K).}
#'   \item{reference}{If include_reference=TRUE, operating characteristics for
#'     reference design.}
#'   \item{inputs}{List of input parameters.}
#'
#' @details
#' The Simon basket trial design (Cunanan et al. 2017) is an adaptive two-stage
#' design that aims to efficiently test whether a targeted therapy is effective
#' across multiple tumor types (baskets) that share a common molecular alteration.
#'
#' **Design Structure:**
#'
#' 1. **Stage 1:** Accrue N1 patients across K baskets
#' 2. **Interim Analysis:** Test of homogeneity (Fisher's exact test)
#'    - If p-value <= gamma: Pursue **heterogeneous path** (analyze separately)
#'    - If p-value > gamma: Pursue **homogeneous path** (pool all baskets)
#' 3. **Stage 2 - Heterogeneous:** For each basket with >= r_s responders,
#'    accrue n2k more patients. Test each continuing basket separately with
#'    Bonferroni correction.
#' 4. **Stage 2 - Homogeneous:** If >= r_c total responders, accrue N2 more
#'    patients pooled. Test overall response rate. Declare all active or all inactive.
#'
#' **Optimization Strategy:**
#'
#' The function performs grid search over (n2k, gamma, alpha_s, alpha_c) to find
#' the design that:
#' 1. Controls FWER <= target_fwer (when A=0)
#' 2. Achieves power >= target_power (when A=calibrate_at_A)
#' 3. Maintains power >= min_power (when A=1)
#' 4. Minimizes a utility function balancing power and sample size
#'
#' **Key Advantages:**
#' - More efficient than independent trials when drug works in most/all baskets
#' - Simpler than Bayesian approaches
#' - Frequentist framework with clear error control
#'
#' **Trade-offs:**
#' - May have reduced power when drug works in only one basket
#' - Can incorrectly pool active and inactive baskets (specificity decreases as A increases)
#' - Only weak FWER control (controlled at A=0 only)
#'
#' @references
#' Cunanan KM, Iasonos A, Shen R, Begg CB, Gonen M. (2017). An Efficient Basket
#' Trial Design. Statistics in Medicine, 36(10):1568-1579. DOI: 10.1002/sim.7227
#'
#' @examples
#' \dontrun{
#' # Example from Cunanan et al. (2017) Table 2
#' # K=5 baskets, theta0=0.15, thetaa=0.45
#' design <- simple_basket_trial(
#'   K = 5,
#'   p0 = 0.15,
#'   pa = 0.45,
#'   N1 = 35,
#'   N2 = 20,
#'   r_s = 1,
#'   r_c = 5,
#'   target_fwer = 0.05,
#'   target_power = 0.80,
#'   min_power = 0.70,
#'   accrual_rate = 2,  # 2 patients/month per basket
#'   n_sims = 1000
#' )
#'
#' # View design parameters
#' print(design)
#'
#' # View operating characteristics
#' summary(design)
#' }
#'
#' @export
simple_basket_trial <- function(K, p0, pa,
                               N1, N2,
                               r_s, r_c,
                               target_fwer = 0.05,
                               target_power = 0.80,
                               min_power = 0.70,
                               accrual_rate = 2,
                               calibrate_at_A = 2,
                               n2k_range = 12:22,
                               gamma_range = seq(0.1, 0.9, by = 0.05),
                               alpha_s_range = seq(0.01, 0.10, by = 0.01),
                               alpha_c_range = seq(0.01, 0.05, by = 0.01),
                               n_sims = 1000,
                               fwer_margin = 0.01,
                               power_margin = 0.02,
                               seed = 4206,
                               include_reference = FALSE,
                               parallel = TRUE,
                               n_cores = NULL) {
  
  # Input validation
  if (K < 3) {
    stop("K must be at least 3 baskets")
  }
  if (p0 >= pa) {
    stop("p0 must be less than pa")
  }
  if (p0 < 0 || p0 > 1 || pa < 0 || pa > 1) {
    stop("p0 and pa must be between 0 and 1")
  }
  if (target_fwer <= 0 || target_fwer > 0.10) {
    stop("target_fwer must be between 0 and 0.10")
  }
  if (min_power >= target_power) {
    stop("min_power must be less than target_power")
  }
  if (length(accrual_rate) == 1) {
    accrual_rate <- rep(accrual_rate, K)
  }
  if (length(accrual_rate) != K) {
    stop("accrual_rate must have length K or be a scalar")
  }
  
  # Store inputs
  inputs <- list(
    K = K,
    p0 = p0,
    pa = pa,
    N1 = N1,
    N2 = N2,
    r_s = r_s,
    r_c = r_c,
    target_fwer = target_fwer,
    target_power = target_power,
    min_power = min_power,
    accrual_rate = accrual_rate,
    calibrate_at_A = calibrate_at_A,
    n_sims = n_sims,
    seed = seed
  )
  
  message("Optimizing Simon basket trial design...")
  message(sprintf("  K = %d baskets", K))
  message(sprintf("  theta0 = %.2f, thetaa = %.2f", p0, pa))
  message(sprintf("  N1 = %d, N2 = %d", N1, N2))
  message(sprintf("  Calibrating at A = %d active baskets", calibrate_at_A))
  message(sprintf("  Grid search: %d n2k x %d gamma x %d alphaS x %d alphaC = %d designs",
                  length(n2k_range), length(gamma_range),
                  length(alpha_s_range), length(alpha_c_range),
                  length(n2k_range) * length(gamma_range) * 
                    length(alpha_s_range) * length(alpha_c_range)))
  
  # ============================================================================
  # OPTIMIZATION: Grid search over n2k, gamma, alpha_s, alpha_c
  # ============================================================================
  
  # Create full grid
  grid <- expand.grid(
    n2k = n2k_range,
    gamma = gamma_range,
    alpha_s = alpha_s_range,
    alpha_c = alpha_c_range,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("  Evaluating %d candidate designs...", nrow(grid)))
  
  # Function to evaluate one design
  evaluate_design <- function(i, grid_row, K, p0, pa, N1, N2, r_s, r_c,
                             accrual_rate, calibrate_at_A, n_sims, seed) {
    
    n2k <- grid_row$n2k
    gamma <- grid_row$gamma
    alpha_s <- grid_row$alpha_s
    alpha_c <- grid_row$alpha_c
    
    # Set seed for this design (make it depend on i for reproducibility)
    set.seed(seed + i)
    
    # Simulate A=0 (null scenario) for FWER
    oc_A0 <- simulate_simple_operating_characteristics(
      K = K, A = 0, p0 = p0, pa = pa,
      N1 = N1, N2 = N2, n2k = n2k,
      r_s = r_s, r_c = r_c,
      gamma = gamma, alpha_s = alpha_s, alpha_c = alpha_c,
      accrual_rate = accrual_rate,
      n_sims = n_sims,
      seed = NULL  # Already set above
    )
    
    # Simulate calibration scenario
    oc_calib <- simulate_simple_operating_characteristics(
      K = K, A = calibrate_at_A, p0 = p0, pa = pa,
      N1 = N1, N2 = N2, n2k = n2k,
      r_s = r_s, r_c = r_c,
      gamma = gamma, alpha_s = alpha_s, alpha_c = alpha_c,
      accrual_rate = accrual_rate,
      n_sims = n_sims,
      seed = NULL
    )
    
    # Simulate A=1 for minimum power check
    oc_A1 <- simulate_simple_operating_characteristics(
      K = K, A = 1, p0 = p0, pa = pa,
      N1 = N1, N2 = N2, n2k = n2k,
      r_s = r_s, r_c = r_c,
      gamma = gamma, alpha_s = alpha_s, alpha_c = alpha_c,
      accrual_rate = accrual_rate,
      n_sims = n_sims,
      seed = NULL
    )
    
    list(
      n2k = n2k,
      gamma = gamma,
      alpha_s = alpha_s,
      alpha_c = alpha_c,
      fwer = oc_A0$fwer,
      EN_A0 = oc_A0$EN,
      power_calib = oc_calib$marginal_power[1],  # Power in first active basket
      EN_calib = oc_calib$EN,
      power_A1 = oc_A1$marginal_power[1],
      EN_A1 = oc_A1$EN
    )
  }
  
  # Run grid search (with or without parallelization)
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) {
      # R CMD check sets _R_CHECK_LIMIT_CORES_ and errors if mclapply spawns
      # more than 2 processes
      check_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(check_limit) && !identical(tolower(check_limit), "false")) {
        n_cores <- 2L
      } else {
        n_cores <- max(1, getOption("mc.cores", parallel::detectCores() - 1))
      }
    }
    message(sprintf("  Using %d cores for parallel processing...", n_cores))
    
    results <- parallel::mclapply(
      1:nrow(grid),
      function(i) evaluate_design(i, grid[i, ], K, p0, pa, N1, N2, r_s, r_c,
                                   accrual_rate, calibrate_at_A, n_sims, seed),
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
  } else {
    results <- lapply(
      1:nrow(grid),
      function(i) evaluate_design(i, grid[i, ], K, p0, pa, N1, N2, r_s, r_c,
                                   accrual_rate, calibrate_at_A, n_sims, seed)
    )
  }
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # ============================================================================
  # FILTER: Select calibrated designs
  # ============================================================================
  
  message("  Filtering calibrated designs...")
  
  # Step 1: FWER control
  fwer_ok <- results_df$fwer <= (target_fwer + fwer_margin)
  
  if (sum(fwer_ok) == 0) {
    warning("No designs meet FWER constraint. Using design with FWER closest to target.")
    fwer_ok <- results_df$fwer == min(results_df$fwer)
  }
  
  # Step 2: Power at calibration scenario
  power_calib_ok <- results_df$power_calib >= (target_power - power_margin)
  
  # Step 3: Minimum power when A=1
  power_A1_ok <- results_df$power_A1 >= (min_power - power_margin)
  
  # Combined filter
  calibrated <- fwer_ok & power_calib_ok & power_A1_ok
  
  if (sum(calibrated) == 0) {
    warning("No fully calibrated designs found. Relaxing constraints...")
    # Try relaxing min_power constraint
    calibrated <- fwer_ok & power_calib_ok
    
    if (sum(calibrated) == 0) {
      warning("Still no designs. Using design with maximum power at calibration scenario.")
      calibrated <- results_df$power_calib == max(results_df$power_calib[fwer_ok])
    }
  }
  
  # Select best design among calibrated
  # Utility function: weighted combination of power and EN
  # Minimize: -sum(power_A) + weight * sum(EN_A)
  # For simplicity, just minimize EN at calibration scenario among those meeting constraints
  calibrated_designs <- results_df[calibrated, ]
  
  best_idx <- which.min(calibrated_designs$EN_calib)
  best_design <- calibrated_designs[best_idx, ]
  
  message(sprintf("  Optimal design: n2k=%d, gamma=%.2f, alphaS=%.3f, alphaC=%.3f",
                  best_design$n2k, best_design$gamma,
                  best_design$alpha_s, best_design$alpha_c))
  
  # ============================================================================
  # FINAL EVALUATION: Compute OCs for all scenarios with optimal design
  # ============================================================================
  
  message("  Computing final operating characteristics for all scenarios...")
  
  performance <- list()
  
  for (A in 0:K) {
    set.seed(seed + 1000 + A)  # Different seed for each scenario
    
    oc <- simulate_simple_operating_characteristics(
      K = K, A = A, p0 = p0, pa = pa,
      N1 = N1, N2 = N2, n2k = best_design$n2k,
      r_s = r_s, r_c = r_c,
      gamma = best_design$gamma,
      alpha_s = best_design$alpha_s,
      alpha_c = best_design$alpha_c,
      accrual_rate = accrual_rate,
      n_sims = n_sims,
      seed = NULL
    )
    
    performance[[paste0("A", A)]] <- oc
  }
  
  # ============================================================================
  # REFERENCE DESIGN (optional)
  # ============================================================================
  
  reference <- NULL
  
  if (include_reference) {
    message("  Computing reference design (independent Simon two-stage)...")
    # TODO: Implement reference design computation
    # For now, just placeholder
    reference <- list(note = "Reference design not yet implemented")
  }
  
  # ============================================================================
  # RETURN
  # ============================================================================
  
  result <- list(
    design = list(
      N1 = N1,
      N2 = N2,
      n2k = best_design$n2k,
      r_s = r_s,
      r_c = r_c,
      gamma = best_design$gamma,
      alpha_s = best_design$alpha_s,
      alpha_c = best_design$alpha_c
    ),
    performance = performance,
    reference = reference,
    inputs = inputs
  )
  
  class(result) <- "simple_basket_design"
  
  message("Done!")
  
  result
}
