#' Optimal Two-Stage Basket Trial Design with Aggregated Futility Analysis
#'
#' @description
#' Calculate optimal two-stage basket trial design with aggregated futility 
#' analysis for multiple tumor indications. This implementation is based on 
#' Jing et al. (2022), which proposes pooling all indications for an early 
#' futility decision while using the pruning and pooling method for final 
#' analysis.
#'
#' @param K Integer. Number of tumor indications (cohorts) in the basket trial.
#' @param p0 Numeric vector of length K. Response rates under the null hypothesis 
#'   for each indication. For homogeneous setting, provide a single value or a 
#'   vector with identical elements.
#' @param pa Numeric vector of length K. Response rates under the alternative 
#'   hypothesis for each indication. Should be > p0.
#' @param alpha Numeric. Global type I error rate (default: 0.05).
#' @param beta Numeric. Expected type II error rate (default: 0.2, i.e., 80% 
#'   expected power).
#' @param constraint_alpha Numeric. Probability of not pruning an inactive 
#'   indication (default: NULL). If specified along with \code{constraint_beta}, 
#'   reduces computational burden by constraining individual indication pruning.
#' @param constraint_beta Numeric. Probability of pruning a truly active 
#'   indication (default: NULL). If specified along with \code{constraint_alpha}, 
#'   uses exact binomial test to determine N and r for each indication.
#' @param n_sims_design Integer. Number of simulations used during design 
#'   optimization to calculate empirical Type I error and power (default: 1000).
#'   Larger values give more accurate optimization but take longer. For initial
#'   exploration, 500-1000 is reasonable; for final designs, use 5000-10000.
#' @param seed Integer. Random seed for reproducibility in design optimization 
#'   (optional).
#'
#' @details
#' This function implements the optimal two-stage basket trial design with 
#' \strong{aggregated futility analysis} as described in Jing et al. (2022). 
#' The key innovation is conducting one pooled futility analysis across all 
#' indications at Stage I, rather than individual futility analyses per indication.
#' 
#' \strong{Design Structure:}
#' 
#' \itemize{
#'   \item \strong{Stage I (Aggregated Futility Analysis):} After S patients 
#'     total are enrolled across all K indications (note: S is the \emph{total} 
#'     across all indications, not per indication), pool all data together and 
#'     test at significance level α₁. If ≥ R₁ responders observed across all 
#'     indications combined, continue to Stage II. Otherwise, stop entire trial 
#'     for futility. The sample size per indication (n₁ₖ) is \emph{not} 
#'     pre-specified and adapts to enrollment rates, with only the total S 
#'     constrained (∑n₁ₖ = S).
#'     
#'   \item \strong{Stage II (Pruning and Pooling):} Enroll additional patients 
#'     so each indication k has Nₖ total patients (note: Nₖ is per indication). 
#'     For final analysis:
#'     \enumerate{
#'       \item \strong{Prune:} Exclude indication k if it has < rₖ responders at 
#'         end of Stage II
#'       \item \strong{Pool:} Test pooled data from remaining (non-pruned) 
#'         indications at significance level α₂
#'       \item \strong{Claim:} Drug is effective in pooled indications if ≥ R₂ 
#'         total responders observed
#'     }
#' }
#' 
#' \strong{Notation Clarification:}
#' \itemize{
#'   \item \strong{S}: Total Stage I sample size \emph{across all K indications}
#'   \item \strong{Nₖ}: Total sample size for indication k (Stage I + Stage II combined)
#'   \item \strong{n₁ₖ}: Stage I sample size for indication k (adapts to enrollment; ∑n₁ₖ = S)
#' }
#' 
#' \strong{Hypotheses:}
#' 
#' \itemize{
#'   \item \strong{Global Null (H₀):} Drug does not work in \emph{any} indication 
#'     (pₖ = p₀ₖ for all k)
#'   \item \strong{Global Alternative (H₁):} Drug works in \emph{at least one} 
#'     indication (pₖ = paₖ for some k, where paₖ > p₀ₖ)
#' }
#' 
#' The design controls:
#' \itemize{
#'   \item \strong{Global Type I error:} P(claim drug works in ≥1 indication | H₀) ≤ α
#'   \item \strong{Expected Type II error:} Average P(fail to detect effect | H₁) ≤ β
#' }
#' 
#' \strong{Optimization:}
#' 
#' Design parameters (S, α₁*, α₂*, r, N) are optimized by minimizing the 
#' expected total sample size under H₀:
#' 
#' \deqn{E[N] = S + \sum_{k=1}^{K}(N_k - S) \cdot P(\sum_{k=1}^{K} X_{1k} \geq R_1 | H_0)}
#' 
#' subject to global Type I error ≤ α and expected Type II error ≤ β.
#' 
#' \strong{Homogeneous vs. Heterogeneous Settings:}
#' 
#' \itemize{
#'   \item \strong{Homogeneous:} All indications have same p₀ and pa. Simpler 
#'     calculations with Nₖ = N and rₖ = r for all k.
#'   \item \strong{Heterogeneous:} Indications may have different p₀ₖ and paₖ 
#'     (e.g., different baseline response rates). More realistic but 
#'     computationally intensive.
#' }
#' 
#' \strong{Computational Constraints:}
#' 
#' For heterogeneous settings with K ≥ 4, the number of possible Stage I 
#' allocations (n₁ₖ) becomes very large. To reduce computational burden, the 
#' function uses:
#' 
#' \enumerate{
#'   \item \strong{Empirical Type I/II error:} 10,000 simulated trials with 
#'     multinomial allocation of S patients to indications
#'   \item \strong{Pruning constraints:} If \code{constraint_alpha} and 
#'     \code{constraint_beta} specified, determines Nₖ and rₖ via exact binomial 
#'     test, reducing unknowns from 3+2K to 5 parameters
#' }
#' 
#' \strong{Advantages over Individual Futility Analysis:}
#' 
#' \itemize{
#'   \item \strong{Earlier decision:} Futility analysis occurs when total S 
#'     patients enrolled (typically ~2-3 months), not waiting for slow-enrolling 
#'     indications (which may take >5 months)
#'   \item \strong{Resource reallocation:} If No-Go decision, resources freed up 
#'     for other compounds sooner
#'   \item \strong{Better TP identification:} Pruning at Stage II with full 
#'     sample size reduces false pruning of active indications
#'   \item \strong{Smaller total N:} Though expected N under H₀ may be slightly 
#'     higher, total sample size is typically lower
#' }
#' 
#' \strong{Trade-offs:}
#' 
#' \itemize{
#'   \item Aggregated analysis may suffer from effect dilution if only a few 
#'     indications are active
#'   \item Slow-enrolling indications may have minimal data at futility analysis
#'   \item Expected sample size under H₀ is typically higher than individual 
#'     futility analysis (price for earlier decision)
#' }
#'
#' @return A list with class "basket_trial_design" containing:
#'   \describe{
#'     \item{param}{Data frame of input parameters: K, p0, pa, alpha, beta}
#'     \item{design}{Data frame with optimal design parameters:
#'       \itemize{
#'         \item \code{S}: Total sample size at Stage I (aggregated futility)
#'         \item \code{alpha1}: Significance level for Stage I test
#'         \item \code{R1}: Critical value for Stage I (continue if ≥ R1 responders)
#'         \item \code{N}: Vector of total sample sizes per indication
#'         \item \code{r}: Vector of pruning thresholds per indication
#'         \item \code{alpha2}: Significance level for Stage II test
#'         \item \code{EN_H0}: Expected total sample size under H₀
#'         \item \code{total_N}: Sum of all Nₖ (maximum possible sample size)
#'       }}
#'     \item{performance}{Data frame with operating characteristics:
#'       \itemize{
#'         \item \code{type1_error}: Empirical global Type I error
#'         \item \code{expected_power}: Empirical expected power
#'         \item \code{prob_continue_h0}: P(continue to Stage II | H₀)
#'       }}
#'   }
#'
#' @references
#' Jing, N., Liu, F., Wu, C., Zhou, H., & Chen, C. (2022). An optimal two-stage 
#' exploratory basket trial design with aggregated futility analysis. 
#' \emph{Contemporary Clinical Trials}, 114, 106296. 
#' \doi{10.1016/j.cct.2022.106296}
#' 
#' Chen, C., Li, X., Yuan, S., Antonijevic, Z., Kalamegham, R., & Beckman, R. A. 
#' (2016). Statistical design and considerations of a phase 3 basket trial for 
#' simultaneous investigation of multiple tumor types in one study. 
#' \emph{Statistics in Biopharmaceutical Research}, 8(3), 248-257.
#'
#' @examples
#' \dontrun{
#' # Example 1: Homogeneous setting (K=4, p0=5%, pa=20%)
#' design_homo <- two_stage_basket_trial(
#'   K = 4,
#'   p0 = c(0.05, 0.05, 0.05, 0.05),
#'   pa = c(0.20, 0.20, 0.20, 0.20),
#'   alpha = 0.05,
#'   beta = 0.2,
#'   constraint_alpha = 0.1,
#'   constraint_beta = 0.2,
#'   seed = 123
#' )
#' print(design_homo)
#' 
#' # Example 2: Heterogeneous setting (K=4, varying p0)
#' # From Table 1 in Jing et al. (2022)
#' design_hetero <- two_stage_basket_trial(
#'   K = 4,
#'   p0 = c(0.05, 0.05, 0.20, 0.20),
#'   pa = c(0.20, 0.20, 0.35, 0.35),
#'   alpha = 0.05,
#'   beta = 0.2,
#'   constraint_alpha = 0.1,
#'   constraint_beta = 0.2,
#'   seed = 123
#' )
#' print(design_hetero)
#' 
#' # Example 3: Small basket (K=3)
#' design_small <- two_stage_basket_trial(
#'   K = 3,
#'   p0 = c(0.01, 0.05, 0.20),
#'   pa = c(0.16, 0.20, 0.35),
#'   alpha = 0.05,
#'   beta = 0.2,
#'   seed = 456
#' )
#' summary(design_small)
#' }
#'
#' @export
two_stage_basket_trial <- function(K,
                                   p0,
                                   pa,
                                   alpha = 0.05,
                                   beta = 0.2,
                                   constraint_alpha = NULL,
                                   constraint_beta = NULL,
                                   n_sims_design = 1000,
                                   seed = NULL) {
  
  # Input validation
  if (K < 2) {
    stop("K must be at least 2 (need multiple indications for basket trial)")
  }
  
  if (length(p0) == 1) {
    p0 <- rep(p0, K)
  }
  if (length(pa) == 1) {
    pa <- rep(pa, K)
  }
  
  if (length(p0) != K || length(pa) != K) {
    stop("p0 and pa must have length K or length 1")
  }
  
  if (any(p0 <= 0) || any(p0 >= 1)) {
    stop("p0 must be between 0 and 1")
  }
  
  if (any(pa <= 0) || any(pa >= 1)) {
    stop("pa must be between 0 and 1")
  }
  
  if (any(pa <= p0)) {
    stop("pa must be greater than p0 for all indications")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  
  if (beta <= 0 || beta >= 1) {
    stop("beta must be between 0 and 1")
  }
  
  # Check if homogeneous setting
  is_homogeneous <- length(unique(p0)) == 1 && length(unique(pa)) == 1
  
  # Use constraints if provided
  use_constraints <- !is.null(constraint_alpha) && !is.null(constraint_beta)
  
  if (use_constraints) {
    if (constraint_alpha <= 0 || constraint_alpha >= 1) {
      stop("constraint_alpha must be between 0 and 1")
    }
    if (constraint_beta <= 0 || constraint_beta >= 1) {
      stop("constraint_beta must be between 0 and 1")
    }
  }
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Optimize design parameters
  cat("Optimizing basket trial design...\n")
  cat(sprintf("  K = %d indications\n", K))
  cat(sprintf("  Setting: %s\n", ifelse(is_homogeneous, "Homogeneous", "Heterogeneous")))
  cat(sprintf("  Target: alpha = %.3f, beta = %.3f\n", alpha, beta))
  
  optimal_design <- .optimize_basket_design(
    K = K,
    p0 = p0,
    pa = pa,
    alpha = alpha,
    beta = beta,
    constraint_alpha = constraint_alpha,
    constraint_beta = constraint_beta,
    n_sims = n_sims_design,
    is_homogeneous = is_homogeneous
  )
  
  # Prepare output
  param <- data.frame(
    K = K,
    p0 = I(list(p0)),
    pa = I(list(pa)),
    alpha = alpha,
    beta = beta,
    homogeneous = is_homogeneous
  )
  
  design <- data.frame(
    S = optimal_design$S,
    alpha1 = optimal_design$alpha1,
    R1 = optimal_design$R1,
    N = I(list(optimal_design$N)),
    r = I(list(optimal_design$r)),
    alpha2 = optimal_design$alpha2,
    EN_H0 = optimal_design$EN_H0,
    total_N = sum(optimal_design$N)
  )
  
  performance <- data.frame(
    type1_error = optimal_design$type1_error,
    expected_power = optimal_design$expected_power,
    prob_continue_h0 = optimal_design$prob_continue_h0
  )
  
  result <- list(
    param = param,
    design = design,
    performance = performance
  )
  
  class(result) <- "basket_trial_design"
  
  return(result)
}


#' Optimize Basket Trial Design Parameters
#'
#' @description
#' Internal function to optimize design parameters using grid search.
#'
#' @keywords internal
#' @noRd
.optimize_basket_design <- function(K, p0, pa, alpha, beta,
                                    constraint_alpha, constraint_beta,
                                    n_sims, is_homogeneous) {
  
  # Define search grids
  # S: total sample size at interim (sum across all indications)
  # Start with rough estimate based on pooled parameters
  p0_pooled <- mean(p0)
  pa_pooled <- mean(pa)
  
  # Single-stage sample size estimate for pooled problem
  n_pooled_est <- ceiling((stats::qnorm(1 - alpha) + stats::qnorm(1 - beta))^2 * 
                          p0_pooled * (1 - p0_pooled) / (pa_pooled - p0_pooled)^2)
  
  # S should be less than total N, roughly 30-60% of total
  # Use coarser grid for speed - can refine around optimal later
  S_min <- max(K * 5, ceiling(0.25 * n_pooled_est * K))
  S_max <- ceiling(0.75 * n_pooled_est * K)
  S_range <- unique(round(seq(S_min, S_max, length.out = 15)))
  
  # alpha1 and alpha2: significance levels
  # Use coarser grid for speed
  alpha_range <- seq(0.05, 0.95, by = 0.05)
  
  # If using constraints, calculate N and r from exact binomial test
  if (!is.null(constraint_alpha) && !is.null(constraint_beta)) {
    # For each indication, find N and r using exact binomial test
    N_vec <- numeric(K)
    r_vec <- numeric(K)
    
    for (k in 1:K) {
      # Find minimum N such that:
      # P(X >= r | p = p0k) <= constraint_alpha (type I error)
      # P(X < r | p = pak) <= constraint_beta (type II error)
      
      # Start with simple sample size estimate
      n_start <- max(10, ceiling((stats::qnorm(1 - constraint_alpha/2) + 
                                   stats::qnorm(1 - constraint_beta))^2 * 
                                  pa[k] * (1 - pa[k]) / (pa[k] - p0[k])^2))
      
      found <- FALSE
      for (N_try in n_start:(n_start + 100)) {
        for (r_try in 0:N_try) {
          alpha_k <- 1 - stats::pbinom(r_try - 1, N_try, p0[k])
          beta_k <- stats::pbinom(r_try - 1, N_try, pa[k])
          
          if (alpha_k <= constraint_alpha && beta_k <= constraint_beta) {
            N_vec[k] <- N_try
            r_vec[k] <- r_try
            found <- TRUE
            break
          }
        }
        if (found) break
      }
      
      if (!found) {
        stop(sprintf("Could not find valid (N, r) for indication %d with given constraints", k))
      }
    }
    
    # Now optimize over S, alpha1, alpha2
    best_EN <- Inf
    best_params <- NULL
    
    cat("  Grid search over S, alpha1, alpha2...\n")
    total_iterations <- length(S_range) * length(alpha_range)^2
    cat(sprintf("  Total combinations to evaluate: %d\n", total_iterations))
    
    iter_count <- 0
    for (S in S_range) {
      if (S >= sum(N_vec)) next  # S must be less than total N
      
      for (alpha1 in alpha_range) {
        for (alpha2 in alpha_range) {
          iter_count <- iter_count + 1
          
          # Progress indicator every 50 iterations
          if (iter_count %% 50 == 0) {
            cat(sprintf("  Progress: %d/%d (%.1f%%), Best EN=%.1f\n", 
                        iter_count, total_iterations, 
                        100 * iter_count / total_iterations, best_EN))
          }
          
          # Calculate R1 and check operating characteristics
          result <- .evaluate_basket_design(
            S = S, alpha1 = alpha1, alpha2 = alpha2,
            N = N_vec, r = r_vec,
            K = K, p0 = p0, pa = pa,
            n_sims = n_sims
          )
          
          # Check if constraints satisfied
          if (result$type1_error <= alpha && result$expected_beta <= beta) {
            if (result$EN_H0 < best_EN) {
              best_EN <- result$EN_H0
              best_params <- list(
                S = S,
                alpha1 = alpha1,
                alpha2 = alpha2,
                R1 = result$R1,
                N = N_vec,
                r = r_vec,
                EN_H0 = result$EN_H0,
                type1_error = result$type1_error,
                expected_power = 1 - result$expected_beta,
                prob_continue_h0 = result$prob_continue_h0
              )
              cat(sprintf("  *** New best found: S=%d, alpha1=%.2f, alpha2=%.2f, EN=%.1f\n",
                          S, alpha1, alpha2, best_EN))
            }
          }
        }
      }
    }
    
  } else {
    # Full optimization without constraints (search over N and r as well)
    # Use adaptive grid search with early pruning for efficiency
    
    cat("  Adaptive grid search with early pruning...\n")
    
    # For simplicity, assume homogeneous N and r if is_homogeneous
    if (is_homogeneous) {
      # PHASE 1: Coarse grid search to find promising region
      cat("  Phase 1: Coarse search...\n")
      
      # Balanced coarsening - not too aggressive
      N_range_coarse <- unique(round(seq(
        from = max(10, ceiling(n_pooled_est * 0.4)),
        to = ceiling(n_pooled_est * 1.3),
        length.out = 4  # Compromise between 3 and 5
      )))
      
      S_range_coarse <- unique(round(seq(S_min, S_max, length.out = 4)))
      alpha_range_coarse <- c(0.20, 0.40, 0.60, 0.80)  # 4 values
      
      best_EN <- Inf
      best_params <- NULL
      best_N <- NULL
      best_r <- NULL
      best_S <- NULL
      
      iter_count <- 0
      
      # Use fewer sims in Phase 1 for speed (but not too few)
      n_sims_phase1 <- max(25, floor(n_sims / 2))
      
      for (N in N_range_coarse) {
        # For each N, try a few r values
        r_min <- max(1, floor(N * p0_pooled * 0.5))
        r_max <- min(N - 1, ceiling(N * pa_pooled * 1.5))
        r_range_coarse <- unique(round(seq(r_min, r_max, length.out = 3)))
        
        for (r in r_range_coarse) {
          for (S in S_range_coarse) {
            if (S >= N * K) next
            
            for (alpha1 in alpha_range_coarse) {
              for (alpha2 in alpha_range_coarse) {
                iter_count <- iter_count + 1
                
                # Use reduced simulations in Phase 1
                result <- .evaluate_basket_design(
                  S = S, alpha1 = alpha1, alpha2 = alpha2,
                  N = rep(N, K), r = rep(r, K),
                  K = K, p0 = p0, pa = pa,
                  n_sims = n_sims_phase1
                )
                
                if (result$type1_error <= alpha && result$expected_beta <= beta) {
                  if (result$EN_H0 < best_EN) {
                    best_EN <- result$EN_H0
                    best_N <- N
                    best_r <- r
                    best_S <- S
                    best_params <- list(
                      S = S,
                      alpha1 = alpha1,
                      alpha2 = alpha2,
                      R1 = result$R1,
                      N = rep(N, K),
                      r = rep(r, K),
                      EN_H0 = result$EN_H0,
                      type1_error = result$type1_error,
                      expected_power = 1 - result$expected_beta,
                      prob_continue_h0 = result$prob_continue_h0
                    )
                    cat(sprintf("  *** Phase 1 best: S=%d, N=%d, r=%d, alpha1=%.2f, alpha2=%.2f, EN=%.1f\n",
                                S, N, r, alpha1, alpha2, best_EN))
                  }
                }
              }
            }
          }
        }
      }
      
      cat(sprintf("  Phase 1 complete: %d evaluations\n", iter_count))
      
      if (is.null(best_params)) {
        stop("Could not find feasible design in coarse search. Try relaxing alpha/beta constraints.")
      }
      
      # PHASE 2: Fine search around best solution
      cat("  Phase 2: Refining around best solution...\n")
      
      # Define narrow search ranges around best solution
      N_range_fine <- unique(round(seq(
        max(10, best_N - 3),
        best_N + 3,
        by = 1
      )))
      
      S_range_fine <- unique(round(seq(
        max(K * 5, best_S - 5),
        min(best_N * K - 1, best_S + 5),
        by = 1
      )))
      
      alpha_range_fine <- seq(
        max(0.05, best_params$alpha1 - 0.15),
        min(0.95, best_params$alpha1 + 0.15),
        by = 0.05
      )
      
      for (N in N_range_fine) {
        r_min <- max(1, floor(N * p0_pooled * 0.5))
        r_max <- min(N - 1, ceiling(N * pa_pooled * 1.5))
        r_range_fine <- unique(round(seq(
          max(r_min, best_r - 2),
          min(r_max, best_r + 2),
          by = 1
        )))
        
        for (r in r_range_fine) {
          for (S in S_range_fine) {
            if (S >= N * K) next
            
            for (alpha1 in alpha_range_fine) {
              for (alpha2 in alpha_range_fine) {
                iter_count <- iter_count + 1
                
                # Early stopping: if this combination can't beat best_EN, skip
                # Heuristic: EN_H0 is at least S (always enroll Stage I)
                if (S > best_EN) next
                
                result <- .evaluate_basket_design(
                  S = S, alpha1 = alpha1, alpha2 = alpha2,
                  N = rep(N, K), r = rep(r, K),
                  K = K, p0 = p0, pa = pa,
                  n_sims = n_sims
                )
                
                if (result$type1_error <= alpha && result$expected_beta <= beta) {
                  if (result$EN_H0 < best_EN) {
                    best_EN <- result$EN_H0
                    best_params <- list(
                      S = S,
                      alpha1 = alpha1,
                      alpha2 = alpha2,
                      R1 = result$R1,
                      N = rep(N, K),
                      r = rep(r, K),
                      EN_H0 = result$EN_H0,
                      type1_error = result$type1_error,
                      expected_power = 1 - result$expected_beta,
                      prob_continue_h0 = result$prob_continue_h0
                    )
                    cat(sprintf("  *** Phase 2 improved: S=%d, N=%d, r=%d, alpha1=%.2f, alpha2=%.2f, EN=%.1f\n",
                                S, N, r, alpha1, alpha2, best_EN))
                  }
                }
              }
            }
          }
        }
      }
      
      cat(sprintf("  Phase 2 complete: Total evaluations=%d\n", iter_count))
      
    } else {
      stop("Full optimization for heterogeneous setting without constraints is not implemented. Please provide constraint_alpha and constraint_beta.")
    }
  }
  
  if (is.null(best_params)) {
    stop("Could not find design satisfying alpha and beta constraints. Try relaxing constraints or adjusting search ranges.")
  }
  
  cat(sprintf("  Optimal design found: S=%d, EN(H0)=%.1f\n", 
              best_params$S, best_params$EN_H0))
  
  return(best_params)
}


#' Evaluate Basket Trial Design
#'
#' @description
#' Internal function to calculate operating characteristics for a given set of
#' design parameters using simulation.
#'
#' @keywords internal
#' @noRd
.evaluate_basket_design <- function(S, alpha1, alpha2, N, r, K, p0, pa, n_sims) {
  
  # Vectorized simulation for speed
  # Pre-generate all random allocations and responses at once
  
  # Simulate Stage I allocations (n_sims x K matrix)
  n1_sims <- t(stats::rmultinom(n_sims, S, rep(1/K, K)))
  
  # Under H0: calculate total responses for all simulations at once
  # Generate all responses in one vectorized call
  responses_h0_matrix <- matrix(0, nrow = n_sims, ncol = K)
  for (k in 1:K) {
    responses_h0_matrix[, k] <- stats::rbinom(n_sims, n1_sims[, k], p0[k])
  }
  responses_h0_total <- rowSums(responses_h0_matrix)
  
  # R1: critical value for Stage I
  R1 <- ceiling(stats::quantile(responses_h0_total, 1 - alpha1))
  
  # Vectorized Type I error calculation
  continue_to_stage2_h0 <- responses_h0_total >= R1
  n_continue_h0 <- sum(continue_to_stage2_h0)
  
  if (n_continue_h0 == 0) {
    # All trials stop early under H0 - perfect Type I error control
    type1_error <- 0
    prob_continue_h0 <- 0
  } else {
    # For trials that continue, simulate Stage II
    idx_continue <- which(continue_to_stage2_h0)
    n1_continue <- n1_sims[idx_continue, , drop = FALSE]
    X1_continue <- responses_h0_matrix[idx_continue, , drop = FALSE]
    
    # Generate Stage II responses only for continuers
    X2_matrix <- matrix(0, nrow = length(idx_continue), ncol = K)
    for (k in 1:K) {
      n_additional <- pmax(0, N[k] - n1_continue[, k])
      X2_matrix[, k] <- ifelse(n_additional > 0, 
                               stats::rbinom(length(idx_continue), n_additional, p0[k]),
                               0)
    }
    
    X_total <- pmin(X1_continue + X2_matrix, matrix(N, nrow = length(idx_continue), ncol = K, byrow = TRUE))
    
    # Vectorized pruning and pooled test
    type1_count <- 0
    for (i in 1:length(idx_continue)) {
      pooled_idx <- which(X_total[i, ] >= r)
      if (length(pooled_idx) > 0) {
        total_pooled_resp <- sum(X_total[i, pooled_idx])
        total_pooled_n <- sum(N[pooled_idx])
        p0_pooled <- sum(N[pooled_idx] * p0[pooled_idx]) / total_pooled_n
        R2 <- stats::qbinom(1 - alpha2, total_pooled_n, p0_pooled)
        if (total_pooled_resp > R2) {
          type1_count <- type1_count + 1
        }
      }
    }
    
    type1_error <- type1_count / n_sims
    prob_continue_h0 <- n_continue_h0 / n_sims
  }
  
  # Calculate expected Type II error (simplified for speed)
  # Sample configurations instead of exhaustive
  max_configs_per_G <- min(5, choose(K, ceiling(K/2)))  # Reduced from 10
  
  total_beta <- 0
  count_configs <- 0
  
  for (G in 1:K) {
    # Sample configurations
    if (choose(K, G) <= max_configs_per_G) {
      active_configs <- utils::combn(K, G)
      n_configs <- ncol(active_configs)
    } else {
      n_configs <- max_configs_per_G
      active_configs <- replicate(n_configs, sample(K, G))
      if (G == 1) {
        active_configs <- matrix(active_configs, nrow = 1)
      }
    }
    
    for (cfg in 1:n_configs) {
      active_idx <- active_configs[, cfg]
      g <- rep(0, K)
      g[active_idx] <- 1
      p_true <- ifelse(g == 1, pa, p0)
      
      # Vectorized H1 simulation
      X1_h1_matrix <- matrix(0, nrow = n_sims, ncol = K)
      for (k in 1:K) {
        X1_h1_matrix[, k] <- stats::rbinom(n_sims, pmin(n1_sims[, k], N[k]), p_true[k])
      }
      
      total_resp_stage1 <- rowSums(X1_h1_matrix)
      continue_to_stage2 <- total_resp_stage1 >= R1
      
      reject_count <- 0
      idx_continue_h1 <- which(continue_to_stage2)
      
      if (length(idx_continue_h1) > 0) {
        n1_continue_h1 <- n1_sims[idx_continue_h1, , drop = FALSE]
        X1_continue_h1 <- X1_h1_matrix[idx_continue_h1, , drop = FALSE]
        
        X2_h1_matrix <- matrix(0, nrow = length(idx_continue_h1), ncol = K)
        for (k in 1:K) {
          n_additional <- pmax(0, N[k] - n1_continue_h1[, k])
          X2_h1_matrix[, k] <- ifelse(n_additional > 0,
                                      stats::rbinom(length(idx_continue_h1), n_additional, p_true[k]),
                                      0)
        }
        
        X_total_h1 <- X1_continue_h1 + X2_h1_matrix
        
        for (i in 1:length(idx_continue_h1)) {
          pooled_idx <- which(X_total_h1[i, ] >= r)
          if (length(pooled_idx) > 0) {
            total_pooled_resp <- sum(X_total_h1[i, pooled_idx])
            total_pooled_n <- sum(N[pooled_idx])
            p0_pooled <- sum(N[pooled_idx] * p0[pooled_idx]) / total_pooled_n
            R2 <- stats::qbinom(1 - alpha2, total_pooled_n, p0_pooled)
            if (total_pooled_resp > R2) {
              reject_count <- reject_count + 1
            }
          }
        }
      }
      
      power_this_config <- reject_count / n_sims
      beta_this_config <- 1 - power_this_config
      total_beta <- total_beta + beta_this_config
      count_configs <- count_configs + 1
    }
  }
  
  expected_beta <- total_beta / count_configs
  
  # Calculate expected sample size under H0
  # E[N|H0] = S + Pr(continue|H0) * E[additional patients|continue]
  additional_under_h0 <- sum(N) - S
  EN_H0 <- S + prob_continue_h0 * additional_under_h0
  
  return(list(
    R1 = R1,
    type1_error = type1_error,
    expected_beta = expected_beta,
    EN_H0 = EN_H0,
    prob_continue_h0 = prob_continue_h0
  ))
}


#' Simulate a Basket Trial from a Design Object
#'
#' @description
#' Simulates a single realization of a basket trial based on a design object 
#' created by `two_stage_basket_trial()`.
#'
#' @param design A basket trial design object from `two_stage_basket_trial()`.
#' @param p_true Numeric vector of length K. True response rates for each 
#'   indication under simulation. Use p0 for null scenario or pa for 
#'   alternative scenario, or mix for partial response.
#' @param allocation_probs Numeric vector of length K. Probabilities for 
#'   allocating patients to each indication at Stage I (default: equal 
#'   allocation). Must sum to 1.
#'
#' @return A list containing:
#'   \item{reject_h0}{Logical; TRUE if null hypothesis is rejected}
#'   \item{stopped_early}{Logical; TRUE if trial stopped at Stage I for futility}
#'   \item{stage1_n}{Vector of length K; patients enrolled per indication at Stage I}
#'   \item{stage1_responses}{Vector of length K; responses per indication at Stage I}
#'   \item{stage1_total_resp}{Integer; total responses across all indications at Stage I}
#'   \item{stage1_decision}{Character; "Continue" or "Stop for futility"}
#'   \item{stage2_total_n}{Vector of length K; total patients per indication (NA if stopped early)}
#'   \item{stage2_total_resp}{Vector of length K; total responses per indication (NA if stopped early)}
#'   \item{pruned_indications}{Vector of indices of pruned indications (NA if stopped early)}
#'   \item{pooled_indications}{Vector of indices of pooled indications (NA if stopped early)}
#'   \item{final_decision}{Character; "Reject H0", "Fail to reject H0", or "Stopped early"}
#'
#' @examples
#' \dontrun{
#' # Create design
#' design <- two_stage_basket_trial(
#'   K = 4,
#'   p0 = c(0.05, 0.05, 0.10, 0.10),
#'   pa = c(0.20, 0.20, 0.25, 0.25),
#'   alpha = 0.05,
#'   beta = 0.2,
#'   constraint_alpha = 0.1,
#'   constraint_beta = 0.2,
#'   seed = 123
#' )
#' 
#' # Simulate under H0 (all null)
#' trial_h0 <- simulate_basket_trial(
#'   design = design,
#'   p_true = design$param$p0[[1]]
#' )
#' 
#' # Simulate under H1 (all alternative)
#' trial_h1 <- simulate_basket_trial(
#'   design = design,
#'   p_true = design$param$pa[[1]]
#' )
#' 
#' # Simulate under partial response (2 active, 2 inactive)
#' p_partial <- c(0.20, 0.20, 0.10, 0.10)
#' trial_partial <- simulate_basket_trial(
#'   design = design,
#'   p_true = p_partial
#' )
#' }
#'
#' @export
simulate_basket_trial <- function(design, p_true, allocation_probs = NULL) {
  
  # Validate input
  if (!inherits(design, "basket_trial_design")) {
    stop("design must be a basket_trial_design object from two_stage_basket_trial()")
  }
  
  K <- design$param$K
  p0 <- design$param$p0[[1]]
  
  if (length(p_true) != K) {
    stop("p_true must have length K")
  }
  
  if (is.null(allocation_probs)) {
    allocation_probs <- rep(1/K, K)
  }
  
  if (length(allocation_probs) != K) {
    stop("allocation_probs must have length K")
  }
  
  if (abs(sum(allocation_probs) - 1) > 1e-6) {
    stop("allocation_probs must sum to 1")
  }
  
  # Extract design parameters
  S <- design$design$S
  R1 <- design$design$R1
  alpha2 <- design$design$alpha2
  N <- design$design$N[[1]]
  r <- design$design$r[[1]]
  
  # === STAGE I: Aggregated Futility Analysis ===
  
  # Allocate S patients to K indications
  stage1_n <- as.vector(stats::rmultinom(1, S, allocation_probs))
  
  # Generate responses
  stage1_responses <- stats::rbinom(K, stage1_n, p_true)
  stage1_total_resp <- sum(stage1_responses)
  
  # Decision at Stage I
  if (stage1_total_resp < R1) {
    # Stop for futility
    return(list(
      reject_h0 = FALSE,
      stopped_early = TRUE,
      stage1_n = stage1_n,
      stage1_responses = stage1_responses,
      stage1_total_resp = stage1_total_resp,
      stage1_decision = "Stop for futility",
      stage2_total_n = NA,
      stage2_total_resp = NA,
      pruned_indications = NA,
      pooled_indications = NA,
      final_decision = "Stopped early"
    ))
  }
  
  # === STAGE II: Pruning and Pooling ===
  
  # Enroll additional patients
  stage2_additional_n <- N - stage1_n
  stage2_additional_resp <- stats::rbinom(K, stage2_additional_n, p_true)
  
  # Total across both stages
  stage2_total_n <- N
  stage2_total_resp <- stage1_responses + stage2_additional_resp
  
  # Pruning step
  pruned_indications <- which(stage2_total_resp < r)
  pooled_indications <- which(stage2_total_resp >= r)
  
  # Check if any indications remain for pooling
  if (length(pooled_indications) == 0) {
    # All pruned, fail to reject H0
    return(list(
      reject_h0 = FALSE,
      stopped_early = FALSE,
      stage1_n = stage1_n,
      stage1_responses = stage1_responses,
      stage1_total_resp = stage1_total_resp,
      stage1_decision = "Continue",
      stage2_total_n = stage2_total_n,
      stage2_total_resp = stage2_total_resp,
      pruned_indications = pruned_indications,
      pooled_indications = pooled_indications,
      final_decision = "Fail to reject H0 (all pruned)"
    ))
  }
  
  # Pooled test
  pooled_total_resp <- sum(stage2_total_resp[pooled_indications])
  pooled_total_n <- sum(N[pooled_indications])
  
  # Calculate R2: critical value for pooled test
  # Weighted average of p0 for pooled indications
  p0_pooled <- sum(N[pooled_indications] * p0[pooled_indications]) / pooled_total_n
  R2 <- stats::qbinom(1 - alpha2, pooled_total_n, p0_pooled)
  
  # Final decision
  reject <- pooled_total_resp > R2
  
  return(list(
    reject_h0 = reject,
    stopped_early = FALSE,
    stage1_n = stage1_n,
    stage1_responses = stage1_responses,
    stage1_total_resp = stage1_total_resp,
    stage1_decision = "Continue",
    stage2_total_n = stage2_total_n,
    stage2_total_resp = stage2_total_resp,
    pruned_indications = pruned_indications,
    pooled_indications = pooled_indications,
    pooled_total_resp = pooled_total_resp,
    pooled_R2 = R2,
    final_decision = ifelse(reject, "Reject H0", "Fail to reject H0")
  ))
}


#' Simulate Operating Characteristics of Basket Trial Design
#'
#' @description
#' Runs multiple simulations to estimate operating characteristics for a 
#' basket trial design, including Type I error, expected power, probability 
#' of identifying true positives, and expected sample sizes.
#'
#' @param design A basket trial design object from `two_stage_basket_trial()`.
#' @param n_sims Integer. Number of simulations to run (default: 1000).
#' @param allocation_probs Numeric vector of length K. Probabilities for 
#'   allocating patients to each indication at Stage I (default: equal). 
#'   Must sum to 1.
#' @param seed Integer. Random seed for reproducibility (optional).
#'
#' @return A list containing:
#'   \item{type1_error}{Estimated global Type I error rate}
#'   \item{expected_power}{Estimated expected power (averaged over all 
#'     alternative configurations)}
#'   \item{power_by_n_active}{Data frame with power for each number of active 
#'     indications (1 to K)}
#'   \item{expected_n_h0}{Expected total sample size under H₀}
#'   \item{expected_n_h1}{Expected total sample size under H₁ (averaged)}
#'   \item{prob_early_stop_h0}{Probability of early stopping under H₀}
#'   \item{prob_early_stop_h1}{Probability of early stopping under H₁ (averaged)}
#'   \item{prob_identify_1plus_tp}{Probability of identifying ≥1 true positive}
#'   \item{prob_identify_2plus_tp}{Probability of identifying ≥2 true positives}
#'   \item{prob_identify_all_tp}{Probability of identifying all true positives}
#'   \item{expected_n_true_positives}{Expected number of true positives identified}
#'   \item{expected_n_false_positives}{Expected number of false positives identified}
#'   \item{summary_by_config}{Data frame with detailed results by configuration}
#'
#' @details
#' This function comprehensively evaluates the basket trial design by:
#' 
#' \enumerate{
#'   \item Simulating under H₀ (all indications inactive) to estimate Type I error
#'   \item Simulating under all 2^K - 1 possible alternative configurations to 
#'     estimate expected power, true positive rates, etc.
#'   \item Averaging metrics across configurations using uniform weighting
#' }
#' 
#' The function provides both aggregate metrics (e.g., expected power) and 
#' configuration-specific results for detailed analysis.
#'
#' @examples
#' \dontrun{
#' # Create design
#' design <- two_stage_basket_trial(
#'   K = 4,
#'   p0 = c(0.05, 0.05, 0.10, 0.10),
#'   pa = c(0.20, 0.20, 0.25, 0.25),
#'   alpha = 0.05,
#'   beta = 0.2,
#'   constraint_alpha = 0.1,
#'   constraint_beta = 0.2,
#'   seed = 123
#' )
#' 
#' # Evaluate operating characteristics
#' oc <- simulate_basket_operating_characteristics(
#'   design = design,
#'   n_sims = 5000,
#'   seed = 456
#' )
#' 
#' cat("Type I Error:", oc$type1_error, "\n")
#' cat("Expected Power:", oc$expected_power, "\n")
#' cat("E[N | H0]:", oc$expected_n_h0, "\n")
#' cat("P(≥1 TP):", oc$prob_identify_1plus_tp, "\n")
#' 
#' # View power by number of active indications
#' print(oc$power_by_n_active)
#' 
#' # View detailed results by configuration
#' head(oc$summary_by_config)
#' }
#'
#' @export
simulate_basket_operating_characteristics <- function(design,
                                                       n_sims = 1000,
                                                       allocation_probs = NULL,
                                                       seed = NULL) {
  
  if (!inherits(design, "basket_trial_design")) {
    stop("design must be a basket_trial_design object from two_stage_basket_trial()")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  K <- design$param$K
  p0 <- design$param$p0[[1]]
  pa <- design$param$pa[[1]]
  
  # === Simulate under H0 (all indications inactive) ===
  cat("Simulating under H0...\n")
  
  results_h0 <- replicate(n_sims, {
    simulate_basket_trial(design, p_true = p0, allocation_probs = allocation_probs)
  }, simplify = FALSE)
  
  type1_error <- mean(sapply(results_h0, function(x) x$reject_h0))
  expected_n_h0 <- mean(sapply(results_h0, function(x) sum(x$stage1_n) + 
                                  ifelse(x$stopped_early, 0, sum(x$stage2_total_n) - sum(x$stage1_n))))
  prob_early_stop_h0 <- mean(sapply(results_h0, function(x) x$stopped_early))
  
  # === Simulate under all alternative configurations ===
  cat("Simulating under H1 (all alternative configurations)...\n")
  
  # Store results for each configuration
  config_results <- list()
  config_index <- 1
  
  for (G in 1:K) {  # Number of active indications
    active_configs <- utils::combn(K, G)
    n_configs <- ncol(active_configs)
    
    for (cfg in 1:n_configs) {
      active_idx <- active_configs[, cfg]
      g <- rep(0, K)
      g[active_idx] <- 1
      
      # True response rates for this configuration
      p_true <- ifelse(g == 1, pa, p0)
      
      # Simulate n_sims trials
      results <- replicate(n_sims, {
        simulate_basket_trial(design, p_true = p_true, 
                            allocation_probs = allocation_probs)
      }, simplify = FALSE)
      
      # Calculate metrics for this configuration
      power <- mean(sapply(results, function(x) x$reject_h0))
      expected_n <- mean(sapply(results, function(x) {
        sum(x$stage1_n) + ifelse(x$stopped_early, 0, 
                                  sum(x$stage2_total_n) - sum(x$stage1_n))
      }))
      prob_early_stop <- mean(sapply(results, function(x) x$stopped_early))
      
      # True positive and false positive analysis
      # TP: correctly identified active indication
      # FP: incorrectly identified inactive indication
      
      tp_counts <- sapply(results, function(x) {
        if (x$stopped_early || length(x$pooled_indications) == 0) return(0)
        sum(g[x$pooled_indications] == 1)
      })
      
      fp_counts <- sapply(results, function(x) {
        if (x$stopped_early || length(x$pooled_indications) == 0) return(0)
        sum(g[x$pooled_indications] == 0)
      })
      
      prob_1plus_tp <- mean(tp_counts >= 1)
      prob_2plus_tp <- mean(tp_counts >= 2)
      prob_all_tp <- mean(tp_counts == G)
      expected_tp <- mean(tp_counts)
      expected_fp <- mean(fp_counts)
      
      # Store results
      config_results[[config_index]] <- list(
        n_active = G,
        active_indications = paste(active_idx, collapse = ","),
        power = power,
        expected_n = expected_n,
        prob_early_stop = prob_early_stop,
        prob_1plus_tp = prob_1plus_tp,
        prob_2plus_tp = prob_2plus_tp,
        prob_all_tp = prob_all_tp,
        expected_tp = expected_tp,
        expected_fp = expected_fp
      )
      
      config_index <- config_index + 1
    }
  }
  
  # Convert to data frame
  summary_by_config <- do.call(rbind, lapply(config_results, as.data.frame))
  
  # Calculate overall expected metrics (uniform average across configs)
  expected_power <- mean(summary_by_config$power)
  expected_n_h1 <- mean(summary_by_config$expected_n)
  prob_early_stop_h1 <- mean(summary_by_config$prob_early_stop)
  prob_identify_1plus_tp <- mean(summary_by_config$prob_1plus_tp)
  prob_identify_2plus_tp <- mean(summary_by_config$prob_2plus_tp)
  prob_identify_all_tp <- mean(summary_by_config$prob_all_tp)
  expected_n_true_positives <- mean(summary_by_config$expected_tp)
  expected_n_false_positives <- mean(summary_by_config$expected_fp)
  
  # Aggregate by number of active indications
  power_by_n_active <- stats::aggregate(
    power ~ n_active,
    data = summary_by_config,
    FUN = mean
  )
  names(power_by_n_active) <- c("n_active_indications", "power")
  
  return(list(
    type1_error = type1_error,
    expected_power = expected_power,
    power_by_n_active = power_by_n_active,
    expected_n_h0 = expected_n_h0,
    expected_n_h1 = expected_n_h1,
    prob_early_stop_h0 = prob_early_stop_h0,
    prob_early_stop_h1 = prob_early_stop_h1,
    prob_identify_1plus_tp = prob_identify_1plus_tp,
    prob_identify_2plus_tp = prob_identify_2plus_tp,
    prob_identify_all_tp = prob_identify_all_tp,
    expected_n_true_positives = expected_n_true_positives,
    expected_n_false_positives = expected_n_false_positives,
    summary_by_config = summary_by_config
  ))
}

