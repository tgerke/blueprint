#' Print method for two_stage_single_arm_tte results
#'
#' @param x Output from \code{two_stage_single_arm_tte}
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.two_stage_design <- function(x, ...) {
  cat("Optimal Two-Stage Design for Time-to-Event Endpoint\n")
  cat("====================================================\n\n")
  
  cat("Input Parameters:\n")
  print(x$param)
  cat("\n")
  
  cat("Single-Stage Design:\n")
  cat(sprintf("  Sample size: %d\n", x$single_stage$nsingle))
  cat(sprintf("  Accrual time: %.2f\n", x$single_stage$tasingle))
  cat(sprintf("  Critical value: %.4f\n", x$single_stage$csingle))
  cat("\n")
  
  cat("Two-Stage Design:\n")
  cat(sprintf("  Stage 1 sample size: %d\n", x$two_stage$n1))
  cat(sprintf("  Stage 1 critical value: %.4f (stop if Z1 <= %.4f)\n", 
              x$two_stage$c1, x$two_stage$c1))
  cat(sprintf("  Total sample size: %d\n", x$two_stage$n))
  cat(sprintf("  Final critical value: %.4f (reject H0 if Z > %.4f)\n", 
              x$two_stage$c, x$two_stage$c))
  cat(sprintf("  Interim analysis time: %.2f\n", x$two_stage$t1))
  cat(sprintf("  Maximum total study length: %.2f\n", x$two_stage$MTSL))
  cat(sprintf("  Expected sample size under H0: %.2f\n", x$two_stage$ES))
  cat(sprintf("  Probability of early stopping under H0: %.4f\n", x$two_stage$PS))
  
  invisible(x)
}


#' Summary method for two_stage_single_arm_tte results
#'
#' @param object Output from \code{two_stage_single_arm_tte}
#' @param ... Additional arguments (not used)
#'
#' @return A summary object
#' @export
summary.two_stage_design <- function(object, ...) {
  structure(
    list(
      param = object$param,
      single_stage = object$single_stage,
      two_stage = object$two_stage,
      efficiency = object$two_stage$ES / object$single_stage$nsingle
    ),
    class = "summary.two_stage_design"
  )
}


#' Print summary for two-stage design
#'
#' @param x Summary object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.summary.two_stage_design <- function(x, ...) {
  cat("Summary of Two-Stage Design\n")
  cat("===========================\n\n")
  
  cat("Design Parameters:\n")
  if ("shape" %in% names(x$param)) {
    cat(sprintf("  Distribution shape parameter: %.5f\n", x$param$shape))
    cat(sprintf("  Survival prob at time %.2f: %.2f\n", x$param$x0, x$param$S0))
  }
  cat(sprintf("  Hazard ratio (alternative): %.4f\n", x$param$hr))
  cat(sprintf("  Restricted follow-up: %.2f\n", x$param$x))
  cat(sprintf("  Accrual rate: %.2f per time unit\n", x$param$rate))
  cat(sprintf("  Type I error: %.3f\n", x$param$alpha))
  cat(sprintf("  Power: %.1f%%\n", (1 - x$param$beta) * 100))
  cat("\n")
  
  cat("Sample Size Comparison:\n")
  cat(sprintf("  Single-stage: %d patients\n", x$single_stage$nsingle))
  cat(sprintf("  Two-stage (maximum): %d patients\n", x$two_stage$n))
  cat(sprintf("  Two-stage (expected under H0): %.1f patients\n", x$two_stage$ES))
  cat(sprintf("  Efficiency (ES/n_single): %.2f%%\n", x$efficiency * 100))
  cat("\n")
  
  cat("Operating Characteristics:\n")
  cat(sprintf("  Probability of early stop under H0: %.1f%%\n", 
              x$two_stage$PS * 100))
  cat(sprintf("  Stage 1 sample size: %d (%.1f%% of total)\n", 
              x$two_stage$n1, x$two_stage$n1 / x$two_stage$n * 100))
  
  invisible(x)
}


#' Print method for basket_trial_design results
#'
#' @param x Output from \code{two_stage_basket_trial}
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.basket_trial_design <- function(x, ...) {
  cat("Optimal Two-Stage Basket Trial Design (Aggregated Futility Analysis)\n")
  cat("===================================================================\n\n")
  
  K <- x$param$K
  p0 <- x$param$p0[[1]]
  pa <- x$param$pa[[1]]
  
  cat("Trial Configuration:\n")
  cat(sprintf("  Number of indications: %d\n", K))
  cat(sprintf("  Setting: %s\n", ifelse(x$param$homogeneous, "Homogeneous", "Heterogeneous")))
  cat(sprintf("  Global Type I error: %.3f\n", x$param$alpha))
  cat(sprintf("  Expected Type II error: %.3f (power: %.1f%%)\n", 
              x$param$beta, (1 - x$param$beta) * 100))
  cat("\n")
  
  cat("Response Rates:\n")
  cat("  Indication  p0     pa    \n")
  cat("  ----------  -----  -----\n")
  for (i in 1:K) {
    cat(sprintf("  %10d  %.3f  %.3f\n", i, p0[i], pa[i]))
  }
  cat("\n")
  
  N <- x$design$N[[1]]
  r <- x$design$r[[1]]
  
  cat("Design Parameters:\n")
  cat("  Stage I (Aggregated Futility Analysis):\n")
  cat(sprintf("    Total sample size across all indications: S = %d patients\n", x$design$S))
  cat(sprintf("    Significance level: alpha1 = %.3f\n", x$design$alpha1))
  cat(sprintf("    Critical value: R1 = %d (continue if >=%d total responders)\n", 
              x$design$R1, x$design$R1))
  cat("    Note: Allocation per indication adapts to enrollment rates\n")
  cat("\n")
  
  cat("  Stage II (Pruning and Pooling):\n")
  cat("    Sample sizes per indication (total including Stage I):\n")
  for (i in 1:K) {
    cat(sprintf("      Indication %d: N=%d, r=%d (prune if <%d responders at Stage II)\n", 
                i, N[i], r[i], r[i]))
  }
  cat(sprintf("    Total maximum sample size: %d patients\n", x$design$total_N))
  cat(sprintf("    Significance level: alpha2 = %.3f\n", x$design$alpha2))
  cat("\n")
  
  cat("Operating Characteristics:\n")
  cat(sprintf("  Expected sample size under H0: %.1f\n", x$design$EN_H0))
  cat(sprintf("  Actual Type I error: %.4f\n", x$performance$type1_error))
  cat(sprintf("  Actual expected power: %.4f\n", x$performance$expected_power))
  cat(sprintf("  P(continue to Stage II | H0): %.4f\n", 
              x$performance$prob_continue_h0))
  
  invisible(x)
}


#' Summary method for basket_trial_design results
#'
#' @param object Output from \code{two_stage_basket_trial}
#' @param ... Additional arguments (not used)
#'
#' @return A summary object
#' @export
summary.basket_trial_design <- function(object, ...) {
  structure(
    list(
      param = object$param,
      design = object$design,
      performance = object$performance,
      efficiency = object$design$EN_H0 / object$design$total_N,
      avg_pruning_threshold = mean(object$design$r[[1]]),
      avg_indication_n = mean(object$design$N[[1]])
    ),
    class = "summary.basket_trial_design"
  )
}


#' Print summary for basket trial design
#'
#' @param x Summary object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.summary.basket_trial_design <- function(x, ...) {
  cat("Summary of Two-Stage Basket Trial Design\n")
  cat("=========================================\n\n")
  
  K <- x$param$K
  p0 <- x$param$p0[[1]]
  pa <- x$param$pa[[1]]
  
  cat("Design Summary:\n")
  cat(sprintf("  Number of indications: %d\n", K))
  cat(sprintf("  Setting: %s\n", ifelse(x$param$homogeneous, "Homogeneous", "Heterogeneous")))
  cat(sprintf("  Target Type I error: %.3f (actual: %.4f)\n", 
              x$param$alpha, x$performance$type1_error))
  cat(sprintf("  Target power: %.1f%% (actual: %.1f%%)\n", 
              (1 - x$param$beta) * 100, x$performance$expected_power * 100))
  cat("\n")
  
  cat("Sample Size Summary:\n")
  cat(sprintf("  Stage I total: %d patients\n", x$design$S))
  cat(sprintf("  Stage II total (maximum): %d patients\n", x$design$total_N))
  cat(sprintf("  Expected under H0: %.1f patients\n", x$design$EN_H0))
  cat(sprintf("  Efficiency (E[N|H0] / total N): %.1f%%\n", x$efficiency * 100))
  cat(sprintf("  Average per indication: %.1f patients\n", x$avg_indication_n))
  cat("\n")
  
  cat("Decision Thresholds:\n")
  cat(sprintf("  Stage I: Continue if >=%d total responders\n", x$design$R1))
  cat(sprintf("  Stage II: Average pruning threshold r = %.1f responders\n", 
              x$avg_pruning_threshold))
  cat("\n")
  
  cat("Response Rate Ranges:\n")
  cat(sprintf("  Null (p0): %.3f to %.3f\n", min(p0), max(p0)))
  cat(sprintf("  Alternative (pa): %.3f to %.3f\n", min(pa), max(pa)))
  cat(sprintf("  Improvement: %.1f%% to %.1f%% absolute\n", 
              min((pa - p0) * 100), max((pa - p0) * 100)))
  
  invisible(x)
}


#' Print method for simple_basket_design results
#'
#' @param x Output from \code{simon_basket_trial}
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.simple_basket_design <- function(x, ...) {
  cat("Simple Basket Trial Design\n")
  cat("=========================\n\n")
  
  cat("Design Type: Adaptive Two-Stage with Heterogeneity Assessment\n")
  cat("Reference: Cunanan et al. (2017) Statistics in Medicine\n\n")
  
  cat("Input Parameters:\n")
  cat(sprintf("  Number of baskets (K): %d\n", x$inputs$K))
  cat(sprintf("  Null response rate (theta0): %.2f\n", x$inputs$p0))
  cat(sprintf("  Alternative response rate (thetaa): %.2f\n", x$inputs$pa))
  cat(sprintf("  Target FWER: %.3f\n", x$inputs$target_fwer))
  cat(sprintf("  Target power: %.2f\n", x$inputs$target_power))
  cat(sprintf("  Minimum power (A=1): %.2f\n", x$inputs$min_power))
  cat(sprintf("  Calibrated at: A=%d active baskets\n", x$inputs$calibrate_at_A))
  cat("\n")
  
  cat("Optimized Design Parameters:\n")
  cat(sprintf("  Stage 1 total sample size (N1): %d patients\n", x$design$N1))
  cat(sprintf("  Stage 2 sample size per basket (n2k, heterogeneous): %d patients\n", 
              x$design$n2k))
  cat(sprintf("  Stage 2 total sample size (N2, homogeneous): %d patients\n", 
              x$design$N2))
  cat(sprintf("  Heterogeneity tuning parameter (gamma): %.2f\n", x$design$gamma))
  cat(sprintf("  Min responders per basket to continue (r_s): %d\n", x$design$r_s))
  cat(sprintf("  Min total responders to continue (r_c): %d\n", x$design$r_c))
  cat(sprintf("  Significance level - separate (alphaS): %.3f\n", x$design$alpha_s))
  cat(sprintf("  Significance level - combined (alphaC): %.3f\n", x$design$alpha_c))
  cat("\n")
  
  cat("Operating Characteristics at Key Scenarios:\n")
  
  # A=0 (null scenario)
  cat(sprintf("  A=0 (all null): FWER=%.3f, EN=%.0f, ET=%.1f months\n",
              x$performance$A0$fwer,
              x$performance$A0$EN,
              x$performance$A0$ET))
  
  # A=1
  cat(sprintf("  A=1: Power=%.2f, EN=%.0f, ET=%.1f months\n",
              x$performance$A1$marginal_power[1],
              x$performance$A1$EN,
              x$performance$A1$ET))
  
  # Calibration scenario
  calib_A <- x$inputs$calibrate_at_A
  calib_key <- paste0("A", calib_A)
  cat(sprintf("  A=%d (calibration): Power=%.2f, EN=%.0f, ET=%.1f months\n",
              calib_A,
              x$performance[[calib_key]]$marginal_power[1],
              x$performance[[calib_key]]$EN,
              x$performance[[calib_key]]$ET))
  
  # A=K (all active)
  K <- x$inputs$K
  all_active_key <- paste0("A", K)
  cat(sprintf("  A=%d (all active): Power=%.2f, EN=%.0f, ET=%.1f months\n",
              K,
              mean(x$performance[[all_active_key]]$marginal_power),
              x$performance[[all_active_key]]$EN,
              x$performance[[all_active_key]]$ET))
  
  cat("\n")
  cat("Use summary() for detailed operating characteristics across all scenarios.\n")
  
  invisible(x)
}


#' Summary method for simple_basket_design results
#'
#' @param object Output from \code{simon_basket_trial}
#' @param ... Additional arguments (not used)
#'
#' @return A summary object
#' @export
summary.simple_basket_design <- function(object, ...) {
  structure(
    list(
      design = object$design,
      performance = object$performance,
      inputs = object$inputs,
      reference = object$reference
    ),
    class = "summary.simple_basket_design"
  )
}


#' Print summary for Simon basket design
#'
#' @param x Summary object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.summary.simple_basket_design <- function(x, ...) {
  cat("Summary of Simple Basket Trial Design\n")
  cat("=====================================\n\n")
  
  cat("Design Parameters:\n")
  cat(sprintf("  K = %d baskets\n", x$inputs$K))
  cat(sprintf("  theta0 = %.2f, thetaa = %.2f\n", x$inputs$p0, x$inputs$pa))
  cat(sprintf("  N1 = %d, N2 = %d, n2k = %d\n", 
              x$design$N1, x$design$N2, x$design$n2k))
  cat(sprintf("  r_s = %d, r_c = %d\n", x$design$r_s, x$design$r_c))
  cat(sprintf("  gamma = %.2f, alphaS = %.3f, alphaC = %.3f\n",
              x$design$gamma, x$design$alpha_s, x$design$alpha_c))
  cat("\n")
  
  cat("Operating Characteristics by Scenario\n")
  cat("======================================\n\n")
  
  # Create summary table
  K <- x$inputs$K
  scenarios <- 0:K
  
  cat("Scenario | FWER/Power (Basket 1) | EN   | ET   | Sens | Spec | Track1%\n")
  cat("---------|------------------------|------|------|------|------|--------\n")
  
  for (A in scenarios) {
    key <- paste0("A", A)
    oc <- x$performance[[key]]
    
    # Format output depending on scenario
    if (A == 0) {
      power_str <- sprintf("%.3f*", oc$fwer)
    } else {
      power_str <- sprintf("%.2f ", oc$marginal_power[1])
    }
    
    sens_str <- if (is.na(oc$sensitivity)) "  -  " else sprintf("%.2f", oc$sensitivity)
    spec_str <- if (is.na(oc$specificity)) "  -  " else sprintf("%.2f", oc$specificity)
    
    cat(sprintf("  A=%-2d   |  %6s              | %4.0f | %4.1f | %4s | %4s | %5.1f%%\n",
                A, power_str, oc$EN, oc$ET, sens_str, spec_str, 
                oc$prop_track1 * 100))
  }
  
  cat("\n* FWER (Family-Wise Error Rate) for A=0; Power for A>0\n")
  cat("Sens = Sensitivity (true positive rate)\n")
  cat("Spec = Specificity (true negative rate)\n")
  cat("Track1% = Proportion using heterogeneous path\n")
  cat("\n")
  
  # Power by basket for heterogeneous scenarios
  cat("Marginal Power by Basket (Active Baskets Only)\n")
  cat("===============================================\n\n")
  
  for (A in 1:K) {
    key <- paste0("A", A)
    oc <- x$performance[[key]]
    
    active_powers <- oc$marginal_power[1:A]
    cat(sprintf("A=%d: ", A))
    cat(sprintf("%.2f", active_powers), sep = ", ")
    cat(sprintf(" (avg: %.2f)\n", mean(active_powers)))
  }
  
  cat("\n")
  
  # Reference design comparison if available
  if (!is.null(x$reference) && !is.null(x$reference$A0)) {
    cat("Comparison to Reference Design (Independent Simon Two-Stage)\n")
    cat("============================================================\n\n")
    
    cat("Scenario | Proposed EN | Reference EN | Reduction | Proposed ET | Reference ET\n")
    cat("---------|-------------|--------------|-----------|-------------|-------------\n")
    
    for (A in scenarios) {
      key <- paste0("A", A)
      prop_en <- x$performance[[key]]$EN
      ref_en <- x$reference[[key]]$EN
      reduction <- (ref_en - prop_en) / ref_en * 100
      prop_et <- x$performance[[key]]$ET
      ref_et <- x$reference[[key]]$ET
      
      cat(sprintf("  A=%-2d   |    %4.0f     |     %4.0f     |  %5.1f%%  |    %4.1f    |    %4.1f\n",
                  A, prop_en, ref_en, reduction, prop_et, ref_et))
    }
    
    cat("\n")
  }
  
  # Interpretation
  cat("Interpretation\n")
  cat("==============\n")
  cat("- Design controls FWER at", sprintf("%.1f%%", x$inputs$target_fwer * 100), 
      "when all baskets are null (A=0)\n")
  cat("- Achieves", sprintf("%.1f%%", x$inputs$target_power * 100), 
      "power when A=", x$inputs$calibrate_at_A, "baskets are active\n")
  cat("- Maintains >=", sprintf("%.1f%%", x$inputs$min_power * 100), 
      "power when only A=1 basket is active\n")
  cat("- Most efficient when drug works in most/all baskets (lower EN as A increases)\n")
  cat("- Trade-off: Specificity decreases as A increases (risk of false positives in\n")
  cat("  inactive baskets due to pooling in homogeneous path)\n")
  
  invisible(x)
}
