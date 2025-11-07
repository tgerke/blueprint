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
  cat(sprintf("    Total sample size: %d patients\n", x$design$S))
  cat(sprintf("    Significance level: %.3f\n", x$design$alpha1))
  cat(sprintf("    Critical value R1: %d (continue if >=%d responders)\n", 
              x$design$R1, x$design$R1))
  cat("\n")
  
  cat("  Stage II (Pruning and Pooling):\n")
  cat("    Sample sizes per indication:\n")
  for (i in 1:K) {
    cat(sprintf("      Indication %d: N=%d, r=%d (prune if <%d responders)\n", 
                i, N[i], r[i], r[i]))
  }
  cat(sprintf("    Total sample size: %d patients\n", x$design$total_N))
  cat(sprintf("    Significance level: %.3f\n", x$design$alpha2))
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
