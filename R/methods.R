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
  cat(sprintf("  Sample size: %d\n", x$Single_stage$nsingle))
  cat(sprintf("  Accrual time: %.2f\n", x$Single_stage$tasingle))
  cat(sprintf("  Critical value: %.4f\n", x$Single_stage$csingle))
  cat("\n")
  
  cat("Two-Stage Design:\n")
  cat(sprintf("  Stage 1 sample size: %d\n", x$Two_stage$n1))
  cat(sprintf("  Stage 1 critical value: %.4f (stop if Z1 <= %.4f)\n", 
              x$Two_stage$c1, x$Two_stage$c1))
  cat(sprintf("  Total sample size: %d\n", x$Two_stage$n))
  cat(sprintf("  Final critical value: %.4f (reject H0 if Z > %.4f)\n", 
              x$Two_stage$c, x$Two_stage$c))
  cat(sprintf("  Interim analysis time: %.2f\n", x$Two_stage$t1))
  cat(sprintf("  Maximum total study length: %.2f\n", x$Two_stage$MTSL))
  cat(sprintf("  Expected sample size under H0: %.2f\n", x$Two_stage$ES))
  cat(sprintf("  Probability of early stopping under H0: %.4f\n", x$Two_stage$PS))
  
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
      single_stage = object$Single_stage,
      two_stage = object$Two_stage,
      efficiency = object$Two_stage$ES / object$Single_stage$nsingle
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
