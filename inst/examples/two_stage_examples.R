# Example 1: Small-Cell Lung Cancer Trial (from Wu et al. 2020)
# ===============================================================

library(blueprint)

# This example reproduces the design from Wu et al. (2020) for a 
# single-arm phase II trial evaluating immunotherapy as second-line 
# treatment for small-cell lung cancer patients.

# Design parameters based on historical data:
# - Null hypothesis: median PFS = 3.5 months (Weibull distribution)
# - Alternative hypothesis: median PFS = 5 months
# - Hazard ratio = (3.5/5)^1.47327 = 0.5913
# - Weibull shape parameter: 1.47327 (fitted from historical data)
# - Restricted follow-up: 5 months (practical constraint)
# - Accrual rate: 2 patients per month
# - Type I error: 5% (one-sided)
# - Power: 80%

design_5mo <- two_stage_single_arm_tte(
  shape = 1.47327,    # Weibull shape parameter
  S0 = 0.5,           # Survival probability at 3.5 months under H0
  x0 = 3.5,           # Time point (median under H0)
  hr = 0.5913,        # Hazard ratio under H1
  x = 5,              # Restricted follow-up time (months)
  rate = 2,           # Accrual rate (patients/month)
  alpha = 0.05,       # Type I error (one-sided)
  beta = 0.2,         # Type II error (80% power)
  dist = "WB"         # Weibull distribution
)

# View the design
print(design_5mo)

# Interpretation:
# - Enroll 9 patients in Stage 1
# - After all Stage 1 patients are enrolled (~4.2 months), conduct interim
# - Each patient followed for up to 5 months or until event
# - If Z1 <= 1.25, stop for futility
# - Otherwise, continue to enroll total of 26 patients
# - Final analysis after last patient followed for 5 months
# - Reject H0 if Z > 0.69


# Example 2: Extended Follow-up
# ==============================

# Same design but with 10 months of follow-up
design_10mo <- two_stage_single_arm_tte(
  shape = 1.47327,
  S0 = 0.5,
  x0 = 3.5,
  hr = 0.5913,
  x = 10,             # Extended follow-up
  rate = 2,
  alpha = 0.05,
  beta = 0.2,
  dist = "WB"
)

print(design_10mo)

# Note: Longer follow-up typically reduces required sample size
# because more information is collected per patient


# Example 3: Different Distributions
# ===================================

# Log-normal distribution example
design_ln <- two_stage_single_arm_tte(
  shape = 0.5,        # Log-normal shape (SD of log-survival)
  S0 = 0.3,           # Survival probability at x0
  x0 = 1,             # Time point
  hr = 0.65,          # Hazard ratio
  x = 1,              # Follow-up time
  rate = 10,          # Accrual rate
  alpha = 0.05,
  beta = 0.2,
  dist = "LN"         # Log-normal distribution
)

print(design_ln)


# Gamma distribution example
design_gm <- two_stage_single_arm_tte(
  shape = 1,          # Gamma shape parameter
  S0 = 0.3,
  x0 = 1,
  hr = 0.65,
  x = 1,
  rate = 10,
  alpha = 0.05,
  beta = 0.2,
  dist = "GM"         # Gamma distribution
)

print(design_gm)


# Log-logistic distribution example
design_lg <- two_stage_single_arm_tte(
  shape = 1,          # Log-logistic shape parameter
  S0 = 0.3,
  x0 = 1,
  hr = 0.65,
  x = 1,
  rate = 10,
  alpha = 0.05,
  beta = 0.2,
  dist = "LG"         # Log-logistic distribution
)

print(design_lg)


# Example 4: Sensitivity Analysis
# ================================

# Evaluate how sample size changes with different hazard ratios
hr_values <- seq(0.5, 0.8, by = 0.05)
results <- data.frame(
  hr = hr_values,
  n_single = NA,
  n_two_stage = NA,
  n1 = NA,
  ES = NA
)

for (i in seq_along(hr_values)) {
  design <- two_stage_single_arm_tte(
    shape = 1.47327,
    S0 = 0.5,
    x0 = 3.5,
    hr = hr_values[i],
    x = 5,
    rate = 2,
    alpha = 0.05,
    beta = 0.2,
    dist = "WB"
  )
  
  results$n_single[i] <- design$Single_stage$nsingle
  results$n_two_stage[i] <- design$Two_stage$n
  results$n1[i] <- design$Two_stage$n1
  results$ES[i] <- design$Two_stage$ES
}

print(results)

# Plot the relationship
plot(results$hr, results$n_single, type = "b", col = "blue",
     xlab = "Hazard Ratio", ylab = "Sample Size",
     main = "Sample Size vs Hazard Ratio",
     ylim = c(0, max(results$n_single)))
lines(results$hr, results$n_two_stage, type = "b", col = "red")
lines(results$hr, results$ES, type = "b", col = "green")
legend("topright", 
       legend = c("Single-stage", "Two-stage (max)", "Two-stage (ES)"),
       col = c("blue", "red", "green"), lty = 1, pch = 1)
