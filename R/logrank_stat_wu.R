#' Calculate One-Sample Log-Rank Test Statistic (Wu et al. 2020)
#'
#' @description
#' Calculates the one-sample log-rank test statistic as described in Wu et al. 
#' (2020). This statistic compares observed events against expected events under 
#' a null hypothesis survival distribution.
#' 
#' The test statistic is defined as:
#' 
#' \deqn{Z = \frac{E - O}{\sqrt{E}}}
#' 
#' where \eqn{E} is the expected number of events under H₀ (sum of cumulative 
#' hazards at observation times) and \eqn{O} is the observed number of events. 
#' The variance approximation uses \eqn{Var(O) \approx E} (Poisson approximation).
#'
#' @param obs_times Numeric vector of observation times for each patient (either 
#'   event time or censoring time).
#' @param event_ind Logical vector indicating whether each observation is an 
#'   event (TRUE) or censored (FALSE).
#' @param dist Character string specifying the null hypothesis distribution. 
#'   Options are:
#'   \itemize{
#'     \item \code{"WB"}: Weibull distribution
#'     \item \code{"LN"}: Log-normal distribution
#'     \item \code{"LG"}: Log-logistic distribution
#'     \item \code{"GM"}: Gamma distribution
#'   }
#' @param shape_h0 Numeric. Shape parameter of the null hypothesis distribution.
#' @param scale_h0 Numeric. Scale parameter of the null hypothesis distribution.
#'
#' @return Numeric. The one-sample log-rank test statistic. Positive values 
#'   indicate fewer observed events than expected under H₀ (favoring H₁, better 
#'   survival). Negative values indicate more observed events than expected 
#'   (favoring H₀).
#'
#' @details
#' The function implements the one-sample log-rank test as described in Wu et al. 
#' (2020), matching the implementation in the \code{OneArm2stage} package.
#' 
#' The expected number of events is calculated as the sum of cumulative hazards 
#' under H₀ at each observation time:
#' 
#' \deqn{E = \sum_{i=1}^{n} H_0(t_i)}
#' 
#' where \eqn{H_0(t) = -\log(S_0(t))} is the cumulative hazard function under 
#' the null hypothesis.
#' 
#' **Sign Convention**: The statistic uses the convention \eqn{Z = (E - O) / \sqrt{E}}, 
#' so positive values favor the alternative hypothesis (fewer events = better survival).
#'
#' @references
#' Wu, J., Chen, L., & Wei, L. J. (2020). Two‐stage phase II survival trial design. 
#' \emph{Pharmaceutical Statistics}, 19(3), 214-229. \doi{10.1002/pst.1983}
#'
#' @keywords internal
#' @export
logrank_stat_wu <- function(obs_times, event_ind, dist, shape_h0, scale_h0) {
  
  # Create cumulative hazard function under H0
  # H(t) = -log(S(t))
  H0 <- switch(dist,
    "WB" = function(t) {
      # H(t) = (t/scale)^shape for Weibull
      (t / scale_h0)^shape_h0
    },
    "LN" = function(t) {
      # H(t) = -log(1 - plnorm(...))
      S <- 1 - stats::plnorm(t, meanlog = scale_h0, sdlog = shape_h0)
      -log(pmax(S, 1e-10))  # Avoid log(0)
    },
    "LG" = function(t) {
      # H(t) = -log(S(t)) where S(t) = 1/(1+(t/scale)^shape)
      -log(1 / (1 + (t/scale_h0)^shape_h0))
    },
    "GM" = function(t) {
      S <- 1 - stats::pgamma(t, shape = shape_h0, scale = scale_h0)
      -log(pmax(S, 1e-10))
    },
    stop("Unsupported distribution: ", dist)
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
    # Degenerate case - no expected events
    Z <- 0
  }
  
  return(Z)
}
