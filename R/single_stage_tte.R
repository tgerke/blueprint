#' Single-Stage Single-Arm Design for Time-to-Event Endpoints
#'
#' @description
#' Calculate single-stage single-arm phase II trial design with a 
#' time-to-event endpoint using the one-sample log-rank test. This implementation
#' is based on the methodology described in Wu (2015).
#'
#' @param shape Numeric. Shape parameter for the parametric survival distribution 
#'   under the null hypothesis.
#' @param S0 Numeric. Survival probability at fixed time point \code{x0} under 
#'   the null hypothesis.
#' @param x0 Numeric. Fixed time point at which \code{S0} is specified.
#' @param hr Numeric. Hazard ratio under the alternative hypothesis (proportional 
#'   hazards model). Values < 1 indicate improvement.
#' @param tf Numeric. Follow-up time. Time from enrollment of last patient to 
#'   end of study. Earlier patients have longer follow-up (unrestricted follow-up).
#' @param rate Numeric. Constant accrual rate (patients per time unit).
#' @param alpha Numeric. Type I error rate. Default is 0.05.
#' @param beta Numeric. Type II error rate. Default is 0.2 (80% power).
#' @param two_sided Logical. Whether to use a two-sided test (default \code{TRUE}).
#'   Wu (2015) describes two-sided tests, while Wu et al. (2020) uses one-sided.
#' @param dist Character. Distribution under the null hypothesis. Options are:
#'   \itemize{
#'     \item \code{"WB"}: Weibull distribution
#'     \item \code{"LN"}: Log-normal distribution
#'     \item \code{"GM"}: Gamma distribution
#'     \item \code{"LG"}: Log-logistic distribution
#'     \item \code{"SP"}: Logspline distribution (requires \code{data} argument)
#'   }
#' @param restricted Logical. Follow-up model (default \code{FALSE}):
#'   \itemize{
#'     \item \code{FALSE}: Unrestricted follow-up (default). 
#'       Patients enrolled earlier have longer follow-up. Total study time = 
#'       accrual time + \code{tf}.
#'     \item \code{TRUE}: Restricted follow-up. Each patient followed for exactly 
#'       \code{tf} time units regardless of enrollment time.
#'   }
#' @param data Data frame. Required only when \code{dist = "SP"}. Must contain 
#'   columns \code{time} and \code{status} for fitting logspline distribution.
#'
#' @details
#' The function implements the single-stage design algorithm described in 
#' Wu (2015). Unlike the two-stage design in Wu et al. (2020), this supports
#' both one-sided and two-sided tests.
#' 
#' \strong{Follow-up Models:}
#' 
#' The \code{restricted} parameter determines how patient follow-up is handled:
#' 
#' \itemize{
#'   \item \strong{Unrestricted} (\code{restricted = FALSE}, default): The last 
#'     patient enrolled is followed for \code{tf} time units. Patients enrolled 
#'     earlier have progressively longer follow-up, up to \code{ta + tf} for the 
#'     first patient (where \code{ta} is the accrual time).
#'     
#'   \item \strong{Restricted} (\code{restricted = TRUE}): Every patient is 
#'     followed for exactly \code{tf} time units after their enrollment 
#'     (administrative censoring at \code{tf}).
#' }
#' 
#' The design parameters are:
#' \itemize{
#'   \item \code{n}: Total sample size
#'   \item \code{ta}: Accrual time
#'   \item \code{c}: Critical value (reject H0 if Z > c for one-sided, or |Z| > c for two-sided)
#' }
#' 
#' Four parametric distributions are supported:
#' \itemize{
#'   \item Weibull: S(t) = exp(-(t/b)^a)
#'   \item Log-normal: S(t) = 1 - Φ((log(t) - μ)/σ)
#'   \item Gamma: S(t) = 1 - I_a(t/b)
#'   \item Log-logistic: S(t) = 1/(1 + (t/b)^a)
#' }
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{param}{Data frame of input parameters}
#'     \item{single_stage}{Data frame with single-stage design parameters:
#'       \code{n}, \code{ta}, \code{c}, \code{MTSL} (maximum total study length),
#'       \code{E} (expected events), \code{P1} (probability of event)}
#'   }
#'
#' @references
#' Wu, J. (2015). Sample size calculation for the one-sample log-rank test. 
#' \emph{Pharmaceutical Statistics}, 14(1), 26-33. \doi{10.1002/pst.1654}
#'
#' @examples
#' \dontrun{
#' # Example: Exponential distribution (shape = 1)
#' # M0 = 0.5, M1 = 0.75, HR = 0.667
#' # Two-sided test (Wu 2015 default)
#' design_two_sided <- single_stage_single_arm_tte(
#'   shape = 1.0,
#'   S0 = 0.5,
#'   x0 = 0.5,
#'   hr = 0.667,
#'   tf = 1.0,
#'   rate = 32.4,  # 81 patients / 2.5 time units
#'   alpha = 0.05,
#'   beta = 0.10,
#'   two_sided = TRUE,
#'   dist = "WB",
#'   restricted = FALSE
#' )
#' 
#' # One-sided test
#' design_one_sided <- single_stage_single_arm_tte(
#'   shape = 1.47327,
#'   S0 = 0.5,
#'   x0 = 3.5,
#'   hr = 0.5913,
#'   tf = 5,
#'   rate = 2,
#'   alpha = 0.05,
#'   beta = 0.2,
#'   two_sided = FALSE,
#'   dist = "WB",
#'   restricted = FALSE
#' )
#' 
#' print(design_two_sided)
#' }
#'
#' @export
single_stage_single_arm_tte <- function(shape = NULL, S0 = NULL, x0 = NULL, hr, tf, 
                                        rate, alpha = 0.05, beta = 0.2, 
                                        two_sided = TRUE,
                                        dist, restricted = FALSE,
                                        data = NULL) {
  
  # Input validation
  if (!dist %in% c("WB", "LN", "GM", "LG", "SP")) {
    stop("dist must be one of: 'WB', 'LN', 'GM', 'LG', 'SP'")
  }
  
  if (dist != "SP" && (is.null(shape) || is.null(S0) || is.null(x0))) {
    stop("shape, S0, and x0 are required for parametric distributions")
  }
  
  if (dist == "SP" && is.null(data)) {
    stop("data is required when dist = 'SP'")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  
  if (beta <= 0 || beta >= 1) {
    stop("beta must be between 0 and 1")
  }
  
  if (hr <= 0) {
    stop("hr must be positive")
  }
  
  if (rate <= 0) {
    stop("rate must be positive")
  }
  
  if (tf <= 0) {
    stop("tf must be positive")
  }
  
  # Initialize survival functions
  surv_funs <- .define_survival_functions(dist, shape, S0, x0, hr, data)
  s <- surv_funs$s
  h <- surv_funs$h
  s0 <- surv_funs$s0
  h0 <- surv_funs$h0
  H0 <- surv_funs$H0
  scale1 <- hr
  
  # Determine critical value based on one-sided vs two-sided
  # For two-sided: alpha/2 in each tail
  # For one-sided: alpha in upper tail
  if (two_sided) {
    z_alpha <- stats::qnorm(1 - alpha/2)
  } else {
    z_alpha <- stats::qnorm(1 - alpha)
  }
  
  z_beta <- stats::qnorm(1 - beta)
  
  # Initial estimate without G(t) - just integrate over tf
  # This gives us a starting point for ta
  g0_init <- function(t) {
    val <- s(scale1, t) * h0(t)
    val[!is.finite(val)] <- 0
    val
  }
  g1_init <- function(t) {
    val <- s(scale1, t) * h(scale1, t)
    val[!is.finite(val)] <- 0
    val
  }
  g00_init <- function(t) {
    s_val <- s(scale1, t)
    val <- s_val * H0(t) * h0(t)
    val[s_val < 1e-10 | !is.finite(val)] <- 0
    val
  }
  g01_init <- function(t) {
    s_val <- s(scale1, t)
    val <- s_val * H0(t) * h(scale1, t)
    val[s_val < 1e-10 | !is.finite(val)] <- 0
    val
  }
  
  p0_init <- stats::integrate(g0_init, 0, tf, rel.tol = 1e-6)$value
  p1_init <- stats::integrate(g1_init, 0, tf, rel.tol = 1e-6)$value
  p00_init <- stats::integrate(g00_init, 0, tf, rel.tol = 1e-6)$value
  p01_init <- stats::integrate(g01_init, 0, tf, rel.tol = 1e-6)$value
  
  s1_init <- sqrt(p1_init - p1_init^2 + 2 * p00_init - p0_init^2 - 
                    2 * p01_init + 2 * p0_init * p1_init)
  s0_init <- sqrt(p0_init)
  om_init <- p0_init - p1_init
  
  n_init <- ceiling((s0_init * z_alpha + s1_init * z_beta)^2 / om_init^2)
  ta_init <- n_init / rate
  
  # Now refine with actual G(t) calculation
  # Iterate to converge on the right ta
  ta <- ta_init
  for (i in 1:10) {
    tau <- ta + tf
    
    if (restricted) {
      G <- function(t) ifelse(t <= tf, 1, 0)
    } else {
      G <- function(t) 1 - stats::punif(t, tf, tau)
    }
    
    g0 <- function(t) {
      val <- s(scale1, t) * h0(t) * G(t)
      val[!is.finite(val)] <- 0
      val
    }
    g1 <- function(t) {
      val <- s(scale1, t) * h(scale1, t) * G(t)
      val[!is.finite(val)] <- 0
      val
    }
    g00 <- function(t) {
      s_val <- s(scale1, t)
      # Set to 0 where survival is too small to avoid H0(t) explosion
      val <- s_val * H0(t) * h0(t) * G(t)
      val[s_val < 1e-10 | !is.finite(val)] <- 0
      val
    }
    g01 <- function(t) {
      s_val <- s(scale1, t)
      # Set to 0 where survival is too small to avoid H0(t) explosion
      val <- s_val * H0(t) * h(scale1, t) * G(t)
      val[s_val < 1e-10 | !is.finite(val)] <- 0
      val
    }
    
    p0 <- stats::integrate(g0, 0, tau, rel.tol = 1e-6)$value
    p1 <- stats::integrate(g1, 0, tau, rel.tol = 1e-6)$value
    p00 <- stats::integrate(g00, 0, tau, rel.tol = 1e-6)$value
    p01 <- stats::integrate(g01, 0, tau, rel.tol = 1e-6)$value
    
    s1 <- sqrt(p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1)
    s0 <- sqrt(p0)
    om <- p0 - p1
    
    n_required <- (s0 * z_alpha + s1 * z_beta)^2 / om^2
    ta_new <- n_required / rate
    
    # Check for convergence
    if (abs(ta_new - ta) < 0.01) break
    ta <- ta_new
  }
  
  n <- ceiling(ta * rate)
  
  # Calculate expected events and probability of event under H1
  # Using final converged values
  expected_events <- n * p1
  prob_event <- p1
  
  # Prepare output
  if (dist == "SP") {
    param <- data.frame(
      hr = hr, alpha = alpha, beta = beta, rate = rate, tf = tf,
      two_sided = two_sided, restricted = restricted, dist = dist
    )
  } else {
    param <- data.frame(
      shape = shape, S0 = S0, x0 = x0, hr = hr, alpha = alpha, beta = beta, 
      rate = rate, tf = tf, two_sided = two_sided, restricted = restricted, dist = dist
    )
  }
  
  single_stage <- data.frame(
    n = n,
    ta = ta,
    c = z_alpha,
    MTSL = ta + tf,
    E = round(expected_events, 1),
    P1 = round(prob_event, 4)
  )
  
  design <- list(
    param = param,
    single_stage = single_stage
  )
  
  class(design) <- "single_stage_design"
  
  return(design)
}
