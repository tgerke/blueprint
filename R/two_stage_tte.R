#' Two-Stage Single-Arm Design for Time-to-Event Endpoints
#'
#' @description
#' Calculate optimal two-stage single-arm phase II trial design with a 
#' time-to-event endpoint using the one-sample log-rank test. This implementation
#' is based on the methodology described in Wu et al. (2020).
#'
#' @param shape Numeric. Shape parameter for the parametric survival distribution 
#'   under the null hypothesis.
#' @param S0 Numeric. Survival probability at fixed time point \code{x0} under 
#'   the null hypothesis.
#' @param x0 Numeric. Fixed time point at which \code{S0} is specified.
#' @param hr Numeric. Hazard ratio under the alternative hypothesis (proportional 
#'   hazards model). Values < 1 indicate improvement.
#' @param tf Numeric. Follow-up time. Interpretation depends on \code{restricted}:
#'   \itemize{
#'     \item If \code{restricted = TRUE}: Each patient followed for exactly \code{tf} 
#'       time units after enrollment (administrative censoring).
#'     \item If \code{restricted = FALSE}: Time from enrollment of last patient to 
#'       end of study. Earlier patients have longer follow-up (unrestricted follow-up).
#'   }
#' @param rate Numeric. Constant accrual rate (patients per time unit).
#' @param alpha Numeric. Type I error rate (one-sided). Default is 0.05.
#' @param beta Numeric. Type II error rate. Default is 0.2 (80% power).
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
#'     \item \code{FALSE}: Unrestricted follow-up (Wu et al. 2020 default). 
#'       Patients enrolled earlier have longer follow-up. Total study time = 
#'       accrual time + \code{tf}.
#'     \item \code{TRUE}: Restricted follow-up. Each patient followed for exactly 
#'       \code{tf} time units regardless of enrollment time.
#'   }
#' @param data Data frame. Required only when \code{dist = "SP"}. Must contain 
#'   columns \code{time} and \code{status} for fitting logspline distribution.
#'
#' @details
#' The function implements the optimal two-stage design algorithm described in 
#' Wu, Chen, Wei, Weiss, and Chauhan (2020). The design allows for early stopping 
#' for futility at the first stage based on a one-sample log-rank test.
#' 
#' \strong{Follow-up Models:}
#' 
#' The \code{restricted} parameter determines how patient follow-up is handled:
#' 
#' \itemize{
#'   \item \strong{Unrestricted} (\code{restricted = FALSE}, default): The last 
#'     patient enrolled is followed for \code{tf} time units. Patients enrolled 
#'     earlier have progressively longer follow-up, up to \code{ta + tf} for the 
#'     first patient (where \code{ta} is the accrual time). This matches the 
#'     approach in Wu et al. (2020) and the \code{OneArm2stage} package.
#'     
#'   \item \strong{Restricted} (\code{restricted = TRUE}): Every patient is 
#'     followed for exactly \code{tf} time units after their enrollment 
#'     (administrative censoring at \code{tf}). This is simpler but may require 
#'     larger sample sizes.
#' }
#' 
#' The design parameters are:
#' \itemize{
#'   \item \code{n1}: Sample size for stage 1
#'   \item \code{c1}: Critical value for stage 1 (stop if Z1 ≤ c1)
#'   \item \code{n}: Total sample size
#'   \item \code{c}: Critical value for final analysis (reject H0 if Z > c)
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
#' @return A list with three elements:
#'   \describe{
#'     \item{param}{Data frame of input parameters}
#'     \item{single_stage}{Data frame with single-stage design parameters:
#'       \code{nsingle}, \code{tasingle}, \code{csingle}}
#'     \item{two_stage}{Data frame with two-stage design parameters:
#'       \code{n1}, \code{c1}, \code{n}, \code{c}, \code{t1}, \code{MTSL}, 
#'       \code{ES}, \code{PS}}
#'   }
#'
#' @references
#' Wu, J., Chen, L., Wei, J., Weiss, H., & Chauhan, A. (2020). Optimal 
#' two-stage phase II survival trial design. \emph{Pharmaceutical Statistics}, 
#' 19(3), 214-229. \doi{10.1002/pst.1983}
#'
#' @examples
#' \dontrun{
#' # Example from Wu et al. (2020) - Small-cell lung cancer trial
#' # Unrestricted follow-up (default)
#' design_unrestricted <- two_stage_single_arm_tte(
#'   shape = 1.47327,
#'   S0 = 0.5,
#'   x0 = 3.5,
#'   hr = 0.5913,
#'   tf = 5,
#'   rate = 2,
#'   alpha = 0.05,
#'   beta = 0.2,
#'   dist = "WB",
#'   restricted = FALSE  # Unrestricted (Wu et al. default)
#' )
#' 
#' # Restricted follow-up (each patient followed exactly 5 time units)
#' design_restricted <- two_stage_single_arm_tte(
#'   shape = 1.47327,
#'   S0 = 0.5,
#'   x0 = 3.5,
#'   hr = 0.5913,
#'   tf = 5,
#'   rate = 2,
#'   alpha = 0.05,
#'   beta = 0.2,
#'   dist = "WB",
#'   restricted = TRUE
#' )
#' 
#' print(design_unrestricted)
#' print(design_restricted)
#' }
#'
#' @export
two_stage_single_arm_tte <- function(shape = NULL, S0 = NULL, x0 = NULL, hr, tf, 
                                   rate, alpha = 0.05, beta = 0.2, dist,
                                   restricted = FALSE,
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
  
  # Core function for computing design parameters
  fct <- function(zi, ceps = 0.0001, alphaeps = 0.0001, nbmaxiter = 100, dist) {
    ta <- as.numeric(zi[1])
    t1 <- as.numeric(zi[2])
    c1 <- as.numeric(zi[3])
    
    # Define survival functions based on distribution
    surv_funs <- .define_survival_functions(dist, shape, S0, x0, hr, data)
    s0 <- surv_funs$s0
    h0 <- surv_funs$h0
    H0 <- surv_funs$H0
    s <- surv_funs$s
    h <- surv_funs$h
    scale1 <- hr
    
    # Define G(t) and tau based on restricted flag
    if (restricted) {
      # Restricted: each patient followed for exactly tf
      G <- function(t) ifelse(t <= tf, 1, 0)
      tau <- tf
    } else {
      # Unrestricted: enrollment distribution matters
      # Last patient gets tf, first patient gets ta + tf
      G <- function(t) 1 - stats::punif(t, tf, ta + tf)
      tau <- ta + tf
    }
    
    # Calculate variance components
    g0 <- function(t) s(scale1, t) * h0(t) * G(t)
    g1 <- function(t) s(scale1, t) * h(scale1, t) * G(t)
    g00 <- function(t) s(scale1, t) * H0(t) * h0(t) * G(t)
    g01 <- function(t) s(scale1, t) * H0(t) * h(scale1, t) * G(t)
    
    p0 <- stats::integrate(g0, 0, tau)$value
    p1 <- stats::integrate(g1, 0, tau)$value
    p00 <- stats::integrate(g00, 0, tau)$value
    p01 <- stats::integrate(g01, 0, tau)$value
    
    sigma2.1 <- p1 - p1^2 + 2 * p00 - p0^2 - 2 * p01 + 2 * p0 * p1
    sigma2.0 <- p0
    om <- p0 - p1
    
    # Stage 1 variance components  
    # G1(t) accounts for enrollment pattern at interim analysis
    # Following Wu et al. (2020) and OneArm2stage implementation:
    # For unrestricted follow-up, G1(t) = 1 - punif(t, t1-ta, t1)
    # This represents the censoring distribution at interim
    if (restricted) {
      # For restricted follow-up, everyone has exactly tf follow-up
      G1 <- function(t) ifelse(t <= tf, 1, 0)
    } else {
      # For unrestricted, use OneArm2stage formula
      G1 <- function(t) 1 - stats::punif(t, t1 - ta, t1)
    }
    g0_1 <- function(t) s(scale1, t) * h0(t) * G1(t)
    g1_1 <- function(t) s(scale1, t) * h(scale1, t) * G1(t)
    g00_1 <- function(t) s(scale1, t) * H0(t) * h0(t) * G1(t)
    g01_1 <- function(t) s(scale1, t) * H0(t) * h(scale1, t) * G1(t)
    
    # For stage 1, integrate to ta + tf (same as stage 2)
    # This matches OneArm2stage implementation
    p0_1 <- stats::integrate(g0_1, 0, tau)$value
    p1_1 <- stats::integrate(g1_1, 0, tau)$value
    p00_1 <- stats::integrate(g00_1, 0, tau)$value
    p01_1 <- stats::integrate(g01_1, 0, tau)$value
    
    sigma2.11 <- p1_1 - p1_1^2 + 2 * p00_1 - p0_1^2 - 2 * p01_1 + 2 * p0_1 * p1_1
    sigma2.01 <- p0_1
    om1 <- p0_1 - p1_1  # CRITICAL: interim effect size (must use om1 not om in cb1 calculation)
    
    q1 <- function(t) s0(t) * h0(t) * G1(t)
    q <- function(t) s0(t) * h0(t) * G(t)
    v1 <- stats::integrate(q1, 0, tau)$value
    v <- stats::integrate(q, 0, tau)$value
    
    rho0 <- sqrt(v1 / v)
    rho1 <- sqrt(sigma2.11 / sigma2.1)
    
    # Iterate to find c
    cL <- -10
    cU <- 10
    alphac <- -100
    iter <- 0
    
    while ((abs(alphac - alpha) > alphaeps || cU - cL > ceps) && iter < nbmaxiter) {
      iter <- iter + 1
      c <- (cL + cU) / 2
      alphac <- .calculate_alpha(c, c1, rho0)
      
      if (alphac > alpha) {
        cL <- c
      } else {
        cU <- c
      }
    }
    
    # Calculate power
    # CRITICAL: Use om1 (interim effect size) for cb1, not om (final effect size)
    cb1 <- sqrt(sigma2.01 / sigma2.11) * (c1 - (om1 * sqrt(rate * t1) / sqrt(sigma2.01)))
    cb <- sqrt(sigma2.0 / sigma2.1) * (c - (om * sqrt(rate * ta)) / sqrt(sigma2.0))
    pwrc <- .calculate_power(cb = cb, cb1 = cb1, rho1 = rho1)
    
    res <- c(cL, cU, alphac, 1 - pwrc, rho0, rho1, cb1, cb)
    return(res)
  }
  
  # Initialize survival functions for single-stage calculation
  surv_funs <- .define_survival_functions(dist, shape, S0, x0, hr, data)
  s <- surv_funs$s
  h <- surv_funs$h
  s0 <- surv_funs$s0
  h0 <- surv_funs$h0
  H0 <- surv_funs$H0
  scale1 <- hr
  
  # For single-stage, we need to find ta such that:
  # n = rate * ta = (s0 * z_alpha + s1 * z_beta)^2 / om^2
  # where s0, s1, om depend on ta through tau = ta + tf and G(t)
  #
  # Following OneArm2stage, we solve this iteratively using uniroot
  # But to avoid convergence issues, we'll use a reasonable approximation:
  # Start with a guess based on tf-only integration (no G function)
  # This gives us a starting point for ta
  
  # Initial estimate without G(t) - just integrate over tf
  g0_init <- function(t) s(scale1, t) * h0(t)
  g1_init <- function(t) s(scale1, t) * h(scale1, t)
  g00_init <- function(t) s(scale1, t) * H0(t) * h0(t)
  g01_init <- function(t) s(scale1, t) * H0(t) * h(scale1, t)
  
  p0_init <- stats::integrate(g0_init, 0, tf, rel.tol = 1e-6)$value
  p1_init <- stats::integrate(g1_init, 0, tf, rel.tol = 1e-6)$value
  p00_init <- stats::integrate(g00_init, 0, tf, rel.tol = 1e-6)$value
  p01_init <- stats::integrate(g01_init, 0, tf, rel.tol = 1e-6)$value
  
  s1_init <- sqrt(p1_init - p1_init^2 + 2 * p00_init - p0_init^2 - 
                    2 * p01_init + 2 * p0_init * p1_init)
  s0_init <- sqrt(p0_init)
  om_init <- p0_init - p1_init
  
  n_init <- ceiling((s0_init * stats::qnorm(1 - alpha) + s1_init * stats::qnorm(1 - beta))^2 / om_init^2)
  ta_init <- n_init / rate
  
  # Now refine with actual G(t) calculation
  # We'll iterate a few times to converge on the right ta
  ta <- ta_init
  for (i in 1:10) {
    tau_single <- ta + tf
    
    if (restricted) {
      G_single <- function(t) ifelse(t <= tf, 1, 0)
    } else {
      G_single <- function(t) 1 - stats::punif(t, tf, tau_single)
    }
    
    g0_single <- function(t) s(scale1, t) * h0(t) * G_single(t)
    g1_single <- function(t) s(scale1, t) * h(scale1, t) * G_single(t)
    g00_single <- function(t) s(scale1, t) * H0(t) * h0(t) * G_single(t)
    g01_single <- function(t) s(scale1, t) * H0(t) * h(scale1, t) * G_single(t)
    
    p0_single <- stats::integrate(g0_single, 0, tau_single, rel.tol = 1e-6)$value
    p1_single <- stats::integrate(g1_single, 0, tau_single, rel.tol = 1e-6)$value
    p00_single <- stats::integrate(g00_single, 0, tau_single, rel.tol = 1e-6)$value
    p01_single <- stats::integrate(g01_single, 0, tau_single, rel.tol = 1e-6)$value
    
    s1_single <- sqrt(p1_single - p1_single^2 + 2 * p00_single - p0_single^2 - 
                        2 * p01_single + 2 * p0_single * p1_single)
    s0_single <- sqrt(p0_single)
    om_single <- p0_single - p1_single
    
    n_required <- (s0_single * stats::qnorm(1 - alpha) + s1_single * stats::qnorm(1 - beta))^2 / om_single^2
    ta_new <- n_required / rate
    
    # Check for convergence
    if (abs(ta_new - ta) < 0.01) break
    ta <- ta_new
  }
  
  nsingle <- ceiling(ta * rate)
  tasingle <- ceiling(ta)
  
  single_stage <- data.frame(
    nsingle = nsingle,
    tasingle = tasingle,
    csingle = stats::qnorm(1 - alpha)
  )
  
  # Optimize two-stage design
  # Note: Grid search initialization affects convergence
  # Starting point based on single-stage design
  atc0 <- data.frame(n = nsingle, t1 = tasingle, c1 = 0.25)
  nbpt <- 11
  pascote <- 1.26
  cote <- 1 * pascote
  # CRITICAL: First iteration uses wide c1 range to explore broadly
  # Subsequent iterations narrow around optimal c1 (updated below)
  c1.lim <- nbpt * c(-1, 1)  # Initial: c(-11, 11) - matches OneArm2stage
  EnH0 <- 10000  # Permissive initial bound on expected sample size under H0
  iter <- 0
  nbmaxiter <- 100
  ceps <- 0.001
  alphaeps <- 0.001
  
  atcs <- NULL
  
  while (iter < nbmaxiter && diff(c1.lim) / nbpt > 0.001) {
    iter <- iter + 1
    
    if (iter %% 2 == 0) nbpt <- nbpt + 1
    cote <- cote / pascote
    
    n.lim <- atc0$n + c(-1, 1) * nsingle * cote
    t1.lim <- atc0$t1 + c(-1, 1) * tasingle * cote
    # After first iteration, narrow c1 range around optimal value
    c1.lim <- atc0$c1 + c(-1, 1) * cote
    ta.lim <- n.lim / rate
    t1.lim <- pmax(0, t1.lim)
    
    n <- seq(n.lim[1], n.lim[2], length.out = nbpt)
    n <- ceiling(n)
    n <- unique(n)
    ta <- n / rate
    t1 <- seq(t1.lim[1], t1.lim[2], length.out = nbpt)
    c1 <- seq(c1.lim[1], c1.lim[2], length.out = nbpt)
    ta <- ta[ta > 0]
    t1 <- t1[t1 >= 0.2 * tasingle & t1 <= 1.2 * tasingle]
    
    z <- expand.grid(list(ta = ta, t1 = t1, c1 = c1))
    z <- z[z$ta > z$t1, ]
    z$pap <- stats::pnorm(z$c1)
    z$eta <- z$ta - pmax(0, z$ta - z$t1) * z$pap
    z$enh0 <- z$eta * rate
    z <- z[z$enh0 <= EnH0, ]
    
    if (nrow(z) == 0) break
    
    resz <- t(apply(z, 1, fct, ceps = ceps, alphaeps = alphaeps, 
                    nbmaxiter = nbmaxiter, dist = dist))
    resz <- as.data.frame(resz)
    names(resz) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
    
    r <- cbind(z, resz)
    r$pap <- stats::pnorm(r$c1)
    r$eta <- r$ta - pmax(0, r$ta - r$t1) * r$pap
    r$etar <- r$eta * rate
    r$tar <- r$ta * rate
    r$c <- r$cL + (r$cU - r$cL) / 2
    r$diffc <- r$cU - r$cL
    r$diffc <- ifelse(r$diffc <= ceps, 1, 0)
    
    r <- r[1 - r$betac >= 1 - beta, ]
    r <- r[order(r$enh0), ]
    r$n <- r$ta * rate
    
    if (nrow(r) > 0) {
      atc <- r[, c("ta", "t1", "c1", "n")][1, ]
      r1 <- r[1, ]
      if (r1$enh0 < EnH0) {
        EnH0 <- r1$enh0
        atc0 <- atc
      }
    } else {
      atc <- data.frame(ta = NA, t1 = NA, c1 = NA, n = NA)
    }
    
    atc$iter <- iter
    atc$enh0 <- ifelse(nrow(r) > 0, r$enh0[1], NA)
    atc$EnH0 <- EnH0
    atc$tai <- ta.lim[1]
    atc$tas <- ta.lim[2]
    atc$ti <- t1.lim[1]
    atc$ts <- t1.lim[2]
    atc$ci <- c1.lim[1]
    atc$cs <- c1.lim[2]
    atc$cote <- cote
    
    if (iter == 1) {
      atcs <- atc
    } else {
      atcs <- rbind(atcs, atc)
    }
  }
  
  # Final design
  a <- atcs[!is.na(atcs$n), ]
  
  if (nrow(a) == 0) {
    stop("Optimization failed to converge")
  }
  
  p <- t(apply(a, 1, fct, ceps = ceps, alphaeps = alphaeps, 
               nbmaxiter = nbmaxiter, dist = dist))
  p <- as.data.frame(p)
  names(p) <- c("cL", "cU", "alphac", "betac", "rho0", "rho1", "cb1", "cb")
  p$i <- seq_len(nrow(p))
  a$i <- seq_len(nrow(a))
  res <- merge(a, p, by = "i", all = TRUE)
  
  res$pap <- round(stats::pnorm(res$c1), 4)
  res$ta <- res$n / rate
  res$Enh0 <- (res$ta - pmax(0, res$ta - res$t1) * res$pap) * rate
  res$diffc <- res$cU - res$cL
  res$c <- round(res$cL + (res$cU - res$cL) / 2, 4)
  res$diffc <- ifelse(res$diffc <= ceps, 1, 0)
  res <- res[order(res$enh0), ]
  des <- res[1, ]
  
  # Prepare output
  if (dist == "SP") {
    param <- data.frame(
      hr = hr, alpha = alpha, beta = beta, rate = rate, tf = tf,
      restricted = restricted, dist = dist
    )
  } else {
    param <- data.frame(
      shape = shape, S0 = S0, hr = hr, alpha = alpha, beta = beta, 
      rate = rate, x0 = x0, tf = tf, restricted = restricted, dist = dist
    )
  }
  
  two_stage <- data.frame(
    n1 = ceiling(des$t1 * rate),
    c1 = des$c1,
    n = ceiling(des$ta * rate),
    c = des$c,
    t1 = des$t1,
    MTSL = des$ta + tf,
    ES = des$Enh0,
    PS = des$pap
  )
  
  design <- list(
    param = param,
    single_stage = single_stage,
    two_stage = two_stage
  )
  
  class(design) <- "two_stage_design"
  
  return(design)
}


#' Define Survival Functions for Different Distributions
#'
#' @description
#' Internal helper function to define survival, hazard, and cumulative hazard 
#' functions for different parametric distributions.
#'
#' @keywords internal
#' @noRd
.define_survival_functions <- function(dist, shape, S0, x0, hr, data) {
  
  if (dist == "WB") {
    # Weibull distribution
    scale0 <- x0 / (-log(S0))^(1/shape)
    
    list(
      s0 = function(u) 1 - stats::pweibull(u, shape, scale0),
      f0 = function(u) stats::dweibull(u, shape, scale0),
      h0 = function(u) {
        f <- stats::dweibull(u, shape, scale0)
        s <- 1 - stats::pweibull(u, shape, scale0)
        f / s
      },
      H0 = function(u) -log(1 - stats::pweibull(u, shape, scale0)),
      s = function(b, u) (1 - stats::pweibull(u, shape, scale0))^b,
      h = function(b, u) {
        f <- stats::dweibull(u, shape, scale0)
        s <- 1 - stats::pweibull(u, shape, scale0)
        b * f / s
      },
      H = function(b, u) -b * log(1 - stats::pweibull(u, shape, scale0))
    )
    
  } else if (dist == "LN") {
    # Log-normal distribution
    scale0 <- log(x0) - shape * stats::qnorm(1 - S0)
    
    list(
      s0 = function(u) 1 - stats::plnorm(u, scale0, shape),
      f0 = function(u) stats::dlnorm(u, scale0, shape),
      h0 = function(u) {
        f <- stats::dlnorm(u, scale0, shape)
        s <- 1 - stats::plnorm(u, scale0, shape)
        f / s
      },
      H0 = function(u) -log(1 - stats::plnorm(u, scale0, shape)),
      s = function(b, u) (1 - stats::plnorm(u, scale0, shape))^b,
      h = function(b, u) {
        f <- stats::dlnorm(u, scale0, shape)
        s <- 1 - stats::plnorm(u, scale0, shape)
        b * f / s
      },
      H = function(b, u) -b * log(1 - stats::plnorm(u, scale0, shape))
    )
    
  } else if (dist == "LG") {
    # Log-logistic distribution
    scale0 <- x0 / (1/S0 - 1)^(1/shape)
    
    list(
      s0 = function(u) 1 / (1 + (u/scale0)^shape),
      f0 = function(u) {
        (shape/scale0) * (u/scale0)^(shape - 1) / (1 + (u/scale0)^shape)^2
      },
      h0 = function(u) {
        num <- (shape/scale0) * (u/scale0)^(shape - 1)
        denom <- (1 + (u/scale0)^shape)^2
        s <- 1 / (1 + (u/scale0)^shape)
        (num / denom) / s
      },
      H0 = function(u) -log(1 / (1 + (u/scale0)^shape)),
      s = function(b, u) (1 / (1 + (u/scale0)^shape))^b,
      h = function(b, u) {
        f <- (shape/scale0) * (u/scale0)^(shape - 1) / (1 + (u/scale0)^shape)^2
        s <- 1 / (1 + (u/scale0)^shape)
        b * f / s
      },
      H = function(b, u) -b * log(1 / (1 + (u/scale0)^shape))
    )
    
  } else if (dist == "GM") {
    # Gamma distribution
    root0 <- function(t) 1 - stats::pgamma(x0, shape, t) - S0
    scale0 <- stats::uniroot(root0, c(0, 10))$root
    
    list(
      s0 = function(u) 1 - stats::pgamma(u, shape, scale0),
      f0 = function(u) stats::dgamma(u, shape, scale0),
      h0 = function(u) {
        f <- stats::dgamma(u, shape, scale0)
        s <- 1 - stats::pgamma(u, shape, scale0)
        f / s
      },
      H0 = function(u) -log(1 - stats::pgamma(u, shape, scale0)),
      s = function(b, u) (1 - stats::pgamma(u, shape, scale0))^b,
      h = function(b, u) {
        f <- stats::dgamma(u, shape, scale0)
        s <- 1 - stats::pgamma(u, shape, scale0)
        b * f / s
      },
      H = function(b, u) -b * log(1 - stats::pgamma(u, shape, scale0))
    )
    
  } else if (dist == "SP") {
    # Logspline distribution
    if (!requireNamespace("logspline", quietly = TRUE)) {
      stop("Package 'logspline' is required for dist = 'SP'")
    }
    
    time <- data$time
    status <- data$status
    fitSP <- logspline::oldlogspline(time[status == 1], time[status == 0], 
                                     lbound = 0)
    
    list(
      s0 = function(u) 1 - logspline::poldlogspline(u, fitSP),
      f0 = function(u) logspline::doldlogspline(u, fitSP),
      h0 = function(u) {
        f <- logspline::doldlogspline(u, fitSP)
        s <- 1 - logspline::poldlogspline(u, fitSP)
        f / s
      },
      H0 = function(u) -log(1 - logspline::poldlogspline(u, fitSP)),
      s = function(b, u) (1 - logspline::poldlogspline(u, fitSP))^b,
      h = function(b, u) {
        f <- logspline::doldlogspline(u, fitSP)
        s <- 1 - logspline::poldlogspline(u, fitSP)
        b * f / s
      },
      H = function(b, u) -b * log(1 - logspline::poldlogspline(u, fitSP))
    )
  }
}


#' Calculate Type I Error for Two-Stage Design
#'
#' @description
#' Internal helper function to calculate the Type I error probability for a 
#' two-stage design given critical values and correlation.
#'
#' @param c2 Numeric. Final stage critical value.
#' @param c1 Numeric. Stage 1 critical value.
#' @param rho0 Numeric. Correlation between stage 1 and final test statistics 
#'   under H0.
#'
#' @return Numeric. Type I error probability.
#'
#' @details
#' Computes P(Z1 > c1 and Z2 > c2 | H0) where (Z1, Z2) are bivariate normal
#' with correlation rho0. Uses numerical integration of the bivariate normal
#' density.
#'
#' @keywords internal
#' @noRd
.calculate_alpha <- function(c2, c1, rho0) {
  fun1 <- function(z, c1, rho0) {
    stats::dnorm(z) * stats::pnorm((rho0 * z - c1) / sqrt(1 - rho0^2))
  }
  alpha_val <- stats::integrate(fun1, lower = c2, upper = Inf, c1, rho0)$value
  return(alpha_val)
}


#' Calculate Power for Two-Stage Design
#'
#' @description
#' Internal helper function to calculate the power for a two-stage design 
#' given standardized critical values and correlation.
#'
#' @param cb Numeric. Standardized final stage critical value under H1.
#' @param cb1 Numeric. Standardized stage 1 critical value under H1.
#' @param rho1 Numeric. Correlation between stage 1 and final test statistics 
#'   under H1.
#'
#' @return Numeric. Power (probability of rejecting H0 when H1 is true).
#'
#' @details
#' Computes P(Z1 > cb1 and Z2 > cb | H1) where (Z1, Z2) are bivariate normal
#' with correlation rho1. The critical values cb and cb1 have been standardized
#' to account for the effect size under H1. Uses numerical integration of the 
#' bivariate normal density.
#'
#' @keywords internal
#' @noRd
.calculate_power <- function(cb, cb1, rho1) {
  fun2 <- function(z, cb1, rho1) {
    stats::dnorm(z) * stats::pnorm((rho1 * z - cb1) / sqrt(1 - rho1^2))
  }
  pwr <- stats::integrate(fun2, lower = cb, upper = Inf, cb1 = cb1, 
                         rho1 = rho1)$value
  return(pwr)
}
