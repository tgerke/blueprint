#' Construct a Simon two-stage design
#'
#' @description
#' Constructs a Simon (1989) two-stage design for a single-arm trial with a
#' binary endpoint from its four structural parameters. Use this when the
#' design is already known (e.g., reproducing a design from a protocol);
#' use [search_simon_designs()] to find designs meeting target error rates.
#'
#' The design: enroll `n1` patients in stage 1 and stop for futility if `r1`
#' or fewer respond. Otherwise enroll `n - n1` more patients and reject the
#' null hypothesis if more than `r` of the `n` total patients respond.
#'
#' @param n1 Integer. Stage 1 sample size.
#' @param r1 Integer. Stage 1 futility boundary: stop if responses <= `r1`.
#' @param n Integer. Total sample size across both stages.
#' @param r Integer. Final boundary: reject the null if responses > `r`.
#' @param p0 Numeric. Response rate under the null hypothesis.
#' @param pa Numeric. Response rate under the alternative hypothesis
#'   (`pa > p0`).
#' @param alpha Optional numeric. The type I error target the design was
#'   selected to meet (metadata; the attained error is computed exactly).
#' @param power Optional numeric. The power target the design was selected
#'   to meet (metadata; the attained power is computed exactly).
#' @param criterion Optional character, `"optimal"` or `"minimax"`, recording
#'   how the design was selected.
#'
#' @return An object of class `c("simon_design", "trial_design")`.
#'
#' @references
#' Simon, R. (1989). Optimal two-stage designs for phase II clinical trials.
#' *Controlled Clinical Trials*, 10(1), 1-10.
#' \doi{10.1016/0197-2456(89)90015-9}
#'
#' @examples
#' # The classic p0 = 0.2 vs. pa = 0.4 optimal design
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' design
#' evaluate(design)
#'
#' @seealso [search_simon_designs()], [evaluate()], [verify()],
#'   [simulate_trial()], [draft_protocol_text()]
#'
#' @export
design_simon <- function(n1, r1, n, r, p0, pa,
                         alpha = NULL, power = NULL, criterion = NULL) {
  check_count <- function(x, name, min) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != trunc(x) ||
        x < min) {
      stop("`", name, "` must be a single whole number >= ", min, ".",
           call. = FALSE)
    }
    as.integer(x)
  }
  n1 <- check_count(n1, "n1", 1)
  r1 <- check_count(r1, "r1", 0)
  n <- check_count(n, "n", 2)
  r <- check_count(r, "r", 0)

  if (r1 >= n1) {
    stop("`r1` must be less than `n1`; otherwise the trial always stops at stage 1.",
         call. = FALSE)
  }
  if (n <= n1) {
    stop("`n` must be greater than `n1`.", call. = FALSE)
  }
  if (r < r1 || r >= n) {
    stop("`r` must satisfy `r1 <= r < n`.", call. = FALSE)
  }
  check_prob <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0 || x >= 1) {
      stop("`", name, "` must be a single number strictly between 0 and 1.",
           call. = FALSE)
    }
    x
  }
  p0 <- check_prob(p0, "p0")
  pa <- check_prob(pa, "pa")
  if (p0 >= pa) {
    stop("`pa` must be greater than `p0`.", call. = FALSE)
  }
  if (!is.null(criterion)) {
    criterion <- match.arg(criterion, c("optimal", "minimax"))
  }

  out <- new_trial_design(
    subclass = "simon_design",
    design = list(n1 = n1, r1 = r1, n = n, r = r),
    hypotheses = list(p0 = p0, pa = pa, alpha = alpha, power = power),
    endpoint = "binary",
    citation = paste(
      "Simon, R. (1989). Optimal two-stage designs for phase II clinical",
      "trials. Controlled Clinical Trials, 10(1), 1-10."
    )
  )
  out$criterion <- criterion
  out
}


# Exact operating characteristics of a Simon two-stage design at true rate p.
# Success means rejecting H0 (responses > r after both stages).
#' @noRd
.simon_exact_ocs <- function(n1, r1, n, r, p) {
  n2 <- n - n1
  pet <- stats::pbinom(r1, n1, p)
  x1 <- seq.int(r1 + 1, n1)
  prob_success <- sum(
    stats::dbinom(x1, n1, p) * (1 - stats::pbinom(r - x1, n2, p))
  )
  list(
    prob_success = prob_success,
    prob_early_stop = pet,
    expected_n = n1 + (1 - pet) * n2
  )
}


#' Search for Simon two-stage designs meeting target error rates
#'
#' @description
#' Searches over stage sizes and boundaries for Simon (1989) two-stage
#' designs with type I error at most `alpha` and power at least `power`.
#' For each total sample size `n`, the design minimizing the expected sample
#' size under the null is retained; the globally best design is flagged
#' `"optimal"` and the design with the smallest maximum sample size is
#' flagged `"minimax"`.
#'
#' Search and construction are separate steps: this function returns a
#' tibble of candidates, and [pick_design()] (or [design_simon()] with any
#' row's parameters) turns a candidate into a design object.
#'
#' @param p0 Numeric. Response rate under the null hypothesis.
#' @param pa Numeric. Response rate under the alternative (`pa > p0`).
#' @param alpha Numeric. Maximum type I error (one-sided). Default 0.05.
#' @param power Numeric. Minimum power. Default 0.80.
#' @param n_max Integer. Largest total sample size to search. Default 100.
#'
#' @return A tibble of class `simon_search` with one row per candidate total
#'   sample size and columns `n1`, `r1`, `n`, `r`, `type1_error`, `power`,
#'   `pet_null`, `expected_n_null`, and `criterion` (`"optimal"`,
#'   `"minimax"`, both, or `NA`).
#'
#' @references
#' Simon, R. (1989). Optimal two-stage designs for phase II clinical trials.
#' *Controlled Clinical Trials*, 10(1), 1-10.
#' \doi{10.1016/0197-2456(89)90015-9}
#'
#' @examples
#' candidates <- search_simon_designs(p0 = 0.05, pa = 0.25, n_max = 30)
#' candidates
#' pick_design(candidates, "optimal")
#'
#' @seealso [pick_design()], [design_simon()]
#'
#' @export
search_simon_designs <- function(p0, pa, alpha = 0.05, power = 0.80,
                                 n_max = 100) {
  if (!is.numeric(p0) || !is.numeric(pa) || p0 <= 0 || pa >= 1 || p0 >= pa) {
    stop("Require 0 < `p0` < `pa` < 1.", call. = FALSE)
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 ||
      !is.numeric(power) || power <= 0 || power >= 1) {
    stop("`alpha` and `power` must be strictly between 0 and 1.",
         call. = FALSE)
  }

  rows <- list()
  for (n in 2:n_max) {
    best <- NULL
    for (n1 in 1:(n - 1)) {
      n2 <- n - n1
      x1 <- 0:n1
      d0 <- stats::dbinom(x1, n1, p0)
      d1 <- stats::dbinom(x1, n1, pa)
      rs <- 0:(n - 1)
      # W[i, j] = P(x1 responses in stage 1) * P(stage 2 pushes total past r)
      # for r = rs[i], x1 = j - 1
      diff <- outer(rs, x1, "-")
      W0 <- sweep(1 - stats::pbinom(diff, n2, p0), 2, d0, "*")
      W1 <- sweep(1 - stats::pbinom(diff, n2, pa), 2, d1, "*")
      # Reverse cumulative sums over x1 so column k holds the rejection
      # probability when continuing requires x1 >= k - 1
      rev_cols <- rev(seq_len(n1 + 1))
      S0 <- t(apply(W0[, rev_cols, drop = FALSE], 1, cumsum))[, rev_cols,
                                                              drop = FALSE]
      S1 <- t(apply(W1[, rev_cols, drop = FALSE], 1, cumsum))[, rev_cols,
                                                              drop = FALSE]

      for (r1 in 0:(n1 - 1)) {
        pet0 <- stats::pbinom(r1, n1, p0)
        en0 <- n1 + (1 - pet0) * n2
        if (!is.null(best) && en0 >= best$expected_n_null) next

        a_vec <- S0[, r1 + 2]
        # smallest r >= r1 meeting the alpha constraint (a_vec decreases in r)
        ok <- which(a_vec <= alpha & rs >= r1)
        if (length(ok) == 0) next
        r_idx <- ok[1]
        pow <- S1[r_idx, r1 + 2]
        if (pow < power) next

        best <- list(
          n1 = n1, r1 = r1, n = n, r = rs[r_idx],
          type1_error = a_vec[r_idx], power = pow,
          pet_null = pet0, expected_n_null = en0
        )
      }
    }
    if (!is.null(best)) rows[[length(rows) + 1]] <- best
  }

  if (length(rows) == 0) {
    stop("No design with n <= ", n_max, " meets the constraints. ",
         "Increase `n_max`.", call. = FALSE)
  }

  out <- tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame)))
  out$criterion <- NA_character_
  out$criterion[which.min(out$expected_n_null)] <- "optimal"
  i_minimax <- which(out$n == min(out$n))
  out$criterion[i_minimax] <- ifelse(
    is.na(out$criterion[i_minimax]), "minimax",
    paste(out$criterion[i_minimax], "minimax", sep = ", ")
  )

  if (out$n[which(out$criterion %in% c("optimal", "optimal, minimax"))] == n_max) {
    warning("The optimal design sits at the search boundary (n = ", n_max,
            "); a larger `n_max` may find a design with smaller expected ",
            "sample size.", call. = FALSE)
  }

  attr(out, "p0") <- p0
  attr(out, "pa") <- pa
  attr(out, "alpha") <- alpha
  attr(out, "power") <- power
  class(out) <- c("simon_search", class(out))
  out
}


#' @export
print.simon_search <- function(x, ...) {
  cat("<Simon two-stage design search>\n")
  cat(sprintf("  p0 = %g, pa = %g | type I error <= %g, power >= %g\n",
              attr(x, "p0"), attr(x, "pa"),
              attr(x, "alpha"), attr(x, "power")))
  cat("  One candidate per total n (minimum expected null sample size);\n")
  cat("  pass to pick_design() to construct the optimal or minimax design.\n\n")
  df <- x
  class(df) <- setdiff(class(df), "simon_search")
  print(df, ...)
  invisible(x)
}


#' Pick a design from a search result
#'
#' @description
#' Turns one row of a design search (e.g., [search_simon_designs()]) into a
#' design object by the named selection criterion.
#'
#' @param candidates A search result, e.g. from [search_simon_designs()].
#' @param criterion Character. `"optimal"` (minimum expected sample size
#'   under the null) or `"minimax"` (minimum maximum sample size).
#' @param ... Passed to methods.
#'
#' @return A `trial_design` object.
#'
#' @examples
#' candidates <- search_simon_designs(p0 = 0.05, pa = 0.25, n_max = 30)
#' pick_design(candidates, "minimax")
#'
#' @export
pick_design <- function(candidates, criterion = c("optimal", "minimax"), ...) {
  UseMethod("pick_design")
}

#' @export
pick_design.default <- function(candidates, criterion = c("optimal", "minimax"),
                                ...) {
  stop("No `pick_design()` method for class <", class(candidates)[1], ">.",
       call. = FALSE)
}

#' @rdname pick_design
#' @export
pick_design.simon_search <- function(candidates,
                                     criterion = c("optimal", "minimax"),
                                     ...) {
  criterion <- match.arg(criterion)
  i <- grep(criterion, candidates$criterion)
  row <- candidates[i, ]
  design_simon(
    n1 = row$n1, r1 = row$r1, n = row$n, r = row$r,
    p0 = attr(candidates, "p0"), pa = attr(candidates, "pa"),
    alpha = attr(candidates, "alpha"), power = attr(candidates, "power"),
    criterion = criterion
  )
}


#' @export
assume_null.simon_design <- function(design, ...) {
  assume_response(design$hypotheses$p0)
}

#' @export
assume_alternative.simon_design <- function(design, ...) {
  assume_response(design$hypotheses$pa)
}


# Validate that a scenario works for a single-arm binary design
#' @noRd
.check_response_scenario <- function(under) {
  if (missing(under) || is.null(under)) {
    stop("Supply a scenario via `under`, e.g. `under = assume_null(design)` ",
         "or `under = assume_response(p = 0.3)`.", call. = FALSE)
  }
  if (!inherits(under, "response_scenario")) {
    stop("`under` must be a response scenario; see `assume_response()`.",
         call. = FALSE)
  }
  if (length(under$p) != 1) {
    stop("Simon designs are single-arm: `under` must have a single `p`.",
         call. = FALSE)
  }
  under
}


#' Simulate trials from a Simon two-stage design
#'
#' @param design A `simon_design` object.
#' @param under A `scenario` from [assume_response()], [assume_null()], or
#'   [assume_alternative()].
#' @param n_sims Integer. Number of trials to simulate.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Not used.
#'
#' @return A tibble with one row per simulated trial: `sim`,
#'   `responses_stage1`, `stopped_early`, `n_enrolled`, `responses_total`,
#'   and `success` (whether the null was rejected).
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' simulate_trial(design, under = assume_alternative(design), n_sims = 10)
#'
#' @export
simulate_trial.simon_design <- function(design, under, n_sims = 1,
                                        seed = NULL, ...) {
  under <- .check_response_scenario(under)
  d <- design$design
  n2 <- d$n - d$n1
  if (!is.null(seed)) set.seed(seed)

  x1 <- stats::rbinom(n_sims, d$n1, under$p)
  continue <- x1 > d$r1
  x2 <- integer(n_sims)
  x2[continue] <- stats::rbinom(sum(continue), n2, under$p)
  total <- x1 + x2

  tibble::tibble(
    sim = seq_len(n_sims),
    responses_stage1 = x1,
    stopped_early = !continue,
    n_enrolled = d$n1 + n2 * as.integer(continue),
    responses_total = total,
    success = continue & total > d$r
  )
}


#' Evaluate operating characteristics of a Simon two-stage design
#'
#' @param design A `simon_design` object.
#' @param under A `scenario`, or a named list of scenarios. Defaults to the
#'   design's null and alternative hypotheses.
#' @param method `"exact"` for closed-form operating characteristics or
#'   `"simulation"` for Monte Carlo estimates.
#' @param n_sims Integer. Number of simulated trials per scenario when
#'   `method = "simulation"`.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Not used.
#'
#' @return A tibble with one row per scenario: `scenario`, `p_true`,
#'   `prob_success` (probability of rejecting the null; type I error under
#'   the null scenario, power under the alternative), `prob_early_stop`,
#'   `expected_n`, and `method`.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' evaluate(design)
#' evaluate(design, under = assume_response(0.3))
#'
#' @export
evaluate.simon_design <- function(design, under = NULL,
                                  method = c("exact", "simulation"),
                                  n_sims = 10000, seed = NULL, ...) {
  method <- match.arg(method)

  if (is.null(under)) {
    scenarios <- list(
      null = assume_null(design),
      alternative = assume_alternative(design)
    )
  } else if (inherits(under, "scenario")) {
    scenarios <- list(custom = under)
  } else if (is.list(under)) {
    scenarios <- under
    if (is.null(names(scenarios))) {
      names(scenarios) <- paste0("scenario_", seq_along(scenarios))
    }
  } else {
    stop("`under` must be a scenario or a list of scenarios.", call. = FALSE)
  }
  scenarios <- lapply(scenarios, .check_response_scenario)

  d <- design$design
  if (method == "simulation" && !is.null(seed)) set.seed(seed)

  ocs <- lapply(scenarios, function(sc) {
    if (method == "exact") {
      .simon_exact_ocs(d$n1, d$r1, d$n, d$r, sc$p)
    } else {
      sims <- simulate_trial(design, under = sc, n_sims = n_sims)
      list(
        prob_success = mean(sims$success),
        prob_early_stop = mean(sims$stopped_early),
        expected_n = mean(sims$n_enrolled)
      )
    }
  })

  tibble::tibble(
    scenario = names(scenarios),
    p_true = unname(vapply(scenarios, function(sc) sc$p, numeric(1))),
    prob_success = unname(vapply(ocs, `[[`, numeric(1), "prob_success")),
    prob_early_stop = unname(vapply(ocs, `[[`, numeric(1), "prob_early_stop")),
    expected_n = unname(vapply(ocs, `[[`, numeric(1), "expected_n")),
    method = method
  )
}


#' @export
print.simon_design <- function(x, ...) {
  d <- x$design
  h <- x$hypotheses
  oc0 <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$p0)
  oca <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$pa)

  label <- if (is.null(x$criterion)) {
    "<Simon two-stage design>"
  } else {
    sprintf("<Simon two-stage design: %s>", x$criterion)
  }
  cat(label, "\n")
  cat(sprintf("  Hypotheses: H0 p = %g vs. H1 p = %g (one-sided)\n",
              h$p0, h$pa))
  cat(sprintf("  Stage 1: enroll %d; stop for futility if <= %d responses\n",
              d$n1, d$r1))
  cat(sprintf("  Stage 2: enroll %d more (%d total); reject H0 if >= %d responses\n",
              d$n - d$n1, d$n, d$r + 1))
  cat("  Exact operating characteristics:\n")
  cat(sprintf("    type I error %.4f | power %.4f\n",
              oc0$prob_success, oca$prob_success))
  cat(sprintf("    PET(H0) %.3f | E[N | H0] %.1f\n",
              oc0$prob_early_stop, oc0$expected_n))
  invisible(x)
}


#' Summarize a Simon two-stage design
#'
#' @param object A `simon_design` object.
#' @param ... Not used.
#'
#' @return A `summary.simon_design` object: the design plus its exact
#'   operating characteristics table.
#'
#' @export
summary.simon_design <- function(object, ...) {
  structure(
    list(design = object, ocs = evaluate(object, method = "exact")),
    class = "summary.simon_design"
  )
}

#' @export
print.summary.simon_design <- function(x, ...) {
  print(x$design)
  cat("\nOperating characteristics (exact):\n")
  print(x$ocs)
  invisible(x)
}


#' Tidy a Simon two-stage design
#'
#' @description
#' `tidy()` returns the design's exact operating characteristics under its
#' null and alternative hypotheses as a tibble; `glance()` returns a one-row
#' summary of the design.
#'
#' @param x A `simon_design` object.
#' @param ... Not used.
#'
#' @return A tibble.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' tidy(design)
#' glance(design)
#'
#' @export
tidy.simon_design <- function(x, ...) {
  evaluate(x, method = "exact")
}

#' @rdname tidy.simon_design
#' @export
glance.simon_design <- function(x, ...) {
  d <- x$design
  h <- x$hypotheses
  oc0 <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$p0)
  oca <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$pa)
  tibble::tibble(
    n1 = d$n1, r1 = d$r1, n = d$n, r = d$r,
    p0 = h$p0, pa = h$pa,
    type1_error = oc0$prob_success,
    power = oca$prob_success,
    pet_null = oc0$prob_early_stop,
    expected_n_null = oc0$expected_n,
    criterion = ifelse(is.null(x$criterion), NA_character_, x$criterion)
  )
}


#' @export
draft_protocol_text.simon_design <- function(design, ...) {
  d <- design$design
  h <- design$hypotheses
  oc0 <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$p0)
  oca <- .simon_exact_ocs(d$n1, d$r1, d$n, d$r, h$pa)

  design_name <- if (is.null(design$criterion)) {
    "Simon two-stage design"
  } else {
    sprintf("Simon %s two-stage design", design$criterion)
  }
  pct <- function(p) sprintf("%g%%", 100 * p)

  paste0(
    "A ", design_name, " (Simon, 1989) will be used. ",
    "The null hypothesis that the true response rate is ", pct(h$p0),
    " will be tested against a one-sided alternative. ",
    "In the first stage, ", d$n1, " patients will be accrued. ",
    "If there are ", d$r1, " or fewer responses in these ", d$n1,
    " patients, the study will be stopped. ",
    "Otherwise, ", d$n - d$n1, " additional patients will be accrued ",
    "for a total of ", d$n, ". ",
    "The null hypothesis will be rejected if ", d$r + 1,
    " or more responses are observed in ", d$n, " patients. ",
    "This design yields a type I error rate of ",
    sprintf("%.3f", oc0$prob_success), " and power of ",
    sprintf("%.3f", oca$prob_success), " when the true response rate is ",
    pct(h$pa), ". ",
    "Under the null hypothesis, the probability of early termination is ",
    round(oc0$prob_early_stop, 2), " and the expected sample size is ",
    round(oc0$expected_n, 1), " patients."
  )
}
