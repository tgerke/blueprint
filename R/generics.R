#' Simulate trials from a design
#'
#' @description
#' Simulates trial realizations from a design object under a data-generating
#' scenario. Methods exist for each blueprint design class; new design
#' families plug in by providing a method.
#'
#' For the modern design classes (e.g., [design_simon()]), pass a `scenario`
#' created with [assume_response()], [assume_null()], or
#' [assume_alternative()] via the `under` argument, and the result is a tidy
#' tibble with one row per simulated trial. Legacy classes
#' (`two_stage_design`, `single_stage_design`, `basket_trial_design`) retain
#' their original arguments (e.g., `shape`, `scale`, `p_true`) and return
#' structure.
#'
#' @param design A `trial_design` object (or a legacy blueprint design).
#' @param ... Method-specific arguments; see the method documentation.
#'
#' @return Method-specific; tidy design classes return a tibble with one row
#'   per simulated trial.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' simulate_trial(design, under = assume_alternative(design), n_sims = 10)
#'
#' @seealso [evaluate()] for operating characteristics, [verify()] to check
#'   simulated against exact operating characteristics.
#'
#' @export
simulate_trial <- function(design, ...) {
  UseMethod("simulate_trial")
}

#' @export
simulate_trial.default <- function(design, ...) {
  stop(
    "No `simulate_trial()` method for class <", class(design)[1], ">.",
    call. = FALSE
  )
}


#' Evaluate operating characteristics of a design
#'
#' @description
#' Computes the operating characteristics of a design under one or more
#' data-generating scenarios, returning a tidy tibble. Where a design family
#' supports it (e.g., Simon two-stage), characteristics are available in
#' closed form (`method = "exact"`); simulation is always available
#' (`method = "simulation"`).
#'
#' @param design A `trial_design` object.
#' @param ... Method-specific arguments. Common ones:
#'   \describe{
#'     \item{under}{A `scenario` (or named list of scenarios). Defaults to
#'       the design's own null and alternative hypotheses.}
#'     \item{method}{`"exact"` or `"simulation"`.}
#'     \item{n_sims}{Number of simulated trials when `method = "simulation"`.}
#'   }
#'
#' @return A tibble with one row per scenario and columns for the operating
#'   characteristics (e.g., `prob_success`, `prob_early_stop`, `expected_n`).
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' evaluate(design)
#' evaluate(design, under = assume_response(0.3))
#'
#' @seealso [verify()] to compare exact and simulated results.
#'
#' @export
evaluate <- function(design, ...) {
  UseMethod("evaluate")
}

#' @export
evaluate.default <- function(design, ...) {
  stop(
    "No `evaluate()` method for class <", class(design)[1], ">.",
    call. = FALSE
  )
}


#' Verify a design's operating characteristics by simulation
#'
#' @description
#' The package's core philosophy made executable: after designing to target
#' operating characteristics, simulate to confirm they are met. `verify()`
#' computes exact and simulated operating characteristics side by side under
#' the design's null and alternative hypotheses and reports their agreement.
#'
#' @param design A `trial_design` object whose [evaluate()] method supports
#'   both `method = "exact"` and `method = "simulation"`.
#' @param n_sims Integer. Number of simulated trials per scenario.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Passed to methods.
#'
#' @return A tibble of class `trial_verification` with columns `scenario`,
#'   `metric`, `exact`, `simulated`, `difference`, and `mc_se` (the Monte
#'   Carlo standard error, for probability metrics).
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' verify(design, n_sims = 1000, seed = 123)
#'
#' @export
verify <- function(design, ...) {
  UseMethod("verify")
}

#' @export
verify.default <- function(design, ...) {
  stop(
    "No `verify()` method for class <", class(design)[1], ">.",
    call. = FALSE
  )
}

#' @rdname verify
#' @export
verify.trial_design <- function(design, n_sims = 10000, seed = NULL, ...) {
  exact <- evaluate(design, method = "exact")
  simulated <- evaluate(design, method = "simulation", n_sims = n_sims,
                        seed = seed)

  metrics <- setdiff(
    intersect(names(exact), names(simulated)),
    c("scenario", "p_true", "method")
  )

  rows <- lapply(seq_len(nrow(exact)), function(i) {
    tibble::tibble(
      scenario = exact$scenario[i],
      metric = metrics,
      exact = unlist(exact[i, metrics], use.names = FALSE),
      simulated = unlist(simulated[i, metrics], use.names = FALSE)
    )
  })
  out <- do.call(rbind, rows)
  out$difference <- out$simulated - out$exact
  out$mc_se <- NA_real_
  prob <- startsWith(out$metric, "prob_")
  out$mc_se[prob] <- sqrt(out$exact[prob] * (1 - out$exact[prob]) / n_sims)

  structure(out, class = c("trial_verification", class(out)), n_sims = n_sims)
}

#' @export
print.trial_verification <- function(x, ...) {
  n_sims <- attr(x, "n_sims")
  cat(sprintf("<verification: exact vs. simulated (%s trials/scenario)>\n",
              format(n_sims, big.mark = ",")))
  df <- x
  class(df) <- setdiff(class(df), "trial_verification")
  print(df, ...)

  prob <- !is.na(x$mc_se)
  flagged <- prob & abs(x$difference) > 3.5 * x$mc_se
  if (any(flagged)) {
    cat("\nNote: some probability metrics differ from exact values by more",
        "than 3.5 Monte Carlo standard errors:\n")
    cat(paste0("  - ", x$scenario[flagged], " / ", x$metric[flagged],
               collapse = "\n"), "\n")
  } else {
    cat("\nAll probability metrics are within 3.5 Monte Carlo standard",
        "errors of their exact values.\n")
  }
  invisible(x)
}


#' Draft protocol text describing a design
#'
#' @description
#' Drafts the sample size justification paragraph for a design: the decision
#' structure, hypotheses, operating characteristics, and citation, in the
#' phrasing conventional for clinical trial protocols. Returns a character
#' string so the text drops directly into a protocol document or a Quarto
#' report.
#'
#' @param design A `trial_design` object.
#' @param ... Passed to methods.
#'
#' @return A character string.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' cat(draft_protocol_text(design))
#'
#' @export
draft_protocol_text <- function(design, ...) {
  UseMethod("draft_protocol_text")
}

#' @export
draft_protocol_text.default <- function(design, ...) {
  stop(
    "No `draft_protocol_text()` method for class <", class(design)[1], ">.",
    call. = FALSE
  )
}


#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance
