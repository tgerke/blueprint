#' Construct a trial design object
#'
#' @description
#' Low-level constructor establishing the common contract for all blueprint
#' design classes. Every design family constructor (e.g., [design_simon()])
#' returns an object built by `new_trial_design()`, so downstream generics
#' ([simulate_trial()], [evaluate()], [verify()]) can rely on a shared
#' structure.
#'
#' The contract: a `trial_design` is a list with elements
#' \describe{
#'   \item{design}{Named list of structural parameters (sample sizes,
#'     decision boundaries).}
#'   \item{hypotheses}{Named list of the hypotheses the design targets
#'     (e.g., `p0`, `pa`, and optionally the `alpha`/`power` targets used
#'     in a search).}
#'   \item{endpoint}{Either `"binary"` or `"tte"`.}
#'   \item{citation}{Optional character string citing the methodology.}
#' }
#'
#' @param subclass Character. The specific design class, e.g. `"simon_design"`.
#' @param design Named list of structural design parameters.
#' @param hypotheses Named list of design hypotheses.
#' @param endpoint Character. Endpoint type, `"binary"` or `"tte"`.
#' @param citation Optional character string citing the methodology.
#'
#' @return An object of class `c(subclass, "trial_design")`.
#'
#' @noRd
new_trial_design <- function(subclass,
                             design,
                             hypotheses,
                             endpoint = c("binary", "tte"),
                             citation = NULL) {
  endpoint <- match.arg(endpoint)

  if (!is.list(design) || is.null(names(design)) || any(names(design) == "")) {
    stop("`design` must be a fully named list.", call. = FALSE)
  }
  if (!is.list(hypotheses) || is.null(names(hypotheses)) ||
      any(names(hypotheses) == "")) {
    stop("`hypotheses` must be a fully named list.", call. = FALSE)
  }

  structure(
    list(
      design = design,
      hypotheses = hypotheses,
      endpoint = endpoint,
      citation = citation
    ),
    class = c(subclass, "trial_design")
  )
}


#' Test if an object is a trial design
#'
#' @param x An object.
#'
#' @return `TRUE` if `x` inherits from class `trial_design`, else `FALSE`.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' is_trial_design(design)
#'
#' @export
is_trial_design <- function(x) {
  inherits(x, "trial_design")
}
