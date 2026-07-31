#' Assume a true response rate scenario
#'
#' @description
#' Creates a scenario object describing a data-generating truth for binary
#' (response) endpoints. Scenarios are deliberately separate from designs:
#' a design encodes the decision structure, while a scenario encodes what
#' the world is actually like. Pass scenarios to [simulate_trial()] and
#' [evaluate()] to study a design under any truth, not just the hypotheses
#' it was built for.
#'
#' @param p Numeric vector of true response probabilities in (0, 1). Length 1
#'   for single-arm designs; length K for basket designs.
#'
#' @return An object of class `c("response_scenario", "scenario")`.
#'
#' @examples
#' # Truth halfway between a null of 0.2 and an alternative of 0.4
#' assume_response(p = 0.3)
#'
#' @seealso [assume_null()] and [assume_alternative()] to derive scenarios
#'   directly from a design's hypotheses; [assume_survival()] for
#'   time-to-event endpoints.
#'
#' @export
assume_response <- function(p) {
  if (!is.numeric(p) || length(p) < 1 || any(p <= 0) || any(p >= 1)) {
    stop("`p` must be numeric with all values strictly between 0 and 1.",
         call. = FALSE)
  }
  structure(
    list(p = p),
    class = c("response_scenario", "scenario")
  )
}


#' Assume a true survival distribution scenario
#'
#' @description
#' Creates a scenario object describing a data-generating truth for
#' time-to-event endpoints.
#'
#' @param dist Character. Survival distribution: one of `"weibull"`,
#'   `"lognormal"`, `"gamma"`, or `"loglogistic"`.
#' @param shape Numeric. Shape parameter of the distribution.
#' @param scale Numeric. Scale parameter of the distribution.
#'
#' @return An object of class `c("survival_scenario", "scenario")`.
#'
#' @examples
#' assume_survival(dist = "weibull", shape = 1.5, scale = 4.2)
#'
#' @seealso [assume_response()] for binary endpoints.
#'
#' @export
assume_survival <- function(dist = c("weibull", "lognormal", "gamma", "loglogistic"),
                            shape,
                            scale) {
  dist <- match.arg(dist)
  if (!is.numeric(shape) || length(shape) != 1 || shape <= 0) {
    stop("`shape` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("`scale` must be a single positive number.", call. = FALSE)
  }
  structure(
    list(dist = dist, shape = shape, scale = scale),
    class = c("survival_scenario", "scenario")
  )
}


#' Derive the null scenario from a design
#'
#' @description
#' Returns the scenario in which the design's null hypothesis is true,
#' constructed from the hypotheses stored in the design object. Users never
#' need to re-derive data-generating parameters by hand.
#'
#' @param design A `trial_design` object.
#' @param ... Passed to methods.
#'
#' @return A `scenario` object.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' assume_null(design)
#'
#' @export
assume_null <- function(design, ...) {
  UseMethod("assume_null")
}

#' @export
assume_null.default <- function(design, ...) {
  stop("No `assume_null()` method for class <", class(design)[1], ">.",
       call. = FALSE)
}


#' Derive the alternative scenario from a design
#'
#' @description
#' Returns the scenario in which the design's alternative hypothesis is true,
#' constructed from the hypotheses stored in the design object.
#'
#' @param design A `trial_design` object.
#' @param ... Passed to methods.
#'
#' @return A `scenario` object.
#'
#' @examples
#' design <- design_simon(n1 = 13, r1 = 3, n = 43, r = 12, p0 = 0.2, pa = 0.4)
#' assume_alternative(design)
#'
#' @export
assume_alternative <- function(design, ...) {
  UseMethod("assume_alternative")
}

#' @export
assume_alternative.default <- function(design, ...) {
  stop("No `assume_alternative()` method for class <", class(design)[1], ">.",
       call. = FALSE)
}


#' @export
print.response_scenario <- function(x, ...) {
  cat("<scenario: response>\n")
  cat("  true response rate(s):", paste(format(x$p), collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.survival_scenario <- function(x, ...) {
  cat("<scenario: survival>\n")
  cat(sprintf("  distribution: %s (shape = %s, scale = %s)\n",
              x$dist, format(x$shape), format(x$scale)))
  invisible(x)
}
