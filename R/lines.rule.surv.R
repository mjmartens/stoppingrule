#' @title Add Stopping Rule Curve to Current Plot (Survival Data)
#' @description
#' Add a survival stopping rule graphically as a curve on current plot for time-to-event data
#'
#'
#' @param x A 'rule.surv' object calculated by \code{calc.rule.surv()} function
#' @param ... Other options to be passed to generic \code{lines()} function
#'
#' @return No return value, function solely modifies current plot
#' @export
#'
#' @examples
#' pocock.rule <- calc.rule.surv(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' OBF.rule <- calc.rule.surv(n = 30, tau = 100, p0 = 0.1, type = "OBF", alpha = 0.05)
#' plot(pocock.rule)
#' lines(OBF.rule, col = "blue")
#'
lines.rule.surv = function(x,...) {
  x = x$Rule
  NextMethod("lines",type='l',...)
}
