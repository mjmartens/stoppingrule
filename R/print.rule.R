#' @title Print Stopping Rule
#' @description Display a stopping rule in tabular form
#'
#' @param x A rule object, containing a matrix describing a stopping rule
#' @param ... Other options to be passed to generic \code{plot} function
#'
#' @return A matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc_rule = calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
#'
#' # Print stopping rule in table
#' print(poc_rule)
print.rule = function(x,...) {
  NextMethod("print",...)
}
