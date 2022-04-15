#' @title Add Stopping Rule Curve to Current Plot
#' @description Add a stopping rule graphically as a curve on current plot
#'
#' @param x A rule object, being a matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param ... Other options to be passed to generic \code{lines} function
#'
#' @return No return value; function solely modifies current plot
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc_rule = calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
#' # Bayesian beta-binomial monitoring in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' bb_rule = calc.rule(ns=1:50,p0=0.20,type="BB",alpha=0.10,param=c(1,7))
#'
#' # Plot stopping boundaries for stopping rules
#' plot(poc_rule,col="blue")
#' lines(bb_rule,col="red")
lines.rule = function(x,...) {
  x = smooth.bnd(x)
    NextMethod("lines",type='l',...)
}
