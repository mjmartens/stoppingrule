#' @title Plot Stopping Rule
#' @description Display a stopping rule graphically as a curve
#'
#' @param x A rule object, being a matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param ... Other options to be passed to generic \code{plot} function
#'
#' @return No return value; function solely generates a plot
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc_rule = calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
#'
#' # Plot stopping boundary with smoothing
#' plot(poc_rule,col="blue")
plot.rule = function(x,...) {
  x = smooth.bnd(x)
    NextMethod("plot",type='l',xlim=c(0,max(x[,1])),ylim=c(0,max(x[,2])+1),xlab="# Evaluable",ylab="# Events",...)
}
