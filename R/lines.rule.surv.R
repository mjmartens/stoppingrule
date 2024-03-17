#' @title Add Stopping Rule Curve to Current Plot (Survival Data)
#' @description
#' Add a survival stopping rule graphically as a curve on current plot for time-to-event data
#'
#'
#' @param x A \code{rule.surv} object calculated by \code{calc.rule.surv()} function
#' @param ... Other options to be passed to generic \code{lines()} function
#'
#' @return No return value, function solely modifies current plot
#' @export
#'
#' @examples
#' poc_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="Pocock",tau=100)
#' gp_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="GP",tau=100,param=c(1,1000))
#' plot(poc_rule)
#' lines(gp_rule,col="red")
#'

lines.rule.surv = function(x,...) {
  f = bdryfcn.surv(x$n,x$p0,x$type,x$tau,x$cval,x$param)
  curve(f,add=TRUE,...)
}
