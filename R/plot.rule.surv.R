#' @title Plot Stopping Rule (Survival Data)
#' @description
#' Display a stopping rule graphically as a curve for time-to-event data
#'
#' @param x A \code{rule.surv} object calculated by \code{calc.rule.surv()} function
#' @param xlim The x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a ‘reversed axis’.
#' @param ylim The y limits of the plot.
#' @param xlab The title for the x axis
#' @param ylab The title for the y axis
#' @param ... Other parameters passed to the \code{plot} function.
#'
#' @return No return value; function solely generates a plot
#' @export
#'
#' @examples
#' poc_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="Pocock",tau=100)
#' gp_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="GP",tau=100,param=c(1,1000))
#' plot(poc_rule)
#' lines(gp_rule,col="red")

plot.rule.surv <- function(x,xlim=c(0,max(x$Rule[,1])),ylim=c(0,max(x$Rule[,2])+1),
                           xlab="Total Exposure Time",ylab="# Events",...) {
  f = bdryfcn.surv(x$n,x$p0,x$type,x$tau,x$cval,x$param)
  curve(f,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
}
