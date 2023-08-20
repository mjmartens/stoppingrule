#' @title Plot Stopping Rule (Survival Data)
#' @description
#' Display a stopping rule graphically as a curve for time-to-event data
#'
#' @param x A 'rule.surv' object calculated by \code{calc.rule.surv()} function
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
#' pocock.rule <- calc.rule.surv(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' plot(pocock.rule, col = "red")

plot.rule.surv <- function(x,xlim=c(0,max(x$Rule[,1])),ylim=c(0,max(x$Rule[,2])+1),
                           xlab="Total Exposure Time",ylab="# Events",...) {
  f = bdryfcn.surv(x$n,x$p0,x$cval,x$tau,x$type,x$param)
  curve(f,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
}
