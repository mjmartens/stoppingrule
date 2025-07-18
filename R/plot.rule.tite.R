#' @title Plot Stopping Rule (TITE Method)
#' @description Display a stopping rule graphically as a curve
#'
#' @param x A \code{rule.tite} object calculated by \code{calc.rule.tite()} function
#' @param xlim The x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a ‘reversed axis’.
#' @param ylim The y limits of the plot.
#' @param xlab The title for the x axis
#' @param ylab The title for the y axis
#' @param ... Other options to be passed to generic \code{plot} function
#'
#' @return No return value; function solely generates a plot
#' @export
#'
#' @examples
#'\dontrun{# Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc_rule = calc.rule.tite(n=50,p0=0.20,alpha=0.10,type="Pocock")
#'
#' # Bayesian beta-extended binomial method in 50 patient cohort at 10% level,
#' # expected toxicity probability of 20%
#' bb_rule = calc.rule.tite(n=50,p0=0.20,alpha=0.10,type="BB",param=c(2,8))
#'
#' # Plot stopping boundary
#' plot(poc_rule,col="blue")
#' lines(bb_rule,col="red")
#'}

plot.rule.tite = function(x,xlim=c(0,x$n),
                         ylim=c(0,max(x$Rule[,2])+1),xlab=" Effective Sample Size",
                         ylab="# Events",...) {
  f = bdryfcn.bin(max(x$ns),x$p0,x$type,x$cval,x$param)
  plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,las=1)
  curve(f,xlim=c(which(x$Rule[,1]>=x$Rule[,2])[1],x$n),add=TRUE,...)
}
