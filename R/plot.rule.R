#' @title Plot Stopping Rule
#' @description Display a stopping rule graphically as a curve
#'
#' @param x A rule object, being a matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param smooth Binary indicator of whether stopping rule boundary should be smoothed by linear interpolation between evaluation points
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
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc_rule = calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
#'
#' # Plot stopping boundary with smoothing
#' plot(poc_rule,col="blue")
plot.rule = function(x,smooth=TRUE,xlim=c(0,max(x[,1])),
                     ylim=c(0,max(x[,2])+1),xlab="# Evaluable",
                     ylab="# Events",...) {
  if(smooth==TRUE)  {x = smooth.bnd(x); ltype = 'l';}
  else {x = unclass(x); x = x[x[,1]>=x[,2],];ltype = 's';}
  plot(x,type=ltype,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
}
