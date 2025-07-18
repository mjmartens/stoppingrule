#' @title Plot Stopping Rule (Binary Data)
#' @description Display a stopping rule graphically as a curve
#'
#' @param x A \code{rule.bin} object calculated by \code{calc.rule.bin()} function
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
#'\dontrun{# Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="Pocock")
#'
#' # Bayesian beta-binomial method of Geller et al. in 50 patient cohort at 10% level,
#' # expected toxicity probability of 20%
#' bb_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="BB",param=c(2,8))
#'
#' # Plot stopping boundary with smoothing
#' plot(poc_rule,col="blue")
#' lines(bb_rule,col="red")
#'}

plot.rule.bin = function(x,smooth=TRUE,xlim=c(0,max(x$ns)),
                     ylim=c(0,max(x$Rule[,2])+1),xlab="# Evaluable",
                     ylab="# Events",...) {
  if(smooth==TRUE){
    f = bdryfcn.bin(max(x$ns),x$p0,x$type,x$cval,x$param)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,las=1)
    curve(f,xlim=c(which(x$Rule[,1]>=x$Rule[,2])[1],max(x$ns)),add=TRUE,...)
  }
  else if(smooth==FALSE){
    tab = x$Rule; tab = tab[tab[,1]>=tab[,2],]; ltype = 's';
    plot(tab,type=ltype,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,las=1,...)
  }
}
