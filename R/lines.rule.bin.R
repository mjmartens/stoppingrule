#' @title Add Stopping Rule Curve to Current Plot (Binary Data)
#' @description Add a binary stopping rule graphically as a curve on current plot
#'
#' @param x A 'rule.bin' object calculated by \code{calc.rule.bin()} function
#' @param smooth Binary indicator of whether stopping rule boundary should be smoothed by linear interpolation between evaluation points
#' @param ... Other options to be passed to generic \code{lines} function
#'
#' @return No return value; function solely modifies current plot
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="Pocock")
#' # Bayesian beta-binomial monitoring in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' bb_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="BB",param=c(1,7))
#'
#' # Plot stopping boundaries for stopping rules
#' plot(poc_rule,col="blue")
#' lines(bb_rule,col="red")

lines.rule.bin = function(x,smooth=TRUE,...) {
  if(smooth==TRUE)  {
    f = bdryfcn.bin(x$n,x$p0,x$cval,x$type,x$param)
    curve(f,add=TRUE,...)
  }
  else if(smooth==FALSE){
    tab = x$Rule; tab = tab[tab[,1]>=tab[,2],]; ltype = 's';
    lines(tab,type=ltype,...)
  }
}