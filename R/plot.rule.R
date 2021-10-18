#' @title Plot Stopping Rule
#' @description Display a stopping rule graphically as a curve
#'
#' @param rule A rule object, being a matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param smooth Indicator (TRUE/FALSE) of whether stopping rule curve should be smoothed
#' @param add Indicator (TRUE/FALSE) of whether stopping rule curve should be added to existing plot
#' @param ...
#'
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc_rule = calc.rule(ns=1:30,p0=0.15,type="Pocock",alpha=0.10)
#'
#' # Plot stopping boundary with smoothing
#' plot(poc_rule,col="blue")
plot.rule = function(rule,smooth=TRUE,add=FALSE,...) {
  bdry_s = smooth.bnd(rule)
  if(add==FALSE) {
    plot(bdry_s[,1],bdry_s[,2],type='l',xlim=c(0,max(bdry_s[,1])),ylim=c(0,max(bdry_s[,2])+1),xlab="# Evaluable",ylab="# Events",...)
  }
  else {
    lines(bdry_s[,1],bdry_s[,2],type='l',...)
  }
}
