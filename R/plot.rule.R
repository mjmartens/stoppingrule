#' @title Plot Stopping Rule
#' @description Display a stopping rule graphically as a curve in a plot
#'
#' @param rule A rule object, being a matrix with two columns: the vector of sample sizes at which the stopping rule is evaluated, and a corresponding vector of rejection boundaries at these sizes
#' @param add Indicator (TRUE/FALSE) of whether stopping rule curve should be added to existing plot
#' @param ...
#'
#' @export
#'
#' @examples
#' poc_rule = calc.rule(ns=1:30,p0=0.15,type="Pocock",alpha=0.10)
#' plot(poc_rule,col="blue")
plot.rule = function(rule,add=FALSE,...) {
  bdry_s = smooth.bnd(rule)
  if(add==FALSE) {
    plot(bdry_s[,1],bdry_s[,2],type='l',xlim=c(0,max(bdry_s[,1])),ylim=c(0,max(bdry_s[,2])+1),xlab="# Evaluable",ylab="# Events",...)
  }
  else {
    lines(bdry_s[,1],bdry_s[,2],type='l',...)
  }
}
