#' @title Search for Calibration Value (Survival Data)
#' @description
#' Internal workhorse function to calculate the calibration constant value that attains level alpha for given method for time-to-event data
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Length of observation period
#' @param p0 The probability of a toxicity occurring in \code{tau} units of time under the null hypothesis
#' @param alpha The nominal type I error/false positive rate for the stopping rule,
#' under an assumption that the cumulative number of events follows a Poisson process
#' over the study duration.
#' @param type The method used for constructing the stopping rule
#' @param maxInf Specification of the maximum information (maximum exposure time) used for designing the
#' stopping rule. Options include the expected exposure time for n patients used H0 ("expected") and the
#' maximum possible exposure time ("maximum"). Default is "expected" (expected exposure time in cohort).
#' @param param A vector of the extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For truncated SPRT, this is the targeted alternative toxicity probability p1. For Bayesian Gamma-Poisson model, this is the vector of hyperparameters (shape,rate) for the gamma prior on the toxicity event rate.
#'
#' @return The calibration constant used for subsequent stopping boundary calculation

findconst.surv = function(n, p0, alpha, type, tau, maxInf="expected", param = NULL){
  inner <- function(n, tau, p0, type, cval, param){
    bnd = calc.bnd.surv(n=n,p0=p0,type=type,tau=tau,cval=cval,maxInf=maxInf,param=param)
    lambda0 <- -log(1 - p0)/tau
    if(maxInf=="maximum") {Umax = n*tau}
    else if(maxInf=="expected") {Umax = n*p0/lambda0}
    bnd$Rule[which(bnd$Rule[,1]>Umax),1] = Umax

    return(stopping.prob.surv(bnd, p = p0))
  }

  if (type == 'Pocock'|type == 'OBF'|type == 'WT'){
    l = qnorm(1-alpha)
    u = qnorm(1 - alpha/(n*p0))
  } else if (type == "GP"){
    lambda0 <- -log(1 - p0)/tau
    if(maxInf=="maximum") {Umax = n*tau}
    else if(maxInf=="expected") {Umax = n*p0/lambda0}

    if (length(param) == 1){
      k <- param
      s <- k/lambda0
    } else {
      k <- param[1] # shape
      s <- param[2] # rate
    }
    l = (1 - pgamma(lambda0, shape = (k + 1), rate = (s + Umax))) + 0.01
    u = 0.999
  } else {
    l = 0.5*qchisq(1-alpha,1)
    u = 0.5*qchisq(1-alpha/(n*p0),1)
  }

  cval = uniroot(function(cval) {
    inner(n = n, tau = tau, p0 = p0, type = type, param = param, cval = cval)$Stop.prob - alpha
  }, lower = l, upper = u, extendInt = "yes")$root

  return(cval)
}
