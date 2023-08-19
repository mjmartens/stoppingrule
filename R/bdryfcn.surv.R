#' @title Stopping Rule Boundary Function (Survival Data)
#' @description Calculate the boundary for a given stopping rule
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The probability of a toxicity occurring in \code{tau} units of time under the null hypothesis
#' @param cval Critical for the stopping rule. For Wang-Tsiatis tests, this is the Delta parameter. For the Bayesian Gamma-Poisson method, this is the threshold on the posterior probability. For the truncated SPRT, this is the threshold on the log likelihood ratio. For the MaxSPRT, this is the threshold on the log generalized likelihood ratio.
#' @param type The method used for constructing the stopping rule. Choices including a Pocock test ("Pocock"),
#' a O'Brein-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian Gamma-Poisson method ("Bayesian"),
#' a modified sequential probability ratio test ("SPRT"), and a maximized SPRT ("MaxSPRT")
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity probability p1. For Bayesian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity event rate.
#'
#' @return A univariate function that defines the rejection boundary at any number of evaluable patients
#' @export

bdryfcn.surv = function(n,p0,cval,tau,type,param=NULL) {
  lambda0 <- -log(1 - p0)/tau
  Umax <- n*tau
  if(type=="Pocock") {
    val = function(x) {
      lambda0*x + cval*sqrt(lambda0 * x)
    }
  }
  else if(type=="OBF") {
    val = function(x) {
      lambda0*x + cval*sqrt(lambda0*Umax)
    }
  }
  else if(type=="WT") {
    val = function(x) {
      lambda0*x + cval*sqrt(lambda0)*Umax^(0.5 - param)*x^param
    }
  }
  else if(type=="Bayesian") {
    val = function(x) {
      f = function(y) {pgamma(lambda0,param[1]+y,param[2]+x-y,lower.tail=FALSE) - cval}
      return(uniroot(f,c(0,x+param[2]))$root)
    }
    val = Vectorize(val)
  }
  else if(type=="SPRT") {
    p1 <- param
    lambda1 <- -log(1 - p1)/tau
    val = function(x) {
      (cval + (lambda1 - lambda0)*x)/log(lambda1/lambda0)
    }
  }
  else if(type=="MaxSPRT") {
    val = function(x) {
      if (x==0){return(0)}
      else {return((cval - lambda0*x)*(lambertWp(1/exp(1)*(cval/(lambda0*x) - 1)))^(-1))}
    }
    val = Vectorize(val)
  }

  return(val)
}
