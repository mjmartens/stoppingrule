#' @title Search for Calibration Value
#' @description Internal workhorse function to calculate the calibration constant value that attains level alpha for given method
#'
#' @param ns A vector of sample sizes at which sequential testing is performed
#' @param p0 The toxicity rate under the null hypothesis
#' @param type The method used for constructing the stopping rule
#' @param alpha The desired type I error / false positive rate for the stopping rule
#' @param l Lower starting value of bracket for calibration constant
#' @param u Upper starting value of bracket for calibration constant
#' @param iter The number of iterations used to search for the boundary
#' @param param Extra parameter(s) needed for certain stopping rule methods. For binomial Wang-Tsiatis tests, this is the Delta parameter. For Bayesian beta-binomial model, this is the pair of hyperparameters for the beta prior on the toxicity rate. For modified SPRT, this is the targeted alternative toxicity rate p1.
#'
#' @return The calibration constant used for subsequent stopping boundary calculation
findconst = function(ns,p0,type,alpha,l,u,iter=50,param) {
  k = length(ns)
  n = tail(ns,1)
  low = l
  upp = u
  for(i in 1:iter) {
    const = upp
    c_curr = (low+upp)/2
    bs = calc.bnd(n,p0,c_curr,type,param)
    bs = bs[ns]
    typeIrate = sum(opchars(cbind(ns,bs),p0)$reject.prob)
    if(typeIrate <= alpha) {
      const = c_curr
      upp = c_curr
    }
    else {low = c_curr}
  }
  return(const)
}
