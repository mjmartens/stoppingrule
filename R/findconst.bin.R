#' @title Search for Calibration Value (Binary Data)
#' @description Internal workhorse function to calculate the calibration constant value that attains level alpha for given method
#'
#' @param ns A vector of sample sizes at which sequential testing is performed
#' @param p0 The toxicity probability under the null hypothesis
#' @param alpha The desired type I error / false positive rate for the stopping rule
#' @param type The method used for constructing the stopping rule
#' @param l Lower starting value of bracket for calibration constant
#' @param u Upper starting value of bracket for calibration constant
#' @param iter The number of iterations used to search for the boundary
#' @param param A vector of the extra parameter(s) needed for certain stopping rule methods. For binomial Wang-Tsiatis tests, this is the Delta parameter. For the Geller et al. method, this is the vector of hyperparameters (a,b) for the beta prior on the toxicity probability. For Chen and Chaloner's method, this is the vector (a,b,p1,nu), containing the hyperparameters (a,b) for the beta prior on the toxicity probability, the targeted alternative toxicity probability p1, and the threshold nu for the posterior probability that the true toxicity probability p > p1. For truncated SPRT, this is the targeted alternative toxicity probability p1.
#'
#' @return The calibration constant used for subsequent stopping boundary calculation

findconst.bin = function(ns,p0,alpha,type,l,u,iter=50,param) {
  k = length(ns)
  n = tail(ns,1)
  low = l
  upp = u
  for(i in 1:iter) {
    const = upp
    c_curr = (low+upp)/2
    bs = calc.bnd.bin(n,p0,type,c_curr,param)
    bs = bs[ns]
    typeIrate = sum(opchars.bin(list(Rule=cbind(ns,bs)),p0)$power)
    if(typeIrate <= alpha) {
      const = c_curr
      upp = c_curr
    }
    else {low = c_curr}
  }
  return(const)
}
