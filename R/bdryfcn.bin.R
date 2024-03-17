#' @title Stopping Rule Boundary Function (Binary Data)
#' @description Calculate the boundary for a given stopping rule
#'
#' @param n Maximum sample size for safety monitoring
#' @param p0 The toxicity probability under the null hypothesis
#' @param type The method used for constructing the stopping rule. Choices include a Pocock test ("Pocock"), an O'Brien-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian beta-binomial method ("BB") proposed by Geller et al. 2003, the Bayesian beta-binomial method ("CC") proposed by Chen and Chaloner 2006, a truncated SPRT ("SPRT"), and a maximized SPRT ("MaxSPRT").
#' @param cval Critical value for stopping rule method. For Wang-Tsiatis tests, this is the Delta parameter. For the Bayesian Beta-Binomial method, this is the threshold on the posterior probability. For the truncated SPRT, this is the threshold on the log likelihood ratio. For the MaxSPRT, this is the threshold on the log generalized likelihood ratio.
#' @param param A vector of the extra parameter(s) needed for certain stopping rule methods. For binomial Wang-Tsiatis tests, this is the Delta parameter. For the Geller et al. method, this is the vector of hyperparameters (a,b) for the beta prior on the toxicity probability. For Chen and Chaloner's method, this is the vector (a,b,p1,nu), containing the hyperparameters (a,b) for the beta prior on the toxicity probability, the targeted alternative toxicity probability p1, and the threshold nu for the posterior probability that the true toxicity probability p > p1. For truncated SPRT, this is the targeted alternative toxicity probability p1.
#'
#' @return A univariate function that defines the rejection boundary at any number of evaluable patients
#' @export

bdryfcn.bin = function(n,p0,type,cval,param=NULL) {
  if(type=="Pocock") {
    val = function(x) {
      if(pbeta(p0,x,1) >= cval) {return(x+0.01)}
      else{
      f = function(y) {pbeta(p0,y,x-y+1) - cval}
      return(uniroot(f,c(0.0001,x))$root)
      }
    }
    val = Vectorize(val)
  }
  else if(type=="OBF") {
    val = function(x) {
      cval * sqrt(n*p0*(1 - p0)) + x*p0
    }
  }
  else if(type=="WT") {
    val = function(x) {
      cval * (x/n)^(param[1]-0.5) * sqrt(x*p0*(1 - p0)) + x*p0
    }
  }
  else if(type=="BB") {
    val = function(x) {
      f = function(y) {pbeta(p0,param[1]+y,param[2]+x-y,lower.tail=FALSE) - cval}
      return(uniroot(f,c(0,x+param[2]))$root)
    }
    val = Vectorize(val)
  }
  else if(type=="SPRT") {
    val = function(x) {
      (cval - x*log((1-param[1])/(1-p0)))/log(param[1]*(1-p0)/((1-param[1])*p0))
    }
  }
  else if(type=="MaxSPRT") {
    val = function(x) {
      x_min = -cval/log(p0)
      if(x<x_min) {return(x+0.01)}
      else {
        f = function(y) {y*log(y/(x*p0)) + ifelse(y<x,(x-y)*log((x-y)/(x-x*p0)),0) - cval}
        return(uniroot(f,c(max(x_min,x*p0),x))$root)
      }
    }
    val = Vectorize(val)
  }

  return(val)
}
