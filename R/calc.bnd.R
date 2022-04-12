#' @title Stopping Boundary Calculation
#' @description Internal workhorse function to calculate stopping boundary for a given method
#'
#' @param n Maximum sample size for safety monitoring
#' @param p0 The toxicity rate under the null hypothesis
#' @param cval Critical value for stopping rule method
#' @param type The method used for constructing the stopping rule
#' @param param Extra parameter(s) needed for certain stopping rule methods. For binomial Wang-Tsiatis tests, this is the Delta parameter. For Bayesian beta-binomial model, this is the pair of hyperparameters for the beta prior on the toxicity rate. For modified SPRT, this is the targeted alternative toxicity rate p1.
#'
#' @return A vector of stopping boundaries at the sample sizes 1, 2, ..., n
calc.bnd = function(n,p0,cval,type,param) {
  bs = NULL
  if(type=="Pocock") {
    bs = qbinom(cval,1:n,p0)+1
  }
  else if(type=="WT") {
    del = param
    bs = ceiling(cval*((1:n)/n)^(del-0.5) * sqrt((1:n)*p0*(1-p0)) + (1:n)*p0)
  }
  else if(type=="BB") {
    a = param[1]
    b = param[2]
    postprob = matrix(0,nrow=n,ncol=n+1)
    for(j in 1:n) {
      postprob[j,1:(j+1)] = pbeta(p0,a+0:j,b+j:0,lower.tail=FALSE)
      if(max(postprob[j,]) > cval) {
        bs[j] = which(postprob[j,] > cval)[1]-1
      }
      else {bs[j] = j+1}
    }
  }
  else if(type=="SPRT") {
    p1 = param
    bs = ceiling((cval - (1:n)*log((1-p1)/(1-p0)))/log(p1*(1-p0)/((1-p1)*p0)))
  }
  else if(type=="MaxSPRT") {
    GLR = matrix(1,nrow=n,ncol=n+1)
    for(j in 1:n) {
      phat = (0:j)/j
      GLR[j,1:(j+1)] = ifelse(phat >= p0,(phat/p0)^(0:j) * ((1-phat)/(1-p0))^(j:0),1)
      if(max(GLR[j,]) > exp(cval)) {
        bs[j] = which(GLR[j,] > exp(cval))[1]-1
      }
      else {bs[j] = j+1}
    }
  }
  return(bs)
}
