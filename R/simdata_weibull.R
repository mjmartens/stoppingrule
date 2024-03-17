#' @title Simulate survival data for safety monitoring under Weibull distribution
#' @description
#' Internal function to simulate survival data from Weibull distribution for evaluating safety
#' monitoring rules. A random sample of size \code{n} is generated from a Weibull
#' distribution with shape parameter \code{s} to attain a toxicity rate of \code{p} at
#' survival time \code{tau}. Enrollment times are also simulated over an accrual period
#' of duration \code{A} under a uniform (0,\code{A}) distribution.
#'
#' @param n Maximum sample size for safety monitoring
#' @param p The probability of a toxicity occurring in \code{tau} units of time under the null hypothesis
#' @param tau Length of observation period
#' @param A Length of accrual period
#' @param s Shape parameter for the Weibull distribution; default value is 1
#'          (exponential distribution)
#'
#' @return A matrix with two columns: patient enrollment time and event time

simdata_weibull = function(n,p,tau,A,s=1){
  lambda = -log(1-p)/tau^s    # rate parameter
  a = runif(n,0,A)            # calendar time of enrollment
  u = runif(n)
  t = (-log(u)/lambda)^(1/s)  # toxicity event time
  return(data.frame(enrtime=a,survtime=t))
}
