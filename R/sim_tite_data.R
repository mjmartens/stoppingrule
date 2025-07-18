#' @title Simulate time-to-event data for safety monitoring under Weibull distribution or power family (TITE Method)
#'
#' @describeIn Internal function to simulate time-to-event data from either Weibull or power family distribution for evaluating safety
#' monitoring rules. A random sample of size \code{n} is generated from a Weibull
#' distribution with shape parameter \code{s} or power family distribution with power parameter \code{s} to attain a toxicity rate matching the first element of \code{p} at
#' event time \code{tau}. Enrollment times are also simulated over an accrual period
#' of duration \code{A} under a uniform (0,\code{A}) distribution.
#'
#' @param n Sample size
#' @param p Vector of cumulative incidence probabilities at time \code{tau} for all causes (length 1 for survival, 2 for competing risks)
#' @param tau Length of observation period
#' @param A Length of accrual period
#' @param family Event time distribution, choices including Weibull distribution ('weibull'), or power family ('power')
#' @param s Shape parameter for Weibull distribution or power parameter for power family
#'
#' @returns A matrix with three columns: patient enrollment time, event time, and event indicator

sim_tite_data = function(n, p, tau, A, family, s) {

  # Simulate recruitment time
  a = runif(n, 0, A)

  # Simulate event time
  if (length(p) == 1) {p = c(p,0)}
  else if (length(p)!=2){stop("Error: p must be length 1 or 2.")}

  if(family=="weibull") {
    lambda = -log(1-p)/tau^s  # Rates that produce target cumulative incidences

    ## Simulate data from truncated distribution given event family
    e = sample(0:2,size=n,replace=TRUE,prob=c(1-sum(p),p))    # Determine event/censoring indicator
    t = ifelse(e>0,(-log(1-p[e]*runif(sum(e>0)))/lambda[e])^(1/s),tau)  # Compute event/censoring time
  } else if (family == "power") {
    ## Simulate data from truncated distribution given event family
    e = sample(0:2, size=n, replace=TRUE, prob=c(1-sum(p), p))  # Determine event/censoring indicator
    t = ifelse(e>0, tau*runif(sum(e>0), 0, 1)^(1/s), tau + 1e-6)       # Compute event/censoring time
  }

  return(data.frame(enrtime = a ,survtime=t, delta = e))
}
