#' @title Operating Characteristics Function (TITE Method)
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
#'
#' @param rule A \code{rule.tite} object with the safety stopping rule for evaluation
#' @param p Vector of cumulative incidence probabilities at time \code{tau} for all causes (length 1 for survival, 2 for competing risks)
#' @param tau Length of observation period
#' @param MC Number of Monte Carlo replicated datasets to simulate
#' @param A Length of accrual period
#' @param family Event time distribution, choices including Weibull distribution ('weibull') and power family ('power')
#' @param s Shape parameter for Weibull distribution or power parameter for power family
#'
#' @return A list containing the toxicity probability \code{p}, and the corresponding rejection probability and expected number of events. If \code{tau} and \code{A} are also specified, the expected number of enrolled patients and the expected calendar time at the point of stopping/study end are also included.

opchars.tite = function (rule, p, tau, MC, A, family, s) {
  n = rule$n
  sims = simtrials.tite(rule = rule, p = p, tau = tau, MC = MC, A=A, family = family, s=s)
  Reject.prob = mean(sims$stopped)
  Exp.tox = mean(sims$n.Toxicity)
  Exp.n = mean(sims$n.Enrolled)
  Exp.caltime = mean(sims$Calendar.Time)
  val = list(p = p, Reject.prob = Reject.prob, Exp.tox = Exp.tox,
             Exp.n = Exp.n, Exp.caltime = Exp.caltime)
  return(val)
}
