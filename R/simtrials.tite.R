#' @title Simulate trials with safety monitoring by the TITE method
#' @description
#' Internal workhorse function used to simulate trials with safety monitoring by
#' the TITE method. The provided stopping \code{rule} is used for
#' monitoring of \code{MC} simulated trials. For each trial, a random sample is
#' generated from either a Weibull distribution with shape parameter \code{s} or
#' power distribution with power parameter to attain a toxicity rate of
#' \code{p}. Enrollment times are simulated over an accrual
#' period of duration \code{A} under a uniform (0,\code{A}) distribution.
#'
#' @param rule A \code{rule.tite} object with the safety stopping rule for evaluation
#' @param p Vector of cumulative incidence probabilities at time \code{tau} for all causes (length 1 for survival, 2 for competing risks)
#' @param tau Length of observation period
#' @param MC Number of Monte Carlo replicated datasets to simulate
#' @param A Length of accrual period
#' @param family Event time distribution, choices including Weibull distribution ('weibull') and power family ('power')
#' @param s Shape parameter for Weibull distribution or power parameter for power family
#'
#'
#' @return A matrix with \code{MC} rows and 14 columns, one row per simulated
#'         trial. Columns include the stopping rule type and design parameters,
#'         the numbers of events and enrolled patients, the total follow-up time
#'         in the cohort, the calendar time when the study ends, the reject/no reject
#'         decision, and the last stage of monitoring reached when the study ends.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(2025)
# Pocock = calc.rule.tite(n=50, p0=0.20, alpha = 0.05, type = "Pocock")
# simdata = simtrials.tite(rule = Pocock, p = 0.25, tau = 100, MC = 1000, A = 365, family = "power", s=1)
# simdata.compt= simtrials.tite(rule = Pocock, p = c(0.2, 0.1), tau = 100, MC = 1000, A = 365, family = "power", s=1)
#'}
simtrials.tite = function (rule, p, tau, MC, A, family, s) {

  p0 = rule$p0
  alpha = rule$alpha
  n = rule$n
  alldata = sim_tite_data(n = n * MC, p = p, tau = tau, A = A, family = family, s = s)

  # alldata$d = 1 * I(alldata$survtime <= tau) # Indicator for an event before tau
  alldata$d = 1 * I(alldata$survtime <= tau)*I(alldata$delta == 1) # Indicator for an event of interest before tau
  alldata$t = pmin(alldata$survtime, tau) # Minimum time between tau and event time
  observed = matrix(NA, n, 4 * MC)
  for (m in 1:MC) {
    # m <- 1
    data = alldata[n * (m - 1) + 1:n, ]
    monitoring.times = sort(data$enrtime[data$d == 1] + data$t[data$d == 1])
    n_events = length(monitoring.times)
    if (n_events > 0) {
      for (i in 1:n_events) {
        # i = 10

        if (length(p) > 1){
          event.of.interest = 1*(monitoring.times[i] >= data$enrtime + data$t)*(data$d == 1)
          # Number of events of interest happened on or before this monitoring time

          event.competing = 1*(monitoring.times[i] >= data$enrtime + data$t)*(data$delta == 2)
          other.ess = (pmax(pmin(monitoring.times[i] - data$enrtime, tau),0)/tau)[(event.of.interest + event.competing) == 0]
          # Partial numbers for pending evaluation or not yet entered patients

          observed[i, 4 * (m - 1) + 1:4] = c(sum(other.ess, event.of.interest, event.competing),
                                             i, monitoring.times[i], sum(other.ess > 0) + sum(event.of.interest, event.competing))

        } else {
          event = 1*(monitoring.times[i] >= data$enrtime + data$t)*(data$d == 1)
          # Number of events of interest happened on or before this monitoring time

          other.ess = (pmax(pmin(monitoring.times[i] - data$enrtime, tau),0)/tau)[event == 0]
          # Partial numbers for pending evaluation or not yet entered patients

          observed[i, 4 * (m - 1) + 1:4] = c(sum(other.ess, event), i, monitoring.times[i], sum(other.ess > 0) + sum(event))
        }


      }
    }
    observed[n, 4 * (m - 1) + 1:4] = c(n, n_events, A + tau, n)
  }
  colnames(observed) <- rep(c("Total.ess", "n.Toxicity",
                              "Calendar.Time", "n.Enrolled"), MC)
  rm(alldata)
  monitor = data.frame(matrix(NA, MC, 15))
  Rule = rule$Rule
  bmin = min(Rule[, 2])
  bmax = max(Rule[, 2])
  if (bmin > 1) {
    Rule = rbind(matrix(c(rep(NA, bmin - 1), 1:(bmin - 1)),
                        ncol = 2), Rule, matrix(c(rep(NA, n - bmax), (bmax + 1):n), ncol = 2))
  }
  else {
    Rule = rbind(Rule, matrix(c(rep(NA, n - bmax), (bmax + 1):n), ncol = 2))
  }
  # colnames(Rule) = c("Total.followup", "Reject bdry")
  colnames(Rule) = c("Effective Sample Size", "Reject bdry")
  # decisions = (Rule[, 1] >= observed[, 4 * 1:MC - 3])
  decisions = (Rule[, 1] >= observed[, 4 * 1:MC - 3])
  rejects = colMaxs(decisions + 0, na.rm = TRUE)
  rejects[rejects == -Inf] = 0
  f = function(x) {
    val = which(x == TRUE)[1]
    if (length(val) == 0 | is.na(val)) {
      val = bmax
    }
    return(val)
  }
  laststage = apply(decisions, 2, f)
  indices = decisions # Flag for an early stopping happened
  indices[which(indices == 0)] = Inf
  indices[n, ] = 1
  Calendar.time = colMins(indices * observed[, 4 * 1:MC - 1],
                          na.rm = TRUE)
  n.Toxicity = colMins(indices * observed[, 4 * 1:MC - 2],
                       na.rm = TRUE)
  # Total.followup = colMins(indices * observed[, 4 * 1:MC - 3], na.rm = TRUE)
  Total.ess = colMins(indices * observed[, 4 * 1:MC - 3], na.rm = TRUE)
  # n.Enrolled = colMins(indices * matrix(1:n, nrow = n, ncol = MC),
  #                      na.rm = TRUE)
  n.Enrolled = colMins(indices * observed[, 4 * 1:MC - 0],
                       na.rm = TRUE)
  if (length(p)==1){
    p = c(p,0)
  }
  monitor[1:MC, ] = cbind(n = n, p0 = p0, p = p[1], pc = p[2],
                          tau = tau,
                          A = A, alpha = alpha, type = rule$type, cval = rule$cval,
                          # Total.followup = Total.followup,
                          Total.ess = Total.ess,
                          n.Toxicity = n.Toxicity,
                          Calendar.Time = Calendar.time, n.Enrolled = n.Enrolled,
                          rejects, laststage)
  colnames(monitor) <- c("n", "p0", "p","p.compt" ,"tau", "A", "alpha",
                         "rule", "cval", "Total.ess", "n.Toxicity", "Calendar.Time",
                         "n.Enrolled", "stopped", "laststage")
  idx = c(1:7, 9:15)
  monitor[, idx] = apply(monitor[, idx], 2, function(x) as.numeric(as.character(x)))
  return(monitor)
}
