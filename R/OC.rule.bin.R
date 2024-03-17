#' @title Operating Characteristics Function (Binary Data)
#' @description Compute operating characteristics for a stopping rule at a set of toxicity rates. Characteristics calculated include the overall rejection probability, the expected number of patients evaluated, and the expected number of events.
#'
#' @param rule A \code{rule.bin} object calculated by \code{calc.rule.bin()} function
#' @param ps A vector of toxicity probabilities at which the operating characteristics will be computed
#' @param tau Length of observation period
#' @param A Length of the enrollment period.
#'
#' @return A matrix with columns containing the toxicity probabilities \code{ps}, the corresponding rejection probabilities, and the corresponding expected number of events. If \code{tau} and \code{A} are also specified, the expected numbers of enrolled patients and the expected calendar time at the point of stopping/study end are also included.
#'
#' @details
#' If \code{tau} and \code{A} are specified, the expected number of events includes events among patients who are still pending evaluation at the time of early stopping, computed under an assumption of a random uniform accrual distribution. Otherwise, only events that occurred prior to stopping are included, as the number of events occurring in pending patients depends on \code{tau} and \code{A}.
#'
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="Pocock")
#'
#' # Bayesian beta-binomial method of Geller et al. in 50 patient cohort at 10% level,
#' # expected toxicity probability of 20%
#' bb_rule = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="BB",param=c(2,8))
#'
#' # Compute operating characteristics at toxicity probabilities of 20%, 25%, 30%, 35%, and 40%
#' OC.rule.bin(rule=poc_rule,ps=seq(0.2,0.4,0.05))
#' OC.rule.bin(rule=bb_rule,ps=seq(0.2,0.4,0.05),tau=30,A=730)

OC.rule.bin = function(rule, ps, tau = NULL, A = NULL) {
  ns = rule$Rule[,1]
  bs = rule$Rule[,2]
  k = length(ns)

  if (!is.null(A) & !is.null(tau)){
    tab = matrix(0,nrow=length(ps),ncol=5)
    for(i in 1:length(ps)) {
      op = opchars.bin(rule = rule,p = ps[i], A = A, tau = tau)
      tab[i,] = c(ps[i],op$power, op$exp.toxicities, op$exp.enrolled, op$exp.calendar)
    }
    colnames(tab) = c("p","Reject Prob","E(events)","E(Enrolled)", "E(CalendarTime)")
    return(tab)
  } else {
    tab = matrix(0,nrow=length(ps),ncol=3)
    for(i in 1:length(ps)) {
      op = opchars.bin(rule = rule,p = ps[i])
      tab[i,] = c(ps[i],op$power, op$exp.toxicities)
    }
    colnames(tab) = c("p","Reject Prob","E(events)")
    return(tab)
  }
}
