#' @title Operating Characteristics Function (TITE Method)
#' @description Compute operating characteristics for a stopping rule at a set of toxicity rates. Characteristics calculated include the overall rejection probability, the expected number of patients evaluated, and the expected number of events.
#'
#' @param rule A \code{rule.tite} object with the safety stopping rule for evaluation
#' @param ps Vector of cumulative incidence probabilities for toxicity at time \code{tau}
#' @param ps.compt Vector of cumulative incidence probabilities for competing risks at time \code{tau} if competing risks are involved. Set to NULL if no competing risks exist.
#' @param tau Length of observation period
#' @param MC Number of Monte Carlo replicated datasets to simulate
#' @param A Length of accrual period
#' @param family Event time distribution, choices including Weibull distribution ('weibull') and power family ('power')
#' @param s Shape parameter for Weibull distribution or power parameter for power family
#'
#'
#' @return A matrix with columns containing the toxicity probabilities \code{ps}, competing risk probability (0 for survival outcome), the corresponding rejection probabilities, and the corresponding expected number of events. If \code{tau} and \code{A} are also specified, the expected numbers of enrolled patients and the expected calendar time at the point of stopping/study end are also included.
#'
#' @details
#' The Weibull family has cumulative distribution functions of the form
#' \deqn{1- \exp(- \lambda t^s),}
#' where \eqn{\lambda} is the rate parameter and \code{s} is the shape parameter.
#' The power family has cumulative distribution functions of the form
#' \deqn{1- t^s}
#' where \code{s} is the power parameter.
#' If \code{tau} and \code{A} are specified, the expected number of events includes events among patients who are still pending evaluation at the time of early stopping, computed under an assumption of a random uniform accrual distribution. Otherwise, only events that occurred prior to stopping are included, as the number of events occurring in pending patients depends on \code{tau} and \code{A}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Bayesian beta-extended binomial method in 50 patient cohort at 10% level,
#' # expected toxicity probability of 20%
#' bb_rule = calc.rule.tite(n=50,p0=0.20,alpha=0.10,type="BB",param=c(2,8))
#'
#' # Compute operating characteristics at toxicity probabilities of 20%, 25%, 30%, 35%, and 40%
#' OC.rule.tite(rule=bb_rule,ps=seq(0.2,0.4,0.05), MC =1000, tau=30,A=730, family = 'weibull', s = 2)
#'}

OC.rule.tite = function (rule, ps, ps.compt=NULL, MC, tau, A, family="power", s=1) {
  if (MC >= 0) {
    ps = cbind(ps,ps.compt)
    if (is.null(ps.compt)){
      tab = matrix(0, nrow = length(ps), ncol = 6)
      for (i in 1:length(ps)){
        op = opchars.tite(rule = rule, p = ps[i,], tau = tau, MC = MC, A = A, family = family, s = s)
        tab[i, ] = c(ps[i], 0, op$Reject.prob, op$Exp.tox, op$Exp.n, op$Exp.caltime)
      }
    } else {
      tab = matrix(0, nrow = nrow(ps), ncol = 6)
      for (i in 1:nrow(ps)){
        op = opchars.tite(rule = rule, p = ps[i,], tau = tau, MC = MC, A = A, family = family, s = s)
        tab[i, ] = c(ps[i,1], ps[i,2], op$Reject.prob, op$Exp.tox, op$Exp.n, op$Exp.caltime)
      }
    }
    colnames(tab) = c("p", "p.compt","Reject Prob", "E(Events)",
                      "E(Enrolled)", "E(Calendar time)")
    return(tab)
  }
  else {
    print("Error: MC must be a nonnegative number")
  }
}
