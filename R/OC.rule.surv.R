#' @title Operating Characteristics Function (Survival Data)
#' @description
#' Compute operating characteristics for a stopping rule at a set of toxicity rates.
#' Characteristics calculated include the overall rejection probability, the expected
#' number of patients evaluated, and the expected number of events for time-to-event data.
#'
#'
#' @param rule A \code{rule.surv} object calculated by \code{calc.rule.surv()} function
#' @param ps A vector of toxicity probabilities at which the operating characteristics will be computed
#' @param MC Number of Monte Carlo replicates to simulate for estimating operating characteristics. If \code{MC} = 0, a Poisson process assumption on the event process is used to compute operating characteristics.
#' @param A Length of the enrollment period. Only required if \code{MC} > 0.
#' @param s Shape parameter for the Weibull distribution used to simulate event times. Default is \code{s} = 1 (exponential). Only required if \code{MC} > 0.
#'
#' @return A matrix with columns containing the toxicity probabilities \code{ps},
#' the corresponding rejection probabilities, and the corresponding expected number
#' of events. If \code{MC} is not NULL, the expected number of enrolled patients and total
#' follow up time are also included.
#'
#' @details
#' Operating characteristics are generated either by Monte Carlo estimation or computed
#' directly under a Poisson process assumption for the event process over time. The
#' Monte Carlo approach assumes a random uniform accrual distribution and a Weibull
#' event time distribution with distribution function \eqn{exp(-\lambda * t^s)}, so
#' it requires specification of the enrollment period length and shape parameter
#' of the event distribution.
#' @export
#'
#' @examples
#' poc_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="Pocock",tau=100)
#' gp_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="GP",tau=60,param=c(1,1000))
#' OC.rule.surv(rule=poc_rule,ps=seq(0.2,0.4,0.05),MC=0)
#' OC.rule.surv(rule=gp_rule,ps=seq(0.2,0.4,0.05),MC=0)
#'
#' set.seed(82426499)
#' ps = seq(0.15,0.35,0.05)
#' wt_rule = calc.rule.surv(n=46,p0=0.15,alpha=0.10,type="WT",tau=100,param=0.25)
#' OC.rule.surv(rule=wt_rule,ps=ps,MC=1000,A=1095)
#'
#' p1h = 0.3418071
#' sp_rule = calc.rule.surv(n=46,p0=0.15,alpha=0.10,type="SPRT",tau=100,param=p1h)
#' OC.rule.surv(rule=sp_rule,ps=ps,MC=1000,A=1095)
#'
#' gp_rule = calc.rule.surv(n=46,p0=0.15,alpha=0.10,type="GP",tau=100,
#'                           param=11.5*c(-log(1-0.15),100))
#' OC.rule.surv(rule=gp_rule,ps=ps,MC=1000,A=1095)
#'

OC.rule.surv = function(rule,ps,MC,A,s=1){
  if(MC>=0) {
    tab = matrix(0,nrow=length(ps),ncol=5-2*I(MC==0))
    for(i in 1:length(ps)) {
    op = opchars.surv(rule,ps[i],MC,A,s=1)
    if(MC==0) {
      tab[i,] = c(ps[i],op$power,op$ED)
      colnames(tab) = c('p',"Reject Prob", "E(Events)")
    }
    else if (MC > 0) {
      tab[i,] = c(ps[i],op$power,op$ED,op$EN,op$EFU)
      colnames(tab) = c('p',"Reject Prob", "E(Events)","E(Enrolled)","E(Total Follow up time)")
    }
    }
    return(tab)
  }
  else {
    print("Error: MC must be a nonnegative number")
  }
}
