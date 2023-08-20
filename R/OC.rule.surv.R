#' @title Operating Characteristics Function (Survival Data)
#' @description
#' Compute operating characteristics for a stopping rule at a set of toxicity rates.
#' Characteristics calculated include the overall rejection probability, the expected
#' number of patients evaluated, and the expected number of events for time-to-event data
#'
#'
#' @param rule A 'rule.surv' object calculated by \code{calc.rule.surv()} function
#' @param ps A vector of toxicity probabilities at which the operating characteristics will be computed
#'
#' @return A matrix with four columns: the toxicity probability \code{ps}, the corresponding
#' rejection probabilities, and the corresponding expected total follow up time and number of events at the point of stopping/study end
#' @export
#'
#' @examples
#' pocock.rule <- calc.rule.surv(n = 30, p0 = 0.1, alpha = 0.05, tau = 100, type = "Pocock")
#' OC.rule.surv(rule = pocock.rule, ps = seq(0.1, 0.5, 0.1))

OC.rule.surv <- function(rule, ps){
  tab = matrix(0,nrow=length(ps),ncol=4)
  for(i in 1:length(ps)) {
    op = opchars.surv(rule,ps[i])
    tab[i,] = c(ps[i],op$power,op$EFU,op$ED)
  }

  colnames(tab) <- c('p',"Reject Prob","E(Total Follow up time)", "E(Events)")
  return(tab)
}
