#' @title Operating Characteristics Function (Binary Data)
#' @description Compute operating characteristics for a stopping rule at a set of toxicity rates. Characteristics calculated include the overall rejection probability, the expected number of patients evaluated, and the expected number of events.
#'
#' @param rule A 'rule.bin' object calculated by \code{calc.rule.bin()} function
#' @param ps A vector of toxicity probabilities at which the operating characteristics will be computed
#'
#' @return A matrix with four columns: the toxicity probabilities \code{ps}, the corresponding rejection probabilities, the corresponding expected numbers of evaluated patients and events at the point of stopping/study end
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc = calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="Pocock")
#'
#' # Compute operating characteristics at toxicity probabilities of 20%, 25%, 30%, 35%, and 40%
#' OC.rule.bin(rule=poc,ps=seq(0.2,0.4,0.05))

OC.rule.bin = function(rule,ps) {
  ns = rule$Rule[,1]
  bs = rule$Rule[,2]
  k = length(ns)
  tab = matrix(0,nrow=length(ps),ncol=4)
  for(i in 1:length(ps)) {
    op = opchars.bin(rule,ps[i])
    tab[i,] = c(ps[i],op$power,op$Eeval,op$ED)
  }
  colnames(tab) = c("p","Reject Prob","E(evalauted)","E(events)")
  return(tab)
}
