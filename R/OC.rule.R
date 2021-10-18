#' @title Operating Characteristics Function
#' @description Compute operating characteristics for a stopping rule at a set of toxicity rates. Characteristics calculated include the overall rejection probability, the expected number of patients evaluated, and the expected number of events.
#'
#' @param rule A matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param ps A vector of toxicity rates at which the operating characteristics will be computed
#'
#' @return A matrix with four columns: the toxicity rates \code{ps}, the corresponding rejection probabilities, the corresponding expected numbers of evaluated patients, and the corresponding expected numbers of events
#' @export
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' poc = calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
#'
#' # Compute operating characteristics at toxicity rates of 20%, 25%, 30%, 35%, and 40%
#' OC.rule(rule=poc,ps=seq(0.2,0.4,0.05))
OC.rule = function(rule,ps) {
  ns = rule[,1]
  bs = rule[,2]
  k = length(ns)
  tab = matrix(0,nrow=length(ps),ncol=4)
  for(i in 1:length(ps)) {
    op = opchars(rule,ps[i])
    tab[i,] = c(ps[i],op$power,op$Eeval,op$ED)
  }
  colnames(tab) = c("p","Reject Prob","E(evalauted)","E(events)")
  return(tab)
}
