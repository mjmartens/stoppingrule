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
  tab = matrix(0, nrow = length(ps), ncol = 4)
  for (i in 1:length(ps)){
    tab[i,1] <- ps[i]
    bnd = list(tau = rule$tau, S = rule$Rule[,2], ud = rule$Rule[,1])
    probs <- stopping.prob(bnd = bnd, p = ps[i])
    tab[i,2] <- probs$Stop.prob
    tab[i,3] <- sum(probs$stage.stop.prob*rule$Rule[,1]) +
      (1 - tab[i,2])*rule$Rule[nrow(rule$Rule),1]
    dmin <- rule$Rule[1,2]
    dmax <- rule$Rule[nrow(rule$Rule),2]
    m <- dmax - dmin + 1
    s <- 0
    for (k in 1:m){
      q <- (dmin + k -1)*probs$stage.stop.prob[k]
      s <- s + q
    }
    t <- 0
    for (k in 0:(dmax - 1)){
      q <- k*probs$last.stage[k+1]
      t <- q + t
    }
    tab[i,4] <- s + t
  }

  colnames(tab) <- c('p',"Reject Prob","E(Total Follow up time)", "E(Events)")
  return(tab)
}
