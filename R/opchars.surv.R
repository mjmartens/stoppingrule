#' @title Operating Characteristics Function (Survival Data)
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
#'
#' @param rule A 'rule.surv' object calculated by \code{calc.rule.surv()} function
#' @param p The toxicity probability
#'
#' @return A list with the following elements: p, the corresponding rejection probability, and the corresponding expected total follow up time and number of events at the point of stopping/study end

opchars.surv = function(rule,p) {
  bnd = list(tau = rule$tau, S = rule$Rule[,2], ud = rule$Rule[,1])
  probs <- stopping.prob(bnd=bnd,p=p)
  power <- probs$Stop.prob
  EFU <- sum(probs$stage.stop.prob*rule$Rule[,1]) +
    (1-power)*rule$Rule[nrow(rule$Rule),1]
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
  ED <- s + t

  return(list(p=p,power=power,EFU=EFU,ED=ED))
}
