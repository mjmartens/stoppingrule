#' @title Operating Characteristics Function (Survival Data)
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
#'
#' @param rule A \code{rule.surv} object calculated by \code{calc.rule.surv()} function
#' @param p The toxicity probability
#' @param MC Number of Monte Carlo replicates to simulate for estimating operating characteristics. If \code{MC} = 0, a Poisson process assumption on the event process is used to compute operating characteristics.
#' @param A Length of the enrollment period. Only required if \code{MC} > 0.
#' @param s Shape parameter for the Weibull distribution used to simulate event times. Only required if \code{MC} > 0.
#'
#' @return A list containing the rejection probability \code{p}, and the corresponding
#' rejection probability and number of events. If \code{MC} is not NULL, the expected
#' number of enrolled patients and total follow up time are also included.

opchars.surv = function(rule,p,MC,A,s=1) {
  if(MC==0) {
  bnd = list(tau = rule$tau, Rule = rule$Rule)
  probs <- stopping.prob.surv(bnd=bnd,p=p)
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
  val = list(p=p,power=power,ED=ED)
  }
  else {
    sims = simtrials.surv(rule,p,MC,A,s)
    power = mean(sims$stopped)
    ED = mean(sims$n.Toxicity)
    EN = mean(sims$n.Enrolled)
    EFU = mean(sims$Calendar.Time)
    val = list(p=p,power=power,ED=ED,EN=EN,EFU=EFU)
  }
  return(val)
}
