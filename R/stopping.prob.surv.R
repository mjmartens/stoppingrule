#' @title Calculating the stopping probability given a rejection boundary (Survival Data)
#' @description
#' Internal workhouse function to calculate the stopping probability given a rejection boundary for time-to-event data
#'
#' @param bnd A list object calculated by \code{calc.bnd.surv} function
#' @param p True toxicity probability
#'
#' @return A list of three: stopping probabilities at each stage, total stopping probability, and non-stopping probabilities of each possible number of events at the last stage.
#' @export

stopping.prob.surv <- function(bnd, p){
  tau <- bnd$tau
  lambda <- -log(1 - p)/tau

  S <- bnd$Rule[,2]
  ud <- bnd$Rule[,1]
  dmin <- S[1]
  dmax <- S[length(S)]
  m <- dmax - dmin + 1

  # compute stopping probabilities
  notStop.prob <- matrix(NA, ncol = m, nrow = dmax)
  for (j in 0:nrow(notStop.prob)-1){
    if (j <= dmin - 1){notStop.prob[j+1,1] <- dpois(j, lambda = lambda*ud[1])}
    else {notStop.prob[j+1,1] <- 0}
  }

  for (k in 2:m){
    for (j in 0:nrow(notStop.prob) -1){
      Delta <- ud[which(S == dmin + k -1)] - ud[which(S == dmin + k -2)]
      if (j <= dmin + k -2){
        p <- 0
        for (i in 0:min(j, dmin + k -3)){
          temp <- notStop.prob[i+1,k-1]*dpois(j - i, lambda = lambda*Delta)
          p <- p + temp
        }
        notStop.prob[j+1, k] <- p
      }
      else {notStop.prob[j+1,k] <- 0}
    }
  }

  P <- diff(c(1, colSums(notStop.prob)))*(-1)
  return(list(stage.stop.prob=P,Stop.prob=sum(P),last.stage=notStop.prob[,m]))
}
