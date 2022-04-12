#' @title Operating Characteristics Function
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and event rate
#'
#' @param rule A matrix with two columns: the sample sizes at which sequential testing is performed, and their corresponding rejection boundaries
#' @param p The event probability
#'
#' @return A list with four objects: the toxicity rate \code{p}, the corresponding rejection probability, the corresponding expected number of evaluated patients, and the corresponding expected number of events
opchars = function(rule,p) {
  ns = rule[,1]
  bs = rule[,2]
  k = length(ns)
  ms = diff(c(0,ns))
  C = matrix(0,nrow=k,ncol=max(ns)+1)
  reject.prob = rep(0,k)
  exp.events = rep(0,k)
  C[1,1:(ms[1]+1)] = dbinom(0:ms[1],ms[1],p)
  if(ms[1]>=bs[1]) {
    reject.prob[1] = sum(C[1,-1:-bs[1]])
    exp.events[1] = sum(bs[1]:ns[1]*C[1,(bs[1]+1):(ns[1]+1)])
  }
  for(j in 2:k) {
    for(y in 0:ns[j]) {
      Ek = max(0,y-ms[j])
      Fk = min(bs[j-1]-1,y)
      if(Ek <= Fk) {
        C[j,y+1] = sum(C[j-1,(Ek+1):(Fk+1)] * dbinom(y-Ek:Fk,ms[j],p))
      }
    }
    reject.prob[j] = sum(C[j,-1:-bs[j]])
    if(j<k) {exp.events[j] = sum(bs[j]:ns[j]*C[j,(bs[j]+1):(ns[j]+1)])}
    else    {exp.events[j] = sum(0:ns[j]*C[j,])}
  }
  Eeval = sum(ns[-k]*reject.prob[-k]) + ns[k]*(1-sum(reject.prob[-k]))
  ED = sum(exp.events)
  power = sum(reject.prob)

  return(list(reject.prob=reject.prob,power=power,Eeval=Eeval,ED=ED))
}
