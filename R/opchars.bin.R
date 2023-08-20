#' @title Operating Characteristics Function (Binary Data)
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
#'
#' @param rule A 'rule.bin' object calculated by \code{calc.rule.bin()} function
#' @param p The toxicity probability
#'
#' @return A list with the following elements: p, the corresponding rejection probability, and the corresponding expected numbers of evaluated patients and events at the point of stopping/study end

opchars.bin = function(rule,p) {
  tab = rule$Rule
  ns = tab[,1]
  bs = tab[,2]
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

  return(list(p=p,power=power,Eeval=Eeval,ED=ED))
}
