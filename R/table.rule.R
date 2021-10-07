#' @title Tabulate Stopping Rule
#' @description Summarize a stopping rule in a condensed tabular format
#'
#' @param rule A matrix with two columns: the vector of sample sizes at which the stopping rule is evaluated, and a corresponding vector of rejection boundaries at these sizes
#'
#' @return A matrix with two columns: first column lists ranges of evaluable patients; the second column lists the corresponding rejection boundaries for these ranges
#' @export
#'
#' @examples
#' poc_rule = calc.rule(ns=1:30,p0=0.15,type="Pocock",alpha=0.10)
#' table.rule(poc_rule)
table.rule = function(rule) {
  n = max(rule[,1])
  idx = NULL
  rule = rule[rule[,2]<=rule[,1],] # Find rows where rejection can happen
  brange = range(rule[,2])
  for(i in 1:(brange[2]-brange[1]+1)){
    idx[i] = which(rule[,2]==i+(brange[1]-1))[1]
  }
  rule = rule[idx,]

  k = nrow(rule)
  n_eval = rep(0,k)
  n_eval[-k] = paste(rule[-k,1],"-",rule[-1,1]-1)
  if((rule[2,1] - rule[1,1]) == 1) {n_eval[1] = rule[1,1]}
  if(rule[k,1] < n) {
    n_eval[k] = paste(rule[k,1],"-",n)
  }
  else {n_eval[k] = n}

  val = cbind(n_eval,rule[,2])
  colnames(val) = c("N evaluable","Reject If N >=")

  return(val)
}
