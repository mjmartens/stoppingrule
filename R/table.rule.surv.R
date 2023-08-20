#' @title Tabulate Stopping Rule (Survival data)
#' @description Summarize a stopping rule in a condensed tabular format
#'
#' @param rule A 'rule.surv' object calculated by \code{calc.rule.surv()} function
#' @param dec Number of decimal places to which the stagewise total follow-up times should be rounded
#'
#' @return A matrix with two columns: total follow up time and their corresponding rejection boundary
#' @export
#'
#' @examples
#' pocock.rule <- calc.rule.surv(n = 30, p0 = 0.1, alpha = 0.05, tau = 100, type = "Pocock")
#' table.rule.surv(pocock.rule)

table.rule.surv <- function(rule,dec=0){
  rule <- round(rule$Rule,digits=dec)
  TFT <- rep(0, nrow(rule))
  TFT[1] <- paste(0,"-",rule[1,1])
  for (k in 2:nrow(rule)){
    TFT[k] <- paste(rule[k-1,1]+1,"-",rule[k,1])
  }
  val <- cbind(TFT, rule[,2])
  colnames(val) <- c("Total follow up time","Reject bdry")
  return(val)
}
