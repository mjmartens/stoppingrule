#' @title Tabulate Stopping Rule (Survival data)
#' @description Summarize a stopping rule in a condensed tabular format
#'
#' @param rule A \code{rule.surv} object calculated by \code{calc.rule.surv()} function
#' @param dec Number of decimal places to which the stagewise total follow-up times should be rounded
#'
#' @return A matrix with two columns: total follow up time and their corresponding rejection boundary
#' @export
#'
#' @examples
#'\dontrun{gp_rule = calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="GP",tau=100,param=c(1,1000))
#' table.rule.surv(gp_rule,2)
#'}

table.rule.surv <- function(rule,dec=0){
  rule <- round(rule$Rule,digits=dec)
  TFT <- rep(0, nrow(rule))
  TFT[1] <- paste(0,"-",sprintf(paste0("%.",dec,"f"),rule[1,1]))
  for (k in 2:nrow(rule)){
    TFT[k] <- paste(sprintf(paste0("%.",dec,"f"),rule[k-1,1]+10^-dec),"-",
                    sprintf(paste0("%.",dec,"f"),rule[k,1]))
  }
  val <- cbind(TFT, rule[,2])
  colnames(val) <- c("Total Follow Up Time","Reject Bdry")
  return(val)
}
