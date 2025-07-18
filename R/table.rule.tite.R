#' @title Tabulate Stopping Rule (TITE Method)
#' @description Summarize a stopping rule in a condensed tabular format
#'
#' @param x A \code{rule.tite} object calculated by \code{calc.rule.tite()} function
#' @param dec Number of decimal places to which the stagewise effective sample sizes should be rounded
#'
#' @return A matrix with two columns: the ranges of effective sample sizes, and corresponding rejection boundaries for these ranges
#' @export
#'
#' @examples
#'
#'\dontrun{
#'# Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' poc_rule = calc.rule.tite(n=50,p0=0.20,alpha=0.10,type="Pocock")
#'
#' # Tabulate stopping boundary
#' table.rule.tite(poc_rule)
#'}

table.rule.tite = function(x,dec=3) {
  rule = x$Rule
  n = x$n
  lower = c(min(rule[,2]), head(rule[,1], -1))
  upper = rule[,1]
  ess_range = NULL

  for (k in 1:(nrow(rule))){
    if(k==1) {
    ess_range[k] <- paste(sprintf(paste0("%.",dec,"f"),lower[k]),"-",
                    sprintf(paste0("%.",dec,"f"),upper[k]))
    }
    else if(k==nrow(rule)){
      ess_range[k] <- paste(sprintf(paste0("%.",dec,"f"),lower[k]+10^-dec),"-",
                            sprintf(paste0("%.",dec,"f"),n))

    }
    else {
      ess_range[k] <- paste(sprintf(paste0("%.",dec,"f"),lower[k]+10^-dec),"-",
                                sprintf(paste0("%.",dec,"f"),upper[k]))
    }
  }
  val = cbind(ess_range, rule[,2])
  colnames(val) = c("Effective Sample Size","Toxicity")

  return(val)
}
