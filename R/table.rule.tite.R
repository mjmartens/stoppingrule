#' @title Tabulate Stopping Rule (TITE Method)
#' @description Summarize a stopping rule in a condensed tabular format
#'
#' @param x A \code{rule.tite} object calculated by \code{calc.rule.tite()} function
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

table.rule.tite = function(x) {
  rule = x$Rule
  n = max(rule[,1])
  lower = c(min(rule[,2]), head(rule[,1], -1))
  upper = rule[,1]
  ess_range = sprintf("%.4f - %.4f", lower, upper)

  # If the last range doesn't end at n, update it to go to n
  if (max(upper) < n) {
    ess_range[length(ess_range)] = sprintf("%.4f - %.4f", upper[length(upper) - 1], n)
  }

  val = cbind(ess_range, rule[,2])
  colnames(val) = c("Effective Sample Size","Toxicity")

  return(val)
}
