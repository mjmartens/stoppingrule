#' @title Stopping Rule Calculation (TITE method)
#' @description Calculate a stopping rule for safety monitoring of time-to-event data using the TITE approach.
#'
#' @param n Maximum sample size
#' @param p0 The toxicity probability under the null hypothesis
#' @param alpha The desired type I error / false positive rate for the stopping rule
#' @param type The method used for constructing the TITE stopping rule. Choices include a Pocock test ("Pocock"), an O'Brien-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian beta-binomial method ("BB") proposed by Geller et al. 2003, a truncated SPRT ("SPRT"), and a maximized SPRT ("MaxSPRT").
#' @param param A vector of the extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For the Geller et al. method, this is the vector of hyperparameters (a,b) for the beta prior on the toxicity probability. For truncated SPRT, this is the targeted alternative toxicity probability p1.
#' @param iter The number of iterations used to search for the boundary
#'
#' @return A rule.tite object, which is a list with the following elements: Rule, a two-column matrix with the effective sample sizes and their corresponding rejection boundaries; \code{n}; \code{p0}; \code{alpha}; \code{type}; \code{param}; and cval, the boundary parameter for the rule
#' @export
#'
#' @references Geller, N.L., Follman, D., Leifer, E.S. and Carter, S.L. (2003). Design of early trials in stem cell transplantation: a hybrid frequentist-Bayesian approach. \emph{Advances in Clinical Trial Biostatistics}.
#' @references Goldman, A.I. (1987). Issues in designing sequential stopping rules for monitoring side effects in clinical trials. \emph{Controlled Clinical Trials}  \strong{8(4)}, 327-37.
#' @references Ivanova, A., Qaqish, B.F. and Schell, M.J. (2005). Continuous toxicity monitoring in phase II trials in oncology. \emph{Biometrics} \strong{61(2)}, 540-545.
#' @references Kulldorff, M., Davis, R.L., Kolczak, M., Lewis, E., Lieu, T. and Platt, R. (2011). A maximized sequential probability ratio test for drug and vaccine safety surveillance. \emph{Sequential Analysis} \strong{30(1)}, 58-78.
#' @references Martens, M.J. and Logan, B.R. (2024). Statistical Rules for Safety Monitoring in Clinical Trials. \emph{Clinical Trials} \strong{21(2)}, 152-161.
#' @references Wang, S.K. and Tsiatis, A.A. (1987). Approximately optimal one-parameter boundaries for group sequential trials. \emph{Biometrics} \strong{193-199}.
#'
#' @examples
#'\dontrun{# Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
# calc.rule.tite(n=50, p0=0.2, alpha = 0.10, type = "Pocock")
#
# # Wang-Tsiatis test with Delta = 0.25 in 50 patient cohort at 10% level, expected toxicity probability of 20%
# calc.rule.tite(n=50,p0=0.20,alpha=0.10,type="WT",param=0.25)
#
# # O'Brien-Fleming test with p0=0.2 in 50 patients cohort at 10% level.
# calc.rule.tite(n=50, p0=0.2,alpha=0.1, type = "OBF")
#
# # Beta-binomial test of Geller et al. 2003 with hyperparameters (1, 9) in 100
# # patient cohort at 5% level, expected toxicity probability of 10%
# calc.rule.tite(n=100,p0=0.10,alpha=0.05,type="BB",param=c(1,9))
#
# # Truncated SPRT test with p0=0.1 and targeted alternative toxicity rate of 0.3 in
# # 100 patients cohort at 5% level.
# calc.rule.tite(n=100,p0=0.10,alpha=0.05,type="SPRT",param=0.3)
#
# # Max SPRT test with p0=0.1 in a 50 patients cohort at 10% level
# calc.rule.tite(n=50, p0=0.1, alpha = 0.1, type = "MaxSPRT")
#'}


calc.rule.tite <- function(n, p0, alpha, type, param = NULL, iter = 50){
  rule <- calc.rule.bin(ns = 1:n, p0 = p0, alpha = alpha, type = type, param = param, iter = iter)

  # Get the minimum and maximum toxicity
  b.min <- min(which(rule$Rule[,1] >= rule$Rule[,2]))
  b.max <- max(rule$Rule[,2]) -1

  # Calculate the effective sample size (ess) for each integer valued toxicity
  bdry.inverse <- bdryfcn.bin.inverse(b = b.min, n = n, p0=p0, type = type, cval = rule$cval, param = param)
  ess = bdry.inverse(b.min:b.max)

  tox <- c(b.min:(b.max+1))
  ess <- c(ess, n)

  tab = cbind(ess, tox)
  colnames(tab) = c("Effective Sample Size","Toxicity")
  val = list(Rule=tab, n=n,p0=p0,type=type,alpha=alpha,param=param,cval=rule$cval)

  class(val) = 'rule.tite'
  return(val)
}
