#' @title Stopping Rule Calculation (Survival Data)
#' @description Calculate a stopping rule for safety monitoring for time-to-event data
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Length of observation period
#' @param p0 The probability of a toxicity occurring in \code{tau} units of time under the null hypothesis
#' @param alpha The nominal type I error/false positive rate for the stopping rule,
#' under an assumption that the cumulative number of events follows a Poisson process
#' over the study duration.
#' @param type The method used for constructing the stopping rule. Choices including a Pocock test ("Pocock"),
#' a O'Brein-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian Gamma-Poisson method ("GP"),
#' a truncated sequential probability ratio test ("SPRT"), and a maximized SPRT ("MaxSPRT")
#' @param maxInf Specification of the maximum information (maximum exposure time) used for designing the
#' stopping rule. Options include the expected exposure time for n patients used H0 ("expected") and the
#' maximum possible exposure time ("maximum"). Default is "expected" (expected exposure time in cohort).
#' @param param A vector of the extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For truncated SPRT, this is the targeted alternative toxicity probability p1. For Bayesian Gamma-Poisson model, this is the vector of hyperparameters (shape,rate) for the gamma prior on the toxicity event rate.
#'
#' @return A rule.surv object, which is a list with the following elements: Rule, a two-column matrix with total follow-up times for each stage and their corresponding rejection boundaries; \code{n}; \code{p0}; \code{alpha}; \code{type}; \code{tau}; \code{param}; and cval, the boundary parameter for the rule
#' @export
#'
#' @references Martens, M. J., Lian, Q., Geller, N. L., Leifer, E. S., and Logan, B. L. (2024). Sequential monitoring of time-to-event safety endpoints in clinical trials. \emph{Clinical Trials}, ePub ahead of print.
#' @references Kulldorff, M., Davis, R. L., Kolczak, M., Lewis, E., Lieu, T., and Platt, R. (2011). A maximized sequential probability ratio test for drug and vaccine safety surveillance. \emph{Sequential Analysis}, \strong{30(1)}, 58–78.
#' @references Zacks, S. and Mukhopadhyay, N. (2006). Exact risks of sequential point estimators of the exponential parameter. \emph{Sequential Analysis}, \strong{25(2)}, 203–226.
#'
#' @examples
#'\dontrun{# Survival Pocock test in 50 patient cohort at 10% level, expected toxicity
#' # probability of 20%, 100 day observation period
#' calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="Pocock",tau=100)
#'
#' # Survival Wang-Tsiatis test with Delta = 0.25 in 50 patient cohort at 10% level,
#' # expected toxicity probability of 20%, 100 day observation period
#' calc.rule.surv(n=50,p0=0.20,alpha=0.10,type="WT",tau=100,param=0.25)
#'
#' # Gamma-Poisson test with hyperparameters (1, 1000) in 100 patient cohort at 5% level,
#' # expected toxicity probability of 10%, 60 day observation period
#' calc.rule.surv(n=100,p0=0.10,alpha=0.05,type="GP",tau=60,param=c(1,1000))
#'
#' # Truncated exponential SPRT with p1 = 0.3 in 100 patient cohort at 5% level,
#' # expected toxicity probability of 10%, 60 day observation period
#' calc.rule.surv(n=100,p0=0.10,alpha=0.05,type="SPRT",tau=60,param=0.3)
#'}

calc.rule.surv = function(n, p0, alpha, type, tau, maxInf="expected", param=NULL){
  cval = findconst.surv(n,p0,alpha,type,tau,param=param,maxInf=maxInf)

  bdry = calc.bnd.surv(n=n,p0=p0,type=type,tau=tau,cval=cval,maxInf=maxInf,
                       param=param)
  val = bdry$Rule
  colnames(val) = c("Total follow up time","Reject bdry")
  val2 = list(Rule=val,n=n,p0=p0,alpha=alpha,type=type,tau=tau,maxInf=maxInf,
              param=param,cval=cval)
  return(structure(val2, class = "rule.surv"))
}
