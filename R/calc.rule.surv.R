#' @title Stopping Rule Calculation (Survival Data)
#' @description Calculate a stopping rule for safety monitoring for time-to-event data
#'
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The probability of a toxicity occurring in \code{tau}] units of time under the null hypothesis
#' @param alpha The desired type I error/false positive rate for the stopping rule
#' @param type The method used for constructing the stopping rule. Choices including a Pocock test ("Pocock"),
#' a O'Brein-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian Gamma-Poisson method ("Bayesian"),
#' a modified sequential probability ratio test ("SPRT"), and a maximized SPRT ("MaxSPRT")
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity probability p1. For Bayesian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity event rate.
#'
#' @return A rule.surv object, which is a list with the following elements: Rule, a two-column matrix with total follow-up times for each stage and their corresponding rejection boundaries; n; p0; alpha; type; tau; param; and cval
#' @export
#'
#' @references Kulldorff, M., Davis, R. L., Kolczak, M., Lewis, E., Lieu, T., and Platt, R. (2011). A maximized sequential probability ratio test for drug and vaccine safety surveillance. \emph{Sequential Analysis}, \strong{30(1)}, 58–78.
#' @references Zacks, S. and Mukhopadhyay, N. (2006). Exact risks of sequential point estimators of the exponential parameter. \emph{Sequential Analysis}, \strong{25(2)}, 203–226.
#'
#' @examples
#' calc.rule.surv(n = 30, p0 = 0.1, alpha = 0.05, type = "Pocock", tau = 100)

calc.rule.surv <- function(n, p0, alpha, type, tau, param=NULL){
  cval <- findconst.surv(n = n, tau = tau, p0 = p0, type = type, param = param, alpha = alpha)

  bdry <- calc.bnd.surv(n = n, p0 = p0, cval = cval, tau = tau, type = type, param = param)
  val <- cbind(bdry$ud, bdry$S)
  colnames(val) <- c("Total follow up time","Reject bdry")
  val2 <- list(Rule=val,n=n,p0=p0,alpha=alpha,type=type,tau=tau,param=param,cval=cval)
  return(structure(val2, class = "rule.surv"))
}
