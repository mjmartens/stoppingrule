#' @title Stopping Rule Calculation (Binary Data)
#' @description Calculate a stopping rule for safety monitoring, treating toxicities as binary data
#'
#' @param ns A vector of sample sizes at which sequential testing is performed
#' @param p0 The toxicity probability under the null hypothesis
#' @param alpha The desired type I error / false positive rate for the stopping rule
#' @param type The method used for constructing the stopping rule. Choices include a Pocock test ("Pocock"), an O'Brien-Fleming test ("OBF"), a Wang-Tsiatis test ("WT"), the Bayesian beta-binomial method ("BB") proposed by Geller et al. 2003, the Bayesian beta-binomial method ("CC") proposed by Chen and Chaloner 2006, a truncated SPRT ("SPRT"), and a maximized SPRT ("MaxSPRT").
#' @param param Extra parameter(s) needed for certain stopping rule methods. For binomial Wang-Tsiatis tests, this is the Delta parameter. For the Geller et al. method, this is the pair of hyperparameters (a,b) for the beta prior on the toxicity probability. For Chen and Chaloner's method, this is the pair of hyperparameters (a,b) for the beta prior on the toxicity probability, the targeted alternative toxicity probability p1, and the threshold nu for the posterior probability that p > p1. For modified SPRT, this is the targeted alternative toxicity probability p1.
#' @param iter The number of iterations used to search for the boundary
#'
#' @return A rule.bin object, which is a list with the following elements: Rule, a two-column matrix with the sample sizes \code{ns} and their corresponding rejection boundaries; ns; p0; alpha; type; param; and cval
#' @export
#'
#' @references Chen, C. and Chaloner, K. (2006). A Bayesian stopping rule for a single arm study: With a case study of stem cell transplantation. \emph{Statistics in Medicine} \strong{25(17)}, 2956-66.
#' @references Geller, N.L., Follman, D., Leifer, E.S. and Carter, S.L. (2003). Design of early trials in stem cell transplantation: a hybrid frequentist-Bayesian approach. \emph{Advances in Clinical Trial Biostatistics}.
#' @references Goldman, A.I. (1987). Issues in designing sequential stopping rules for monitoring side effects in clinical trials. \emph{Controlled Clinical Trials}  \strong{8(4)}, 327-37.
#' @references Ivanova, A., Qaqish, B.F. and Schell, M.J. (2005). Continuous toxicity monitoring in phase II trials in oncology. \emph{Biometrics} \strong{61(2)}, 540-545.
#' @references Kulldorff, M., Davis, R.L., Kolczak, M., Lewis, E., Lieu, T. and Platt, R. (2011). A maximized sequential probability ratio test for drug and vaccine safety surveillance. \emph{Sequential Analysis} \strong{30(1)}, 58-78.
#' @references Martens, M.J. and Logan, B.R. (2023). Statistical Rules for Safety Monitoring in Clinical Trials. \emph{Clinical Trials} \strong{Article in press}.
#' @references Pocock, S.J. (1977). Group sequential methods in the design and analysis of clinical trials. \emph{Biometrika} \strong{64(2)}, 191-199.
#' @references Wang, S.K. and Tsiatis, A.A. (1987). Approximately optimal one-parameter boundaries for group sequential trials. \emph{Biometrics} \strong{193-199}.
#'
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity probability of 20%
#' calc.rule.bin(ns=1:50,p0=0.20,alpha=0.10,type="Pocock")

calc.rule.bin = function(ns,p0,alpha,type,param=NULL,iter=50) {
  k = length(ns)
  n = tail(ns,1)
  bs = 2:(n+1)
  if(type=="Pocock") {
    l = 0
    u = 1
    const = findconst.bin(ns,p0,alpha,"Pocock",l,u,iter,NA)
  }
  else if(type=="OBF") {
    l = qnorm(1-alpha)
    u = qnorm(1-alpha/n)
    const = findconst.bin(ns,p0,alpha,"OBF",l,u,iter)
  }
  else if(type=="WT") {
    l = qnorm(1-alpha)
    u = qnorm(1-alpha/n)
    const = findconst.bin(ns,p0,alpha,"WT",l,u,iter,param)
  }
  else if(type=="BB") {
    const = findconst.bin(ns,p0,alpha,"BB",0,1,iter,param)
  }
  else if(type=="CC") {
    const = findconst.bin(ns,p0,alpha,"CC",0,1,iter,param)
  }
  else if(type=="SPRT") {
    l = 0.5*qchisq(1-alpha,1)
    u = 0.5*qchisq(1-alpha/n,1)
    const = findconst.bin(ns,p0,alpha,"SPRT",l,u,iter,param)
  }
  else if(type=="MaxSPRT") {
    l = 0.5*qchisq(1-alpha,1)
    u = 0.5*qchisq(1-alpha/n,1)
    const = findconst.bin(ns,p0,alpha,"MaxSPRT",l,u,iter,NA)
  }
  bs = calc.bnd.bin(n,p0,type,const,param)

  # Truncate boundary to match stagewise sample sizes requested
  if(type=="Pocock") {const = 1-const}
  bdry = bs[ns]
  tab = cbind(ns,bdry)
  colnames(tab) = c("N evaluable","Reject bdry")
  val = list(Rule=tab,ns=ns,p0=p0,type=type,alpha=alpha,param=param,cval=const)
  class(val) = "rule.bin"
  return(val)
}
