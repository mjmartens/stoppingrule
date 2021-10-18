#' @title Stopping Rule Calculation
#' @description Calculate a stopping rule for safety monitoring
#'
#' @param ns A vector of sample sizes at which sequential testing is performed
#' @param p0 The toxicity rate under the null hypothesis
#' @param type The method used for constructing the stopping rule
#' @param param Extra parameter(s) needed for Bayesian and SPRT methods
#' @param alpha The desired type I error / false positive rate for the stopping rule
#' @param iter The number of iterations used to search for the boundary
#'
#' @return A matrix with two columns: the sample sizes \code{ns} and their corresponding rejection boundaries
#' @export
#'
#' @references Goldman, A.I. (1987). Issues in designing sequential stopping rules for monitoring side effects in clinical trials. \emph{Controlled clinical trials}  \strong{8(4)}, 327-37.
#' @references Geller, N.L., Follman, D., Leifer, E.S. and Carter, S.L. (2003). Design of early trials in stem cell transplantation: a hybrid frequentist-Bayesian approach. \emph{Advances in Clinical Trial Biostatistics}.
#' @references Ivanova, A., Qaqish, B.F. and Schell, M.J. (2005). Continuous toxicity monitoring in phase II trials in oncology. \emph{Biometrics} \strong{61(2)}, 540-545.
#' @references Kulldorff, M., Davis, R.L., Kolczak, M., Lewis, E., Lieu, T. and Platt, R. (2011). A maximized sequential probability ratio test for drug and vaccine safety surveillance. \emph{Sequential analysis} \strong{30(1)}, 58-78.
#' @examples
#' # Binomial Pocock test in 50 patient cohort at 10% level, expected toxicity rate of 20%
#' calc.rule(ns=1:50,p0=0.20,type="Pocock",alpha=0.10)
calc.rule = function(ns,p0,type,param=NULL,alpha,iter=50) {
  k = length(ns)
  n = tail(ns,1)
  bs = 2:(n+1)
  if(type=="Pocock") {
    l = 0
    u = 1
    const = findconst(ns,p0,"Pocock",alpha,l,u,iter,NA)
  }
  else if(type=="WT") {
    l = qnorm(1-alpha)
    u = qnorm(1-alpha/n)
    const = findconst(ns,p0,"WT",alpha,l,u,iter,param)
  }
  else if(type=="BB") {
    const = findconst(ns,p0,"BB",alpha,0,1,iter,param)
  }
  else if(type=="SPRT") {
    l = 0.5*qchisq(1-alpha,1)
    u = 0.5*qchisq(1-alpha/n,1)
    const = findconst(ns,p0,"SPRT",alpha,l,u,iter,param)
  }
  else if(type=="MaxSPRT") {
    l = 0.5*qchisq(1-alpha,1)
    u = 0.5*qchisq(1-alpha/n,1)
    const = findconst(ns,p0,"MaxSPRT",alpha,l,u,iter,NA)
  }
  bs = calc.bnd(n,p0,const,type,param)

  # Truncate boundary to match stagewise sample sizes requested
  bdry = bs[ns]
  val = cbind(ns,bdry)
  colnames(val) = c("N evaluable","Reject bdry")

  class(val) <- c("rule", class(val))
  return(val)
}
