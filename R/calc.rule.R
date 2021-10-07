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
#' @return A matrix with two columns: first column is the vector of sample sizes \code{ns}, second column is a corresponding vector of rejection boundaries at these sizes
#' @export
#'
#' @references Blah
#' @examples
#' calc.rule(ns=1:30,p0=0.15,type="Pocock",alpha=0.10)
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
