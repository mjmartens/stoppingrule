#' @title Stopping Boundary Calculation (Survival Data)
#' @description
#' Internal workhorse function to calculate stopping boundary for a given method for time-to-event data
#' @param n Maximum sample size for safety monitoring
#' @param tau Observation period
#' @param p0 The probability of a toxicity occurring in \code{tau} units of time under the null hypothesis
#' @param cval Critical for the stopping rule. For Wang-Tsiatis tests, this is the Delta parameter. For the Bayesian Gamma-Poisson method, this is the threshold on the posterior probability. For the truncated SPRT, this is the threshold on the log likelihood ratio. For the MaxSPRT, this is the threshold on the log generalized likelihood ratio.
#' @param type The method used for constructing the stopping rule
#' @param param Extra parameter(s) needed for certain stopping rule methods. For Wang-Tsiatis tests, this is the Delta parameter. For modified SPRT, this is the targeted alternative toxicity probability p1. For Bayesian Gamma-Poisson model, this is the pair of hyperparameters for the gamma prior on the toxicity event rate.
#'
#' @return A list of three items: tau, number of events that can trigger a stop, and the corresponding total follow up time.

calc.bnd.surv <- function(n, p0, cval, tau, type, param = NULL){
  lambda0 <- -log(1 - p0)/tau
  Umax <- n*tau

  if (type == "Pocock"){
    f <- function(U){
      lambda0*U + cval*sqrt(lambda0*U)
    }
    f.inverse <- function(D){
      (D + 0.5*(cval^2 - sqrt(4*D*cval^2 + cval^4)))/lambda0
    }
  }

  else if (type == "OBF"){
    f <- function(U){
      lambda0*U + cval*sqrt(lambda0*Umax)
    }
    f.inverse <- function(D){
      (D - cval*sqrt(lambda0*Umax))/lambda0
    }
  }

  else if (type == "WT"){
    Delta <- param
    f <- function(U){
      lambda0*U + cval*sqrt(lambda0)*Umax^(0.5 - Delta)*U^Delta
    }
    f.inverse <- function(D){
      uniroot(function(u){f(u) - D}, lower = 0, upper = 2*n*tau)$root
    }
  }

  else if (type == "SPRT"){
    p1 <- param
    lambda1 <- -log(1 - p1)/tau
    f <- function(U){
      (cval + (lambda1 - lambda0)*U)/log(lambda1/lambda0)
    }
    f.inverse <- function(D){
      (D*(log(lambda1) - log(lambda0)) - cval)/(lambda1 - lambda0)
    }
  }

  else if (type == "MaxSPRT"){
    f <- function(U){
      if (U ==0 ){return(0)}
      else {return((cval - lambda0*U)*(lambertWp(1/exp(1)*(cval/(lambda0*U) - 1)))^(-1))}
    }
    f.inverse <- function(D){
      -D/lambda0*lambertWp(-exp(-cval/D - 1))
    }
  }

  if (type != "Bayesian"){
    S <- seq(from = max(ceiling(f(0)),1), to = ceiling(f(Umax)))
    # FirstPositive <- which(S > 0)[1]
    # S <- S[FirstPositive:length(S)]
    dmin <- S[1]
    dmax <- S[length(S)]
    m <- dmax - dmin + 1

    ud <- rep(NA, m)
    for (i in 1:m){
      ud[i] <- f.inverse(S[i])
    }
    ud[length(ud)] <- Umax
  }


  if (type == "Bayesian"){
    # hyperparameters
    if (length(param) == 1){
      k <- param
      l <- k/lambda0
    } else {
      k <- param[1] # shape
      l <- param[2] # rate
    }

    # find dmin and dmax
    post <- NULL
    for (d in 1:n){
      post[d] <- 1 - pgamma(lambda0, shape = (k + d), rate = l)
    }
    dmin = min(which(post >= cval))

    post1 <- NULL
    for (d in 1:n){
      post1[d] <- 1 - pgamma(lambda0, shape = (k + d), rate = (l + Umax))
    }
    dmax = min(which(post1 >= cval))
    S = seq(dmin, dmax)
    ud <- rep(NA, length.out = length(S))
    for (j in 1:(length(S) - 1)){
      inner <- function(ud){
        1 - pgamma(lambda0, shape = (k + S[j]), rate = (l + ud))
      }
      if( j == 1){
        ud[j] <- uniroot(function(ud){inner(ud) - cval}, lower = 0, upper = Umax,
                         extendInt = "no", trace = 2)$root
      } else {
        ud[j] <- uniroot(function(ud){inner(ud) - cval}, lower = ud[j-1], upper = Umax,
                         extendInt = "no", trace = 2)$root
      }
    }
    ud[length(ud)] <- Umax
  }

  return(list(tau = tau, S = S, ud = ud))
}
