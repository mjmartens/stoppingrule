#' @title Operating Characteristics Function (Binary Data)
#' @description Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
#'
#' @param rule A \code{rule.bin} object calculated by \code{calc.rule.bin()} function
#' @param p The toxicity probability
#' @param tau Length of observation period
#' @param A Length of the enrollment period.
#'
#' @return A list containing the toxicity probability \code{p}, and the corresponding rejection probability and expected number of events. If \code{tau} and {A} are also specified, the expected number of enrolled patients and the expected calendar time at the point of stopping/study end are also included.

opchars.bin = function(rule, p, tau = NULL, A = NULL){
  n_k = rule$Rule[,1]
  b_k = rule$Rule[,2]
  M = length(n_k)
  m_k = diff(c(0,n_k))
  D = matrix(0, nrow = M, ncol=max(n_k)+1) # rows are stages, columns are number of events
  reject.prob = rep(0, M)

  ### rejection.prob: stagewise stopping probabilities
  ### D: stagewise outcome probabilities
  for (s in 0:n_k[1]){
    D[1, s+1] = dbinom(s, m_k[1], p)
  }
  for (k in 2:M){
    for (s in 0:n_k[k]){
      lower = max(0, s - m_k[k])
      upper = min(b_k[k - 1] - 1, s)
      if (lower <= upper){
        D[k, s+1] = sum(D[k-1, (lower+1):(upper+1)]*dbinom(s - lower:upper, m_k[k],p))
      }
    }
  }

  for (k in 1:M){
    reject.prob[k] = sum(D[k,-1: - b_k[k]])
  }

  power = sum(reject.prob)


  # Calculate the below two characteritcs only when A and tau were provided
  if (!is.null(A) & !is.null(tau)){
    ## Expected calendar time when study ends
    exp.calendar <- sum(A*n_k*reject.prob/(max(n_k) + 1)) + A*max(n_k)*(1 - power)/(max(n_k)+1) + tau

    ## Expected number of enrolled patients
    temp <- matrix(NA, nrow = M, ncol = max(n_k) - n_k[1])
    for (k in 1:M){
      for (j in (n_k[k] + 1): max(n_k)){
        if(max(n_k) >= (n_k[k] + 1)){
          temp[k, j - 1] <- pbeta(tau/A, j - n_k[k], max(n_k) + 1 - (j - n_k[k]))*reject.prob[k]
        }
      }
    }

    partIII <- sum(temp, na.rm = TRUE)
    exp.enrolled <- sum(n_k*reject.prob) + max(n_k)*(1 - power) + partIII
  }


  ### Expected number of toxicities
  partI <- vector(mode='list', length = M)
  for (k in 1:M){
    temp <- c()
    for (j in b_k[k]:n_k[k]){
      if (n_k[k] >= b_k[k]){
        temp1 <- j*D[k, j + 1] # D matrix's columns are # of events, which starts from zero
        temp <- c(temp, temp1)
      }
    }
    partI[[k]] <- temp
  }
  partI <- sum(unlist(lapply(partI, sum)))

  partII <- c()
  for (j in 0:(b_k[M] - 1)){
    temp <- j*D[M, j + 1]
    partII <- c(partII, temp)
  }
  partII <- sum(partII)

  exp.toxicities.old <- partI + partII
  # Calculate expected number of toxicities in those pending evaluation
  if (!is.null(A) & !is.null(tau)){
    partIII <- vector(mode = "list", length = M)
    for (k in 1:M){
      temp <- c()
      for(j in (n_k[k] + 1):max(n_k)){
        if (max(n_k) >= (n_k[k] + 1)){
          temp1 <- pbeta(tau/A, j - n_k[k], max(n_k) + 1 - (j - n_k[k]))*reject.prob[k]
          temp <- c(temp, temp1)
        }
      }
      partIII[[k]] <- temp
    }
    partIII <- sum(unlist(lapply(partIII, sum)))*p

    exp.toxicities.new <- partI + partII + partIII
  }

  if (!is.null(A) & !is.null(tau)){
    return(list(
      # stage.outcome.prob = D,
      #stage.reject.prob = reject.prob,
      power = power,
      exp.toxicities = exp.toxicities.new,
      exp.enrolled = exp.enrolled,
      exp.calendar = exp.calendar))
  } else {
    return(list(
      power = power,
      exp.toxicities = exp.toxicities.old))
  }
}
