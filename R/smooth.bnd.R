#' @title Smooth Stopping Rule Boundary
#' @description Function to Smooth Stopping RUle Boundary from Step to Piecewise Linear Function to
#'
#' @param rule A matrix with two columns: the vector of sample sizes at which the stopping rule is evaluated, and a corresponding vector of rejection boundaries at these sizes
#'
#' @return A matrix with two columns: the vector of sample sizes at which the stopping rule is evaluated, and a corresponding vector of smoothed rejection boundaries at these sizes
smooth.bnd = function(rule) {
  n = rule[rule[,1]>=rule[,2],1]
  b = rule[rule[,1]>=rule[,2],2]
  m = length(b)
  bmin = min(b)
  bidx = which(b[-1]-b[-m]>0)
  if(bidx[1]>1) {
    b[1:(bidx[1]-1)] = bmin - (bidx[1]-1):1/bidx[1]
  }
  for(i in 1:(length(bidx)-1)) {
    b[(bidx[i]+1):(bidx[i+1])] = b[bidx[i]] + 1:(bidx[i+1] - bidx[i])/(bidx[i+1] - bidx[i])
    if(i==(length(bidx)-1) & tail(bidx,1) != m) {
      b[(bidx[i+1]+1):m] = b[bidx[i+1]] + 1:(m - bidx[i+1])/(bidx[i+1] - bidx[i])
    }
  }
  return(cbind(n,b))
}
