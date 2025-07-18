#' @title Stopping Rule Boundary Function
#' @description A wrapper function to calculate the boundary for a given stopping rule for safety monitoring for time-to-event data or binary data
#'
#' @param data.type Data and monitoring strategy, choices include 'bin' for binary data, 'surv' for time-to-event data using Poisson approach, and 'tite' for time-to-event data using TITE method.
#' @param ... Other options to be passed to the corresponding stopping rule calculation. Please refer to the corresponding data type-specific \code{bdryfcn()} function for more details
#'
#' @return A univariate function that defines the rejection boundary at any number of evaluable patients (binary data), amount of follow-up time (time-to-event data), or effective sample size accounting for partial follow up (TITE method).
#' @export
#'
bdryfcn = function(data.type, ...){
  if (tolower(data.type) == "surv"){
    bdryfcn.surv(...)
  }
  else if (tolower(data.type) == "bin"){
    bdryfcn.bin(...)
  } else if (tolower(data.type == "tite")){
    bdryfcn.bin.inverse(...)
  } else {print("Error: data.type must be 'bin', 'surv', or 'tite'")}
}
