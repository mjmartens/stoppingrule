#' @title Stopping Rule Boundary Function
#' @description A wrapper function to calculate the boundary for a given stopping rule for safety monitoring for time-to-event data or binary data
#'
#' @param data.type Type of data, choices include 'bin' for binary data and 'surv' for time-to-event data
#' @param ... Other options to be passed to the corresponding stopping rule calculation. Please refer to the corresponding data type-specific \code{bdryfcn()} function for more details
#'
#' @return A univariate function that defines the rejection boundary at any number of evaluable patients (binary data) or amount of follow-up time (time-to-event data)
#' @export
#'
bdryfcn = function(data.type, ...){
  if (tolower(data.type) == "surv"){
    bdryfcn.surv(...)
  }
  else if (tolower(data.type) == "bin"){
    bdryfcn.bin(...)
  }
  else {print("Error: data.type must be 'bin' or 'surv'")}
}
