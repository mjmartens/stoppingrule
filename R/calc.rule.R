#' @title Stopping Rule Calculation
#' @description A wrapper function to calculate a stopping rule for safety monitoring for time-to-event data or binary data
#'
#' @param data.type Type of data, choices include 'bin' for binary data and 'surv' for time-to-event data
#' @param ... Other options to be passed to the corresponding stopping rule calculation. Please refer to the corresponding data type-specific \code{calc.rule()} function for more details
#'
#' @return Please refer to the corresponding data type-specific \code{calc.rule()} function for details on its output
#' @export
#'
#' @examples
#' calc.rule(data.type="bin",ns=1:50,p0=0.20,alpha=0.10,type="WT",param=0.25)
#' calc.rule(data.type="surv",n=50,p0=0.20,alpha=0.10,type="WT",tau=100,param=0.25)
#'
calc.rule <- function(data.type, ...){
  if (tolower(data.type) == "surv"){
    calc.rule.surv(...)
  } else if (tolower(data.type) == "bin"){
    calc.rule.bin(...)
  }
  else {print("Error: data.type must be 'bin' or 'surv'")}
}
