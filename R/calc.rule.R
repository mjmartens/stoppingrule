#' @title Stopping Rule Calculation
#' @description A wrapper function to calculate a stopping rule for safety monitoring for time-to-event data or binary data
#'
#' @param data.type Data and monitoring strategy, choices include 'bin' for binary data, 'surv' for time-to-event data using Poisson approach, and 'tite' for time-to-event data using TITE method.
#' @param ... Other options to be passed to the corresponding stopping rule calculation. Please refer to the corresponding data type-specific \code{calc.rule()} function for more details
#'
#' @return Please refer to the corresponding data type-specific \code{calc.rule()} function for details on its output
#' @export
#'
#' @examples
#' \dontrun{
#' calc.rule(data.type="bin",ns=1:50,p0=0.20,alpha=0.10,type="WT",param=0.25)
#' calc.rule(data.type="surv",n=50,p0=0.20,alpha=0.10,type="WT",tau=100,param=0.25)
#' calc.rule(data.type="tite",n=100,p0=0.10,alpha=0.05,type="BB",param=c(1,9))
#'}
calc.rule <- function(data.type, ...){
  if (tolower(data.type) == "surv"){
    calc.rule.surv(...)
  } else if (tolower(data.type) == "bin"){
    calc.rule.bin(...)
  } else if (tolower(data.type) == "tite"){
    calc.rule.tite()
  } else {print("Error: data.type must be 'bin', 'surv', or 'tite'")}
}
