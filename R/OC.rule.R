#' @title Operating Characteristics Function
#' @description A wrapper function to compute operating characteristics for a stopping rule at a set of toxicity rates.
#'
#' @param data.type Type of data, choices include 'bin' for binary data and 'surv' for time-to-event data
#' @param ... Other options to be passed to the corresponding operating characteristics calculation. Please refer to the corresponding \code{OC.rule()} function for more details
#'
#' @return Please refer to the corresponding data type-specific \code{OC.rule()} function for more details
#' @export
#'
#' @examples
#' bb_rule = calc.rule(data.type="bin",ns=1:50,p0=0.20,alpha=0.10,type="BB",param=c(2,8))
#' gp_rule = calc.rule(data.type="surv",n=50,p0=0.20,alpha=0.10,type="GP",tau=60,param=c(1,1000))
#' OC.rule(data.type="bin",rule=bb_rule,ps=seq(0.1, 0.5, 0.1))
#' OC.rule(data.type="bin",rule=bb_rule,ps=seq(0.1, 0.5, 0.1),tau=60,A=730)
#' OC.rule(data.type="surv",rule=gp_rule,ps=seq(0.1, 0.5, 0.1),MC=1000, A=730)

OC.rule <- function(data.type,...){
  if (tolower(data.type) == "surv"){
    OC.rule.surv(...)
  } else if (tolower(data.type) == "bin"){
    OC.rule.bin(...)
  }
  else {print("Error: data.type must be 'bin' or 'surv'")}
}
