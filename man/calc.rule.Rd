% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.rule.R
\name{calc.rule}
\alias{calc.rule}
\title{Stopping Rule Calculation}
\usage{
calc.rule(data.type, ...)
}
\arguments{
\item{data.type}{Type of data, choices include 'bin' for binary data and 'surv' for time-to-event data}

\item{...}{Other options to be passed to the corresponding stopping rule calculation. Please refer to the corresponding data type-specific \code{calc.rule()} function for more details}
}
\value{
Please refer to the corresponding data type-specific \code{calc.rule()} function for details on its output
}
\description{
A wrapper function to calculate a stopping rule for safety monitoring for time-to-event data or binary data
}
\examples{
calc.rule(data.type="bin",ns=1:50,p0=0.20,alpha=0.10,type="WT",param=0.25)
calc.rule(data.type="surv",n=50,p0=0.20,alpha=0.10,type="WT",tau=100,param=0.25)

}
