% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opchars.surv.R
\name{opchars.surv}
\alias{opchars.surv}
\title{Operating Characteristics Function (Survival Data)}
\usage{
opchars.surv(rule, p)
}
\arguments{
\item{rule}{A 'rule.surv' object calculated by \code{calc.rule.surv()} function}

\item{p}{The toxicity probability}
}
\value{
A list with the following elements: p, the corresponding rejection probability, and the corresponding expected total follow up time and number of events at the point of stopping/study end
}
\description{
Internal workhorse function to calculate operating characteristics for a given stopping rule and toxicity probability
}