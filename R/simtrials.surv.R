#' @title Simulate trials with safety monitoring by survival data stopping rules
#' @description
#' Internal workhorse function used to simulate trials with safety monitoring by
#' survival data stopping rules. The provided stopping \code{rule} is used for
#' monitoring of \code{MC} simulated trials. For each trial, a random sample is
#' generated from a Weibull distribution with shape parameter \code{s} to attain
#' a toxicity rate of \code{p}. Enrollment times are simulated over an accrual
#' period of duration \code{A} under a uniform (0,\code{A}) distribution.
#'
#' @param rule A \code{rule.surv} object with the safety stopping rule for evaluation
#' @param p The probability of a toxicity occurring in \code{tau} units of time
#'          under the null hypothesis
#' @param MC Number of Monte Carlo replicated datasets to simulate
#' @param A Length of accrual period
#' @param s Shape parameter for the Weibull distribution; default value is 1
#'          (exponential distribution)
#'
#' @return A matrix with \code{MC} rows and 14 columns, one row per simulated
#'         trial. Columns include the stopping rule type and design parameters,
#'         the numbers of events and enrolled patients, the total follow-up time
#'         in the cohort, the calendar time when the study ends, the reject/no reject
#'         decision, and the last stage of monitoring reached when the study ends.
#' @export
#'
#' @examples
#' set.seed(13)
#' wt_rule = calc.rule.surv(n=46,p0=0.15,alpha=0.10,type="WT",tau=100,param=0.25)
#' sims = simtrials.surv(rule=wt_rule,p=0.15,MC=1000,A=1095)
#' c(mean(sims$stopped),mean(sims$n.Toxicity),mean(sims$n.Enrolled),mean(sims$Calendar.Time))
#' sims = simtrials.surv(rule=wt_rule,p=0.35,MC=1000,A=1095)
#' c(mean(sims$stopped),mean(sims$n.Toxicity),mean(sims$n.Enrolled),mean(sims$Calendar.Time))
#'
#' gp_rule = calc.rule.surv(n=46,p0=0.15,alpha=0.10,type="GP",tau=100,param=11.5*c(-log(1-0.15),100))
#' sims = simtrials.surv(rule=gp_rule,p=0.15,MC=1000,A=1095)
#' c(mean(sims$stopped),mean(sims$n.Toxicity),mean(sims$n.Enrolled),mean(sims$Calendar.Time))
#' sims = simtrials.surv(rule=gp_rule,p=0.35,MC=1000,A=1095)
#' c(mean(sims$stopped),mean(sims$n.Toxicity),mean(sims$n.Enrolled),mean(sims$Calendar.Time))
#'

simtrials.surv = function(rule,p,MC,A,s=1){
  # Generate all data and sample indicator
  n = rule$n
  p0 = rule$p0
  alpha = rule$alpha
  tau = rule$tau
  alldata = simdata_weibull(n*MC,p,rule$tau,A,s)
  alldata$d = 1*I(alldata$survtime<=tau)  # toxicity event indicator
  alldata$t = pmin(alldata$survtime,tau)  # toxicity observation time

  # Summarize accumulated data for each dataset at toxicity times
  observed = matrix(NA,n,4*MC)
  for(m in 1:MC) {
    data = alldata[n*(m-1)+1:n,]
    monitoring.times = sort(data$enrtime[data$d == 1] + data$t[data$d == 1])
    n_events = length(monitoring.times)
    if(n_events>0) {
      for (i in 1:n_events){
        follow.up = pmax(pmin(monitoring.times[i] - data$enrtime,data$t),0)
        observed[i,4*(m-1)+1:4] = c(sum(follow.up),i,monitoring.times[i],sum(follow.up > 0))
      }
    }
    observed[n,4*(m-1)+1:4] = c(sum(data$t),n_events,A+tau,n)
  }
  colnames(observed) <- rep(c("Total.followup", "n.Toxicity","Calendar.Time","n.Enrolled"),MC)
  rm(alldata)

  # Generate each stopping rule, apply it across datasets
  monitor = data.frame(matrix(NA,MC,14))
  Rule = rule$Rule

  # Fill in rows without stopping possible
  bmin = min(Rule[,2])
  bmax = max(Rule[,2])
  if(bmin>1) {Rule = rbind(matrix(c(rep(NA,bmin-1),1:(bmin-1)),ncol=2),Rule,matrix(c(rep(NA,n-bmax),(bmax+1):n),ncol=2))}
  else {Rule = rbind(Rule,matrix(c(rep(NA,n-bmax),(bmax+1):n),ncol=2))}
  colnames(Rule) = c("Total.followup","Reject bdry")

  # Compare stopping rule to observed data
  #suppressWarnings(decisions = (Rule[,1] >= observed[,4*1:MC-3]))
  decisions = (Rule[,1] >= observed[,4*1:MC-3])
  rejects = colMaxs(decisions+0,na.rm=TRUE)
  rejects[rejects==-Inf] = 0
  f = function(x) {
    val = which(x==TRUE)[1]
    if(length(val)==0 | is.na(val)) {val=bmax}
    return(val)
  }
  laststage = apply(decisions,2,f)

  # Use pointwise matrix products to streamline performance metric calculations
  indices = decisions
  indices[which(indices==0)] = Inf
  indices[n,] = 1
  Calendar.time = colMins(indices*observed[,4*1:MC-1],na.rm=TRUE)
  n.Toxicity = colMins(indices*observed[,4*1:MC-2],na.rm=TRUE)
  Total.followup = colMins(indices*observed[,4*1:MC-3],na.rm=TRUE)
  n.Enrolled = colMins(indices*matrix(1:n,nrow=n,ncol=MC),na.rm=TRUE)

  monitor[1:MC,] = cbind(n=n,p0=p0,p=p,tau=tau,A=A,alpha=alpha,type="GP",cval=rule$cval,Total.followup=Total.followup,
                                    n.Toxicity=n.Toxicity,Calendar.Time=Calendar.time,n.Enrolled=n.Enrolled,rejects,laststage)
  colnames(monitor) <- c("n","p0","p","tau","A","alpha","rule","cval","Total.followup", "n.Toxicity","Calendar.Time","n.Enrolled","stopped","laststage")
  idx = c(1:6,8:14)
  monitor[,idx] = apply(monitor[,idx],2,function(x) as.numeric(as.character(x)))

  return(monitor)
}
