#' @title Plot Stopping Rule (Survival Data)
#' @description
#' Display a stopping rule graphically as a curve for time-to-event data
#'
#' @param x A 'rule.surv' object calculated by \code{calc.rule.surv()} function
#' @param ... Other parameters passed to the \code{plot} function.
#'
#' @return
#' @export
#'
#' @examples
#' pocock.rule <- calc.rule.surv(n = 30, tau = 100, p0 = 0.1, type = "Pocock", alpha = 0.05)
#' plot(pocock.rule, col = "red")

plot.rule.surv <- function(x,smooth=TRUE,xlim=c(0,max(x$Rule[,1])),
                           ylim=c(0,max(x$Rule[,2])+1),xlab="Total Exposure Time",
                           ylab="# Events",...){
  if(smooth==TRUE){
    f = bdryfcn.surv(x$n,x$p0,x$cval,x$type,x$param)
    curve(f,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  }
  else if(smooth==FALSE){
    tab = x$Rule; tab = tab[tab[,1]>=tab[,2],]; ltype = 'l';
    plot(tab,type=ltype,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  }

  #x = x$Rule
  #NextMethod('plot', type = "l", xlim = c(0, max(x[,1])),
  #           ylim = c(0, max(x[,2])+1), xlab = "Total follow up time",
  #           ylab = "# Events",... )
}
