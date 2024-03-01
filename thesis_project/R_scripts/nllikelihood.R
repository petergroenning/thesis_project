## The negative loglikelihood
nllikelihood <- function(xm, fit, D, firstorder=TRUE, c=3, n.ahead=1, printit=TRUE){
  if(printit){ print(format(xm,digits=2)) }
  fit$xm <- xm
  ## loglikelihood initialization
  nll <- 0
  ## use predict to compute the likelihood using the parameters in xm.
  Pred <- try(predict(fit, newdata=D, firstorderinputinterpolation=firstorder, n.ahead=n.ahead))
  if(class(Pred) == "try-error"){ return(NA) }
  ## add logligelihood for all outputs individually
  out_pred <- Pred$output$pred 
  for(i in 1:ncol(out_pred)){
    nm <- names(out_pred)[i]
    yhat <- out_pred[ ,nm]
    sd <- Pred$output$sd[ ,nm]
    y <- D[ ,nm]
    ## c is hubers psi. when the normalized residual is larger than c or smaller than
    ## minus c, the loglikelihood continues as a square root (for robustness)
    ressq <- (y-yhat)^2 / sd^2
    ressq[ressq>c^2] <- c*(2*sqrt(ressq[ressq>c^2])-c)
    nll <- nll + 0.5*( sum(log(Pred$output$sd^2)+ressq)  )
  }
  l <- ncol(out_pred) ## dimension of output vector
  N <- nrow(out_pred) ## number of observations
  nll <- nll + 0.5*l*N*log(2*pi)
  if(printit){ paste("Neg. loglikelihood:", print(nll)) }
  return(nll)
}
