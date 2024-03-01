library(ctsmr)

model <- ctsm$new()
model$addObs(X ~ x)
model$addSystem(dx ~ s * dw)
model$setVariance(X ~ sigma_X^2)

model$setParameter(s = c(init = 0.1, lb = 0, ub = 1))
model$setParameter(sigma_X = 0.1)
model$setParameter(x0 = 0.1)

data <- data.frame(X = cumsum(rnorm(100, 0, 1)), t = 1:100)

# predict(fit1, newdata = data, firstorderinputinterpolation=TRUE,)
makefit <- function(model){
    ## Needs to be runned for simulate to work
    model$AnalyseModel()
    ## Make a fit object to pass to simulate.ctsmr
    fit <- list()
    fit$model <- model
    fit$xm <- model$ParameterValues$initial
    names(fit$xm) <- row.names(model$ParameterValues)
    ## Can be nessary such that xm has the right length
    fit$xm <- fit$xm[model$pars]
    ##
    class(fit) <- "ctsmr"
    ##
    return(fit)
}

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
  print(nll)
  if(printit){ paste("Neg. loglikelihood:", print(nll)) }
  return(nll)
}

fit <- makefit(model)

fitnlminb <- nlminb(fit$xm, nllikelihood, fit=fit, D=data[interval,], firstorder=TRUE, c=3, n.ahead=5,
                    control=list(iter.max=400))

