library(ctsmr)
source('thesis_project/models/ctsm_utils.R')
library(ggplot2)
library(tidyverse)
library(splines)
make_model <- function(data, old_fit){

    obj <- ctsm$new()
    if (missing(old_fit)){
        old_fit <- NULL
    }

    # Add observation equations and variances
    obj$addObs(X0 ~ x0m)
    obj$addObs(X1 ~ x1m)
    obj$addObs(X2 ~ x2m)
    obj$addObs(X3 ~ x3m)
    obj$addObs(X4 ~ x4m)
    obj$addObs(X5 ~ x5m)
    obj$addObs(X6 ~ x6m)
    obj$addObs(X7 ~ x7m)
    obj$addObs(X8 ~ x8m)
    obj$addObs(X9 ~ x9m)
    obj$addObs(X10 ~ x10m)
    obj$addObs(X11 ~ x11m)
    obj$addObs(X12 ~ x12m)
    obj$addObs(X13 ~ x13m)
    obj$addObs(X14 ~ x14m)
    obj$addObs(X15 ~ x15m)

    # Set observation equation variances
    obj$setVariance(X0 ~ sigma_X0^2)
    obj$setVariance(X1 ~ sigma_X1^2)
    obj$setVariance(X2 ~ sigma_X2^2)
    obj$setVariance(X3 ~ sigma_X3^2)
    obj$setVariance(X4 ~ sigma_X4^2)
    obj$setVariance(X5 ~ sigma_X5^2)
    obj$setVariance(X6 ~ sigma_X6^2)
    obj$setVariance(X7 ~ sigma_X7^2)
    obj$setVariance(X8 ~ sigma_X8^2)
    obj$setVariance(X9 ~ sigma_X9^2)
    obj$setVariance(X10 ~ sigma_X10^2)
    obj$setVariance(X11 ~ sigma_X11^2)
    obj$setVariance(X12 ~ sigma_X12^2)
    obj$setVariance(X13 ~ sigma_X13^2)
    obj$setVariance(X14 ~ sigma_X14^2)
    obj$setVariance(X15 ~ sigma_X15^2)

    # Add system equations
    obj$addSystem(dx0m ~ 1/(c1*b1x0+c2*b2x0+c3*b3x0+c4*b4x0+c5*b5x0) *  ((x1m - x0m) * (R + FbotVol * Cbot)) * dt + (Tbot - x0m) * FbotIn * Vbot * dt  + sigma_x0 * dw0)
    obj$addSystem(dx1m ~ 1/(c1*b1x1+c2*b2x1+c3*b3x1+c4*b4x1+c5*b5x1) *  ((x0m - x1m) * (R + FbotVol * Cbot) + (x2m - x1m) * (R + FbotVol * Cbot)) * dt + sigma_x1 * dw1)
    obj$addSystem(dx2m ~ 1/(c1*b1x2+c2*b2x2+c3*b3x2+c4*b4x2+c5*b5x2) *  ((x1m - x2m) * (R + FbotVol * Cbot) + (x3m - x2m) * (R + FbotVol * Cbot)) * dt + sigma_x2 * dw2)
    obj$addSystem(dx3m ~ 1/(c1*b1x3+c2*b2x3+c3*b3x3+c4*b4x3+c5*b5x3) *  ((x2m - x3m) * (R + FbotVol * Cbot) + (x4m - x3m) * (R + FbotVol * Cbot)) * dt + sigma_x3 * dw3)
    obj$addSystem(dx4m ~ 1/(c1*b1x4+c2*b2x4+c3*b3x4+c4*b4x4+c5*b5x4) *  ((x3m - x4m) * (R + FbotVol * Cbot) + (x5m - x4m) * (R + FbotVol * Cbot)) * dt + sigma_x4 * dw4)
    obj$addSystem(dx5m ~ 1/(c1*b1x5+c2*b2x5+c3*b3x5+c4*b4x5+c5*b5x5) *  ((x4m - x5m) * (R + FbotVol * Cbot) + (x6m - x5m) * (R + FbotVol * Cbot)) * dt + sigma_x5 * dw5)
    obj$addSystem(dx6m ~ 1/(c1*b1x6+c2*b2x6+c3*b3x6+c4*b4x6+c5*b5x6) *  ((x5m - x6m) * (R + FbotVol * Cbot) + (x7m - x6m) * (R + FbotVol * Cbot)) * dt + sigma_x6 * dw6)
    obj$addSystem(dx7m ~ 1/(c1*b1x7+c2*b2x7+c3*b3x7+c4*b4x7+c5*b5x7) *  ((x6m - x7m) * (R + FbotVol * Cbot) + (x8m - x7m) * (R + FbotVol * Cbot)) * dt + sigma_x7 * dw7)
    obj$addSystem(dx8m ~ 1/(c1*b1x8+c2*b2x8+c3*b3x8+c4*b4x8+c5*b5x8) *  ((x7m - x8m) * (R + FbotVol * Cbot) + (x9m - x8m) * (R + FbotVol * Cbot)) * dt + sigma_x8 * dw8)
    obj$addSystem(dx9m ~ 1/(c1*b1x9+c2*b2x9+c3*b3x9+c4*b4x9+c5*b5x9) *  ((x8m - x9m) * (R + FbotVol * Cbot) + (x10m - x9m) * (R + FbotVol * Cbot)) * dt + sigma_x9 * dw9)
    obj$addSystem(dx10m ~ 1/(c1*b1x10+c2*b2x10+c3*b3x10+c4*b4x10+c5*b5x10) * ((x9m - x10m) * (R + FbotVol * Cbot) + (x11m - x10m) *  (R + FtopVol * Ctop)) * dt + sigma_x10 * dw10)
    obj$addSystem(dx11m ~ 1/(c1*b1x11+c2*b2x11+c3*b3x11+c4*b4x11+c5*b5x11) * ((x10m - x11m) * (R + FtopVol * Ctop) + (x12m - x11m) * (R + FtopVol * Ctop)) * dt + sigma_x11 * dw11)
    obj$addSystem(dx12m ~ 1/(c1*b1x12+c2*b2x12+c3*b3x12+c4*b4x12+c5*b5x12) * ((x11m - x12m) * (R + FtopVol * Ctop) + (x13m - x12m) * (R + FtopVol * Ctop)) * dt + sigma_x12 * dw12)
    obj$addSystem(dx13m ~ 1/(c1*b1x13+c2*b2x13+c3*b3x13+c4*b4x13+c5*b5x13) * ((x12m - x13m) * (R + FtopVol * Ctop) + (x14m - x13m) * (R + FtopVol * Ctop)) * dt + sigma_x13 * dw13)
    obj$addSystem(dx14m ~ 1/(c1*b1x14+c2*b2x14+c3*b3x14+c4*b4x14+c5*b5x14) * ((x13m - x14m) * (R + FtopVol * Ctop) + (x15m - x14m) * (R + FtopVol * Ctop)) * dt + sigma_x14 * dw14)
    obj$addSystem(dx15m ~ 1/(c1*b1x15+c2*b2x15+c3*b3x15+c4*b4x15+c5*b5x15) * ((x14m - x15m) * (R + FtopVol * Ctop)) * dt + (Ttop - x15m) * FtopIn * Vtop * dt + sigma_x15 * dw15)


    obj$addInput("Ttop", "Tmid", "Tbot", "Ftop", "Fbot", "Fmid", "ambientTemp", "FtopIn", "FbotIn", "FtopVol", "FbotVol", "FmidVol")
    obj$addInput("b1x0","b1x1","b1x1","b1x2","b1x3","b1x4","b1x5","b1x6","b1x7","b1x8","b1x9","b1x10","b1x11","b1x12","b1x13","b1x14","b1x15")
    obj$addInput("b2x0","b2x1","b2x1","b2x2","b2x3","b2x4","b2x5","b2x6","b2x7","b2x8","b2x9","b2x10","b2x11","b2x12","b2x13","b2x14","b2x15")
    obj$addInput("b3x0","b3x1","b3x1","b3x2","b3x3","b3x4","b3x5","b3x6","b3x7","b3x8","b3x9","b3x10","b3x11","b3x12","b3x13","b3x14","b3x15")
    obj$addInput("b4x0","b4x1","b4x1","b4x2","b4x3","b4x4","b4x5","b4x6","b4x7","b4x8","b4x9","b4x10","b4x11","b4x12","b4x13","b4x14","b4x15")
    obj$addInput("b5x0","b5x1","b5x1","b5x2","b5x3","b5x4","b5x5","b5x6","b5x7","b5x8","b5x9","b5x10","b5x11","b5x12","b5x13","b5x14","b5x15")
    ###### SET PARAMETERS ######
    # Observation Noise
    obs_std <- 0.15
    obj$setParameter(sigma_X0 = obs_std)
    obj$setParameter(sigma_X1 = obs_std)
    obj$setParameter(sigma_X2 = obs_std)
    obj$setParameter(sigma_X3 = obs_std)
    obj$setParameter(sigma_X4 = obs_std)
    obj$setParameter(sigma_X5 = obs_std)
    obj$setParameter(sigma_X6 = obs_std)
    obj$setParameter(sigma_X7 = obs_std)
    obj$setParameter(sigma_X8 = obs_std)
    obj$setParameter(sigma_X9 = obs_std)
    obj$setParameter(sigma_X10 = obs_std)
    obj$setParameter(sigma_X11 = obs_std)
    obj$setParameter(sigma_X12 = obs_std)
    obj$setParameter(sigma_X13 = obs_std)
    obj$setParameter(sigma_X14 = obs_std)
    obj$setParameter(sigma_X15 = obs_std)

    # System Noise
    lb <- 1e-5
    ub <- 10
    init <- 1e-2
    obj$setParameter(sigma_x0 = c(init =init, lb = lb, ub = ub))
    obj$setParameter(sigma_x1 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x2 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x3 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x4 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x5 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x6 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x7 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x8 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x9 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x10 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x11 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x12 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x13 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x14 = c(init = init, lb = lb, ub = ub))
    obj$setParameter(sigma_x15 = c(init = init, lb = lb, ub = ub))


    # System Parameters
    lb <- 1e-5
    ub <- 1000
    init <- 1
    obj$setParameter(R = c(init=getParam('R', init, old_fit),lb = lb, ub = ub))
    obj$setParameter(C = c(init=getParam('C', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Vbot = c(init=getParam('Vbot', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Vtop =c(init=getParam('Vtop', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Utop = c(init=getParam('Utop', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Cbot =c(init=getParam('Cbot', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Ctop = c(init=getParam('Ctop', init, old_fit), lb = lb, ub = ub))
    lb <- -100
    ub <- 1000
    init <- 1
    obj$setParameter(c1 = c(init=getParam('c1', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(c2 = c(init=getParam('c2', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(c3 = c(init=getParam('c3', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(c4 = c(init=getParam('c4', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(c5 = c(init=getParam('c5', init, old_fit), lb = lb, ub = ub))

    # # System Parameters
    # lb <- 0
    # ub <- 1000
    # init <- 1e-5
    # # obj$setParameter(R = c(init=getParam('R', init, old_fit),lb = lb, ub = ub))
    # # obj$setParameter(C = c(init=getParam('C', init, old_fit), lb = lb, ub = ub))
    # # obj$setParameter(Vbot = c(init=getParam('Vbot', init, old_fit), lb = lb, ub = ub))
    # # obj$setParameter(Vmid = c(init=getParam('Vmid', init, old_fit), lb = lb, ub = ub))
    # # obj$setParameter(Vtop =c(init=getParam('Vtop', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(R =  c(init = 5 ,lb = lb, ub = ub))
    # # obj$setParameter(C = c(init = 10, lb = lb, ub = ub))
    # obj$setParameter(Vbot = c(init = 0, lb = -1, ub = ub))
    # obj$setParameter(Vmid = c(init = 0, lb = -1, ub = ub))
    # obj$setParameter(Vtop = c(init = 0, lb = -1, ub = ub))
    # obj$setParameter(Cbot = c(init=0, lb = -10, ub = 10))
    # obj$setParameter(Ctop = c(init=0, lb = -10, ub = 10))

    # lb <- -1000
    # ub <- 1000
    # init <- 10
    # obj$setParameter(c1 = c(init=init, lb = lb, ub = ub))
    # obj$setParameter(c2 = c(init=init, lb = lb, ub = ub))
    # obj$setParameter(c3 = c(init=init, lb = lb, ub = ub))
    # obj$setParameter(c4 = c(init=init, lb = lb, ub = ub))
    # obj$setParameter(c5 = c(init=init, lb = lb, ub = ub))

    # Initial states
    obj$setParameter(x0m = data$X0[1])
    obj$setParameter(x1m = data$X1[1])
    obj$setParameter(x2m = data$X2[1])
    obj$setParameter(x3m = data$X3[1])
    obj$setParameter(x4m = data$X4[1])
    obj$setParameter(x5m = data$X5[1])
    obj$setParameter(x6m = data$X6[1])
    obj$setParameter(x7m = data$X7[1])
    obj$setParameter(x8m = data$X8[1])
    obj$setParameter(x9m = data$X9[1])
    obj$setParameter(x10m = data$X10[1])
    obj$setParameter(x11m = data$X11[1])
    obj$setParameter(x12m = data$X12[1])
    obj$setParameter(x13m = data$X13[1])
    obj$setParameter(x14m = data$X14[1])
    obj$setParameter(x15m = data$X15[1])

    return(obj)
}


###########################################################
# Make data
############################################################

make_data<-function(){
    water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
    inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
    vars<-read.csv('data/processed/dronninglund/variables.csv')
    ambientTemp<-vars$temp_dry
    # inputs[,8:13] <- inputs[,8:13] # Convert to m3/hr
    # inputs[,2:4] <- inputs[,2:4]  # Convert to m3/hr
    data <- cbind(water_sensors, inputs,ambientTemp)
    data <- data[, !duplicated(colnames(data))]
    # dropna
    data <- data[complete.cases(data),]

    return(data)
}
data<-make_data()
df<-5
bs<-bs(1:16,df=df,intercept = TRUE)
for (i in 1:df){
    for (j in 0:15){
        data[,paste0('b',i,'x',j)]<-bs[(j+1),i]
    }
}

interval<-8500:9000
# Make model object
load('models/oldCTSMfit/model3_splines.RData')
model <- make_model(data = data[interval,],old_fit = fit)


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
#   if(printit){ print(format(xm,digits=2)) }

  fit$xm <- xm
  ## loglikelihood initialization
  nll <- 0
  ## use predict to compute the likelihood using the parameters in xm.
  Pred <- try(predict(fit, newdata=D, firstorderinputinterpolation=firstorder, n.ahead=n.ahead))
  if(class(Pred) == "try-error"){ return(NA) }
  ## add logligelihood for all outputs individually
  out_pred <- Pred$output$pred 
  nm <- names(out_pred)
    yhat <- out_pred
    y <- D[,nm]
return (mean((y - yhat)^2))

#   for(i in 1:ncol(out_pred)){
#     nm <- names(out_pred)[i]
#     yhat <- out_pred[ ,nm]
#     sd <- Pred$output$sd[ ,nm]
#     y <- D[ ,nm]
#     ## c is hubers psi. when the normalized residual is larger than c or smaller than
#     ## minus c, the loglikelihood continues as a square root (for robustness)
#     ressq <- (y-yhat)^2 / sd^2
#     ressq[ressq>c^2] <- c*(2*sqrt(ressq[ressq>c^2])-c)
#     nll <- nll + 0.5*( sum(log(Pred$output$sd^2)+ressq)  )
#   }
# #   print()
#     l <- ncol(out_pred) ## dimension of output vector
#     N <- nrow(out_pred) ## number of observations
# #   nll <- nll + 0.5*l*N*log(2*pi)
#     print(paste("Neg. loglikelihood:", print(nll)))



    # return(nll)
    }

fit <- makefit(model)


nlminb(start = fit$xm, objective = nllikelihood, 
                            #     lower = model$ParameterValues$lower,
                            #  upper = model$ParameterValues$upper,
                             fit=fit,
                             D= data[interval,],
                             firstorder=FALSE,
                             c=3, n.ahead=1, printit=TRUE,
                    control=list(iter.max=10))



  fit$xm <- xm
  ## loglikelihood initialization
  nll <- 0
  ## use predict to compute the likelihood using the parameters in xm.
  Pred <- try(predict(fit, newdata=D, firstorderinputinterpolation=TRUE, n.ahead=1))
  if(class(Pred) == "try-error"){ return(NA) }
  ## add logligelihood for all outputs individually
  out_pred <- Pred$output$pred 
  nm <- names(out_pred)
    yhat <- out_pred
    y <- D[,nm]

which(is.na((yhat-y)^2))
mean((yhat-y)^2)
head(mean(head(y - yhat)))
mean(head(yhat))
mean(yhat$X10)
l <- (yhat-y)^2
mean(l[,])
mean(unname(l))
