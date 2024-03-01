##-----------------------------------------------------------------
sapply(dir("functions", full.names=TRUE), source)
## Read data
Xtherm <- read.table("data/5MinTHERM1.csv",sep=",",stringsAsFactors=FALSE,header=TRUE)
Xtherm$timedate <- as.POSIXct(Xtherm$time)
Xtherm$time <- NULL
Xtherm$t = as.numeric(difftime(Xtherm$timedate, Xtherm$time[1], units="hours"))
Xtherm$Ph <- Xtherm$P/1000
plot(is.na(Xtherm$Ph))
Xtherm <- Xtherm[(max(which(is.na(Xtherm$Ph)))+100):nrow(Xtherm),]
## Average indoor temperature
iTemp <- grep("^T[0-7]",names(Xtherm))
Xtherm$yTi <- apply(Xtherm[,iTemp], 1, mean, na.rm=TRUE)
##
Xtherm$Ps <- Xtherm$G
Xtherm$Te <- Xtherm$Ta
##
iTemp <- grep("yTi|Ph|Ps|Ta",names(Xtherm))
##plot.ts(Xtherm[,iTemp])
##-----------------------------------------------------------------


##-----------------------------------------------------------------
## The control function
ftherm <- rlang::expr((1/(1+exp(-1*(Tset-Ti)*slope))) * flowMax * cw )
## ## Plot the sigmoid function
## Ti <- seq(18,22, by=0.1)
## Tset <- 20
## slope <- 3
## cw <- 0.001118 # kWh/(l C) The specific heat capacity of water
## flowMax <- 180
## plot(Ti, eval(ftherm, list(Tset, Ti, slope, cw, flowMax)))
##-----------------------------------------------------------------


##-----------------------------------------------------------------
## Setup the data
library(ctsmr)
## Resample
X <- Xtherm #resample.data.frame(Xtherm, 5*60, Xtherm$timedate[1], timename="timedate")
X$Tfor <- 60
##X$Tfor[floor(nrow(X)/2):nrow(X)] <- 30
X$Tret <- NA#25
X$Tset <- 20
X$Tset[floor(nrow(X)/2):nrow(X)] <- 23
X$thermvent <- 1
##-----------------------------------------------------------------


##-----------------------------------------------------------------
dev.new()
model <- setprm(TiTh_comp1())
fit <- makefit(model)
sim <- simulate(fit, newdata=X)
##
plotsim(sim, X)

dev.new()
model <- TiTh_comp2()
model <- setprm(model)
fit <- makefit(model)

## Same two calls
simulate(fit, newdata=X)
val <- ctsmr:::simulate.ctsmr(fit, newdata=X)

## Plot
plotsim(val, X)

dev.new()
## In one line
plotsim(simulate(makefit(setprm(TiTh_comp1())), newdata=X), X)
plotsim(simulate(makefit(setprm(TiTh_comp3())), newdata=X), X)
##-----------------------------------------------------------------



## ##-----------------------------------------------------------------
## ## Estimate the sigmoid slope


## ## Fit with nlminb
## model$AnalyseModel()
## ## Make a fit object to pass to simulate.ctsmr
## fit <- list()
## fit$model <- model
## fit$xm <- model$ParameterValues$initial
## names(fit$xm) <- row.names(model$ParameterValues)
## ## Can be nessary such that xm has the right length
## fit$xm <- fit$xm[model$pars]
## class(fit) <- "ctsmr"

## ## Resample
## Xtherm <- resample.data.frame(XthermOr, 30*60, Xtherm$timedate[1], timename="timedate")

## ## The likelihood of the initial parameters set in fit
## nllikelihood(fit$xm[model$pars], fit, D=Xtherm, firstorder=TRUE, c=3, n.ahead=1)

## ## ------------------------------------------------------------------
## ## Fit with different optimizers

## ## Fit with ctsmr
## fitctsmr <- model$estimate(Xtherm, firstorder=TRUE)

## fitnlminb <- nlminb(fitctsmr$xm, nllikelihood, fit=fit, D=Xtherm, firstorder=TRUE, c=3, n.ahead=1,
##                     control=list(iter.max=400))

## fitnlminb

## ## Plot the simulated output from the three estimated models
## insert_xm <- function(fit, xm){
##     fit$xm <- xm
##     names(fit$xm) <- names(xm)
##     fit$xm <- fit$xm[fit$model$pars]
##     return(fit)
## }

## ## ctsmr
## (fitCtsmr <- insert_xm(fit, fitctsmr$xm))
## simCtsmr <- simulate(fitCtsmr, newdata=Xtherm)
## ## nlminb
## (fitNlminb <- insert_xm(fit, fitnlminb$par))
## simNlminb <- simulate(fitNlminb, newdata=Xtherm)


## nllikelihood(fitCtsmr$xm[model$pars], fitCtsmr, D=Xtherm, firstorder=TRUE, c=3, n.ahead=1)
## nllikelihood(fitNlminb$xm[model$pars], fitNlminb, D=Xtherm, firstorder=TRUE, c=3, n.ahead=1)


## ##-----------------------------------------------------------------
## setpar("ts", mfrow=c(5,1))
## X <- Xtherm

## ## One-step predictions
## X$residualCtsmr <- X$Ph - ctsmr:::predict.ctsmr(fitCtsmr, newdata=X)$output$pred$Ph
## X$residualNlminb <- X$Ph - ctsmr:::predict.ctsmr(fitNlminb, newdata=X)$output$pred$Ph

        
## plot(X$timedate, X$Ph, type="l", ylab="Heat power (kW)")
## lines(X$timedate, simCtsmr$output$sim$Ph, type="l", col=2)
## lines(X$timedate, simNlminb$output$sim$Ph, type="l", col=3)
## legend("topright", c("Heat power","Predicted heat power"), lty=1, col=c(1,2), bg="white")
## plot(X$timedate, X$yTi, type="l", ylab="Indoor temperature (C)")
## plot(X$timedate, X$Te, type="l", ylab="Outdoor temperature (C)")
## plot(X$timedate, X$Ps, type="l", ylab="Solar radiation (W/m2)")
## plotTSXAxis(X$timedate, format="%m-%d %H:%M")
## plot(X$timedate, X$residualCtsmr, type="l", ylab="one-step residuals", col=2)
## lines(X$timedate, X$residualNlminb, type="l", ylab="one-step residuals", col=3)

## ##-----------------------------------------------------------------


