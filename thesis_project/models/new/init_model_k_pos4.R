library(ctsmr)

data <- read.csv('data/processed/1m_data.csv')
model <- ctsm$new()

    #### Physical Parameters #####
    # VOLUMES
    model$setParameter(V0 = 808)
    model$setParameter(V1 = 1051)
    model$setParameter(V2 = 1326)
    model$setParameter(V3 = 1633)
    model$setParameter(V4 = 1973)
    model$setParameter(V5 = 2344)
    model$setParameter(V6 = 2747)
    model$setParameter(V7 = 3182)
    model$setParameter(V8 = 3649)
    model$setParameter(V9 = 4149)
    model$setParameter(V10 = 4680)
    model$setParameter(V11 = 5243)
    model$setParameter(V12 = 5838)
    model$setParameter(V13 = 6465)
    model$setParameter(V14 = 7125)
    model$setParameter(V15 = 7816)

    # Add observation equations and variances
    model$addObs(X0 ~ x0m)
    model$addObs(X1 ~ x1m)
    model$addObs(X2 ~ x2m)
    model$addObs(X3 ~ x3m)
    model$addObs(X4 ~ x4m)
    model$addObs(X5 ~ x5m)
    model$addObs(X6 ~ x6m)
    model$addObs(X7 ~ x7m)
    model$addObs(X8 ~ x8m)
    model$addObs(X9 ~ x9m)
    model$addObs(X10 ~ x10m)
    model$addObs(X11 ~ x11m)
    model$addObs(X12 ~ x12m)
    model$addObs(X13 ~ x13m)
    model$addObs(X14 ~ x14m)
    model$addObs(X15 ~ x15m)

    # Set observation equation variances
    model$setVariance(X0 ~ sigma_X^2)
    model$setVariance(X1 ~ sigma_X^2)
    model$setVariance(X2 ~ sigma_X^2)
    model$setVariance(X3 ~ sigma_X^2)
    model$setVariance(X4 ~ sigma_X^2)
    model$setVariance(X5 ~ sigma_X^2)
    model$setVariance(X6 ~ sigma_X^2)
    model$setVariance(X7 ~ sigma_X^2)
    model$setVariance(X8 ~ sigma_X^2)
    model$setVariance(X9 ~ sigma_X^2)
    model$setVariance(X10 ~ sigma_X^2)
    model$setVariance(X11 ~ sigma_X^2)
    model$setVariance(X12 ~ sigma_X^2)
    model$setVariance(X13 ~ sigma_X^2)
    model$setVariance(X14 ~ sigma_X^2)
    model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    obs_std <- 1/sqrt(12)*1.25e-4
    model$setParameter(sigma_X = obs_std)

    model$addSystem(dx0m~dt/V0*((x1m-x0m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(Tbot-x0m)*(FbotIn*exp(vbot)))+exp(sigma_x)*dwx0)
    model$addSystem(dx1m~dt/V1*((x2m-x1m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x0m-x1m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx1)
    model$addSystem(dx2m~dt/V2*((x3m-x2m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x1m-x2m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx2)
    model$addSystem(dx3m~dt/V3*((x4m-x3m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x2m-x3m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx3)
    model$addSystem(dx4m~dt/V4*((x5m-x4m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x3m-x4m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx4)
    model$addSystem(dx5m~dt/V5*((x6m-x5m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x4m-x5m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx5)
    model$addSystem(dx6m~dt/V6*((x7m-x6m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x5m-x6m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx6)
    model$addSystem(dx7m~dt/V7*((x8m-x7m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x6m-x7m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx7)
    model$addSystem(dx8m~dt/V8*((x9m-x8m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x7m-x8m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx8)
    model$addSystem(dx9m~dt/V9*((x10m-x9m)*((FbotOut+exp(l)*FbotOutk1+exp(l)^2*FbotOutk2+exp(l)^3*FbotOutk3)*exp(f)+k)+(x8m-x9m)*((FbotIn+exp(l)*FbotInk1+exp(l)^2*FbotInk2+exp(l)^3*FbotInk3)*exp(f)+k))+exp(sigma_x)*dwx9)
    model$addSystem(dx10m~dt/V10*((x11m-x10m)*((FtopIn+exp(l)*FtopInk1+exp(l)^2*FtopInk2+exp(l)^3*FtopInk3)*exp(f)+k)+(x9m-x10m)*((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3)*exp(f)+k)+(Tmid-x10m)*FmidIn*exp(vmid))+exp(sigma_x)*dwx10)
    model$addSystem(dx11m~dt/V11*((x12m-x11m)*((FtopIn+exp(l)*FtopInk1+exp(l)^2*FtopInk2+exp(l)^3*FtopInk3)*exp(f)+k)+(x10m-x11m)*((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3)*exp(f)+k))+exp(sigma_x)*dwx11)
    model$addSystem(dx12m~dt/V12*((x13m-x12m)*((FtopIn+exp(l)*FtopInk1+exp(l)^2*FtopInk2+exp(l)^3*FtopInk3)*exp(f)+k)+(x11m-x12m)*((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3)*exp(f)+k))+exp(sigma_x)*dwx12)
    model$addSystem(dx13m~dt/V13*((x14m-x13m)*((FtopIn+exp(l)*FtopInk1+exp(l)^2*FtopInk2+exp(l)^3*FtopInk3)*exp(f)+k)+(x12m-x13m)*((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3)*exp(f)+k))+exp(sigma_x)*dwx13)
    model$addSystem(dx14m~dt/V14*((x15m-x14m)*((FtopIn+exp(l)*FtopInk1+exp(l)^2*FtopInk2+exp(l)^3*FtopInk3)*exp(f)+k)+(x13m-x14m)*((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3)*exp(f)+k))+exp(sigma_x)*dwx14)
    model$addSystem(dx15m~dt/V15*((x14m-x15m)*(((FtopOut+exp(l)*FtopOutk1+exp(l)^2*FtopOutk2+exp(l)^3*FtopOutk3))*exp(f)+k)+(Ttop-x15m)*(FtopIn*exp(vtop)))+exp(sigma_x)*dwx15)





    # Inputs
    model$addInput("Ttop", "Tbot", "Tmid", "FtopIn", "FmidIn", "FbotIn", "FtopOut", "FbotOut","FmidOut", "ambientTemp", "Ftop", "Fbot", "Fmid")
    data$FbotInk1 <- c(0, data$FbotIn[1:(nrow(data)-1)])
    data$FbotInk2 <- c(0,0, data$FbotIn[1:(nrow(data)-2)])
    data$FbotInk3 <- c(0,0,0, data$FbotIn[1:(nrow(data)-3)])

    data$FbotOutk1 <- c(0, data$FbotOut[1:(nrow(data)-1)])
    data$FbotOutk2 <- c(0,0, data$FbotOut[1:(nrow(data)-2)])
    data$FbotOutk3 <- c(0,0,0, data$FbotOut[1:(nrow(data)-3)])

    data$FtopInk1 <- c(0, data$FtopIn[1:(nrow(data)-1)])
    data$FtopInk2 <- c(0,0, data$FtopIn[1:(nrow(data)-2)])
    data$FtopInk3 <- c(0,0,0, data$FtopIn[1:(nrow(data)-3)])

    data$FtopOutk1 <- c(0, data$FtopOut[1:(nrow(data)-1)])
    data$FtopOutk2 <- c(0,0, data$FtopOut[1:(nrow(data)-2)])
    data$FtopOutk3 <- c(0,0,0, data$FtopOut[1:(nrow(data)-3)])

    model$addInput('FbotInk1', 'FbotInk2', 'FbotInk3', 'FbotOutk1', 'FbotOutk2', 'FbotOutk3', 'FtopInk1', 'FtopInk2', 'FtopInk3', 'FtopOutk1', 'FtopOutk2', 'FtopOutk3')

    model$setParameter(sigma_x=log(c(init=0.1,lb=1e-12,ub=1)))

    model$setParameter(f=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(k=(c(init=1,lb=-100,ub=100)))
    model$setParameter(Cp=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(v=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(vtop=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(vmid=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(vbot=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(l = log(c(init=0.5,lb=1e-12,ub=1)))

setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)
    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]
    model$ParameterValues$initial <- newparams
    return(model)
}
setInitialState <- function(model, data){
    model$setParameter(x0m = data$X0[1])
    model$setParameter(x1m = data$X1[1])
    model$setParameter(x2m = data$X2[1])
    model$setParameter(x3m = data$X3[1])
    model$setParameter(x4m = data$X4[1])
    model$setParameter(x5m = data$X5[1])
    model$setParameter(x6m = data$X6[1])
    model$setParameter(x7m = data$X7[1])
    model$setParameter(x8m = data$X8[1])
    model$setParameter(x9m = data$X9[1])
    model$setParameter(x10m = data$X10[1])
    model$setParameter(x11m = data$X11[1])
    model$setParameter(x12m = data$X12[1])
    model$setParameter(x13m = data$X13[1])
    model$setParameter(x14m = data$X14[1])
    model$setParameter(x15m = data$X15[1])
    return(model)
}

data <- data[3000:7500,]
model <- setInitialState(model, data)

fit <- model$estimate(data, firstorder = TRUE)

p <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
r <- p$output$pred - data[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
write.csv(r, 'residuals/base_model_k_pos4.csv')
fit$loglik
aic <- -2*fit$loglik + 2*(6)
aic


# p2 <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 24*2*7)


# save(fit, file = 'models/base_model_k_pos.Rdata')
# load('models/base_model_k_pos3.Rdata')
# s <- simulate(fit, newdata = data, firstorderinputinterpolation=TRUE)
# plot(s$output$sim$X15, type = 'l')
# lines(data$X15, col = 'red')
