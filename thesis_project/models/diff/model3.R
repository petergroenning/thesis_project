library(ctsmr)


data <- read.csv('data/processed/data2.csv')
data <- data[3000:4000,]
model <- ctsm$new()

model$addObs(X1 ~ x1) # Input Layer (3m)
model$addObs(X2 ~ x2) # (3m)
model$addObs(X3 ~ x3) # (3m)
model$addObs(X4 ~ x4) # Input Layer (3m)
model$addObs(X5 ~ x5) # (3m)
model$addObs(X6 ~ x6) # Input Layer (1m)

model$setVariance(X1 ~ sigma_X^2)
model$setVariance(X2 ~ sigma_X^2)
model$setVariance(X3 ~ sigma_X^2)
model$setVariance(X4 ~ sigma_X^2)
model$setVariance(X5 ~ sigma_X^2)
model$setVariance(X6 ~ sigma_X^2)
model$setParameter(sigma_X = 2.083e-5)

model$setParameter(V1 = 3185)
model$setParameter(V2 = 5950)
model$setParameter(V3 = 9578)
model$setParameter(V4 = 14072)
model$setParameter(V5 = 19428)
model$setParameter(V6 = 7816)

model$addInput('FbotIn','FbotOut','FmidIn','FtopIn','FtopOut','Ttop','Tmid','Tbot')

model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*exp(vbot)+(x2-x1)*(FbotOut*exp(f1)))+dt*Y1+dw1*sigma_x)
model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f2))+(x3-x2)*(FbotOut*exp(f2)))+dt*Y1+dw2*sigma_x)
model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f3))+(x4-x3)*(FbotOut*exp(f3)))+dt*Y2+dw3*sigma_x)
model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f4))+(x5-x4)*(FbotOut*exp(f4))+(Tmid-x4)*FmidIn*exp(vmid))+dt*Y2+dw4*sigma_x)
model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*exp(f5))+(x6-x5)*(FtopIn*exp(f5)))+dw5*sigma_x)
model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*exp(f6))+(Ttop-x6)*FtopIn*exp(vtop))+dw6*sigma_x)


model$addSystem(dY1~dwY1*sigma_X)
model$addSystem(dY2 ~ dwY2*sigma_X)

# model$addSystem(dYtop ~dt *(FtopOut*exp(fin2)-FtopOut*exp(fout2)-Ytop*a) + dwY2*sigma_X)

model$setParameter(vbot = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(vmid = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(vtop = log(c(init = 0.1, lb = 1e-12, ub = 10)))

model$setParameter(sigma_x = log(c(init = 1, lb = 1e-12, ub = 2)))

model$setParameter(f1 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f2 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f3 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f4 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f5 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f6 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(Y1 = 0)
model$setParameter(Y2 = 0)

model$setParameter(a = log(c(init = 9, lb = 1e-12, ub = 10)))
model$setParameter(fin1 = log(c(init = 0.0000001, lb = 1e-12, ub = 1)))
model$setParameter(fout1 = log(c(init = 0.0000001, lb = 1e-12, ub = 1)))
model$setParameter(fin2 = log(c(init = 0.00001, lb = 1e-12, ub = 1)))
model$setParameter(fout2 = log(c(init = 0.000001, lb = 1e-12, ub = 1)))

model$setParameter(x1 = data$X1[1])
model$setParameter(x2 = data$X2[1])
model$setParameter(x3 = data$X3[1])
model$setParameter(x4 = data$X4[1])
model$setParameter(x5 = data$X5[1])
model$setParameter(x6 = data$X6[1])


setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)

    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]

    model$ParameterValues$initial <- newparams
    return(model)
}

load_model <- function(model_name){
    load(model_name)
    return(fit)
}
old_fit <- load_model('models/diff/model4.RData')
model <- setParams(model, fit)

fit <- model$estimate(data, firstorder=TRUE)


p <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
px1 <- p$output$pred$X1
px2 <- p$output$pred$X2
px3 <- p$output$pred$X3
px4 <- p$output$pred$X4
px5 <- p$output$pred$X5
px6 <- p$output$pred$X6

rx1 <- px1 - data$X1
rx2 <- px2 - data$X2
rx3 <- px3 - data$X3
rx4 <- px4 - data$X4
rx5 <- px5 - data$X5
rx6 <- px6 - data$X6

acf(rx3)
plot(rx2)
dir.create('models/diff', showWarnings = FALSE)
save(fit, file = 'models/diff/model5.RData')
