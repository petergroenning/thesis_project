library(ctsmr)


data <- read.csv('data/processed/allsensors.csv')
data <- data[3000:7000,]
model <- ctsm$new()

model$addObs(X00~x1)
model$addObs(X02~x1+x2)
model$addObs(X16~x1+x2+x3)
model$addObs(X21~x4)
model$addObs(X26~x4+x5)
model$addObs(X30~x6)

model$setVariance(X00 ~ sigma_X)
model$setVariance(X02 ~ sigma_X)
model$setVariance(X16 ~ sigma_X)
model$setVariance(X21 ~ sigma_X)
model$setVariance(X26 ~ sigma_X)
model$setVariance(X30 ~ sigma_X)
model$setParameter(sigma_X = (1.25e-4)/6)



model$addInput('FbotIn','FbotOut','FmidIn','FtopIn','FtopOut','Ttop','Tmid','Tbot')

model$addSystem(dx1~dt/1000*((Tbot-x1)*FbotIn*exp(vbot)+x2*FbotOut*exp(f1))+dwx1*exp(sigma_x))
model$addSystem(dx2~dt/1000*(x3*FbotOut*exp(f2)-x2*FbotOut*exp(f1)-(Tbot-x1)*FbotIn*exp(vbot))+dwx2*exp(sigma_x))
model$addSystem(dx3~dt/1000*(x2*FbotOut*exp(f2)+(Tmid-x4)*FmidIn*exp(vmid))+dwx3*exp(sigma_x))
model$addSystem(dx4~dt/1000*((Tmid-x4)*FmidIn*exp(vmid)-x3*FbotOut*exp(f2)+x5*FtopIn*exp(f4))+dwx4*exp(sigma_x))
model$addSystem(dx5~dt/1000*((Ttop-x6)*FtopIn*exp(vtop)-(Tmid-x4)*FmidIn*exp(vmid))+dwx5*exp(sigma_x))
model$addSystem(dx6~dt/1000*((Ttop-x6)*FtopIn*exp(vtop)-x5*FtopOut*exp(f5))+dwx6*exp(sigma_x))



model$setParameter(vbot = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(vmid = log(c(init = 0.001, lb = 1e-12, ub = 10)))
model$setParameter(vtop = log(c(init = 0.001, lb = 1e-12, ub = 10)))
model$setParameter(sigma_x = log(c(init = 1, lb = 1e-12, ub = 2)))
model$setParameter(f1 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f2 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f4 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f5 = log(c(init = 0.1, lb = 1e-12, ub = 1)))



model$setParameter(x1 = data$X00[1])
model$setParameter(x2 = data$X02[1]-data$X00[1])
model$setParameter(x3 = data$X16[1]-data$X02[1]-data$X00[1])
model$setParameter(x4 = data$X21[1])
model$setParameter(x5 = data$X26[1]-data$X21[1])
model$setParameter(x6 = data$X31[1])
fit <- model$estimate(data, firstorder=TRUE)



p <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
px0 <- p$output$pred$X00
px02 <- p$output$pred$X02
px30 <- p$output$pred$X30
rx0 <- px0 - data$X00
rx1 <- px02 - data$X02
rx2 <- px30 - data$X30
plot(rx0[2:3000])
acf(rx0[2:3000])

plot(px0[200:300])
lines(data$X00[200:300], col='red')
