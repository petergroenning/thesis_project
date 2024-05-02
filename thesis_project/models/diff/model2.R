library(ctsmr)


data <- read.csv('data/processed/diff.csv')
data <- data[3000:3500,]
model <- ctsm$new()

names(data)
layers <- c('X3','X4','X5')
for (i in 3:5){

    model$addObs(as.formula(paste0(layers[i-2],'~x',i)))
    model$setVariance(as.formula(paste0(layers[i-2],'~sigma_X')))
}
model$addInput('FbotIn','FbotOut','FmidIn','FtopIn','FtopOut','Ttop','Tmid','Tbot','X2','X6')
model$setParameter(sigma_X = (1.25e-4)/6)
model$setParameter(sigma_x = log(c(init = 1, lb = 1e-12, ub = 2)))


model$addSystem(dx3~dt/1000*(-x3*FbotIn*exp(f1)+X2*FbotOut*exp(f1))+dwx3*exp(sigma_x))
model$addSystem(dx4~dt/1000*(-x4*FbotIn*exp(f2)+x3*FbotOut*exp(f2))+dwx4*exp(sigma_x))
model$addSystem(dx5~dt/1000*(-x5*FbotIn*exp(f3)+x4*FbotOut*exp(f3))+dwx5*exp(sigma_x))






model$setParameter(vbot = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(vmid = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(vtop = log(c(init = 0.1, lb = 1e-12, ub = 10)))
model$setParameter(f1 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f2 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f3 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f4 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f5 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f6 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f7 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f8 = log(c(init = 0.1, lb = 1e-12, ub = 1)))
model$setParameter(f9 = log(c(init = 0.1, lb = 1e-12, ub = 1)))


# model$setParameter(x1 = data$X1[1])
# model$setParameter(x2 = data$X2[1])
model$setParameter(x3 = data$X3[1])
model$setParameter(x4 = data$X4[1])
model$setParameter(x5 = data$X5[1])
# model$setParameter(x6 = data$X6[1])
# model$setParameter(x7 = data$X7[1])
# model$setParameter(x8 = data$X8[1])
# model$setParameter(x9 = data$X9[1])



fit <- model$estimate(data, firstorder=TRUE)

p <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
px3 <- p$output$pred$X3
px4 <- p$output$pred$X4
px5 <- p$output$pred$X5

rx3 <- px3 - data$X3
rx4 <- px4 - data$X4
rx5 <- px5 - data$X5

acf(rx5)


