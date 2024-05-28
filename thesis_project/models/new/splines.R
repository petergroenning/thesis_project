library(splines)

data <- read.csv('data/processed/1m_2_data.csv')

# idx <- 3500
# d <- t(data[idx,2:33])
# tmin <- 5
# tmax <- 90
# d <- (d-tmin)/(tmax-tmin)
# x <- (1:32)/32


# fit <- lm(d~bs(x, knots = c(0.2,0.3,1), degree = 1, intercept = TRUE))
# fit$effects
# fit$residuals
# plot(fit$residuals)
# plot(x,d)
# lines(x,(fitted(fit)))

# mean(fit$residuals^2)
# b <- bs(x, knots = c(0.2), degree = 3)

# plot(x,b[,1], type = 'l', ylim = c(0,1), lwd = 5)
# for (i in 2:4){
#     lines(x,b[,i], type = 'l', col = i, lw = 5)
# }




# fit$coefficients[4]



data$dX3 <- c(0,diff(data$X3))
data$dX5 <- c(0,diff(data$X5))
data$dX5k1 <- c(0,0,data$dX5[1:(nrow(data)-2)])

data$FbotOutk1 <- c(0,data$FbotOut[1:(nrow(data)-1)])

fit <- lm(dX3~(X4-X3)*Fbot+Fbot*(X2-X3), data = data)
summary(fit)

fit <- lm(dX5~(X4-X5)*(FbotIn)+(FbotOut)*(X6-X5)+dX5k1, data = data)
summary(fit)



plot(data$dX5[2000:2100])
lines(fitted(fit)[2000:2100], col = 'red')

plot(data$X5[1:1500])
lines(data$X5[1]+cumsum(fitted(fit)[1:1500]), col = 'red')


