data <- read.csv('data/processed/data2.csv')


x <- data$FbotOut[1:8000]
y <- data$X1[1:8000]
r1 <- read.csv('residuals/simple3/nonlinear.csv')$X3
r1 <- r[!is.na(r1)]

r2 <- read.csv('residuals/simple3/model3.csv')$X3
r2 <- r2[!is.na(r2)]

modelx <- arima(x, order = c(3,1,3), seasonal = list(order = c(2,0,1), period = 24))

# Filter y with arima model



modely <- arima(r1, order = c(3,1,3), seasonal = list(order = c(2,0,1), period = 24),
                fixed = coef(modelx))


acf(residuals(modely), na.action = na.pass, lag.max = 100)
pacf(residuals(modely), na.action = na.pass, lag.max = 100)


plot((r1))
acf(r2, na.action = na.pass, lag.max = 100)
plot(modely$residuals)
plot(model$residuals)

p <- model$residuals + x

plot(x[1800:2000], type = 'l', ylim = c(-20,200))
lines(p[1800:2000], col = 'red')

acf(p)
