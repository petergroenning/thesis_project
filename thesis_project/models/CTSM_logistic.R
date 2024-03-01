library(ctsmr)
source('thesis_project/models/ctsm_utils.R')
library(ggplot2)


make_model <- function(name, data, old_fit){

    obj <- ctsm$new()
    if (missing(old_fit)){
        old_fit <- NULL
    }

    # Add observation equations and variances
    obj$addObs(X0 ~ T1)
    obj$addObs(X15 ~ T1 + T2 + T3)

    obj$addObs(X1 ~ T1 + T2/(1 + exp((a1 - 1/15)/s1)) + T3/(1 + exp((a2 - 1/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X2 ~ T1 + T2/(1 + exp((a1 - 2/15)/s1)) + T3/(1 + exp((a2 - 2/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X3 ~ T1 + T2/(1 + exp((a1 - 3/15)/s1)) + T3/(1 + exp((a2 - 3/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X4 ~ T1 + T2/(1 + exp((a1 - 4/15)/s1)) + T3/(1 + exp((a2 - 4/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X5 ~ T1 + T2/(1 + exp((a1 - 5/15)/s1)) + T3/(1 + exp((a2 - 5/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X6 ~ T1 + T2/(1 + exp((a1 - 6/15)/s1)) + T3/(1 + exp((a2 - 6/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X7 ~ T1 + T2/(1 + exp((a1 - 7/15)/s1)) + T3/(1 + exp((a2 - 7/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X8 ~ T1 + T2/(1 + exp((a1 - 8/15)/s1)) + T3/(1 + exp((a2 - 8/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X9 ~ T1 + T2/(1 + exp((a1 - 9/15)/s1)) + T3/(1 + exp((a2 - 9/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X10 ~ T1 + T2/(1 + exp((a1 - 10/15)/s1)) + T3/(1 + exp((a2 - 10/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X11 ~ T1 + T2/(1 + exp((a1 - 11/15)/s1)) + T3/(1 + exp((a2 - 11/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X12 ~ T1 + T2/(1 + exp((a1 - 12/15)/s1)) + T3/(1 + exp((a2 - 12/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X13 ~ T1 + T2/(1 + exp((a1 - 13/15)/s1)) + T3/(1 + exp((a2 - 13/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    obj$addObs(X14 ~ T1 + T2/(1 + exp((a1 - 14/15)/s1)) + T3/(1 + exp((a2 - 14/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))
    # obj$addObs(X15 ~ T1 + T2/(1 + exp((a1 - 15/15)/s1)) + T3/(1 + exp((a2 - 15/15)/s2)) - T2/(1 + exp(a1/s1)) - T3/(1 + exp(a2/s2)))


    # Set observation equation variances
    obj$setVariance(X0 ~ sigma_X0^2)
    obj$setVariance(X15 ~ sigma_X15^2)
    obj$setVariance(X10 ~ sigma_X10^2)


    ###### SET PARAMETERS ######
    # Observation Noise
    obs_var <- 1e-3
    obj$setParameter(sigma_X0 = obs_var)
    obj$setParameter(sigma_X15 = obs_var)
    obj$setParameter(sigma_X10 = obs_var)


    # Add system equations
    obj$addSystem(T1 ~ sigma_T1 * dw1)
    obj$addSystem(T2 ~ sigma_T2 * dw2)
    obj$addSystem(T3 ~ sigma_T3 * dw3)
    obj$addSystem(a1 ~ FtopIn * c1 * dt + sigma_a1 * dw4)
    obj$addSystem(s1 ~ FtopIn * c2 * dt + sigma_s1 * dw5)
    obj$addSystem(a2 ~ FtopIn * c3 * dt + sigma_a2 * dw6)
    obj$addSystem(s2 ~ FtopIn * c4 * dt + sigma_s2 * dw7)


    # Add input     
    obj$addInput('FtopIn')

    
    # System Noise
    lb <- 0
    ub <- 5
    init <- 0.0001
    obj$setParameter(sigma_T1 = c(init, lower = lb, upper = ub))
    obj$setParameter(sigma_T2 = c(init, lower = lb, upper = ub))
    obj$setParameter(sigma_T3 = c(init, lower = lb, upper = ub))
    obj$setParameter(sigma_a1 = c(0.001, lower = lb, upper = 2))
    obj$setParameter(sigma_s1 = c(0.001, lower = lb, upper = 2))
    obj$setParameter(sigma_a2 = c(0.001, lower = lb, upper = 2))
    obj$setParameter(sigma_s2 = c(0.001, lower = lb, upper = 2))
    obj$setParameter(c1 = c(0, lower = -1, upper = 1))
    obj$setParameter(c2 = c(0, lower = -1, upper = 1))
    obj$setParameter(c3 = c(0, lower = -1, upper = 1))
    obj$setParameter(c4 = c(0, lower = -1, upper = 1))

    # Set initial state
    obj$setParameter(T1 = data$X0[1])
    obj$setParameter(T2 = data$X15[1])
    obj$setParameter(T3 = data$X15[1] - data$X0[1])
    obj$setParameter(a1 = c(0.5, lower = 0, upper = 1))
    obj$setParameter(s1 = c(5, lower = 0, upper = 10))
    obj$setParameter(a2 = c(0.5, lower = 0, upper = 1))
    obj$setParameter(s2 = c(5, lower = 0, upper = 10))




    return(obj)
}


###########################################################
# Make data
############################################################
water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
inputs[,8:13] <- inputs[,8:13] / 1000 # Convert to m3/hr

data <- cbind(water_sensors, inputs)[7200:7500,]
data <- data[, !duplicated(colnames(data))]
# dropna
data <- data[complete.cases(data),]


# Make model object
obj <- make_model(name = name, data = data)

fit <- obj$estimate(data)

obj$ParameterValues
summary(fit)

preds <- predict(fit, newdata = data, n.ahead = 10)
plot(preds$state$pred$T1, type = 'l', col = 'blue', lwd = 2)
points(data$X0, col = 'red')

ttop <- preds$state$pred$T3 + preds$state$pred$T2 + preds$state$pred$T1
plot(ttop, type = 'l', col = 'blue', lwd = 1)
points(data$X15, col = 'red')

plot(preds$state$pred$a1, type = 'l', col = 'blue', lwd = 1)
plot(preds$state$pred$a2, type = 'l', col = 'blue', lwd = 1)

plot(preds$state$pred$s1, type = 'l', col = 'blue', lwd = 1)
plot(preds$state$pred$s2, type = 'l', col = 'blue', lwd = 1)

a1 <- preds$state$pred$a1
s1 <- preds$state$pred$s1
t1 <- preds$state$pred$T1
t2 <- preds$state$pred$T2
t3 <- preds$state$pred$T3
x <- seq(0, 1, length.out = 16)
i <- 10
y <- c(data[i, 'X0'], data[i, 'X1'], data[i, 'X2'], data[i, 'X3'], data[i, 'X4'], data[i, 'X5'], data[i, 'X6'], data[i, 'X7'], data[i, 'X8'], data[i, 'X9'], data[i, 'X10'], data[i, 'X11'], data[i, 'X12'], data[i, 'X13'], data[i, 'X14'], data[i, 'X15'])

plot(x, t1[i] + t2[i] / (1 + exp((a1[i] - x)/s1[i])) + t2[i] / (1 + exp((a1[i] - x)/s1[i])) - 1/(1 + exp(a1[i]/s1[i])) - 1/(1 + exp(a1[i]/s1[i])), type = 'l', col = 'blue', lwd = 1)
points(seq(0,1, length.out = 16),y)

