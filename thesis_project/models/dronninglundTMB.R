library(ctsmrTMB)
library(ggplot2)
library(patchwork)
rm(list = ls())
source('thesis_project/models/ctsm_utils.R')


###########################################################
# Make data
############################################################
water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
head(water_sensors)
data <- cbind(water_sensors, inputs)[5000:5500,]
data <- data[, !duplicated(colnames(data))]

# data[,2:17] <- data[,2:17] + 10
# data[,22:24] <- data[,22:24] + 10
#### MAKE CTSM OBJECT ####
name <- 'dronninglundTMB'
old_fit <- load_best_fit(name)
obj = ctsmrTMB$new()
# Set name of model (and the created .cpp file)
obj$set_modelname(name)

# Set path where generated C++ files are saved.
# This will create a cppfiles folder in your current working directory if it doesnt exist
obj$set_cppfile_directory(".cppfiles")

# Add observation equations and variances
obj$add_observations(X0 ~ x0)
obj$add_observations(X1 ~ x1)
obj$add_observations(X2 ~ x2)
obj$add_observations(X3 ~ x3)
obj$add_observations(X4 ~ x4)
obj$add_observations(X5 ~ x5)
obj$add_observations(X6 ~ x6)
obj$add_observations(X7 ~ x7)
obj$add_observations(X8 ~ x8)
obj$add_observations(X9 ~ x9)
obj$add_observations(X10 ~ x10)
obj$add_observations(X11 ~ x11)
obj$add_observations(X12 ~ x12)
obj$add_observations(X13 ~ x13)
obj$add_observations(X14 ~ x14)
obj$add_observations(X15 ~ x15)

# Add system equations
obj$add_systems(dx0 ~ ((x1 - x0)    * taux1 + (Tbot - x0) * Fbot * taubot - xloss * x0) * dt * c0 + sigma_x0 * dw0)
obj$add_systems(dx1 ~ ((x0 - x1)    * taux1 + (x2 - x1) * taux2 - x1loss * x1) * dt * c1+ sigma_x1 * dw1)
obj$add_systems(dx2 ~ ((x1 - x2)    * taux2 + (x3 - x2) * taux3 - x2loss * x2) * dt * c2+ sigma_x2 * dw2)
obj$add_systems(dx3 ~ ((x2 - x3)    * taux3 + (x4 - x3) * taux4 - x3loss * x3) * dt * c3+ sigma_x3 * dw3)
obj$add_systems(dx4 ~ ((x3 - x4)    * taux4 + (x5 - x4) * taux5 - x4loss * x4) * dt * c4+ sigma_x4 * dw4)
obj$add_systems(dx5 ~ ((x4 - x5)    * taux5 + (x6 - x5) * taux6 - x5loss * x5) * dt * c5+ sigma_x5 * dw5)
obj$add_systems(dx6 ~ ((x5 - x6)    * taux6 + (x7 - x6) * taux7 - x6loss * x6) * dt * c6+ sigma_x6 * dw6)
obj$add_systems(dx7 ~ ((x6 - x7)    * taux7 + (x8 - x7) * taux8 - x7loss * x7) * dt * c7+ sigma_x7 * dw7)
obj$add_systems(dx8 ~ ((x7 - x8)    * taux8 + (x9 - x8) * taux9 - x8loss * x8) * dt * c8+ sigma_x8 * dw8)
obj$add_systems(dx9 ~ ((x8 - x9)    * taux9 + (x10 - x9) * taux10 +  (Tmid - x9) * Fmid * taumid  - x9loss * x9) * dt* c9 + sigma_x9 * dw9)
obj$add_systems(dx10 ~ ((x9 - x10)  * taux10 + (x11 - x10) * taux11 - x10loss * x10) * dt * c10+ sigma_x10 * dw10)
obj$add_systems(dx11 ~ ((x10 - x11) * taux11 +  (x12 - x11) * taux12 - x11loss * x11) * dt* c11 + sigma_x11 * dw11)
obj$add_systems(dx12 ~ ((x11 - x12) * taux12 + (x13 - x12) * taux13 - x12loss * x12) * dt* c12 + sigma_x12 * dw12)
obj$add_systems(dx13 ~ ((x12 - x13) * taux13 + (x14 - x13) * taux14 - x13loss * x13) * dt* c13 + sigma_x13 * dw13)
obj$add_systems(dx14 ~ ((x13 - x14) * taux14 + (x15 - x14) * taux15 - x14loss * x14) * dt* c14 + sigma_x14 * dw14)
obj$add_systems(dx15 ~ ((x14 - x15) * taux15 + (Ttop - x15) * Ftop * tautop  - x15loss * x15) * dt* c15 + sigma_x15 * dw15)

# Set observation equation variances
obj$add_observation_variances(X0 ~ sigma_X0^2)
obj$add_observation_variances(X1 ~ sigma_X1^2)
obj$add_observation_variances(X2 ~ sigma_X2^2)
obj$add_observation_variances(X3 ~ sigma_X3^2)
obj$add_observation_variances(X4 ~ sigma_X4^2)
obj$add_observation_variances(X5 ~ sigma_X5^2)
obj$add_observation_variances(X6 ~ sigma_X6^2)
obj$add_observation_variances(X7 ~ sigma_X7^2)
obj$add_observation_variances(X8 ~ sigma_X8^2)
obj$add_observation_variances(X9 ~ sigma_X9^2)
obj$add_observation_variances(X10 ~ sigma_X10^2)
obj$add_observation_variances(X11 ~ sigma_X11^2)
obj$add_observation_variances(X12 ~ sigma_X12^2)
obj$add_observation_variances(X13 ~ sigma_X13^2)
obj$add_observation_variances(X14 ~ sigma_X14^2)
obj$add_observation_variances(X15 ~ sigma_X15^2)



# Add temperature parameters
lower <- 0
upper <- 1
init <- 1e-3
obj$add_parameters(taux1 = c(init = get_TMB_parameter_value(old_fit, "taux1", init), lb = lower, ub = upper))
obj$add_parameters(taux2 = c(init = get_TMB_parameter_value(old_fit, "taux2", init), lb = lower, ub = upper))
obj$add_parameters(taux3 = c(init = get_TMB_parameter_value(old_fit, "taux3", init), lb = lower, ub = upper))
obj$add_parameters(taux4 = c(init = get_TMB_parameter_value(old_fit, "taux4", init), lb = lower, ub = upper))
obj$add_parameters(taux5 = c(init = get_TMB_parameter_value(old_fit, "taux5", init), lb = lower, ub = upper))
obj$add_parameters(taux6 = c(init = get_TMB_parameter_value(old_fit, "taux6", init), lb = lower, ub = upper))
obj$add_parameters(taux7 = c(init = get_TMB_parameter_value(old_fit, "taux7", init), lb = lower, ub = upper))
obj$add_parameters(taux8 = c(init = get_TMB_parameter_value(old_fit, "taux8", init), lb = lower, ub = upper))
obj$add_parameters(taux9 = c(init = get_TMB_parameter_value(old_fit, "taux9", init), lb = lower, ub = upper))
obj$add_parameters(taux10 = c(init = get_TMB_parameter_value(old_fit, "taux10", init), lb =lower, ub = upper))
obj$add_parameters(taux11 = c(init = get_TMB_parameter_value(old_fit, "taux11", init), lb = lower, ub = upper))
obj$add_parameters(taux12 = c(init = get_TMB_parameter_value(old_fit, "taux12", init), lb = lower, ub = upper))
obj$add_parameters(taux13 = c(init = get_TMB_parameter_value(old_fit, "taux13", init), lb = lower, ub = upper))
obj$add_parameters(taux14 = c(init = get_TMB_parameter_value(old_fit, "taux14", init), lb = lower, ub = upper))
obj$add_parameters(taux15 = c(init = get_TMB_parameter_value(old_fit, "taux15", init), lb = lower, ub = upper))
obj$add_parameters(tautop = c(init = get_TMB_parameter_value(old_fit, "tautop", init), lb = lower, ub = upper))
obj$add_parameters(taubot = c(init = get_TMB_parameter_value(old_fit, "taubot", init), lb = lower, ub = upper))
obj$add_parameters(taumid = c(init = get_TMB_parameter_value(old_fit, "taumid", init), lb = lower, ub = upper))

obj$add_parameters(c0 = c(init = get_TMB_parameter_value(old_fit, "c0", init), lb = lower, ub = upper))
obj$add_parameters(c1 = c(init = get_TMB_parameter_value(old_fit, "c1", init), lb = lower, ub = upper))
obj$add_parameters(c2 = c(init = get_TMB_parameter_value(old_fit, "c2", init), lb = lower, ub = upper))
obj$add_parameters(c3 = c(init = get_TMB_parameter_value(old_fit, "c3", init), lb = lower, ub = upper))
obj$add_parameters(c4 = c(init = get_TMB_parameter_value(old_fit, "c4", init), lb = lower, ub = upper))
obj$add_parameters(c5 = c(init = get_TMB_parameter_value(old_fit, "c5", init), lb = lower, ub = upper))
obj$add_parameters(c6 = c(init = get_TMB_parameter_value(old_fit, "c6", init), lb = lower, ub = upper))
obj$add_parameters(c7 = c(init = get_TMB_parameter_value(old_fit, "c7", init), lb = lower, ub = upper))
obj$add_parameters(c8 = c(init = get_TMB_parameter_value(old_fit, "c8", init), lb = lower, ub = upper))
obj$add_parameters(c9 = c(init = get_TMB_parameter_value(old_fit, "c9", init), lb = lower, ub = upper))
obj$add_parameters(c10 = c(init = get_TMB_parameter_value(old_fit, "c10", init), lb = lower, ub = upper))
obj$add_parameters(c11 = c(init = get_TMB_parameter_value(old_fit, "c11", init), lb = lower, ub = upper))
obj$add_parameters(c12 = c(init = get_TMB_parameter_value(old_fit, "c12", init), lb = lower, ub = upper))
obj$add_parameters(c13 = c(init = get_TMB_parameter_value(old_fit, "c13", init), lb = lower, ub = upper))
obj$add_parameters(c14 = c(init = get_TMB_parameter_value(old_fit, "c14", init), lb = lower, ub = upper))
obj$add_parameters(c15 = c(init = get_TMB_parameter_value(old_fit, "c15", init), lb = lower, ub = upper))

obj$add_parameters(xloss = c(init = get_TMB_parameter_value(old_fit, "xloss", init), lb = lower, ub = upper))
obj$add_parameters(x1loss = c(init = get_TMB_parameter_value(old_fit, "x1loss", init), lb = lower, ub = upper))
obj$add_parameters(x2loss = c(init = get_TMB_parameter_value(old_fit, "x2loss", init), lb = lower, ub = upper))
obj$add_parameters(x3loss = c(init = get_TMB_parameter_value(old_fit, "x3loss", init), lb = lower, ub = upper))
obj$add_parameters(x4loss = c(init = get_TMB_parameter_value(old_fit, "x4loss", init), lb = lower, ub = upper))
obj$add_parameters(x5loss = c(init = get_TMB_parameter_value(old_fit, "x5loss", init), lb = lower, ub = upper))
obj$add_parameters(x6loss = c(init = get_TMB_parameter_value(old_fit, "x6loss", init), lb = lower, ub = upper))
obj$add_parameters(x7loss = c(init = get_TMB_parameter_value(old_fit, "x7loss", init), lb = lower, ub = upper))
obj$add_parameters(x8loss = c(init = get_TMB_parameter_value(old_fit, "x8loss", init), lb = lower, ub = upper))
obj$add_parameters(x9loss = c(init = get_TMB_parameter_value(old_fit, "x9loss", init), lb = lower, ub = upper))
obj$add_parameters(x10loss = c(init = get_TMB_parameter_value(old_fit, "x10loss", init), lb = lower, ub = upper))
obj$add_parameters(x11loss = c(init = get_TMB_parameter_value(old_fit, "x11loss", init), lb = lower, ub = upper))
obj$add_parameters(x12loss = c(init = get_TMB_parameter_value(old_fit, "x12loss", init), lb = lower, ub = upper))
obj$add_parameters(x13loss = c(init = get_TMB_parameter_value(old_fit, "x13loss", init), lb = lower, ub = upper))
obj$add_parameters(x14loss = c(init = get_TMB_parameter_value(old_fit, "x14loss", init), lb = lower, ub = upper))
obj$add_parameters(x15loss = c(init = get_TMB_parameter_value(old_fit, "x15loss", init), lb = lower, ub = upper))



obj$add_inputs(Ttop, Ftop, Tbot, Fbot, Tmid, Fmid)
# Add algebracis
obj$add_algebraics(sigma_x0 ~ exp(logsigma_x0))
obj$add_algebraics(sigma_x1 ~ exp(logsigma_x1))
obj$add_algebraics(sigma_x2 ~ exp(logsigma_x2))
obj$add_algebraics(sigma_x3 ~ exp(logsigma_x3))
obj$add_algebraics(sigma_x4 ~ exp(logsigma_x4))
obj$add_algebraics(sigma_x5 ~ exp(logsigma_x5))
obj$add_algebraics(sigma_x6 ~ exp(logsigma_x6))
obj$add_algebraics(sigma_x7 ~ exp(logsigma_x7))
obj$add_algebraics(sigma_x8 ~ exp(logsigma_x8))
obj$add_algebraics(sigma_x9 ~ exp(logsigma_x9))
obj$add_algebraics(sigma_x10 ~ exp(logsigma_x10))
obj$add_algebraics(sigma_x11 ~ exp(logsigma_x11))
obj$add_algebraics(sigma_x12 ~ exp(logsigma_x12))
obj$add_algebraics(sigma_x13 ~ exp(logsigma_x13))
obj$add_algebraics(sigma_x14 ~ exp(logsigma_x14))
obj$add_algebraics(sigma_x15 ~ exp(logsigma_x15))
# # Set system noise
lower <- 1e-5
upper <- 1
init <- 1e-3
obj$add_parameters(logsigma_x0 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x0", init), lb = lower, ub = upper)))
obj$add_parameters(logsigma_x1 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x1", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x2 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x2", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x3 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x3", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x4 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x4", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x5 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x5", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x6 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x6", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x7 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x7", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x8 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x8", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x9 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x9", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x10 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x10", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x11 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x11", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x12 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x12", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x13 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x13", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x14 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x14", init), lb = lower, ub =  upper)))
obj$add_parameters(logsigma_x15 = log(c(init = get_TMB_parameter_value(old_fit, "logsigma_x15", init), lb = lower, ub =  upper)))



# Set observation noise
lower <-1e-9
upper <- 1e-3
init <- 1e-4
# obj$add_parameters(sigma_X0 = c(init = get_TMB_parameter_value(old_fit, "sigma_X0", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X1 = c(init = get_TMB_parameter_value(old_fit, "sigma_X1", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X2 = c(init = get_TMB_parameter_value(old_fit, "sigma_X2", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X3 = c(init = get_TMB_parameter_value(old_fit, "sigma_X3", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X4 = c(init = get_TMB_parameter_value(old_fit, "sigma_X4", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X5 = c(init = get_TMB_parameter_value(old_fit, "sigma_X5", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X6 = c(init = get_TMB_parameter_value(old_fit, "sigma_X6", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X7 = c(init = get_TMB_parameter_value(old_fit, "sigma_X7", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X8 = c(init = get_TMB_parameter_value(old_fit, "sigma_X8", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X9 = c(init = get_TMB_parameter_value(old_fit, "sigma_X9", init), lb = lower, ub =  upper))
# obj$add_parameters(sigma_X10 = c(init = get_TMB_parameter_value(old_fit, "sigma_X10", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X11 = c(init = get_TMB_parameter_value(old_fit, "sigma_X11", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X12 = c(init = get_TMB_parameter_value(old_fit, "sigma_X12", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X13 = c(init = get_TMB_parameter_value(old_fit, "sigma_X13", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X14 = c(init = get_TMB_parameter_value(old_fit, "sigma_X14", init), lb = lower, ub = upper))
# obj$add_parameters(sigma_X15 = c(init = get_TMB_parameter_value(old_fit, "sigma_X15", init), lb = lower, ub = upper))

### SETTING IT ALL TO ZERO FOR NOW ###
obj$add_parameters(sigma_X0=lower)
obj$add_parameters(sigma_X1=lower)
obj$add_parameters(sigma_X2=lower)
obj$add_parameters(sigma_X3=lower)
obj$add_parameters(sigma_X4=lower)
obj$add_parameters(sigma_X5=lower)
obj$add_parameters(sigma_X6=lower)
obj$add_parameters(sigma_X7=lower)
obj$add_parameters(sigma_X8=lower)
obj$add_parameters(sigma_X9=lower)
obj$add_parameters(sigma_X10=lower)
obj$add_parameters(sigma_X11=lower)
obj$add_parameters(sigma_X12=lower)
obj$add_parameters(sigma_X13=lower)
obj$add_parameters(sigma_X14=lower)
obj$add_parameters(sigma_X15=lower)


# Set initial states
x0 <- c(data$X0[1], data$X1[1], data$X2[1], data$X3[1], data$X4[1], data$X5[1], data$X6[1], data$X7[1], data$X8[1], data$X9[1], data$X10[1], data$X11[1], data$X12[1], data$X13[1], data$X14[1], data$X15[1])
N <- length(x0)
init_noise <- 1e-7
obj$set_initial_state(x0, init_noise*diag(N))

# Estimate
fit <- obj$estimate(data, 
                method = 'ekf',
                ode.solver = 'rk4',
                compile = FALSE,
                use.hessian = FALSE,
                ode.timestep = 0.5,
                control = list(maxit = 50, trace = 1))

# save_TMB_fit(fit, name)
load('models/dronninglundTMB/fit_ll_19671.RData') # (Best??)

summary(fit)
fit$nll
fit$nll
fit$par.fixed
fit$states$sd$prior

# # Plot
plot(fit$states$mean$prior$x0, type = 'l', col = 'blue', ylim = c(0.1, 0.2))
points(fit$observations$mean$prior$X0, col = 'red')


plot(fit$states$mean$prior$x1, type = 'l', col = 'blue', ylim = c(0.1, 0.2))
points(fit$observations$mean$prior$X1, col = 'red')

plot(fit$states$mean$prior$x2, col = 'green', type = 'l')
points(fit$observations$mean$prior$X2, col = 'red')

plot(fit$states$mean$prior$x3, col = 'black', type = 'l')
points(fit$observations$mean$prior$X3, col = 'red')

plot(fit$states$mean$prior$x4, col = 'blue', type = 'l')
points(fit$observations$mean$prior$X4, col = 'red')

plot(fit$states$mean$prior$x5, col = 'green', type = 'l')
points(fit$observations$mean$prior$X5, col = 'red')

plot(fit$states$mean$prior$x6, col = 'black', type = 'l')
points(fit$observations$mean$prior$X6, col = 'red')


plot(fit$states$mean$prior$x7, col = 'blue', type = 'l')
points(fit$observations$mean$prior$X7, col = 'red')

plot(fit$states$mean$prior$x8, col = 'green', type = 'l')
points(fit$observations$mean$prior$X8, col = 'red')

plot(fit$states$mean$prior$x9, col = 'black', type = 'l')
points(fit$observations$mean$prior$X9, col = 'red')

plot(fit$states$mean$prior$x10, col = 'blue', type = 'l')
points(fit$observations$mean$prior$X10, col = 'red')

plot(fit$states$mean$prior$x11, col = 'green', type = 'l')
points(fit$observations$mean$prior$X11, col = 'red')

plot(fit$states$mean$prior$x12, col = 'black', type = 'l')
points(fit$observations$mean$prior$X12, col = 'red')

plot(fit$states$mean$prior$x13, col = 'blue', type = 'l')
points(fit$observations$mean$prior$X13, col = 'red')

plot(fit$states$mean$prior$x14, col = 'green', type = 'l')
points(fit$observations$mean$prior$X14, col = 'red')

plot(fit$states$mean$prior$x15, col = 'black', type = 'l')
points(fit$observations$mean$prior$X15, col = 'red')


plot(fit$residuals$mean['X0'][,1], col = 'blue')
plot(fit$residuals$mean['X1'][,1], col = 'blue')
plot(fit$residuals$mean['X2'][,1], col = 'blue')
plot(fit$residuals$mean['X3'][,1], col = 'blue')
plot(fit$residuals$mean['X4'][,1], col = 'blue')
plot(fit$residuals$mean['X5'][,1], col = 'blue')
plot(fit$residuals$mean['X6'][,1], col = 'blue')
plot(fit$residuals$mean['X7'][,1], col = 'blue')
plot(fit$residuals$mean['X8'][,1], col = 'blue')
plot(fit$residuals$mean['X9'][,1], col = 'blue')
plot(fit$residuals$mean['X10'][,1], col = 'blue')
plot(fit$residuals$mean['X11'][,1], col = 'blue')
plot(fit$residuals$mean['X12'][,1], col = 'blue')
plot(fit$residuals$mean['X13'][,1], col = 'blue')
plot(fit$residuals$mean['X14'][,1], col = 'blue')
plot(fit$residuals$mean['X15'][,1], col = 'blue')



X <- cbind(water_sensors, inputs)[5500:6000,]

pred <- obj$predict(data, k.ahead = 10)


plot(pred[pred$k.ahead==1,]$x0, type = 'l', ylim = c(0,1))
points(X$X0[10:length(X$X0)])

pred$k.ahead


length(pred$x0)

length(pred[pred$k.ahead==1,]$x0)
length(X$X0)
