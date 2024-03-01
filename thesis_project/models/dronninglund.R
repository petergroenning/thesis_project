library(ctsmr)
library(ggplot2)
library(patchwork)
rm(list = ls())
source('thesis_project/models/ctsm_utils.R')


###########################################################
# Make data
############################################################
water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
inputs  <- read.csv('data/processed/dronninglund/inputs.csv')

data <- cbind(water_sensors, inputs)[3000:3200,]
data <- data[, !duplicated(colnames(data))]


#### MAKE CTSM OBJECT ####
name <- 'dronninglund'
old_fit <- load_best_fit("dronninglund")
obj = ctsm$new()


# Add observation equations and variances
obj$addObs(X0 ~ x0)
obj$addObs(X1 ~ x1)
obj$addObs(X2 ~ x2)
obj$addObs(X3 ~ x3)
# obj$addObs(X4 ~ x4)
# obj$addObs(X5 ~ x5)
# obj$addObs(X6 ~ x6)
# obj$addObs(X7 ~ x7)


# Add system equations
obj$addSystem(dx0 ~ ((x1 - x0) * 1/(tau_x0) + (Tbot - x0) * Fbot * rho_x0 + (Tmid - x0) * Fmid * psi_x0 + (Ttop - x0) * Ftop * phi_x0) * dt + sigma_x0 * dw0)
obj$addSystem(dx1 ~ (-(x1 - x0)/ tau_x0 + (x2 - x1) / tau_x1 + (Tbot - x1) * Fbot * rho_x1 + (Tmid - x1) * Fmid * psi_x1 + (Ttop - x1) * Ftop * phi_x1) * dt + sigma_x1 * dw1)
obj$addSystem(dx2 ~ (-(x2 - x1)/ tau_x1 + (x3 - x2) / tau_x2 + (Tbot - x2) * Fbot * rho_x2 + (Tmid - x2) * Fmid * psi_x2 + (Ttop - x2) * Ftop * phi_x2) * dt + sigma_x2 * dw2)
obj$addSystem(dx3 ~ (-(x3 - x2)/ tau_x2  + (Tbot - x3) * Fbot * rho_x3 + (Tmid - x3) * Fmid * psi_x3 + (Ttop - x3) * Ftop * phi_x3) * dt + sigma_x3 * dw3)
# obj$addSystem(dx4 ~ (-(x4 - x3)/ tau_x3 + (x5 - x4) / tau_x4 + (Tbot - x4) * Fbot * rho_x4 + (Tmid - x4) * Fmid * psi_x4 + (Ttop - x4) * Ftop * phi_x4) * dt + sigma_x4 * dw4)
# obj$addSystem(dx5 ~ (-(x5 - x4)/ tau_x4 + (x6 - x5) / tau_x5 + (Tbot - x5) * Fbot * rho_x5 + (Tmid - x5) * Fmid * psi_x5 + (Ttop - x5) * Ftop * phi_x5) * dt + sigma_x5 * dw5)
# obj$addSystem(dx6 ~ (-(x6 - x5)/ tau_x5 + (x7 - x6) / tau_x6 + (Tbot - x6) * Fbot * rho_x6 + (Tmid - x6) * Fmid * psi_x6 + (Ttop - x6) * Ftop * phi_x6) * dt + sigma_x6 * dw6)
# obj$addSystem(dx7 ~ (-(x7 - x6)/ tau_x6 + (Tbot - x7) * Fbot * rho_x7 + (Tmid - x7) * Fmid * psi_x7 + (Ttop - x7) * Ftop * phi_x7) * dt + sigma_x7 * dw7)

# Set observation equation variances
obj$setVariance(X0 ~ sigma_X0)
obj$setVariance(X1 ~ sigma_X1)
obj$setVariance(X2 ~ sigma_X2)
obj$setVariance(X3 ~ sigma_X3)

# Add temperature parameters
obj$setParameter(tau_x0 = c(init = get_CTSM_parameter_value(old_fit, "tau_x0", 1), lb = 1e-5, ub = 1000))
obj$setParameter(tau_x1 = c(init = get_CTSM_parameter_value(old_fit, "tau_x1", 1), lb = 1e-5, ub = 1000))
obj$setParameter(tau_x2 = c(init = get_CTSM_parameter_value(old_fit, "tau_x2", 1), lb = 1e-5, ub = 1000))
# obj$setParameter(tau_x3 = c(init = get_CTSM_parameter_value(old_fit, "tau_x3", 1), lb = 1e-5, ub = 1000))
# obj$setParameter(tau_x4 = c(init = get_CTSM_parameter_value(old_fit, "tau_x4", 1), lb = 1e-5, ub = 1000))
# obj$setParameter(tau_x5 = c(init = get_CTSM_parameter_value(old_fit, "tau_x5", 1), lb = 1e-5, ub = 1000))
# obj$setParameter(tau_x6 = c(init = get_CTSM_parameter_value(old_fit, "tau_x6", 1), lb = 1e-5, ub = 1000))

# Set system bottom diffuser parameters
obj$setParameter(rho_x0 = c(init = get_CTSM_parameter_value(old_fit, "rho_x0", 1), lb = 0, ub = 1000))
obj$setParameter(rho_x1 = c(init = get_CTSM_parameter_value(old_fit, "rho_x1", 1), lb = 0, ub = 1000))
obj$setParameter(rho_x2 = c(init = get_CTSM_parameter_value(old_fit, "rho_x2", 1), lb = 0, ub = 1000))
obj$setParameter(rho_x3 = c(init = get_CTSM_parameter_value(old_fit, "rho_x3", 1), lb = 0, ub = 1000))
# obj$setParameter(rho_x4 = c(init = get_CTSM_parameter_value(old_fit, "rho_x4", 1), lb = 0, ub = 1000))
# obj$setParameter(rho_x5 = c(init = get_CTSM_parameter_value(old_fit, "rho_x5", 1), lb = 0, ub = 1000))
# obj$setParameter(rho_x6 = c(init = get_CTSM_parameter_value(old_fit, "rho_x6", 1), lb = 0, ub = 1000))
# obj$setParameter(rho_x7 = c(init = get_CTSM_parameter_value(old_fit, "rho_x7", 1), lb = 0, ub = 1000))

# Set system middle diffuser parameters
obj$setParameter(psi_x0 = c(init = get_CTSM_parameter_value(old_fit, "psi_x0", 1), lb = 0, ub = 1000))
obj$setParameter(psi_x1 = c(init = get_CTSM_parameter_value(old_fit, "psi_x1", 1), lb = 0, ub = 1000))
obj$setParameter(psi_x2 = c(init = get_CTSM_parameter_value(old_fit, "psi_x2", 1), lb = 0, ub = 1000))
obj$setParameter(psi_x3 = c(init = get_CTSM_parameter_value(old_fit, "psi_x3", 1), lb = 0, ub = 1000))
# obj$setParameter(psi_x4 = c(init = get_CTSM_parameter_value(old_fit, "psi_x4", 1), lb = 0, ub = 1000))
# obj$setParameter(psi_x5 = c(init = get_CTSM_parameter_value(old_fit, "psi_x5", 1), lb = 0, ub = 1000))
# obj$setParameter(psi_x6 = c(init = get_CTSM_parameter_value(old_fit, "psi_x6", 1), lb = 0, ub = 1000))
# obj$setParameter(psi_x7 = c(init = get_CTSM_parameter_value(old_fit, "psi_x7", 1), lb = 0, ub = 1000))

# Set system top diffuser parameters
obj$setParameter(phi_x0 = c(init = get_CTSM_parameter_value(old_fit, "phi_x0", 1), lb = 0, ub = 1000))
obj$setParameter(phi_x1 = c(init = get_CTSM_parameter_value(old_fit, "phi_x1", 1), lb = 0, ub = 1000))
obj$setParameter(phi_x2 = c(init = get_CTSM_parameter_value(old_fit, "phi_x2", 1), lb = 0, ub = 1000))
obj$setParameter(phi_x3 = c(init = get_CTSM_parameter_value(old_fit, "phi_x3", 1), lb = 0, ub = 1000))
# obj$setParameter(phi_x4 = c(init = get_CTSM_parameter_value(old_fit, "phi_x4", 1), lb = 0, ub = 1000))
# obj$setParameter(phi_x5 = c(init = get_CTSM_parameter_value(old_fit, "phi_x5", 1), lb = 0, ub = 1000))
# obj$setParameter(phi_x6 = c(init = get_CTSM_parameter_value(old_fit, "phi_x6", 1), lb = 0, ub = 1000))
# obj$setParameter(phi_x7 = c(init = get_CTSM_parameter_value(old_fit, "phi_x7", 1), lb = 0, ub = 1000))


# # Set system noise
obj$setParameter(sigma_x0 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x0", 0.1), lb = 1e-5, ub = 5))
obj$setParameter(sigma_x1 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x1", 0.1), lb = 1e-5, ub = 5))
obj$setParameter(sigma_x2 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x2", 0.1), lb = 1e-5, ub = 5))
obj$setParameter(sigma_x3 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x3", 0.1), lb = 1e-5, ub = 5))
# obj$setParameter(sigma_x4 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x4", 0.1), lb = 1e-5, ub = 5))
# obj$setParameter(sigma_x5 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x5", 0.1), lb = 1e-5, ub = 5))
# obj$setParameter(sigma_x6 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x6", 0.1), lb = 1e-5, ub = 5))
# obj$setParameter(sigma_x7 = c(init = get_CTSM_parameter_value(old_fit, "sigma_x7", 0.1), lb = 1e-5, ub = 5))


# obj$setParameter(sigma_x0 = 0.5)
# obj$setParameter(sigma_x1 = 0.5)
# obj$setParameter(sigma_x2 = 0.5)
# obj$setParameter(sigma_x3 = 0.5)
# obj$setParameter(sigma_x4 = 0.5)
# obj$setParameter(sigma_x5 = 0.5)
# obj$setParameter(sigma_x6 = 0.5)
# obj$setParameter(sigma_x7 = 0.5)


# Set observation noise
obj$setParameter(sigma_X0 = 0.1)
obj$setParameter(sigma_X1 = 0.1)
obj$setParameter(sigma_X2 = 0.1)
obj$setParameter(sigma_X3 = 0.1)
# obj$setParameter(sigma_X4 = 0.1)
# obj$setParameter(sigma_X5 = 0.1)
# obj$setParameter(sigma_X6 = 0.1)
# obj$setParameter(sigma_X7 = 0.1)



# Set initial states
obj$setParameter(x0 = data$X0[1])
obj$setParameter(x1 = data$X1[1])
obj$setParameter(x2 = data$X2[1])
obj$setParameter(x3 = data$X3[1])
# obj$setParameter(x4 = data$X4[1])
# obj$setParameter(x5 = data$X5[1])
# obj$setParameter(x6 = data$X6[1])
# obj$setParameter(x7 = data$X7[1])


# Add inputs
obj$addInput("Tbot", "Tmid", "Ttop", "Fbot", "Fmid", "Ftop")


# Estimate
fit <- obj$estimate(data)
save_ctsm_fit(fit, name)

X <- cbind(water_sensors, inputs)[3000:3005,]
X <- X[, !duplicated(colnames(X))]
tmp <- predict(fit, newdata=X)


plot(tmp$state$pred$x0)
lines(X$X1)

summary(fit)





