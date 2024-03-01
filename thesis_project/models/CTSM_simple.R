library(ctsmr)
source('thesis_project/models/ctsm_utils.R')
library(ggplot2)


make_model <- function(data, old_fit){

    obj <- ctsm$new()
    if (missing(old_fit)){
        old_fit <- NULL
    }

    # Add observation equations and variances
    obj$addObs(X0 ~ xa)
    obj$addObs(X1 ~ x1)
    obj$addObs(X2 ~ x2)
    obj$addObs(X3 ~ x3)
    obj$addObs(X4 ~ x4)
    obj$addObs(X5 ~ x5)
    obj$addObs(X6 ~ x6)
    obj$addObs(X7 ~ x7)
    obj$addObs(X8 ~ x8)
    obj$addObs(X9 ~ x9)
    obj$addObs(X10 ~ x1b)
    obj$addObs(X11 ~ x11)
    obj$addObs(X12 ~ x12)
    obj$addObs(X13 ~ x13)
    obj$addObs(X14 ~ x14)
    obj$addObs(X15 ~ x15)

    # Set observation equation variances
    obj$setVariance(X0 ~ sigma_Xa^2)
    obj$setVariance(X1 ~ sigma_X1^2)
    obj$setVariance(X2 ~ sigma_X2^2)
    obj$setVariance(X3 ~ sigma_X3^2)
    obj$setVariance(X4 ~ sigma_X4^2)
    obj$setVariance(X5 ~ sigma_X5^2)
    obj$setVariance(X6 ~ sigma_X6^2)
    obj$setVariance(X7 ~ sigma_X7^2)
    obj$setVariance(X8 ~ sigma_X8^2)
    obj$setVariance(X9 ~ sigma_X9^2)
    obj$setVariance(X10 ~ sigma_X1b^2)
    obj$setVariance(X11 ~ sigma_X11^2)
    obj$setVariance(X12 ~ sigma_X12^2)
    obj$setVariance(X13 ~ sigma_X13^2)
    obj$setVariance(X14 ~ sigma_X14^2)
    obj$setVariance(X15 ~ sigma_X15^2)

    # Add system equations
    obj$addSystem(dxa ~ 1/C0 *  ((x1 - xa) * (R + Cbot*(FbotIn + FbotOut))) * dt + (Tbot - xa) * Vbot * FbotIn * dt +  sigma_xa * dw0)
    obj$addSystem(dx1 ~ 1/C0 *  ((xa - x1) * (R + Cbot*(FbotIn + FbotOut)) + (x2 - x1) * R ) * dt + sigma_x1 * dw1)
    obj$addSystem(dx2 ~ 1/C0 *  ((x1 - x2) * R + (x3 - x2) * R) * dt + sigma_x2 * dw2)
    obj$addSystem(dx3 ~ 1/C0 *  ((x2 - x3) * R + (x4 - x3) * R) * dt + sigma_x3 * dw3)
    obj$addSystem(dx4 ~ 1/C0 *  ((x3 - x4) * R + (x5 - x4) * R) * dt + sigma_x4 * dw4)
    obj$addSystem(dx5 ~ 1/C2 *  ((x4 - x5) * R + (x6 - x5) * R) * dt + sigma_x5 * dw5)
    obj$addSystem(dx6 ~ 1/C2 *  ((x5 - x6) * R + (x7 - x6) * R) * dt + sigma_x6 * dw6)
    obj$addSystem(dx7 ~ 1/C2 *  ((x6 - x7) * R + (x8 - x7) * R) * dt + sigma_x7 * dw7)
    obj$addSystem(dx8 ~ 1/C2 *  ((x7 - x8) * R + (x9 - x8) * R) * dt + sigma_x8 * dw8)
    obj$addSystem(dx9 ~ 1/C2 *  ((x8 - x9) * R + (x1b - x9) * R) * dt + sigma_x9 * dw9)
    obj$addSystem(dx1b ~ 1/C4 * ((x9 - x1b) * R + (x11 - x1b) *  R) * dt + sigma_x1b * dw10)
    obj$addSystem(dx11 ~ 1/C4 * ((x1b - x11) * R + (x12 - x11) * R) * dt + sigma_x11 * dw11)

    obj$addSystem(dx12 ~ 1/C4 * ((x11 - x12) * R + (x13 - x12) * R) * dt + sigma_x12 * dw12)
    obj$addSystem(dx13 ~ 1/C4 * ((x12 - x13) * R + (x14 - x13) * R) * dt + sigma_x13 * dw13)
    obj$addSystem(dx14 ~ 1/C4 * ((x13 - x14) * R + (x15 - x14) * (R + Ctop*(FtopIn + FtopOut))) * dt + sigma_x14 * dw14)
    obj$addSystem(dx15 ~ 1/C5* ((x14 - x15) * (R + Ctop*(FtopIn + FtopOut))) * dt + (Ttop - x15) * Vtop * FtopIn * dt + sigma_x15 * dw15)

    # (R + Cmid*(FmidIn + FmidOut))
    # (Tmid - x1b) * Vmid * FmidIn * dt 
    # Add input 
    obj$addInput("Ttop","Tmid", "Tbot", "FtopIn", "FtopOut", "FbotIn", "FbotOut", "FmidIn", "FmidOut")
    
    ###### SET PARAMETERS ######
    # Observation Noise
    obs_std <- 0.15
    obj$setParameter(sigma_Xa = obs_std)
    obj$setParameter(sigma_X1 = obs_std)
    obj$setParameter(sigma_X2 = obs_std)
    obj$setParameter(sigma_X3 = obs_std)
    obj$setParameter(sigma_X4 = obs_std)
    obj$setParameter(sigma_X5 = obs_std)
    obj$setParameter(sigma_X6 = obs_std)
    obj$setParameter(sigma_X7 = obs_std)
    obj$setParameter(sigma_X8 = obs_std)
    obj$setParameter(sigma_X9 = obs_std)
    obj$setParameter(sigma_X1b = obs_std)
    obj$setParameter(sigma_X11 = obs_std)
    obj$setParameter(sigma_X12 = obs_std)
    obj$setParameter(sigma_X13 = obs_std)
    obj$setParameter(sigma_X14 = obs_std)
    obj$setParameter(sigma_X15 = obs_std)

    # System Noise
    lb <- 1e-5
    ub <- 10
    init <- 1
    obj$setParameter(sigma_xa = c(init = getParam('sigma_x0', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x1 = c(init = getParam('sigma_x1', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x2 = c(init = getParam('sigma_x2', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x3 = c(init = getParam('sigma_x3', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x4 = c(init = getParam('sigma_x4', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x5 = c(init = getParam('sigma_x5', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x6 = c(init = getParam('sigma_x6', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x7 = c(init = getParam('sigma_x7', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x8 = c(init = getParam('sigma_x8', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x9 = c(init = getParam('sigma_x9', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x1b = c(init = getParam('sigma_x10', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x11 = c(init = getParam('sigma_x11', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x12 = c(init = getParam('sigma_x12', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x13 = c(init = getParam('sigma_x13', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x14 = c(init = getParam('sigma_x14', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(sigma_x15 = c(init = getParam('sigma_x15', init, old_fit), lb = lb, ub = ub))

    # obj$setParameter(sigma_xa = getParam('sigma_x0', init, old_fit))
    # obj$setParameter(sigma_x1 = getParam('sigma_x1', init, old_fit))
    # obj$setParameter(sigma_x2 = getParam('sigma_x2', init, old_fit))
    # obj$setParameter(sigma_x3 = getParam('sigma_x3', init, old_fit))
    # obj$setParameter(sigma_x4 = getParam('sigma_x4', init, old_fit))
    # obj$setParameter(sigma_x5 = getParam('sigma_x5', init, old_fit))
    # obj$setParameter(sigma_x6 = getParam('sigma_x6', init, old_fit))
    # obj$setParameter(sigma_x7 = getParam('sigma_x7', init, old_fit))
    # obj$setParameter(sigma_x8 = getParam('sigma_x8', init, old_fit))
    # obj$setParameter(sigma_x9 = getParam('sigma_x9', init, old_fit))
    # obj$setParameter(sigma_x1b = getParam('sigma_x10', init, old_fit))
    # obj$setParameter(sigma_x11 = getParam('sigma_x11', init, old_fit))
    # obj$setParameter(sigma_x12 = getParam('sigma_x12', init, old_fit))
    # obj$setParameter(sigma_x13 = getParam('sigma_x13', init, old_fit))
    # obj$setParameter(sigma_x14 = getParam('sigma_x14', init, old_fit))
    # obj$setParameter(sigma_x15 = getParam('sigma_x15', init, old_fit))



    # System Parameters
    lb <- 1e-5
    ub <- 1000
    init <- 1
    obj$setParameter(R = c(init = getParam('R', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Cbot = c(init = getParam('Cbot', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Cmid = c(init = getParam('Cmid', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Ctop = c(init = getParam('Ctop', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Vbot = c(init = getParam('Vbot', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Vmid = c(init = getParam('Vmid', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(Vtop = c(init = getParam('Vtop', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(C0 = c(init = getParam('C', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C1 = c(init = getParam('C1', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(C2 = c(init = getParam('C2', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C3 = c(init = getParam('C3', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(C4 = c(init = getParam('C1', init, old_fit), lb = lb, ub = ub))
    obj$setParameter(C5 = c(init = getParam('C1', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C6 = c(init = getParam('C1', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C7 = c(init = getParam('C1', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C8 = c(init = getParam('C2', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C9 = c(init = getParam('C2', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C10 = c(init = getParam('C2', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C11 = c(init = getParam('C2', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C12 = c(init = getParam('C3', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C13 = c(init = getParam('C3', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C14 = c(init = getParam('C3', init, old_fit), lb = lb, ub = ub))
    # obj$setParameter(C15 = c(init = getParam('C4', init, old_fit), lb = lb, ub = ub))

    # Initial states
    obj$setParameter(xa = data$X0[1])
    obj$setParameter(x1 = data$X1[1])
    obj$setParameter(x2 = data$X2[1])
    obj$setParameter(x3 = data$X3[1])
    obj$setParameter(x4 = data$X4[1])
    obj$setParameter(x5 = data$X5[1])
    obj$setParameter(x6 = data$X6[1])
    obj$setParameter(x7 = data$X7[1])
    obj$setParameter(x8 = data$X8[1])
    obj$setParameter(x9 = data$X9[1])
    obj$setParameter(x1b = data$X10[1])
    obj$setParameter(x11 = data$X11[1])
    obj$setParameter(x12 = data$X12[1])
    obj$setParameter(x13 = data$X13[1])
    obj$setParameter(x14 = data$X14[1])
    obj$setParameter(x15 = data$X15[1])

    return(obj)
}


###########################################################
# Make data
############################################################
water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
inputs[,8:13] <- inputs[,8:13] / 1000 # Convert to m3/hr

data <- cbind(water_sensors, inputs)[8500:9000,]
data <- data[, !duplicated(colnames(data))]
# dropna
data <- data[complete.cases(data),]


# Load old fit
model <- 'oldCTSMfit'
name <- 'fit_ll_6593'

load(paste0("models/", model, "/", name, ".RData"))

# Make model object
obj <- make_model(data = data, old_fit = fit)
obj$setOptions(list("maxNumberOfEval" = 250, "eps" = 1e-8))
fit <- obj$estimate(data, firstorder = TRUE, threads = 16)

obj$ParameterValues

# save(fit, file = paste0("models/", model, "/", name, ".RData"))
save_ctsm_fit(fit, model)
summary(fit)
fit$loglik

fit

fit$loglik
load("models/oldCTSMfit/fit_ll_17821_c.RData")
load("models/oldCTSMfit/fit_ll_6593.RData")
new_data <- cbind(water_sensors, inputs)[7000:12000,]
new_data <- new_data[, !duplicated(colnames(new_data))]
# # # # dropna
new_data <- new_data[complete.cases(new_data),]
# new_data <- data
preds1 <- predict(fit, newdata = new_data, n.ahead = 24)
se <- preds1$state$sd
# # # rename columns
colnames(se) <- paste0('X', 0:15, '_se')
library(ggplot2)
data_preds1 <- cbind(preds1$state$pred, new_data, se )
ggplot(data_preds1, aes(x = t)) +
  geom_line(aes(y = x1b, color = 'Predicted')) + 
  geom_point(aes(y = X10, color = 'Observed')) + 
  # geom_line(aes(y = FbotIn / 20, color = 'FbotIn')) +
  # geom_line(aes(y = FbotOut / 20, color = 'FbotOut')) +
  # geom_line(aes(y = FtopIn / 20, color = 'FtopIn')) +
  # geom_ribbon(aes(ymin = x9 - X9_se, ymax = x9 + X9_se), fill = 'blue', alpha = 0.2) +
  theme_minimal()



# plot(preds$state$pred$x9, type = 'l', col = 'blue')
# points(new_data$X9, col = 'red')
# # # # # Residuals
ggplot(data_preds1, aes(x = t), ) + 
  geom_point(aes(y = X10 - x1b), color = 'blue') +
  geom_hline(yintercept = 0, color = 'red') +
  ylab('Residuals') +
  theme(text = element_text(size = 20))

fit$model

# # # fit$xm[['C0']]
# # # fit$xm[['C1']]
# # # fit$xm[['C2']]
# # # fit$xm[['C3']]
# # # fit$xm[['C4']]
# # # load(paste0("models/", model, "/", name, ".RData"))
# # new_data <- cbind(water_sensors, inputs)[7200:8000,]
# # new_data <- new_data[, !duplicated(colnames(new_data))]
# # # # dropna
# # new_data <- new_data[complete.cases(new_data),]

new_data <- cbind(water_sensors, inputs)[7200:7224,]
new_data <- new_data[, !duplicated(colnames(new_data))]
# # # # dropna
new_data <- new_data[complete.cases(new_data),]

sim <- simulate(fit, newdata = new_data)
plot(sim$state$sim$xa, type = 'l', col = 'blue')
points(new_data$X0, col = 'red')

# length(sim$state$sim$x3)
# length(new_data$X3)
# summary(fit)


