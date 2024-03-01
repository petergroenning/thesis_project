library(ctsmrTMB)
library(ggplot2)
library(patchwork)
source('thesis_project/models/ctsm_utils.R')


make_model <- function(name, data, old_fit){

    obj <- ctsmrTMB$new()
    # Set name of model (and the created .cpp file)
    obj$set_modelname(name)

    if (missing(old_fit)){
        old_fit <- NULL
    }

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
    obj$add_systems(dx0 ~ 1/C *  ((x1 - x0) * (R + FbotOut * Cbot)) * dt + (Tbot - x0) * FbotIn * Vbot * dt + sigma_x0 * dw0)
    obj$add_systems(dx1 ~ 1/C *  ((x0 - x1) * (R + FbotOut * Cbot) + (x2 - x1) * R) * dt + sigma_x1 * dw1)
    obj$add_systems(dx2 ~ 1/C *  ((x1 - x2) * R + (x3 - x2) * R) * dt + sigma_x2 * dw2)
    obj$add_systems(dx3 ~ 1/C *  ((x2 - x3) * R + (x4 - x3) * R) * dt + sigma_x3 * dw3)
    obj$add_systems(dx4 ~ 1/C *  ((x3 - x4) * R + (x5 - x4) * R) * dt + sigma_x4 * dw4)
    obj$add_systems(dx5 ~ 1/C *  ((x4 - x5) * R + (x6 - x5) * R) * dt + sigma_x5 * dw5)
    obj$add_systems(dx6 ~ 1/C *  ((x5 - x6) * R + (x7 - x6) * R) * dt + sigma_x6 * dw6)
    obj$add_systems(dx7 ~ 1/C *  ((x6 - x7) * R + (x8 - x7) * R) * dt + sigma_x7 * dw7)
    obj$add_systems(dx8 ~ 1/C *  ((x7 - x8) * R + (x9 - x8) * R) * dt + sigma_x8 * dw8)
    obj$add_systems(dx9 ~ 1/C *  ((x8 - x9) * R + (x10 - x9) * (R + FmidOut * Cmid)) * dt + sigma_x9 * dw9)
    obj$add_systems(dx10 ~ 1/C * ((x9 - x10) * (R + FmidOut * Cmid) + (x11 - x10) * (R + FmidOut * Cmid)) * dt + (Tmid - x10) * FmidIn * Vmid * dt + sigma_x10 * dw10)
    obj$add_systems(dx11 ~ 1/C *  ((x10 - x11) * (R + FmidOut * Cmid) + (x12 - x11) * R) * dt + sigma_x11 * dw11)
    obj$add_systems(dx12 ~ 1/C *  ((x11 - x12) * R + (x13 - x12) * R) * dt + sigma_x12 * dw12)
    obj$add_systems(dx13 ~ 1/C *  ((x12 - x13) * R + (x14 - x13) * R) * dt + sigma_x13 * dw13)
    obj$add_systems(dx14 ~ 1/C *  ((x13 - x14) * R + (x15 - x14) * (R + FtopOut * Ctop)) * dt + sigma_x14 * dw14)
    obj$add_systems(dx15 ~ 1/C *  ((x14 - x15) * (R + FtopOut * Ctop)) * dt + (Ttop - x15) * FtopIn * Vtop *dt + sigma_x15 * dw15)

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
    lower <- 1e-10
    upper <- 100
    init <- 1
    obj$add_parameters(Ctop = c(init = get_TMB_parameter_value(old_fit, "Ctop", init), lb = lower, ub = upper))
    obj$add_parameters(Cmid = c(init = get_TMB_parameter_value(old_fit, "Cmid", init), lb = lower, ub = upper))
    obj$add_parameters(Cbot = c(init = get_TMB_parameter_value(old_fit, "Cbot", init), lb = lower, ub = upper))
    obj$add_parameters(R = c(init = get_TMB_parameter_value(old_fit, "R", init), lb = lower, ub = upper))
    obj$add_parameters(C = c(init = get_TMB_parameter_value(old_fit, "C", init), lb = lower, ub = upper))
    obj$add_parameters(Vtop = c(init = get_TMB_parameter_value(old_fit, "Vtop", init), lb = lower, ub = upper))
    obj$add_parameters(Vmid = c(init = get_TMB_parameter_value(old_fit, "Vmid", init), lb = lower, ub = upper))
    obj$add_parameters(Vbot = c(init = get_TMB_parameter_value(old_fit, "Vbot", init), lb = lower, ub = upper))


    obj$add_inputs(Ttop, Tmid, Tbot, FtopIn, FtopOut, FmidIn, FmidOut, FbotIn, FbotOut)

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
    # obj$add_algebraics(Fbot ~ FbotIn + FbotOut)
    # obj$add_algebraics(Fmid ~ FmidIn + FmidOut)
    # obj$add_algebraics(Ftop ~ FtopIn + FtopOut)

    # # Set system noise
    lower <- 1e-7
    upper <- 4
    init <- 1
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
    lower <-1e-2
    upper <- 3
    init <- 1e-5
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

    x0 <- c(data$X0[1], data$X1[1], data$X2[1], data$X3[1], data$X4[1], data$X5[1], data$X6[1], data$X7[1], data$X8[1], data$X9[1], data$X10[1], data$X11[1], data$X12[1], data$X13[1], data$X14[1], data$X15[1])
    N <- length(x0)
    init_noise <- 1e-5
    obj$set_initial_state(x0, init_noise*diag(N))
    return(obj)
}


###########################################################
# Make data
############################################################
water_sensors <- read.csv('data/processed/dronninglund/water_sensors.csv')
inputs  <- read.csv('data/processed/dronninglund/inputs.csv') 
head(water_sensors)
data <- cbind(water_sensors, inputs)[7200:7500,]
data <- data[, !duplicated(colnames(data))]

# Get old fit
name <- 'CTSMsimple'
fit <- load_best_fit(name)

# Make model object
obj <- make_model(name = name, data = data)


# # # Fit model
fit <- obj$estimate(data,
                    method = 'ekf',
                    ode.solver = 'euler',
                    use.hessian = TRUE,
                    ode.timestep = 1,
                    compile = FALSE,
                    control = list(maxit = 10,
                                    tol = 1e-6,
                                    rel.tol = 1e-15,
                                    sing.tol = 1e-20,
                                    trace = 1, 
                                    eval.max = 10))


fit
save_TMB_fit(fit, name)



# x0 <- c(data$X0[1], data$X1[1], data$X2[1], data$X3[1], data$X4[1], data$X5[1], data$X6[1], data$X7[1], data$X8[1], data$X9[1], data$X10[1], data$X11[1], data$X12[1], data$X13[1], data$X14[1], data$X15[1])

# # p <- obj$predict(data, k.ahead = 24, x0 = x0, ode.solver = 'rk4')


# nstep <- 24

# mask <- p$k.ahead == nstep

# true <- data$X14[(1+nstep):(nrow(data))]
# pred <- p[mask, ]$x14
# plot(pred - true)


# plot(pred, type = 'l', col = 'black')
# points(true, col = 'red')




# plot(fit$states$sd$prior$x0)

# summary(fit)
# fit$nll

