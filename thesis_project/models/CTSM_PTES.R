library(ctsmrTMB)
library(ggplot2)
library(patchwork)
rm(list = ls())
############################################################
# Make data
############################################################
data <- read.csv('data/processed/data.csv')
depths <- read.csv('data/processed/depths.csv')



############################################################
# Model creation and estimation
############################################################

# Create model object
obj = ctsmrTMB$new()

# Set name of model (and the created .cpp file)
obj$set_modelname("PTES")

# Set path where generated C++ files are saved.
# This will create a cppfiles folder in your current working directory if it doesnt exist
obj$set_cppfile_directory(".cppfiles")


#### ADD OBSERVATION EQUATIONS AND VARIANCES ####
min_depth = 0.5
obs_noise <- 0.1
step <- 6
depths <- depths[depths$DEPTH > min_depth,]

T1 <- 0
T2 <- 0.8
b <- 0
k1 <- 3
k2 <- 3

for (i in seq(1, nrow(depths), step)){
 
    name <- depths[i, 'POINT_NAME']
    d <- depths[i, 'DEPTH']

    # Add observation equations
    string = paste0(name, ' ~ ', T1, '/(1+exp(', k1,'*(',d,'-0))) + ', T2, '/(1+exp(', k2, '*(',d,'-a2))) + ', b)
    call = as.formula(string)
    obj$add_observations(call)

    # Add observation equation variances
    string = paste0(name, ' ~ sigma_', name)
    call = as.formula(string)
    obj$add_observation_variances(call)
    
  
}
obj
# Add observation variance parameters
rownames<-lapply(depths$POINT_NAME, function(x) paste0('sigma_',x))
initial <- rep(obs_noise, length(rownames))
lower <- rep(NA, length(rownames))
upper <- rep(NA, length(rownames))

params <- matrix(data = c(initial, lower, upper), 
                ncol = 3,
                byrow = FALSE,
                dimnames = list(rownames, c('initial', 'lower', 'upper')))

params <- params[seq(1, nrow(depths), step),]

obj$add_parameters(params)
obj

# Add inputs
obj$add_inputs(
  U_IN
)
obj$add_inputs(
  U_OUT
)




### ADD SYSTEM EQUATIONS ###
obj$add_systems(
  # da1 ~ (phi_a1 * U_IN + psi_a1 * U_OUT) * dt + sigma_a1 * dw,
  da2 ~ (phi_a2 * U_IN + psi_a2 * U_OUT) * dt + sigma_a2 * dw
  # dk1 ~ (phik1 * U_IN - psik1 * U_OUT) * dt + sigmak1 * dw,
  # dk2 ~ (phik2 * U_IN - psik2 * U_OUT) * dt + sigmak2 * dw
)



# Specify algebraic relations
# obj$add_algebraics(
#   phi_a1  ~ exp(logphi_a1),
#   psi_a1  ~ exp(logpsi_a1),

#   phi_a2  ~ exp(logphi_a2),
#   psi_a2  ~ exp(logpsi_a2)

  # phik1  ~ exp(logphik1),
  # psik1  ~ exp(logpsik1),

  # phik2  ~ exp(logphik2),
  # psik2  ~ exp(logpsik2)

# )

system_noise <- 1e-6

# Specify parameter initial values and lower/upper bounds in estimation
obj$add_parameters(
      # T1 = 0.2,
      # T2 = 0.8,
      # b = 0,
      # k1 = 1,
      # k2 = 1,
      # phi_a1 = c(init = 0, lower=-1, upper=1),
      # psi_a1 = c(init=0, lower=-1, upper=1),
      phi_a2 = c(init = 0, lower=-1, upper=1),
      psi_a2 = c(init=0, lower=-1, upper=1),
      # logphik1 = log(c(init = 1, lower=1e-5, upper=50)),
      # logpsik1 = log(c(init=1, lower=1e-5, upper=50)),
      # logphik2 = log(c(init = 1, lower=1e-5, upper=50)),
      # logpsik2 = log(c(init=1, lower=1e-5, upper=50)),
      # sigma_a1 = system_noise,
      sigma_a2 = system_noise
      # sigmak1 = system_noise,
      # sigmak2 = system_noise
)

obj
data_sub <- data[,2:(ncol(data)-3)]
U_IN <- data$U_IN
U_OUT <- data$U_OUT

# U_IN <- (U_IN - min(U_IN)) / (max(U_IN) - min(U_IN))
# U_OUT <- (U_OUT - min(U_OUT)) / (max(U_OUT) - min(U_OUT))


norm_data <- (data_sub - 50) / (90 - 50)
head(norm_data)
new_data <- cbind(data['t'], norm_data, U_IN, U_OUT)

# Set initial state mean and covariance
obj$set_initial_state(c(7), system_noise*diag(1))

# Carry out estimation using extended kalman filter method with stats::nlminb as optimizer
fit <- obj$estimate(data=new_data, method="ekf", ode.solver="rk4", 
                  compile = FALSE,
                  use.hessian=TRUE,
                  control = list(maxit = 1000, trace = 1, reltol = 1e-3, abstol = 1e-3))

fit$par.fixed

fit$states$mean$prior
# pred <- obj$predict(data=new_data, k.ahead=10, method="ekf", ode.solver="rk4")
