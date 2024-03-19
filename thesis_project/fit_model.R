library(ctsmr)
rm(list=ls())
# Load the base model
source('thesis_project/ctsm_models/base_model.R')
# Makefit
source('thesis_project/R_scripts/makefit.R')

# Likelihood function
source('thesis_project/R_scripts/nllikelihood.R')

# makemodel
model_name <- 'model9d'
old_model_name <- 'model9d'
source(paste0('thesis_project/ctsm_models/', model_name, '.R'))
model <- make_model(model)

# Data
source('thesis_project/R_scripts/makedata.R')
data <- make_data(start = '2015-04-01', end = '2015-05-07', splines = 7)
load_model<-function(path){
    load(path)
    return(fit)
}
#fit <- load_model(paste0('models/CTSM/',old_model_name, '.RData'))


#model <- setParams(model, fit)
model <- setInitialState(model, data)
model$options$eps <- 1e-2


fit <- model$estimate(data = data, firstorder = TRUE)
# fit <- makefit(model)
# save(fit, file = paste0('models/CTSM/', model_name, '.RData'))
summary(fit)
states<-c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')
data <- make_data(start = '2015-04-06', end = '2015-05-12', splines = 7)
obs<-data[,states]
n.ahead<-24
initial <- as.numeric(obs[1,])
names(initial) <- names(fit$xm[1:16])
sim <- simulate(fit, newdata = data, firstorderinputinterpolation=TRUE, x0 = initial)
write.csv(as.data.frame(sim$output$sim), 'sim.csv')
write.csv(data, 'data.csv')
