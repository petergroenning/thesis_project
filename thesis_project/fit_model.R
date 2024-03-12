library(ctsmr)
rm(list=ls())
# Load the base model
source('thesis_project/ctsm_models/base_model.R')
# Makefit
source('thesis_project/R_scripts/makefit.R')

# Likelihood function
source('thesis_project/R_scripts/nllikelihood.R')

# makemodel
model_name <- 'model9'
old_model_name <- 'model5'
source(paste0('thesis_project/ctsm_models/', model_name, '.R'))
model <- make_model(model)

# Data
source('thesis_project/R_scripts/makedata.R')
data <- make_data(start = '2015-04-01', end = '2015-05-07', splines = 7)
# data <- read.csv('data/processed/ht/data_grouped.csv')
load_model<-function(path){
    load(path)
    return(fit)
}
fit <- load_model(paste0('models/CTSM/',old_model_name, '.RData'))


model <- setParams(model, fit)
model <- setInitialState(model, data)

model$options$eps <- 1e-3


fit <- model$estimate(data = data, firstorder = TRUE)

save(fit, file = paste0('models/CTSM/', model_name, '.RData'))
summary(fit)
