library(ctsmr)
library(forecast)
rm(list=ls())
# Load the base model
source('thesis_project/models/base_model.R')
# Makefit
source('thesis_project/R_scripts/makefit.R')

# Likelihood function
source('thesis_project/R_scripts/nllikelihood.R')

# makemodel
model_name <- 'simple4'
source(paste0('thesis_project/models/', model_name, '.R'))
model <- make_model(model)
# Data
source('thesis_project/R_scripts/makedata.R')
data <- make_data(start = '2017-04-01', end = '2017-12-01', splines = 5)




fit <- lm(X0~FbotIn*Tbot+FbotIn+FbotOut*X1+t, data = data)
summary(fit)
par(mfrow=c(2,1))
plot(fitted(fit)-data$X0)
acf(fitted(fit)-data$X0)





