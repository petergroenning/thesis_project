library(ctsmr)
library(jsonlite)
rm(list=ls())
# Load the base model
source('thesis_project/models/base_model.R')
# Makefit
source('thesis_project/R_scripts/makefit.R')

# Likelihood function
source('thesis_project/R_scripts/nllikelihood.R')


# Data
source('thesis_project/R_scripts/makedata.R')
load_model <- function(model_name){
    load(model_name)
    return(fit)
}


layer <- 15
types <- c('linear','model','nonlinear2_lag3')
model_name <- paste0('simple', layer)

fits <- list()
preds <- list()
sds <- list()
res <- list()
for (type in types){
    path <- paste0('models/', model_name, '/', type, '.RData')
    fit <- load_model(path)
    fits[[type]] <- fit
    p <- predict(fit, newdata = fit$data[[1]], firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred[,1]
    sd <- p$output$sd[,1]
    r <- pred - fit$data[[1]][paste0('X', layer)]
    preds[[type]] <- pred
    sds[[type]] <- sd
    res[[type]] <- r[,1]
}


res_interval <- function(interval){
    par(mfcol=c(3,length(types)))
    for (type in types){
        r <- res[[type]]
        print(type)
        plot(r[interval], type = 'o', main = paste0('Residuals ', type))
        abline(h=0, col='red')
        acf(r[interval], na.action = na.pass, main = paste0('ACF ', type), ylim=c(-0.2,1))
        abline(h=0, col='red')
        qqnorm(r[interval], main = paste0('QQ plot ', type))
        qqline(r[interval])

        mse <- mean(r[interval]^2, na.rm = TRUE)
        print(paste0('MSE ', type, ': ', mse))

        aic <- -2*fits[[type]]$loglik + 2*(length(fits[[type]]$xm)-2)
        print(paste0('AIC ', type, ': ', aic))

    }
  
}
res_interval(60:7000) 

# for (type in types){
#     print(summary(fits[[type]]))
# }
