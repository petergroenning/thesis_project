library(ctsmr)
rm(list=ls())
# Load the base model
# Load the data
source('thesis_project/R_scripts/makedata.R')
source('thesis_project/models/base_model.R')
source('thesis_project/R_scripts/makefit.R')
source('thesis_project/R_scripts/nllikelihood.R')

model_name <- 'model3'
source(paste0('thesis_project/models/', model_name, '.R'))
model <- make_model(model)


load_model <- function(model_name){
    load(model_name)
    return(fit)
}

data <- make_data(start = '2017-02-01', end = '2017-12-01', fill = c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14'))

# old_fit <- load_model('models/simple0.Rdata')
model <- setInitialState(model, data)
# model$estimate(data = data, firstorder = TRUE)
# model <- setParams(model, old_fit)

set_parameters <- function(fit){
    for (x in 0:15){
        model_name <- paste0('simple', x)
        if (x == 0){
            old_fit <- load_model(paste0('models/', model_name, '/model3.RData'))
        } else if (x == 15){
            old_fit <- load_model(paste0('models/', model_name, '/model3.RData'))
        } else {
            old_fit <- load_model(paste0('models/', model_name, '/model4.RData'))
        }

        names_new <- c()
        names_old <- c()
        if (x == 0){
            names_new <- c(names_new, 'v0','f10','k0')
            names_old <- c(names_old, 'v','f1','k1')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else if (x == 15){
            names_new <- c(names_new,'v15','u','f215','k15')
            names_old <- c(names_old,'v','u','f2','k1')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else if (x==10){
            names_old <- c(names_old,'fmid')
            names_new <- c(names_new,'fmid')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else {
            names_new <- c(paste0('f1',x), paste0('f2',x), paste0('a',x), paste0('b',x),  paste0('sigma_x',x), paste0('sigma_X',x), paste0('sigma_y',x), paste0('Y',x,'m0'),paste0('x',x,'m0'))
            names_old <- c('f1','f2','a','b','sigma_x','sigma_X','sigma_y', paste0('Y',x,'m0'),paste0('x',x,'m0'))
            fit$xm[names_new] <- old_fit$xm[names_old]
        }

        fit$xm[names_new] <- old_fit$xm[names_old]

    # if (x == 0){
    #     names_new <- c('v1bot', 'v2bot')
    #     names_old <- c('v1', 'v2')
    #     fit$xm[names_new] <- old_fit$xm[names_old]
    # } else if (x == 15){
    #     names_new <- c('v1top', 'v2top')
    #     names_old <- c('v1', 'v2')
    #     fit$xm[names_new] <- old_fit$xm[names_old]
    # } else if (x==10){
    #     names <- c('fmid')
    #     fit$xm[names] <- old_fit$xm[names]
    # }
    }


return(fit)
}
fit <- makefit(model)
fit <- set_parameters(fit)

P <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
pred <- P$output$pred


eval <- function(interval = 200:7000, name = 'X15'){
    res <- pred[name] - data[name]
    res <- res[interval,]

    par(mfrow=c(4,1))
    plot(res)
    abline(h=0, col='red')
    acf(res, na.action = na.pass)
    pacf(res, na.action = na.pass)
    qqnorm(res)
    qqline(res)
}
eval(interval= 100:7000, name = 'X4')


par(mfrow = c(1,1))

plot(pred$X15[500:1000], type = 'o')
points(data$X15[500:1000], col='red')




s <- simulate(fit, newdata = data, firstorderinputinterpolation=TRUE)$state$sim
plot(s$x15m)
points(data$X15, col='red')
