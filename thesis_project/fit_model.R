library(ctsmr)
rm(list=ls())
# Load the base model
# Load the data
source('thesis_project/R_scripts/makedata.R')
source('thesis_project/R_scripts/makefit.R')
source('thesis_project/R_scripts/nllikelihood.R')

# makemodel
# layer <- 1
# model_types <- c('linear_lag2', 'linear_lag3', 'linear1_lag1', 'linear1_lag2', 'linear1_lag3', 'nonlinear_lag1', 'nonlinear_lag2', 'nonlinear_lag3', 'nonlinear1_lag1', 'nonlinear1_lag2', 'nonlinear1_lag3')

fit_model <- function(layer, model_type, model_name, fill_state = TRUE){
    source('thesis_project/models/base_model.R')

    state <- paste0('X', layer)
    if (layer == 0){
        inputs <- paste0('X', layer+1)
    } else if (layer < 15){
        inputs <- c(paste0('X', layer-1), paste0('X', layer+1))
    } else {
        inputs <- paste0('X', layer-1)
    }




    print(inputs)
    source(paste0('thesis_project/models/', model_name, '.R'))
    model <- make_model(model, type = model_type)

    if (fill_state){
           inputs <- c(inputs,paste0('X',layer))
    } else{
        model_type <- paste0(model_type, 'missing_obs')
    }    # Data
    data <- make_data(start = '2017-04-01', end = '2017-11-01', fill=inputs)
    model <- setInitialState(model, data)

    model$options$solutionMethod <- 0
    model$options$eps <- 1e-6
    model$options$nIEKF <- 1
    model$options$initialVarianceScaling <- 1e-6
    model$options$odeeps <- 1e-3


    # fit <- makefit(model)
    
    # pars <- fit$xm[c(1,3,4,5,6,7,8,9,10,11)]
    # print(pars)
    # nllikelihood(pars, fit, data, firstorder=TRUE, c=3, n.ahead=1, printit=TRUE)

    # res <- nlminb(start = pars, objective = nllikelihood, fit = fit, D = data, firstorder = TRUE, c = 3, n.ahead = 1, printit = FALSE, lower = 0, upper = 100)
    fit <- model$estimate(data = data, firstorder = TRUE)
    # fit$xm[2:(length(fit$xm)-1)] <- res$par
    # fit$data <- data
    model_path <- paste0('models/', model_name, '/')
    dir.create(model_path, showWarnings = FALSE)

    

    save(fit, file = paste0('models/', model_name, '/', model_type,  '.RData'))


    ## EVALUATION ##
    aic <- -2*fit$loglik + 2*(length(fit$xm)-2)
    aic


    p <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred
    sd <- p$output$sd

    res <- ((pred[,1])-(data[,layer+2]))


    times <- strptime(data$X, format = '%Y-%m-%d %H:%M:%S')
    hours <- format(times, format = '%H')
    hours <- as.numeric(hours)


    res_interval <- function(interval){
        par(mfrow=c(4,1))
        plot(res[interval], type = 'o')
        abline(h=0, col='red')
        acf(res[interval], na.action = na.pass)
        plot(res[interval]~hours[interval])
        abline(h=0, col='red')
        qqnorm(res[interval])
        qqline(res[interval])
    }
    res_interval(100:7000)
    return(fit)
}




layers <- 5
model_types <- c('test')

for (layer in layers){
    print(layer)
    for (model_type in model_types){
        print(model_type)
        model_name <- paste0('simple', layer)

        fit<-fit_model(layer, model_type, model_name = model_name)
        print(summary(fit))
    }
}
fit
summary(fit)
# layer <- 1
# model_types <- c('nonlinear1_lag3')
# for (model_type in model_types){
#     print(model_type)
#     fit<-fit_model(layer, model_type)
#     print(summary(fit))
# }
# # summary(fit)
# inputs <- c('X2','X4')

# pred <- p$output$pred
# sd <- p$output$sd

# res <- ((pred[,1])-(data[,layer+2]))


# times <- strptime(data$X, format = '%Y-%m-%d %H:%M:%S')
# hours <- format(times, format = '%H')
# hours <- as.numeric(hours)


# res_interval <- function(interval){
#     par(mfrow=c(5,1))
#     plot(res[interval], type = 'o')
#     abline(h=0, col='red')
#     plot(data$FbotIn[interval], type = 'o')

#     acf(res[interval], na.action = na.pass)
#     plot(res[interval]~hours[interval])
#     abline(h=0, col='red')
#     qqnorm(res[interval]/sd[interval,])
#     qqline(res[interval]/sd[interval,])
# }
# res_interval(50:7000)
# summary(fit)

