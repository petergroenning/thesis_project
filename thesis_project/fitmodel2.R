library(ctsmr)
rm(list=ls())
# Load the base model
# Load the data
source('thesis_project/R_scripts/makedata.R')
source('thesis_project/models/base_model.R')
source('thesis_project/R_scripts/makefit.R')
source('thesis_project/R_scripts/nllikelihood.R')

model_name <- 'model2'
source(paste0('thesis_project/models/', model_name, '.R'))
model <- make_model(model)


load_model <- function(model_name){
    load(model_name)
    return(fit)
}

data <- make_data(start = '2017-02-01', end = '2017-12-01', fill = c('X10','X11','X12','X13','X14','X15'))

# old_fit <- load_model('models/simple0.Rdata')
model <- setInitialState(model, data)
# model$estimate(data = data, firstorder = TRUE)
# model <- setParams(model, old_fit)

set_parameters <- function(fit){
    for (x in 0:15){
        model_name <- paste0('simple', x)
        old_fit <- load_model(paste0('models/', model_name, '/model.RData'))

        names_new <- c(paste0('f',x), paste0('k',x), paste0('a',x),  paste0('sigma_x',x), paste0('sigma_X',x), paste0('sigma_y',x), paste0('Y',x,'m0'),paste0('x',x,'m0'))
        names_old <- c('f','k','a','sigma_x','sigma_X','sigma_y', paste0('Y',x,'m0'),paste0('x',x,'m0'))
        if (x == 0){
            names_new <- c(names_new, 'v1bot', 'v2bot')
            names_old <- c(names_old, 'v1', 'v2')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else if (x == 15){
            names_new <- c(names_new,'v1top', 'v2top','u')
            names_old <- c(names_old,'v1', 'v2','u')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else if (x==10){
            names_old <- c(names_old,'fmid','v')
            names_new <- c(names_new,'fmid','v10')
            fit$xm[names_new] <- old_fit$xm[names_old]
        } else {
            names_new <- c(names_new,paste0('v',x))
            names_old <- c(names_old,'v')
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


eval <- function(interval = 200:7000, name = 'X0'){
    res <- pred[name] - data[name]
    res <- res[interval,]

    par(mfrow=c(5,1))
    plot(pred[interval,name], type = 'o')
    points(data[interval,name], col='red')
    plot(res)
    abline(h=0, col='red')
    acf(res, na.action = na.pass)
    pacf(res, na.action = na.pass)
    qqnorm(res)
    qqline(res)
}
eval(interval= 100:10000, name = 'X15')


par(mfrow = c(1,1))

plot(pred$X15[500:1000], type = 'o')
points(data$X15[500:1000], col='red')


s <- predict(fit, newdata = data[1:100,], firstorderinputinterpolation=TRUE, n.ahead = 24)

# nllikelihood(params, fit, data, param_names = names, firstorder=TRUE, c=3, n.ahead=1, printit=TRUE)
# res <- nlminb(params, nllikelihood, fit = fit, D=data, param_names = names, firstorder=TRUE, c=3, n.ahead=1, printit=TRUE,
#                 lower = c(0,0,0,0,0,0,0,0,0), upper = c(10,10,10,10,10,10,10,10,10), control = list(trace = 3))


# par <- c(3.965,-2.112, 0.037,-0.598, 0.073,4.404, 0.050, 0.016, 0.069)
# params
# fit$xm[names] <- par
# p<- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)

# r <- p$output$pred$X0 - data$X0
# r <- r[10:7000]
# acf(r, na.action = na.pass, main = 'ACF')
# plot(r)

