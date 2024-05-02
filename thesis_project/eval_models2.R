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


layer <- 3
types <- c('linear','nonlinear_lag3','test1','test2','test3', 'test4')
# types <- c('model3','model4')

model_name <- paste0('simple', layer)

fits <- list()
preds <- list()
sds <- list()
res <- list()
true <- list()
for (type in types){
    path <- paste0('models/', model_name, '/', type, '.RData')
    fit <- load_model(path)
    fits[[type]] <- fit
    p <- predict(fit, newdata = fit$data[[1]], firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred[,1]
    # state <- p$state$pred[,1]
    sd <- p$output$sd[,1]
    r <- pred - fit$data[[1]][paste0('X', layer)]
    preds[[type]] <- pred
    sds[[type]] <- sd
    res[[type]] <- r[,1]
    true[[type]] <- fit$data[[1]][paste0('X', layer)]

    # save residuls and sd
    res_sd <- data.frame(res = r, sd = sd)
    dir.create(paste0('residuals/', model_name), showWarnings = FALSE)
    write.csv(res_sd, paste0('residuals/', model_name, '/', type, '.csv'))
    
}


res_interval <- function(interval){
    par(mfcol=c(4,length(types)))
    for (type in types){

        r <- res[[type]]
        sd <- sds[[type]]
        print(type)
        plot(r[interval], main = paste0('Residuals ', type))
        abline(h=0, col='red')
        acf(r[interval], na.action = na.pass, main = paste0('ACF ', type), ylim=c(-0.5,1))
        pacf(r[interval], na.action = na.pass, main = paste0('PACF ', type), ylim=c(-0.5,1))
        abline(h=0, col='red')
        qqnorm(r[interval]/sd[interval], main = paste0('QQ plot ', type))
        qqline(r[interval]/sd[interval])

        mse <- mean(r[interval]^2, na.rm = TRUE)
        print(paste0('MSE ', type, ': ', mse))

        aic <- -2*fits[[type]]$loglik + 2*(length(fits[[type]]$xm)-2)
        print(paste0('AIC ', type, ': ', aic))

    }
  
}

res_interval(20:3000) 
# par(mfrow = c(1,1))

# summary(fits$model3)
# summary(fits$model4)
# summary(fits$model5)

test_sim <- function(interval = 1:7000){
    par(mfcol = c(2,length(types)))
    
    for (type in types){
        print(type)
        s <- simulate(fits[[type]], newdata = fits[[type]]$data[[1]][interval,], firstorderinputinterpolation=TRUE)$state$sim[interval,2]
        x <- fits[[type]]$data[[1]][interval,paste0('X', layer)]
        plot(s, col = 'red', type = 'l', main = paste0('Simulated ', type))
        points(x, col = 'blue')
        plot(s-x, col = 'blue', type = 'o', main = paste0('Simulated ', type))
        mse <- mean((s-x)^2, na.rm = TRUE)
        print(paste0('MSE ', type, ': ', mse))
    }
}


test_sim(1:1000)


s <- simulate(fits$test3, newdata = fits$test2$data[[1]][1:3000,], firstorderinputinterpolation=TRUE)$state$sim[1:3000,1]
plot(s)
p <- predict(fits$test3, newdata = fits$test2$data[[1]][1:3000,], firstorderinputinterpolation=TRUE, n.ahead = 24*2*7)

r <- p$output$pred$X3 -fits$test3$data[[1]][1:3000,]$X3
r

plot(r)
# s1 <- simulate(fits$nl, newdata = fits$nl$data[[1]][1:3000,], firstorderinputinterpolation=TRUE)$state$sim[1:3000,1]
# par(mfrow = c(2,1))
# plot(s)
# lines(fits$nl2$data[[1]][1:3000,]$X2, col = 'red')
# points(s1)
# sc <- simulate(fits$model1, newdata = fits$model1$data[[1]][1:7000,], firstorderinputinterpolation=TRUE)$state$sim[1:7000,1]
# xc <- fits$modelc$data[[1]][1:7000,]$X2
# plot(sc-xc, col = 'blue',ylim = c(-6,10))
# sb <- simulate(fits$model1, newdata = fits$model1$data[[1]][1:7000,], firstorderinputinterpolation=TRUE)$state$sim[1:7000,1]
# xb <- fits$model2$data[[1]][1:7000,]$X2
# points(sb-xb, col = 'red')

# mean((sc-xc)^2)
# mean((sb-xb)^2)
# s1 <- simulate(fits$model1, newdata = fits$model1$data[[1]][1:7000,], firstorderinputinterpolation=TRUE)$state$sim[1:7000,1]
# s2 <- simulate(fits$model2, newdata = fits$model2$data[[1]][1:7000,], firstorderinputinterpolation=TRUE)$state$sim[1:7000,1]

# x1 <- fits$model1$data[[1]][1:7000,]$X2
# x2 <- fits$model2$data[[1]][1:7000,]$X2

# plot(s1-x1, col = 'blue')
# # points(s2-x2, col = 'red')

# # points(s1-x, col = 'red')
# # s <- simulate(fits$model, newdata = fits$model$data[[1]][1:7000,])
# # par(mfrow = c(2,1))
# # plot(s1, col = 'blue', type = 'o')
# # # points(fits$model$data[[1]][1:1000,]$X2, col = 'red')
# # plot(s1$state$sim$x1m, col = 'blue', type = 'l')
# # points(fits$model1$data[[1]][1:7000,]$X2, col = 'red')


# layers <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
# type <- 'model4'
# f1 <- c()
# f2 <- c()
# k1 <- c()
# a <- c()
# sx <- c()
# sy <- c()
# sX <- c()
# v <- c()
# for (layer in layers){
#     model <- paste0('simple', layer)
#     path <- paste0('models/', model, '/', type, '.RData')
#     fit <- load_model(path)
#     fits[[type]] <- fit
#     print(summary(fit))
#     f1 <- c(f1, unname(fit$xm['f1']))
#     f2 <- c(f2, unname(fit$xm['f2']))
#     k1 <- c(k1, unname(fit$xm['k1']))
#     a <- c(a, unname(fit$xm['a']))
#     # v <- c(v, unname(fit$xm['v']))
#     sx <- c(sx, unname(fit$xm['sigma_x']))
#     sy <- c(sy, unname(fit$xm['sigma_y']))
#     sX <- c(sX, unname(fit$xm['sigma_X']))
# }
# par(mfrow=c(1,1))
# plot(a, type = 'o')
# # par(mfrow=c(1,1))
# plot(f, type = 'o')
# plot(k, type = 'o')
# plot(v, type = 'o')
# plot(sy, type = 'o')
# print as np array



# sx
# par(mfrow = c(1,1))
# plot(0:15,(k1), type = 'o')
# plot(0:15,(f1), type = 'o')
# plot(1:15,(v[!is.na(v)]), type = 'o')
# plot(pred[4900:5100], type = 'o')
# points(state[4900:5100], col = 'blue')
# lines((pred-r[,1])[4900:5100], col = 'red')


