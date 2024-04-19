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


fit <- load_model('models/coupled/model2.Rdata')

p <- predict(fit, newdata = fit$data[[1]], firstorderinputinterpolation=TRUE, n.ahead = 1)
pred <- p$output$pred
state <- p$state$pred
sds <- p$output$sd
res <- pred - fit$data[[1]][c("X0", "X1", "X2",'X3', 'X4', 'X5','X6', 'X7', 'X8', 'X9', 'X10', 'X11', 'X12', 'X13', 'X14', 'X15')]



# # save residuls and sd
# res_sd <- data.frame(res = r, sd = sd)
# dir.create(paste0('residuals/', model_name), showWarnings = FALSE)
# write.csv(res_sd, paste0('residuals/', model_name, '/', type, '.csv'))




res_interval <- function(interval, names = c('X1')){

    par(mfcol=c(3,length(names)))
    for (name in names){
        r <- res[,name]
        sd <- sds[,name]
        print(name)
        plot(r[interval], main = paste0('Residuals ', name))
        abline(h=0, col='red')
        acf(r[interval], na.action = na.pass, main = paste0('ACF ', name), ylim=c(-0.2,1))
        abline(h=0, col='red')
        qqnorm(r[interval]/sd[interval], main = paste0('QQ plot ', name))
        qqline(r[interval]/sd[interval])

        # mse <- mean(r[interval]^2, na.rm = TRUE)
        # print(paste0('MSE ', name, ': ', mse))

        # aic <- -2*fits[[type]]$loglik + 2*(length(fits[[type]]$xm)-2)
        # print(paste0('AIC ', type, ': ', aic))

    }
  
}
res_interval(1:100, c('X1')) 
summary(fit)

par(mfrow = c(1,1))
interval <- 100:150
s <- state[interval,]
plot(s$Y0m, type = 'o', col = 'blue')
lines(s$Y1m, col = 'red', type = 'o')
lines(s$Y2m, type = 'o')

points(s$x0m, col = 'purple')
points(s$x1m, col = 'orange')
points(s$x2m, col = '#04fd25')



plot(s$Y0m-s$x0m, type = 'o')
lines(c(NA,diff(s$x0m))/3, type = 'o',col ='green')
