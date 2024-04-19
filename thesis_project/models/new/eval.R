


types <- c('lowpass')

load_model <- function(model_name){
    load(model_name)
    return(fit)
}
data <- read.csv('data/processed/data2.csv')
D <- data[1:5000,]
fits <- list()
preds <- list()
sds <- list()
res <- list()
true <- list()
for (type in types){
    path <- paste0('models/ctsm/', type, '.RData')
    fit <- load_model(path)
    fits[[type]] <- fit
    p <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred
    # state <- p$state$pred[,1]
    sd <- p$output$sd
    r <- pred - D[c('X1','X2','X3','X4','X5','X6')]
    preds[[type]] <- pred
    sds[[type]] <- sd
    res[[type]] <- r
    true[[type]] <- D[c('X1','X2','X3','X4','X5','X6')]

    # save residuls and sd
    # res_sd <- data.frame(res = r, sd = sd)
    # dir.create(paste0('residuals/', model_name), showWarnings = FALSE)
    # write.csv(res_sd, paste0('residuals/', model_name, '/', type, '.csv'))
    
}


res_interval <- function(interval, name = 'X1'){
    par(mfcol=c(4,length(types)))
    for (type in types){

        r <- res[[type]]
        sd <- sds[[type]]
        print(type)
        plot(r[interval,][,name], main = paste0('Residuals ', type))
        abline(h=0, col='red')
        acf(r[interval,][,name], na.action = na.pass, main = paste0('ACF ', type), ylim=c(-0.5,1))
        pacf(r[interval,][,name], na.action = na.pass, main = paste0('PACF ', type), ylim=c(-0.5,1))
        abline(h=0, col='red')
        qqnorm(r[interval,][,name]/sd[interval,][,name], main = paste0('QQ plot ', type))
        qqline(r[interval,][,name]/sd[interval,][,name])

        # mse <- mean(r[interval,]^2, na.rm = TRUE)
        # print(paste0('MSE ', type, ': ', mse))

        aic <- -2*fits[[type]]$loglik + 2*(length(fits[[type]]$xm)-2)
        print(paste0('AIC ', type, ': ', aic))

    }
  
}
res_interval(1:2000, 'X6') 

names(fits$model4$data[[1]])



# p <- predict(fits$model3, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 24)




# state <- p$state$pred

# D
# par(mfrow = c(1,1))
# plot(t(state[6000,1:6]), type = 'o')
# points(t(D[6000,2:7]), col = 'red', pch = 16)
# # state[200,7:12]

# res <- ((state[1:6])-(D[,2:7]))
# plot(res$x6)
