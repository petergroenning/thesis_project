6


types <- c('basic2','basic4','augmented_3','augmented_6')
load_model <- function(model_name){
    load(model_name)
    return(fit)
}
data <- read.csv('data/processed/data2.csv')
# data$dX1 <- c(0,diff(data$X1))
# data$dX2 <- c(0,diff(data$X2))
# data$dX3 <- c(0,diff(data$X3))
# data$dX4 <- c(0,diff(data$X4))
# data$dX5 <- c(0,diff(data$X5))
# data$dX6 <- c(0,diff(data$X6))
# data$FbotInk1 <- c(0, data$FbotIn[1:(length(data$FbotIn)-1)])
# data$FbotOutk1 <- c(0, data$FbotOut[1:(length(data$FbotOut)-1)])
# data$FmidInk1 <- c(0, data$FmidIn[1:(length(data$FmidIn)-1)])
# data$FtopInk1 <- c(0, data$FtopIn[1:(length(data$FtopIn)-1)])
# data$FtopOutk1 <- c(0, data$FtopOut[1:(length(data$FtopOut)-1)])
# data$Ttopk1 <- c(0, data$Ttop[1:(length(data$Ttop)-1)])
# data$Tmidk1 <- c(0, data$Tmid[1:(length(data$Tmid)-1)])
# data$Tbotk1 <- c(0, data$Tbot[1:(length(data$Tbot)-1)])

D <- data[3000:7500,]
fits <- list()
preds <- list()
sds <- list()
res <- list()
true <- list()
for (type in types){
    path <- paste0('models/ctsm2/', type, '.RData')
    fit <- load_model(path)
    fits[[type]] <- fit
    p <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred
    # state <- p$state$pred[,1]
    sd <- p$output$sd
    # r <- pred - D[c('X1','X2','X3','X4','X5','X6','X7','X8','X9')]
    r <- pred -  D[c('X1','X2','X3','X4','X5','X6')]
    preds[[type]] <- pred
    sds[[type]] <- sd
    res[[type]] <- r
    # true[[type]] <- D[c('X1','X2','X3','X4','X5','X6','X7','X8','X9')]
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

        mse <- mean(t(r[interval,]^2), na.rm = TRUE)
        print(paste0('MSE ', type, ': ', mse))

        aic <- -2*fits[[type]]$loglik + 2*(length(fits[[type]]$xm)-18)
        print(paste0('AIC ', type, ': ', aic))

    }
  

}

res_interval(20:4000, 'X1')

# r <- read.csv('residuals/ctsm2/augmented_3_1.csv')

# acf(r$X2[20:4000])
# par(mfrow = c(3,1))
# interval <- 3400:3475
# plot(D$X3[interval])
# lines(preds$basic5_1$X3[interval], col = 'red')
# lines(preds$model2$X3[interval], col = 'blue')
# lines(preds$model3$X3[interval], col = 'green')

# plot(diff(diff(preds$basic5_1$X3))[interval], type = 'o')
# abline(h=0)

# plot(D$Fbot[interval], type = 'l')
# abline(h=0)
