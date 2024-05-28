6


types <- c('2state_basic','2state_model1','2state_model3','2state_model4')
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
data$FbotInk1 <- c(0, data$FbotIn[1:(length(data$FbotIn)-1)])
# data$FbotOutk1 <- c(0, data$FbotOut[1:(length(data$FbotOut)-1)])
# data$FmidInk1 <- c(0, data$FmidIn[1:(length(data$FmidIn)-1)])
data$FtopInk1 <- c(0, data$FtopIn[1:(length(data$FtopIn)-1)])
data$dFtopIn <- c(0,diff(data$FtopIn))
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
    path <- paste0('models/', type, '.RData')
    fit <- load_model(path)
    fits[[type]] <- fit
    p <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 1)
    pred <- p$output$pred
    # state <- p$state$pred[,1]
    sd <- p$output$sd
    # r <- pred - D[c('X1','X2','X3','X4','X5','X6','X7','X8','X9')]
    r <- pred -  D[c('X1','X2','X3','X4','X5','X6')]
    # r <- pred -  D[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
    preds[[type]] <- pred
    sds[[type]] <- sd
    res[[type]] <- r
    # true[[type]] <- D[c('X1','X2','X3','X4','X5','X6','X7','X8','X9')]
    true[[type]] <- D[c('X1','X2','X3','X4','X5','X6')]
    # true[[type]] <- D[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]

    # save residuls and sd
    # res_sd <- data.frame(res = r, sd = sd)
    # dir.create(paste0('residuals/', model_name), showWarnings = FALSE)
    # write.csv(res_sd, paste0('residuals/', model_name, '/', type, '.csv'))
    
}

res_interval <- function(interval, name = 'X1'){
    par(mfcol=c(5,length(types)))
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

        cpgram(r[interval,][,name], main = paste0('CPgram ', type))


        mse <- mean(t(r[interval,]^2), na.rm = TRUE)
        print(paste0('MSE ', type, ': ', mse))

        aic <- -2*fits[[type]]$loglik + 2*(9)
        print(paste0('AIC ', type, ': ', aic))

        llik <- fits[[type]]$loglik
        print(paste0('Log Likelihood ', type, ': ', llik))
    }
  

}

res_interval(20:4000, 'X6')
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

# write.csv(r, 'residuals/2state_basic2_2.csv')
# 
par(mfrow = c(1,1))
ccf(r$X1, D$Fbot)
# abline(v=0, col='red')
x <- ccf(r$X6, D$FtopIn)
y <- ccf(r$X4, D$FmidIn)
z <- ccf(r$X1, D$FbotIn)
df <- data.frame(x$lag, x$acf, y$acf, z$acf)
write.csv(df, 'ccf.csv')
# x$series
# fits$'2state_basic2'

# sim <- simulate(fits$'2state_model3',  newdata = D, firstorderinputinterpolation=TRUE)
# pred <- predict(fits$'2state_model3', newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 7*24)
# s <- sim$output$sim
# p <- pred$output$pred

# e <- p - D[c('X1','X2','X3','X4','X5','X6')]

# plot(e$X4, type = 'l')
abline(h=0.029, col='red')

cp <- function (ts, taper = 0.1, main = paste("Series: ", deparse1(substitute(ts))), 
    ci.col = "blue") {
    main
    if (NCOL(ts) > 1) 
        stop("only implemented for univariate time series")
    x <- as.vector(ts)
    x <- x[!is.na(x)]
    x <- spec.taper(scale(x, TRUE, FALSE), p = taper)
    y <- Mod(fft(x))^2/length(x)
    y[1L] <- 0
    n <- length(x)
    x <- (0:(n/2)) * frequency(ts)/n
    if (length(x)%%2 == 0) {
        n <- length(x) - 1
        y <- y[1L:n]
        x <- x[1L:n]
    }
    else y <- y[seq_along(x)]
    xm <- frequency(ts)/2
    mp <- length(x) - 1
    crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
    oldpty <- par(pty = "s")
    on.exit(par(oldpty))
    plot(x, cumsum(y)/sum(y), type = "s", xlim = c(0, xm), ylim = c(0, 
        1), xaxs = "i", yaxs = "i", xlab = "frequency", ylab = "")
    lines(c(0, xm * (1 - crit)), c(crit, 1), col = ci.col, lty = 2)
    lines(c(xm * crit, xm), c(0, 1 - crit), col = ci.col, lty = 2)
    title(main = main)
    return(list(pgram = y, cum = cumsum(y)/sum(y), freq = x, 
        crit = crit))
}

x1 <- cp(r$X1)
x2 <- cp(r$X2)
x3 <- cp(r$X3)
x4 <- cp(r$X4)
x5 <- cp(r$X5)
x6 <- cp(r$X6)

df <- data.frame(x1$cum, x2$cum, x3$cum, x4$cum, x5$cum, x6$cum, x1$freq)
write.csv(df, 'cp.csv')
# sim <- simulate(fits$'2state_model3',  newdata = D, firstorderinputinterpolation=TRUE)
# pred <- predict(fits$'2state_model3', newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 7*24)
data <- read.csv('data/processed/data2.csv')


plot(data$X1, type = 'l')
D <- data[19500:26000,]
plot(D$X1)


fit <- fits$'2state_model3'
x0 <- as.numeric(D[1,2:7])
fit$xm
# names(x0) <- names(fit$xm[1:6])
fit$xm[1:6] <- x0
fit$xm[7:12] <- x0
fit$xm
p1d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 1*24)
p2d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 2*24)
p1w <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 7*24)
p2w <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 14*24)
p1m <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 30*24)
sim <- simulate(fit,  newdata = D, firstorderinputinterpolation=TRUE)


write.csv(p1d$output$pred, 'predictions/p1d.csv')
write.csv(p2d$output$pred, 'predictions/p2d.csv')
write.csv(p1w$output$pred, 'predictions/p1w.csv')
write.csv(p2w$output$pred, 'predictions/p2w.csv')
write.csv(p1m$output$pred, 'predictions/p1m.csv')
write.csv(sim$output$sim, 'predictions/sim.csv')

