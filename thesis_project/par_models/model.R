library(ctsmr)
source('thesis_project/R_scripts/makefit.R')
source('thesis_project/R_scripts/nllikelihood.R')
load_data <- function(start = '201', end = '2020'){
    data <- read.csv('data/processed/water_normalised.csv', header = TRUE, sep = ',')
    names(data) <- c('time', paste0('X', 1:(ncol(data)-1)))

    data <- data[which(data$time >= start & data$time <= end),]

    data <- data[,2:ncol(data)]
    data$t <- 1:nrow(data)
    return(data)
}
data <- load_data('2017-01-01', '2017-12-31')
model <- ctsm$new()


# Add Inputs
for (i in 1:32){
    Y <- paste0('X',i)
    x <- i / 32
    Ym <- paste0('b0+b1/(1+exp((m1-',x,')*k1))+b2/(1+exp((m2-',x,')*k2))')
    model$addObs(as.formula(paste0(Y,'~',Ym)))
    model$setVariance(as.formula(paste0(Y,'~sigma_Y^2')))
}
model$setParameter(sigma_Y = c(init=0.1,lb=0,ub=1))

model$addSystem(db0~sigma_b*dw0)
model$addSystem(db1~sigma_b*dw1)
model$addSystem(db2~sigma_b*dw2)
model$addSystem(dm1~sigma_m*dw3)
model$addSystem(dm2~sigma_m*dw4)
model$addSystem(dk1~sigma_k*dw5)
model$addSystem(dk2~sigma_k*dw6)


model$setParameter(b0 = c(init =0.1, lower = 0, upper = 1))
model$setParameter(b1 = c(init =0.1, lower = 0, upper = 1))
model$setParameter(b2 = c(init =0.2, lower = 0, upper = 1))
model$setParameter(m1 = c(init =0.1, lower = 0, upper = 1))
model$setParameter(m2 = c(init =0.4, lower = 0, upper = 1))
model$setParameter(k1 = c(init =10, lower = 0.2, upper = 100))
model$setParameter(k2 = c(init =10, lower = 0.2, upper = 100))

model$setParameter(sigma_b = c(init =1e-2, lower = 0, upper = 1))
model$setParameter(sigma_m = c(init =1e-2, lower = 0, upper = 1))
model$setParameter(sigma_k = c(init =1e-2, lower = 0, upper = 1))


model 
# fit <- model$estimate(data)
fit <- makefit(model)
fit$xm
p <- predict(fit, data, n.ahead = 1, c = 0)

res <- nlminb(start = fit$xm, objective = nllikelihood, fit = fit, D = data, firstorder = TRUE, c = 3, n.ahead = 1, printit = TRUE, lower = 0, upper = 100)

fit$xm <- res$par

p <- predict(fit, data, n.ahead = 1, c = 0)


i <- 3

plot_p <- function(i){
    x <- 1:32/32
    plot(x,t(p$output$pred)[,i], type = 'o', col = 'red') 
    points(x,t(data)[1:32,i], col = 'blue')

}
plot_p(30)
