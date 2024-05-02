library(ctsmr)
source('thesis_project/R_scripts/makefit.R')
data <- read.csv('data/processed/1m_data.csv')

load_model <- function(model_path) {
  load(model_path)
  return(fit)
}

data <- data[3000:7500,]
model1 <- load_model('models/ctsm2/basic2.RData')
model2 <- load_model('models/ctsm2/augmented_3.RData')
model3 <- load_model('models/ctsm2/augmented_5.RData')
model4 <- load_model('models/ctsm2/augmented_6.RData')


model1$Name <- 'basic2_new'
model2$Name <- 'augmented_3_new'
model3$Name <- 'augmented_5_new'
model4$Name <- 'augmented_6'

models <- list(model4)
horizons <- c(1, 24, 48, 7*24, 2*7*24, 4*7*24)




for (model in models) {
  for (horizon in horizons) {
    fit <- model
    # fit <- makefit(model$model)
    p <- predict(fit, newdata = data, firstorderinterpolation = TRUE, n.ahead = horizon)
    # r <- p$output$pred - data[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
    # write.csv(r, paste0('residuals/', model$Name, '_', horizon, '.csv'))
    r <- p$output$pred - data[c('X1','X2','X3','X4','X5','X6')]
    write.csv(r, paste0('residuals/ctsm2/', model$Name, '_', horizon, '.csv'))
  }
  s <- simulate(fit, newdata = data, firstorderinterpolation = TRUE)
  # write.csv(r, paste0('residuals/', model$Name, '_sim.csv'))
  # r <- s$output$sim - data[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
  r <- s$output$sim - data[c('X1','X2','X3','X4','X5','X6')]
write.csv(r, paste0('residuals/ctsm2/', model$Name, '_sim.csv'))
}
# 


#########
# p1 <- predict(model2, data, firstorderinterpolation = TRUE, n.ahead = 1,c=1)
# p2 <- predict(model3, data, firstorderinterpolation = TRUE, n.ahead = 1,c=1)
# r1 <- p1$output$pred - data[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
# r2 <- p2$output$pred - data[c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')]
# plot(r1$X11[interval]^2-r2$X11[interval]^2, type='l')

# plot(p1$output$pred$X10[interval], type='l')
# lines(p2$output$pred$X10[interval], col='red')
# points(data$X10[interval], col='blue')

# mean(t(r1)^2)
# mean(t(r2)^2)

# interval<-3000:4000
# par(mfrow = c(2,1))
# acf(r1$X10[interval])
# acf(r2$X10[interval])


