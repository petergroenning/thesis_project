load_model <- function(model_name){
    load(model_name)
    return(fit)
}

fit <- load_model('models/base_model6.Rdata')

data <- read.csv('data/processed/1m_2_data.csv')

D <- data[19500:26000,]


p1d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 1*24)
write.csv(p1d$state$pred, 'predictions/results_simple/p_1d.csv')
write.csv(p1d$output$sd, 'predictions/results_simple/sd_obs_1d.csv')
write.csv(p1d$state$sd, 'predictions/results_simple/sd_state_1d.csv')

p2d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 2*24)
write.csv(p2d$state$pred, 'predictions/results_simple/p_2d.csv')
write.csv(p2d$output$sd, 'predictions/results_simple/sd_obs_2d.csv')
write.csv(p2d$state$sd, 'predictions/results_simple/sd_state_2d.csv')

p3d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 3*24)
write.csv(p3d$state$pred, 'predictions/results_simple/p_3d.csv')
write.csv(p3d$output$sd, 'predictions/results_simple/sd_obs_3d.csv')
write.csv(p3d$state$sd, 'predictions/results_simple/sd_state_3d.csv')

p4d <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 4*24)
write.csv(p4d$state$pred, 'predictions/results_simple/p_4d.csv')
write.csv(p4d$output$sd, 'predictions/results_simple/sd_obs_4d.csv')
write.csv(p4d$state$sd, 'predictions/results_simple/sd_state_4d.csv')

p1w <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 7*24)
write.csv(p1w$state$pred, 'predictions/results_simple/p_1w.csv')
write.csv(p1w$output$sd, 'predictions/results_simple/sd_obs_1w.csv')
write.csv(p1w$state$sd, 'predictions/results_simple/sd_state_1w.csv')

p2w <- predict(fit, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 14*24)
write.csv(p2w$state$pred, 'predictions/results_simple/p_2w.csv')
write.csv(p2w$output$sd, 'predictions/results_simple/sd_obs_2w.csv')
write.csv(p2w$state$sd, 'predictions/results_simple/sd_state_2w.csv')

