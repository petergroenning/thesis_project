


load_model <- function(model_name){
    load(model_name)
    return(fit)
}

transformp <- function(y){
    yp <- c(rep(y[1],6), rep(y[2],6), rep(y[3],6), rep(y[4],6), rep(y[5],6), rep(y[6],2))
    return(yp)
}
data <- read.csv('data/processed/data3.csv')
D <- read.csv('data/processed/data2.csv')
D$dX1 <- c(0,diff(D$X1))
D$dX2 <- c(0,diff(D$X2))
D$dX3 <- c(0,diff(D$X3))
D$dX4 <- c(0,diff(D$X4))
D$dX5 <- c(0,diff(D$X5))
D$dX6 <- c(0,diff(D$X6))
D$FbotInk1 <- c(0, D$FbotIn[1:(length(D$FbotIn)-1)])
D$FbotOutk1 <- c(0, D$FbotOut[1:(length(D$FbotOut)-1)])
D$FmidInk1 <- c(0, D$FmidIn[1:(length(D$FmidIn)-1)])
D$FtopInk1 <- c(0, D$FtopIn[1:(length(D$FtopIn)-1)])
D$FtopOutk1 <- c(0, D$FtopOut[1:(length(D$FtopOut)-1)])
D$Ttopk1 <- c(0, D$Ttop[1:(length(D$Ttop)-1)])
D$Tmidk1 <- c(0, D$Tmid[1:(length(D$Tmid)-1)])
D$Tbotk1 <- c(0, D$Tbot[1:(length(D$Tbot)-1)])

D <- D[3000:7000,]

fit1 <- load_model('models/ctsm2/basic2.RData')
fit2 <- load_model('models/ctsm2/augmented_3.RData')
p1 <- predict(fit1, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 24*2*7)
p2 <- predict(fit2, newdata = D, firstorderinputinterpolation=TRUE, n.ahead = 24*2*7)

# p1 <- simulate(fit1, newdata = D, firstorderinputinterpolation=TRUE)
# p2 <- simulate(fit2, newdata = D, firstorderinputinterpolation=TRUE)
d <- t(data[3000:7000,2:33])

p1$state$sim

pred1 <- t(p1$output$pred)
# predy1 <- t(p2$state$sim[,7:12])

pred2 <- t(p2$output$pred)


# g <- pred1-predy1

x <- seq(0.5,16,0.5)
xp <- c(seq(1,16,3))

idx <- 3000
par(mfrow=c(1,1))
plot(x,d[,idx], ylim = c(10,90))
points(xp, D[idx,c('X1','X2','X3','X4','X5','X6')], col = 'green', pch=12, bg = 12)

# lines(x,transformp(pred1[,idx]), col = 'blue', lty =2)
lines(xp, pred1[,idx], type = 'o', col = 'blue')

# lines(x,transformp(pred2[,idx]), col = 'red', lty =2)
lines(xp, pred2[,idx], type = 'o', col = 'red')
legend('topleft', c('G1', 'G2'), col = c('blue', 'red'), lty = c(1,1))
# g[,idx]

fit2$loglik


res1 <- t(pred1) - D[,2:7]
res2 <- t(pred2) - D[,2:7]

mean(t(res1^2))
mean(t(res2^2))
plot(x,approx(xp, pred1[,idx], x, rule = 2)$y)
lines(x,d[,idx])


plot(res2$X1)
