library(tidyverse)
library(ggplot2)
source('thesis_project/R_scripts/makedata.R')


load_model<-function(path){
    load(path)
    return(fit)
}

fit1<-load_model('models/CTSM/model3.RData')
fit2 <-load_model('models/CTSM/model8.RData')
fit3 <- load_model('models/CTSM/model9.RData')

summary(fit3)
fit3$model
fit3$loglik
summary(fit2)
states<-c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')
data <- make_data(start = '2015-04-01', end = '2015-05-01', splines = 7)
obs<-data[,states]


n.ahead<-24
initial <- as.numeric(obs[1,])
names(initial) <- names(fit1$xm[1:16])

p1<-predict(fit1, newdata = data, firstorderinputinterpolation=TRUE, n.ahead=n.ahead, x0 = initial)$output$pred
p2<-predict(fit2, newdata = data,  firstorderinputinterpolation=TRUE, n.ahead=n.ahead, x0 = initial)$output$pred
p3<-predict(fit3, newdata = data,  firstorderinputinterpolation=TRUE, n.ahead=n.ahead, x0 = initial)$output$pred

p1<-as.data.frame(p1)
p2<-as.data.frame(p2)
p3<-as.data.frame(p3)

fit3$model

# res1<-p1-obs
# res2<-p2-obs
# res3<-p3-obs

# res1<-res1 %>% gather(key='state',value='residual')
# res2<-res2 %>% gather(key='state',value='residual')
# res3<-res3 %>% gather(key='state',value='residual')

# #renmaestatetonumber
# res1$state<-as.factor(gsub('X','',res1$state))
# res2$state<-as.factor(gsub('X','',res2$state))
# res3$state<-as.factor(gsub('X','',res3$state))

# res1$model<-'model3'
# res2$model<-'model6'
# res3$model<-'model7'

# res<-rbind(res1, res2)
# res$state<-factor(res$state,levels=0:15)

# ggplot(res,aes(x=state,y=residual,fill=model))+
#     geom_boxplot()+
#     theme_minimal()+
#     theme(axis.text.x = element_text(angle = 90,size=20), axis.text.y = element_text(size=20))+
#     labs(title='Residuals for each state',x='State',y='Residual',fontsize=200)
fit2$model
fit2$loglik
summary(fit2)
interval <- 100:500


ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X0[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X0[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X0[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X0',x='Time',y='State', color = "Legend")+
    geom_line(aes(y=X0),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X1[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X1[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X1[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X1',x='Time',y='State')+
    geom_line(aes(y=X1),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X2[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X2[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X2[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X2',x='Time',y='State')+
    geom_line(aes(y=X2),color='red')+
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X3[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X3[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X3[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X3',x='Time',y='State')+
    geom_line(aes(y=X3),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X4[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X4[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X4[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X4',x='Time',y='State')+
    geom_line(aes(y=X4),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X5[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X5[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X5[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X5',x='Time',y='State')+
    geom_line(aes(y=X5),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X6[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X6[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X6[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X6',x='Time',y='State')+
    geom_line(aes(y=X6),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X7[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X7[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X7[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X7',x='Time',y='State')+
    geom_line(aes(y=X7),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X8[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X8[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X8[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X8',x='Time',y='State')+
    geom_line(aes(y=X8),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X9[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X9[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X9[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X9',x='Time',y='State')+
    geom_line(aes(y=X9),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X10[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X10[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X10[interval]),color='green')+
    labs(title='Comparison of model predictions X10',x='Time',y='State')+
    geom_line(aes(y=X10),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X11[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X11[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X11[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X11',x='Time',y='State')+
    geom_line(aes(y=X11),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X12[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X12[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X12[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X12',x='Time',y='State')+
    geom_line(aes(y=X12),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X13[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X13[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X13[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X13',x='Time',y='State')+
    geom_line(aes(y=X13),color='red')+    
    theme_minimal()

ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X14[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X14[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X14[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X14',x='Time',y='State')+
    geom_line(aes(y=X14),color='red')+    
    theme_minimal()


ggplot(data[interval,],aes(x=t))+
    geom_line(aes(y=p1$X15[interval]), color = 'blue', size = 1)+
    geom_line(aes(y=p2$X15[interval]),color='black', size = 1)+
    geom_line(aes(y=p3$X15[interval]),color='green', size = 1)+
    labs(title='Comparison of model predictions X15',x='Time',y='State')+
    geom_line(aes(y=X15),color='red')+    
    theme_minimal()

# fit3 <- load_model('models/CTSM/model3.RData')

states<-c('X0','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15')
data <- make_data(start = '2015-04-01', end = '2016-05-01', splines = 7)
obs<-data[,states]


n.ahead<-24
initial <- as.numeric(obs[1,])
names(initial) <- names(fit1$xm[1:16])

sim <- simulate(fit3, newdata = data, firstorderinputinterpolation=TRUE, x0 = initial)
p3 <- predict(fit3, newdata = data, firstorderinputinterpolation=TRUE, n.ahead=n.ahead, x0 = initial)$output
p3p <- as.data.frame(p3$pred)
p3s <- as.data.frame(p3$sd)

write.csv(p3p, 'p3.csv')
write.csv(p3s, 'p3s.csv')
write.csv(as.data.frame(sim$output$sim), 'sim.csv')
write.csv(data, 'data.csv')


interval <- 300:500
ggplot(data = p3[interval,], aes(x=interval))+
    geom_line(aes(y = X0), size = 1, color = 'green') + 
    geom_line(aes(y = X1), size = 1, color = '#ff00f2') +
    geom_line(aes(y = X2), size = 1, color = '#ff4800') +
    geom_line(aes(y = X3), size = 1, color = '#ff8c00') +
    geom_line(aes(y = data[interval, 'X1']), size = 1, color = '#ff00f2', lty = 2) +
    geom_line(aes(y = data[interval, 'X0']), size = 1, color = 'green', lty = 2) +
    geom_line(aes(y = data[interval, 'X2']), size = 1, color = '#ff4800', lty = 2) +
    geom_line(aes(y = data[interval, 'X3']), size = 1, color = '#ff8c00', lty = 2) 
    # geom_line(aes(y = 20+1/100 * (data[interval, 'FbotOut']*((X1-X0)*9/10 + (X2-X1)*1/10) + data[interval,'FbotIn']*0.15*(X0-X1))), size = 1, color = 'blue')+
    # geom_line(aes(y = 20+1/100 * (data[interval, 'FbotOut']*((X1-X0)*1/2 + (X2-X1)*1/2) + data[interval,'FbotIn']*0.15*(X0-X1))), size = 1, color = '#0084ff')
    # geom_line(aes(y = 20+1/100 * data[interval, 'FbotVol']*((X2-X1)*1/2 + (X3-X2)*1/2) ))

    # geom_line(aes(y = data[interval, 'FbotOut'] / 100))
    # geom_line(aes(y =20+1/100*(data[interval, 'FbotOut']*(X1-X0) + (data[interval,'Tbot']-X0) * data[interval, 'FbotIn'] * 1.5e-1)), size = 1, color = 'blue')
summary(fit3)
fit3$loglik



