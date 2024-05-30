library(ctsmr)

data <- read.csv('data/processed/data2_2.csv')

data$FtopInk1 <- c(0, data$FtopIn[1:(length(data$FtopIn)-1)])
data$dFtopIn <- c(0,diff(data$FtopIn))

setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)
    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]
    model$ParameterValues$initial <- newparams
    return(model)
}
setInitialState <- function(model, data){
    model$setParameter(x1 = data$X1[1])
    model$setParameter(x2 = data$X2[1])
    model$setParameter(x3 = data$X3[1])
    model$setParameter(x4 = data$X4[1])
    model$setParameter(x5 = data$X5[1])
    model$setParameter(x6 = data$X6[1])
    model$setParameter(x7 = data$X7[1])


    model$setParameter(y1 = data$X1[1])
    model$setParameter(y2 = data$X2[1])
    model$setParameter(y3 = data$X3[1])
    model$setParameter(y4 = data$X4[1])
    model$setParameter(y5 = data$X5[1])
    model$setParameter(y6 = data$X6[1])
    model$setParameter(y7 = data$X7[1])
    return(model)
}

makemodel <- function(type = 'basic', data = NULL){
    model <- ctsm$new()
    if (type != 'model5'){
    model$addObs(X1 ~ x1) # Input Layer (3m)
    model$addObs(X2 ~ x2) # (3m)
    model$addObs(X3 ~ x3) # (3m)
    model$addObs(X4 ~ x4) # Input Layer (3m)
    model$addObs(X5 ~ x5) # (3m)
    model$addObs(X6 ~ x6) # (0.5)
    model$addObs(X7 ~ x7) # Input Layer (0.5)


    model$setVariance(X1 ~ sigma_X2^2)
    model$setVariance(X2 ~ sigma_X2^2)
    model$setVariance(X3 ~ sigma_X2^2)
    model$setVariance(X4 ~ sigma_X2^2)
    model$setVariance(X5 ~ sigma_X2^2)
    model$setVariance(X6 ~ sigma_X1^2)
    model$setVariance(X7 ~ sigma_X1^2)
    } else {
    model$addObs(X6 ~ x6) # Input Layer (1m)
    model$setVariance(X6 ~ sigma_X1^2)
    }

    model$setParameter(V1 = 3185)
    model$setParameter(V2 = 5950)
    model$setParameter(V3 = 9578)
    model$setParameter(V4 = 14072)
    model$setParameter(V5 = 19428)
    model$setParameter(V6 = 7816/2)
    model$setParameter(V7 = 7816/2)



    model$setParameter(Abot = 697)
    model$setParameter(A1 = 1475)
    model$setParameter(A2 = 2540)
    model$setParameter(A3 = 3894)
    model$setParameter(A4 = 5535)
    model$setParameter(A5 = 7465)
    model$setParameter(A6 = 7465/4)
    model$setParameter(Atop = 8172)

    model$setParameter(S1 = 114 + 130 + 146)
    model$setParameter(S2 = 162 + 178 + 194)
    model$setParameter(S3 = 210 + 226 + 242)
    model$setParameter(S4 = 258 + 274 + 290)
    model$setParameter(S5 = 306 + 322 + 338)
    model$setParameter(S6 = 354)



    obs_std1<- 1/sqrt(6)*1.25e-4
    obs_std2<- 1/sqrt(36)*1.25e-4
    model$setParameter(sigma_X1 = obs_std1)
    model$setParameter(sigma_X2 = obs_std2)


    if (type == 'simple'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*exp(vbot)+(x2-x1)*(FbotOut*exp(f))-exp(ubot)*(Abot+S1)*x1)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f))+(x3-x2)*(FbotOut*exp(f))-exp(u)*(S2)*x2)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f))+(x4-x3)*(FbotOut*exp(f))-exp(u)*(S3)*x3)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f))+(x5-x4)*(FbotOut*exp(f))-exp(u)*(S4)*x4+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*exp(f))+(x6-x5)*(FtopIn*exp(f))-exp(u)*(S5)*x5)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*exp(f))+(x7-x6)*(FtopIn*exp(f))-exp(u)*(S6)*x6)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*(FtopOut*exp(f))+(Ttop-x7)*FtopIn*exp(vtop)+exp(utop)*(ambientTemp-x7)*Atop)+dw7*exp(sigma_x))

    } else if (type == 'basic'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*exp(vbot)+(x2-x1)*(FbotOut*exp(f))-exp(ubot)*(Abot+S1)*x1)+dt*(y1-x1)*exp(mu_x)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f))+(x3-x2)*(FbotOut*exp(f))-exp(u)*(S2)*x2)+dt*(y2-x2)*exp(mu_x)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f))+(x4-x3)*(FbotOut*exp(f))-exp(u)*(S3)*x3)+dt*(y3-x3)*exp(mu_x)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f))+(x5-x4)*(FbotOut*exp(f))-exp(u)*(S4)*x4+(Tmid-x4)*FmidIn*exp(vmid))+dt*(y4-x4)*exp(mu_x)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*exp(f))+(x6-x5)*(FtopIn*exp(f))-exp(u)*(S5)*x5)+dt*(y5-x5)*exp(mu_x)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*exp(f))+(x7-x6)*(FtopIn*exp(f))-exp(u)*(S6)*x6)+dt*(y6-x6)*exp(mu_x)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*(FtopOut*exp(f))+(Ttop-x7)*FtopIn*exp(vtop)+exp(utop)*(ambientTemp-x7)*Atop)+dt*(y7-x7)*exp(mu_x)+dw7*exp(sigma_x))

    model$addSystem(dy1~dt*(x1-y1)*exp(mu_y)+dwy1*exp(sigma_y))
    model$addSystem(dy2~dt*(x2-y2)*exp(mu_y)+dwy2*exp(sigma_y))
    model$addSystem(dy3~dt*(x3-y3)*exp(mu_y)+dwy3*exp(sigma_y))
    model$addSystem(dy4~dt*(x4-y4)*exp(mu_y)+dwy4*exp(sigma_y))
    model$addSystem(dy5~dt*(x5-y5)*exp(mu_y)+dwy5*exp(sigma_y))
    model$addSystem(dy6~dt*(x6-y6)*exp(mu_y)+dwy6*exp(sigma_y))
    model$addSystem(dy7~dt*(x7-y7)*exp(mu_y)+dwy7*exp(sigma_y))

    }else if (type == 'basic1'){
    
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f))+(x4-x3)*(FbotOut*exp(f))-exp(u)*(S3)*x3)+dt*(y3-x3)*exp(mu_x)+dw3*exp(sigma_x))
    

    }



    ###############
    model$addInput("FtopInk1",'dFtopIn')
    model$addInput("Ttop", "Tbot", "Tmid", "FtopIn", "FmidIn", "FbotIn", "FtopOut", "FbotOut","FmidOut", "ambientTemp", "Ftop", "Fbot", "Fmid")

    model$setParameter(sigma_x=log(c(init=0.1,lb=1e-12,ub=1)))
    model$setParameter(sigma_y=log(c(init=0.1,lb=1e-12,ub=1)))

    model$setParameter(f=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(f2=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(fbot=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(fmid=log(c(init=1,lb=1e-12,ub=10)))


    model$setParameter(utop=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(ubot=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(u=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(vtop=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(vtop1=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(vmid=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(vbot=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(mu_x=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(mu_y=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(mu_xmid=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(mu_ymid=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(mu_xbot=log(c(init=1,lb=1e-12,ub=10)))
    model$setParameter(mu_ybot=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(a = c(init = 0, lb = -1, ub = 1))
    model$setParameter(mu_x1=log(c(init=1,lb=1e-12,ub=10)))
        model$setParameter(mu_x6=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(mu_y1=log(c(init=1,lb=1e-12,ub=10)))
        model$setParameter(mu_y6=log(c(init=1,lb=1e-12,ub=10)))

    model$setParameter(u1=log(c(init=1e-3,lb=1e-12,ub=10)))

    model$setParameter(xtop = c(init = data$X6[1], lb = 0, ub = 100))




    return (model)
    }   


load_model <- function(model_name){
    load(model_name)
    return(fit)
}

# old_fit <- load_model('models/2state_model1.Rdata')

data <- data[3000:7500,] # 1 hour 
type <- 'basic'
model <- makemodel(type = type, data = data)
# model <- setParams(model, old_fit)
model <- setInitialState(model, data)

model$options$initialVarianceScaling <- 1e-5

fit <- model$estimate(data, firstorder = TRUE)

save(fit, file = paste0('models/augmented',type,'.Rdata'))
p  <- predict(fit, newdata = data, firstorderinputinterpolation=TRUE, n.ahead = 1)
r <- p$output$pred - data[c('X1','X2','X3','X4','X5','X6','X7')]
# # r <- p$output$pred - data[c('X6')]
# # write.csv(r, paste0('residuals/2state_',type,'.csv'))

par(mfrow = c(1,1))
plot(r$X7)
# ccf(r$X6, data$FtopIn)
# plot(r$X6)

# ccf(r$X6, data$ambientTemp)
