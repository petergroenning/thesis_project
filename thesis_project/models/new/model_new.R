library(ctsmr)

data <- read.csv('data/processed/data2.csv')


makemodel <- function(type = 'basic', data = NULL){
    model <- ctsm$new()

    model$addObs(X1 ~ x1) # Input Layer (3m)
    model$addObs(X2 ~ x2) # (3m)
    model$addObs(X3 ~ x3) # (3m)
    model$addObs(X4 ~ x4) # Input Layer (3m)
    model$addObs(X5 ~ x5) # (3m)
    model$addObs(X6 ~ x6) # Input Layer (1m)

    model$setVariance(X1 ~ sigma_X)
    model$setVariance(X2 ~ sigma_X)
    model$setVariance(X3 ~ sigma_X)
    model$setVariance(X4 ~ sigma_X)
    model$setVariance(X5 ~ sigma_X)
    model$setVariance(X6 ~ sigma_X)

    model$setParameter(V1 = 3185)
    model$setParameter(V2 = 5950)
    model$setParameter(V3 = 9578)
    model$setParameter(V4 = 14072)
    model$setParameter(V5 = 19428)
    model$setParameter(V6 = 7816)

    model$setParameter(A1 = 1475)
    model$setParameter(A2 = 2540)
    model$setParameter(A3 = 3894)
    model$setParameter(A4 = 5535)
    model$setParameter(A5 = 7465)

    model$setParameter(Cp = c(init = 1000, lb = 1, ub = 10000))
    model$setParameter(lst = 0.63)







    if (type == 'basic1'){
    model$addSystem(dx1~dt/V1*((x2-x1)*FbotIn*exp(f1)+(x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(x5-x6)*FtopIn*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))
    
    } else if (type == 'basic2'){
    model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))
    
    }else if (type == 'basic3'){
    model$addSystem(dx1~dt/V1*((x2-x1)*(FbotOut*exp(f)+exp(k))+(Tbot-x1)*FbotIn*exp(vbot))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f)+exp(k))+(x3-x2)*(FbotOut*exp(f)+exp(k)))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f)+exp(k))+(x4-x3)*(FbotOut*exp(f)+exp(k)))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f)+exp(k))+(x5-x4)*(FbotOut*exp(f)+exp(k))+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*exp(f)+exp(k))+(x6-x5)*(FtopIn*exp(f)+exp(k)))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*exp(f)+exp(k))+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))

    }else if (type == 'augmented_1'){
    model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*(Y1-x1)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*(Y2-x2)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*(Y3-x3)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*(Y4-x4)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*(Y5-x5)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*(Y6-x6)+dw6*exp(sigma_x))

    model$addSystem(dY1~dt*(x1-Y1)+dwy1*exp(sigma_y))
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*exp(sigma_y))
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*exp(sigma_y))
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*exp(sigma_y))
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*exp(sigma_y))
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*exp(sigma_y))

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))

    }else if (type == 'augmented_2'){
      model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*(Y1-x1)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*(Y2-x2)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*(Y3-x3)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*(Y4-x4)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*(Y5-x5)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*(Y6-x6)+dw6*exp(sigma_x))

    model$addSystem(dY1~dt*(x1-Y1)+dwy1*exp(sigma_y))
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*exp(sigma_y))
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*exp(sigma_y))
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*exp(sigma_y))
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*exp(sigma_y))
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*exp(sigma_y))

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))


    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)
    model$setParameter(vbot = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vtop = 0)

    }else if (type == 'augmented_3'){
      model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*(Y1-x1)*exp(a)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*(Y2-x2)*exp(a)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*(Y3-x3)*exp(a)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*(Y4-x4)*exp(a)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*(Y5-x5)*exp(a)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*(Y6-x6)*exp(a)+dw6*exp(sigma_x))

    model$addSystem(dY1~dt*(x1-Y1)*exp(b)+dwy1*exp(sigma_y))
    model$addSystem(dY2~dt*(x2-Y2)*exp(b)+dwy2*exp(sigma_y))
    model$addSystem(dY3~dt*(x3-Y3)*exp(b)+dwy3*exp(sigma_y))
    model$addSystem(dY4~dt*(x4-Y4)*exp(b)+dwy4*exp(sigma_y))
    model$addSystem(dY5~dt*(x5-Y5)*exp(b)+dwy5*exp(sigma_y))
    model$addSystem(dY6~dt*(x6-Y6)*exp(b)+dwy6*exp(sigma_y))

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))
    model$setParameter(a = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(b = log(c(init=1, lower=1e-12, upper=10)))

    }else if (type == 'augmented_4'){
      model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*(Y1-x1)*exp(a)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*(Y2-x2)*exp(a)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*(Y3-x3)*exp(a)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*(Y4-x4)*exp(a)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*(Y5-x5)*exp(a)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*(Y6-x6)*exp(a)+dw6*exp(sigma_x))

    model$addSystem(dY1~dt*((x1-Y1)*exp(b)+(10-Y1)*exp(c))+dwy1*exp(sigma_y))
    model$addSystem(dY2~dt*((x2-Y2)*exp(b)+(10-Y2)*exp(c))+dwy2*exp(sigma_y))
    model$addSystem(dY3~dt*((x3-Y3)*exp(b)+(10-Y3)*exp(c))+dwy3*exp(sigma_y))
    model$addSystem(dY4~dt*((x4-Y4)*exp(b)+(10-Y4)*exp(c))+dwy4*exp(sigma_y))
    model$addSystem(dY5~dt*((x5-Y5)*exp(b)+(10-Y5)*exp(c))+dwy5*exp(sigma_y))
    model$addSystem(dY6~dt*((x6-Y6)*exp(b)+(10-Y6)*exp(c))+dwy6*exp(sigma_y))

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))
    model$setParameter(a = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(b = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(c = log(c(init=1, lower=1e-12, upper=10)))
    }else if (type == 'augmented_5'){
      model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*Y1+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*Y2+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*Y3+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*Y4+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*Y5+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*Y6+dw6*exp(sigma_x))

    model$addSystem(dY1~dwy1*exp(sigma_y))
    model$addSystem(dY2~dwy2*exp(sigma_y))
    model$addSystem(dY3~dwy3*exp(sigma_y))
    model$addSystem(dY4~dwy4*exp(sigma_y))
    model$addSystem(dY5~dwy5*exp(sigma_y))
    model$addSystem(dY6~dwy6*exp(sigma_y))

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)


    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))
    model$setParameter(a = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(b = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(c = log(c(init=1, lower=1e-12, upper=10)))
    }else if (type == 'augmented_6'){
    model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*(Y1-x1)*exp(a)+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2)+(2)*exp(v))+dt*(Y2-x2)*exp(a)+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3)+(2)*exp(v))+dt*(Y3-x3)*exp(a)+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(2)*exp(v)+(Tmid-x4)*FmidIn*exp(vmid))+dt*(Y4-x4)*exp(a)+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5)+(2)*exp(v))+dt*(Y5-x5)*exp(a)+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*(Y6-x6)*exp(a)+dw6*exp(sigma_x))

    model$addSystem(dY1~dt*(x1-Y1)*exp(b)+dwy1*exp(sigma_y))
    model$addSystem(dY2~dt*(x2-Y2)*exp(b)+dwy2*exp(sigma_y))
    model$addSystem(dY3~dt*(x3-Y3)*exp(b)+dwy3*exp(sigma_y))
    model$addSystem(dY4~dt*(x4-Y4)*exp(b)+dwy4*exp(sigma_y))
    model$addSystem(dY5~dt*(x5-Y5)*exp(b)+dwy5*exp(sigma_y))
    model$addSystem(dY6~dt*(x6-Y6)*exp(b)+dwy6*exp(sigma_y))

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))
    model$setParameter(a = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(b = log(c(init=1, lower=1e-12, upper=10)))

    model$setParameter(v = log(c(init=1, lower=1e-12, upper=10)))

    }else if (type == 'basic4'){
      model$addSystem(dx1~dt/V1*((x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2)+(Tbot-x2)*FbotIn/V1*exp(vbot2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))

    # model$addSystem(dY1~dt*(x1-Y1)*exp(b)+dwy1*exp(sigma_y))
    # model$addSystem(dY2~dt*(x2-Y2)*exp(b)+dwy2*exp(sigma_y))
    # model$addSystem(dY3~dt*(x3-Y3)*exp(b)+dwy3*exp(sigma_y))
    # model$addSystem(dY4~dt*(x4-Y4)*exp(b)+dwy4*exp(sigma_y))
    # model$addSystem(dY5~dt*(x5-Y5)*exp(b)+dwy5*exp(sigma_y))
    # model$addSystem(dY6~dt*(x6-Y6)*exp(b)+dwy6*exp(sigma_y))

    # model$setParameter(Y1 = data$X1[1])
    # model$setParameter(Y2 = data$X2[1])
    # model$setParameter(Y3 = data$X3[1])
    # model$setParameter(Y4 = data$X4[1])
    # model$setParameter(Y5 = data$X5[1])
    # model$setParameter(Y6 = data$X6[1])

    model$setParameter(sigma_y = log(c(init=0.1, lower=1e-12, upper=1)))
    model$setParameter(a = log(c(init=1, lower=1e-12, upper=10)))
    model$setParameter(b = log(c(init=1, lower=1e-12, upper=10)))

    model$setParameter(vbot2= c(init = 0, lower = -10, upper = 10))
    model$setParameter(v2= c(init = 0, lower = -10, upper = 10))
    model$setParameter(v3= c(init = 0, lower = -10, upper = 10))
    model$setParameter(v4= c(init = 0, lower = -10, upper = 10))

    }
    
    
    # else if (type == 'basic4'){
    # model$addSystem(dx1~dt/V1*(((x2-x1)*(FbotIn+(x2-x1)+(FbotOut+exp(k1)))/exp(f1)+(Tbot-x1)*FbotIn*exp(vbot)))+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn+exp(k2))+(x3-x2)*(FbotOut+exp(k2)))/exp(f2)+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn+exp(k3))+(x4-x3)*(FbotOut+exp(k3)))/exp(f3)+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*(((x3-x4)*(FbotIn+exp(k4))+(x5-x4)*(FbotOut+exp(k4)))/exp(f5)+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut+exp(k5))+(x6-x5)*(FtopIn+exp(k5)))/exp(f5)+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*(((x5-x6)*(FtopOut+exp(k6))+(x5-x6)*FtopIn)/exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))
    
    # }else if (type == 'model'){
    # model$addSystem(dx1~dt/V1*((x2-x1)*FbotIn*exp(f1)+(x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*Y1+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*Y2+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*Y3+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*Y4+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*Y5+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(x5-x6)*FtopIn*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*Y6+dw6*exp(sigma_x))

    # model$addSystem(dY1~dwy1*exp(sigma_y))
    # model$addSystem(dY2~dwy2*exp(sigma_y))
    # model$addSystem(dY3~dwy3*exp(sigma_y))
    # model$addSystem(dY4~dwy4*exp(sigma_y))
    # model$addSystem(dY5~dwy5*exp(sigma_y))
    # model$addSystem(dY6~dwy6*exp(sigma_y))

    # model$setParameter(Y1 = 0)
    # model$setParameter(Y2 = 0)
    # model$setParameter(Y3 = 0)
    # model$setParameter(Y4 = 0)
    # model$setParameter(Y5 = 0)
    # model$setParameter(Y6 = 0)

    # model$setParameter(f1 = 0)
    # model$setParameter(f2 = 0)
    # model$setParameter(f3 = 0)
    # model$setParameter(f4 = 0)
    # model$setParameter(f5 = 0)
    # model$setParameter(f6 = 0)
    # model$setParameter(vbot = 0)
    # model$setParameter(vmid = 0)
    # model$setParameter(vtop = 0)


    # }  else if (type == 'model2'){
    # model$addSystem(dx1~dt/V1*((x2-x1)*FbotIn*exp(f1)+(x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dt*Y1+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*Y2+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*Y3+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dt*Y4+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))+dt*Y5+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(x5-x6)*FtopIn*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dt*Y6+dw6*exp(sigma_x))
    

    # model$addSystem(dY1~dt*(FbotIn-a*Y1)+dwY1*exp(sigma_y))
    # model$addSystem(dY2~dt*(FbotOut-a*Y2)+dwY2*exp(sigma_y))
    # model$addSystem(dY3~dt*(FtopIn-a*Y3)+dwY3*exp(sigma_y))
    # model$addSystem(dY4~dt*(FtopOut-a*Y4)+dwY4*exp(sigma_y))

    # model$setParameter(Y1 = 0)
    # model$setParameter(Y2 = 0)
    # model$setParameter(Y3 = 0)
    # model$setParameter(Y4 = 0)

    # model$setParameter(a = c(init = 1, lb = 0, ub = 10))


    # }  else if (type == 'model3'){
    # model$addSystem(dx1~dt/V1*((x2-x1)*Y1*exp(f1)+(x2-x1)*Y2*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*Y1*exp(f2)+(x3-x2)*Y2*exp(f2))+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*Y1*exp(f3)+(x4-x3)*Y2*exp(f3))+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*((x3-x4)*Y1*exp(f4)+(x5-x4)*Y2*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*Y4*exp(f5)+(x6-x5)*Y3*exp(f5))+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*((x5-x6)*Y4*exp(f6)+(x5-x6)*Y3*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))+dw6*exp(sigma_x))
    

    # model$addSystem(dY1~dt*(FbotIn+a*Y1)+dwY1*exp(sigma_y))
    # model$addSystem(dY2~dt*(FbotOut+a*Y2)+dwY2*exp(sigma_y))
    # model$addSystem(dY3~dt*(FtopIn+a*Y3)+dwY3*exp(sigma_y))
    # model$addSystem(dY4~dt*(FtopOut+a*Y4)+dwY4*exp(sigma_y))

    # model$setParameter(Y1 = 0)
    # model$setParameter(Y2 = 0)
    # model$setParameter(Y3 = 0)
    # model$setParameter(Y4 = 0)

    # model$setParameter(a = c(init = -1, lb = -10, ub = 0))
    # } else if (type == 'model4'){
    # model$addSystem(dx1~dt/V1*((x2-x1)*FbotIn*exp(f1)+(x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))*Y1+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))*Y2+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))*Y3+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))*Y4+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))*Y5+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(x5-x6)*FtopIn*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))*Y6+dw6*exp(sigma_x))

    # model$addSystem(dY1~dwy1*exp(sigma_y))
    # model$addSystem(dY2~dwy2*exp(sigma_y))
    # model$addSystem(dY3~dwy3*exp(sigma_y))
    # model$addSystem(dY4~dwy4*exp(sigma_y))
    # model$addSystem(dY5~dwy5*exp(sigma_y))
    # model$addSystem(dY6~dwy6*exp(sigma_y))

    # model$setParameter(Y1 = 0)
    # model$setParameter(Y2 = 0)
    # model$setParameter(Y3 = 0)
    # model$setParameter(Y4 = 0)
    # model$setParameter(Y5 = 0)
    # model$setParameter(Y6 = 0)

    # model$setParameter(f1 = 0)
    # model$setParameter(f2 = 0)
    # model$setParameter(f3 = 0)
    # model$setParameter(f4 = 0)
    # model$setParameter(f5 = 0)
    # model$setParameter(f6 = 0)
    # model$setParameter(vbot = 0)
    # model$setParameter(vmid = 0)
    # model$setParameter(vtop = 0)
    # }else if (type == 'model5'){
    # model$addSystem(dx1~dt/V1*((x2-x1)*FbotIn*exp(f1)+(x2-x1)*FbotOut*exp(f1)+(Tbot-x1)*FbotIn*exp(vbot))*Y1+dw1*exp(sigma_x))
    # model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))*Y2+dw2*exp(sigma_x))
    # model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))*Y3+dw3*exp(sigma_x))
    # model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4)+(Tmid-x4)*FmidIn*exp(vmid))*Y4+dw4*exp(sigma_x))
    # model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*exp(f5)+(x6-x5)*FtopIn*exp(f5))*Y5+dw5*exp(sigma_x))
    # model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*exp(f6)+(x5-x6)*FtopIn*exp(f6)+(Ttop-x6)*FtopIn*exp(vtop))*Y6+dw6*exp(sigma_x))

    # model$addSystem(dY1~dwy1*exp(sigma_y))
    # model$addSystem(dY2~dwy2*exp(sigma_y))
    # model$addSystem(dY3~dwy3*exp(sigma_y))
    # model$addSystem(dY4~dwy4*exp(sigma_y))
    # model$addSystem(dY5~dwy5*exp(sigma_y))
    # model$addSystem(dY6~dwy6*exp(sigma_y))

    # model$setParameter(Y1 = 0)
    # model$setParameter(Y2 = 0)
    # model$setParameter(Y3 = 0)
    # model$setParameter(Y4 = 0)
    # model$setParameter(Y5 = 0)
    # model$setParameter(Y6 = 0)

    # model$setParameter(f1 = 0)
    # model$setParameter(f2 = 0)
    # model$setParameter(f3 = 0)
    # model$setParameter(f4 = 0)
    # model$setParameter(f5 = 0)
    # model$setParameter(f6 = 0)
    # model$setParameter(vbot = 0)
    # model$setParameter(vmid = 0)
    # model$setParameter(vtop = 0)
    # }

    model$setParameter(sigma_X = 2.083e-5)
    model$setParameter(sigma_x = log(c(init=0.1, lower=1e-12, upper=1)))
    # model$setParameter(sigma_y = log(c(init=0.001, lower=1e-12, upper=0.5)))
    model$setParameter(f1 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f2 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f3 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f4 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f5 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f6 = log(c(init = 1, lb = 1e-12, ub = 10)))
    # model$setParameter(f = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k1 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k2 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k3 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k4 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k5 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k6 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k7 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k8 = log(c(init = 1, lb = 1e-12, ub = 10)))


    model$setParameter(k = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vbot = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vmid = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vtop = log(c(init = 1, lb = 1e-12, ub = 10)))



    
   



    model$addInput('Ttop','Tmid','Tbot','FtopIn','FtopOut','FmidIn','FmidOut','FbotIn','FbotOut','Ftop','Fmid','Fbot', 'ambientTemp')

    return(model)
}


setInitialState <- function(model, data){
    model$setParameter(x1 = data$X1[1])
    model$setParameter(x2 = data$X2[1])
    model$setParameter(x3 = data$X3[1])
    model$setParameter(x4 = data$X4[1])
    model$setParameter(x5 = data$X5[1])
    model$setParameter(x6 = data$X6[1])

    return(model)
}

setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)

    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]

    model$ParameterValues$initial <- newparams
    return(model)
}

load_fit <- function(path){
    load(path)
    return(fit)
}



basic2<- load_fit('models/ctsm2/basic3.RData')


D <- data[3000:7500,]

type <- 'basic4'
model <- makemodel(type, data = D)
model <- setInitialState(model, D)
model <- setParams(model, basic2)

fit <- model$estimate(data = D, firstorder = TRUE)
model_dir <- 'models/ctsm2/'
dir.create(model_dir, showWarnings = FALSE)

save(fit, file = paste0(model_dir, type ,'.RData'))



