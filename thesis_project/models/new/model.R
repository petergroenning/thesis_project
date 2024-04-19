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

    model$setVariance(X1 ~ sigma_X^2)
    model$setVariance(X2 ~ sigma_X^2)
    model$setVariance(X3 ~ sigma_X^2)
    model$setVariance(X4 ~ sigma_X^2)
    model$setVariance(X5 ~ sigma_X^2)
    model$setVariance(X6 ~ sigma_X^2)

    model$setParameter(V1 = 3185)
    model$setParameter(V2 = 5950)
    model$setParameter(V3 = 9578)
    model$setParameter(V4 = 14072)
    model$setParameter(V5 = 19428)
    model$setParameter(V6 = 7816)


    if (type == 'basic'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)

    } else if (type == 'gradient'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(x2-x1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(x3-x1)/2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(x4-x2)/2)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x3)/2)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(x6-x4)/2)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*(x6-x4))+dw6*sigma_x)
    
    } else if (type == 'random'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*Y6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*(0)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*(0)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*(0)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*(0)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*(0)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*(0)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    } else if (type == 'lowpass'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)+dt*(Y1-x1)*a+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)+dt*(Y2-x2)*a+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)+dt*(Y3-x3)*a+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)+dt*(Y4-x4)*a+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)+dt*(Y5-x5)*a+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u)+dt*(Y6-x6)*a+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*(x1-Y1)*b+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*(x2-Y2)*b+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*(x3-Y3)*b+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*(x4-Y4)*b+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*(x5-Y5)*b+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*(x6-Y6)*b+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    } else if (type == 'gradient_state'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*Y6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut-Y1*FbotIn)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut+(Y1-Y2)*FbotIn)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut+(Y2-Y3)*FbotIn)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut+(Y3-Y4)*FbotIn-Y4*FmidIn)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn+(Y4-Y5)*FtopOut)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut-Y5*FtopIn)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'gradient_state2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*Y6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut*k1-Y1*FbotIn*k1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut*k2+(Y1-Y2)*FbotIn*k2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut*k3+(Y2-Y3)*FbotIn*k3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut*k4+(Y3-Y4)*FbotIn-Y4*FmidIn)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn*k5+(Y4-Y5)*FtopOut*k5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut*k6-Y5*FtopIn)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))



    }else if (type == 'g'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)*Y1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)*Y2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)*Y3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)*Y4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)*Y5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)*Y6+(Ttop-x6)*FtopIn*Y6*vtop+(ambientTemp-x6)*u)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut*k-Y1*FbotIn*k)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut*k+(Y1-Y2)*FbotIn*k)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut*k+(Y2-Y3)*FbotIn*k)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut*k+(Y3-Y4)*FbotIn-Y4*FmidIn)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn*k+(Y4-Y5)*FtopOut*Y5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut*k-Y5*FtopIn)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))
    
    }else if (type == 'g2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u+Ftop*Y6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut*k-Y1*FbotIn*k)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut*k+(Y1-Y2)*FbotIn*k)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut*k+(Y2-Y3)*FbotIn*k)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut*k+(Y3-Y4)*FbotIn-Y4*FmidIn)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn*k+(Y4-Y5)*FtopOut*Y5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut*k-Y5*FtopIn)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))
    
    }else if (type == 'glowpass'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1)*a)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2)*a)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3)*a)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*(Y4-x4)*a)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5)*a)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)*a)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'glowpass2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)*Y1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)*Y2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)*Y3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)*Y4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)*Y5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop)*Y6+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'model'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1)*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2)*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3)*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4-Fbot*(Y4-x4)*a4+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5)*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)*a6+(ambientTemp-x6)*u)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b4)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b6)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'model2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2))+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3))+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4-Fbot*(Y4-x4)+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5))+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)+(ambientTemp-x6)*u)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'model3'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1)*a)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2)*a)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3)*a)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4-Fbot*(Y4-x4)*a+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5)*a)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)*a+(ambientTemp-x6)*u)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b)+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    }else if (type == 'model4'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)*(Y1-x1)*a+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)*(Y2-x2)*a+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)*(Y3-x3)*a+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)*(Y4-x4)*a+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)*(Y5-x5)*a+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u)*(Y6-x6)*a+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b)+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    }else if (type == 'model5'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)*Y1*a1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)*Y2*a2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)*Y3*a3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)*Y4*a4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)*Y5*a5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u)*Y6*a6+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*b1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*b2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*b3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*b4)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*b5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*b6)+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    }
    


    model$setParameter(sigma_X = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_y = c(init=0.1, lower=0, upper=0.5))

    model$setParameter(f1 = c(init=0.1, lower=0, upper=10))
    model$setParameter(f2 = c(init=0.1, lower=0, upper=10))
    model$setParameter(f3 = c(init=0.1, lower=0, upper=10))
    model$setParameter(f4 = c(init=0.1, lower=0, upper=10))
    model$setParameter(f5 = c(init=0.1, lower=0, upper=10))
    model$setParameter(f6 = c(init=0.1, lower=0, upper=10))

    model$setParameter(k1 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k2 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k3 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k4 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k5 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k6 = c(init=0.1, lower=0, upper=10))
    model$setParameter(k = c(init=0.1, lower=0, upper=10))

    model$setParameter(vtop = c(init=0.1, lower=0, upper=10))
    model$setParameter(vmid = c(init=0.1, lower=0, upper=10))
    model$setParameter(vbot = c(init=0.1, lower=0, upper=10))

    model$setParameter(a = c(init=0.1, lower=-10, upper=10))
    model$setParameter(b = c(init=0.1, lower=0, upper=10))
    model$setParameter(u=c(init=0.1, lower=0, upper=10))

    model$setParameter(a1 = c(init=1, lower=0, upper=10))
    model$setParameter(a2 = c(init=1, lower=0, upper=10))
    model$setParameter(a3 = c(init=1, lower=0, upper=10))
    model$setParameter(a4 = c(init=1, lower=0, upper=10))
    model$setParameter(a5 = c(init=1, lower=0, upper=10))
    model$setParameter(a6 = c(init=1, lower=0, upper=10))

    model$setParameter(b1 = c(init=1, lower=0, upper=10))
    model$setParameter(b2 = c(init=1, lower=0, upper=10))
    model$setParameter(b3 = c(init=1, lower=0, upper=10))
    model$setParameter(b4 = c(init=1, lower=0, upper=10))
    model$setParameter(b5 = c(init=1, lower=0, upper=10))
    model$setParameter(b6 = c(init=1, lower=0, upper=10))




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


D <- data[2000:7000,]

type <- 'lowpass'
model <- makemodel(type, data <- D)
model <- setInitialState(model, D)


fit <- model$estimate(data = D, firstorder = TRUE)
model_dir <- 'models/ctsm/'
dir.create(model_dir, showWarnings = FALSE)

save(fit, file = paste0(model_dir, type ,'.RData'))



