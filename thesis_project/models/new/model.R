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
    
    } else if (type == 'basic2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)
    
    }else if (type == 'basic22'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)
    
    }else if (type == 'basic3'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)
    
    
    } else if (type == 'basic4'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u)+dw6*sigma_x)
    
    
    }  else if (type == 'basic5'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*(vbot+a1*(x2-x1))*FbotIn+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3)-Fbot*(x4-x3)*(x3-x2)*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5)+Ftop*(x6-x5)*(x5-x4)*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*(vtop+(x6-x5)*a6)*FtopIn)+dw6*sigma_x)
    
    
    } else if (type == 'basic6'){
    
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn-Fbot*(x5-x4)*(x4-x3)*a3)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn+(ambientTemp-x6)*u)+dw6*sigma_x)
    
    
    } else if (type == 'basic3_2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1))+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2)+(x3-x2)*(FbotOut*f2)-Fbot*(x3-x2)*(x2-x1)*a1)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3)-Fbot*(x4-x3)*(x3-x2)*a2)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5)+Ftop*(x6-x5)*(x5-x4)*a4)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)
    
    }else if (type == 'b7'){
    
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2)+(x3-x2)*(FbotOut*f2)+FbotVol*(x3-2*x2+x1)*a1)+dw2*(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3)+FbotVol*(x4-2*x3+x2)*a2)+dw3*(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4)+(Tmid-x4)*FmidIn*vmid+FbotVol*(x5-2*x4+x3)*a3)+dw4*(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5)+FtopVol*(x6-2*x5+x4)*a4)+dw5*(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*(sigma_x))
    
    
    }  else if (type == 'b8'){
    
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*(sigma_x1))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2)+(x3-x2)*(FbotOut*f2+k2)+FbotVol*(x3-2*x2+x1)*a1)+dw2*(sigma_x2))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3+k3)+FbotVol*(x4-2*x3+x2)*a2)+dw3*(sigma_x3))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid+FbotVol*(x5-2*x4+x3)*a3)+dw4*(sigma_x4))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5+k5)+FtopVol*(x6-2*x5+x4)*a4)+dw5*(sigma_x5))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*(sigma_x6))
    
    
    }else if (type == 'b1_lowpass'){
    # Basic 1 with low pass filter
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1+(Y1-x1)*a)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2+(Y2-x2)*a)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3+(Y3-x3)*a)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid+(Y4-x4)*a)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+(Y5-x5)*a)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+(Y6-x6)*a)+dw6*sigma_x)

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

    }    else if (type == 'b3_lowpass'){
    # Basic 1 with low pass filter
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1+(Y1-x1)*a)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2+(Y2-x2)*a-b1*Fbot*(x3-x2)*(x2-x1))+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3+(Y3-x3)*a-b1*Fbot*(x4-x3)*(x3-x2))+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid+(Y4-x4)*a-b2*Fbot*(x5-x4)*(x4-x3))+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+(Y5-x5)*a+b2*Ftop*(x6-x5)*(x5-x4))+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+(Y6-x6)*a)+dw6*sigma_x)

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

    }  else if (type == 'b1_lowpass2'){
    # Basic 1 with low pass filter
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)+dt*(Y1-x1)*a1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)+dt*(Y2-x2)*a2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)+dt*(Y3-x3)*a3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)+dt*(Y4-x4)*a4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)+dt*(Y5-x5)*a5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop)+dt*(Y6-x6)*a6+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*(x1-Y1)*b1+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*(x2-Y2)*b2+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*(x3-Y3)*b3+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*(x4-Y4)*b4+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*(x5-Y5)*b5+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*(x6-Y6)*b6+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])
    }else if (type == 'gradient_state'){
    # Basic model with gradient as hidden state
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
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn+(x3-x2)*FbotOut-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn+(x4-x3)*FbotOut-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn+(x5-x4)*FbotOut+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut+(x6-x5)*FtopIn+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut+(Ttop-x6)*FtopIn*vtop+Ftop*Y6)+dw6*sigma_x)

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



    }else if (type == 'gradient_state3'){
    # Basic model with gradient as hidden state
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*Y6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut*k1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut*k2+(Y1-Y2)*FbotIn*k2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut*k3+(Y2-Y3)*FbotIn*k3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut*k4+(Y3-Y4)*FbotIn*k4)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn*k5+(Y4-Y5)*FtopOut*k5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))

    }else if (type == 'gradient_state4'){
    # Basic model with gradient as hidden state
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*Y1*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*Y2*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*Y3*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*Y4*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*Y5*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*Y6*a6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Y2-Y1)*FbotOut*k1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((Y3-Y2)*FbotOut*k2+(Y1-Y2)*FbotIn*k2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((Y4-Y3)*FbotOut*k3+(Y2-Y3)*FbotIn*k3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((Y5-Y4)*FbotOut*k4+(Y3-Y4)*FbotIn*k4)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((Y6-Y5)*FtopIn*k5+(Y4-Y5)*FtopOut*k5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((Y5-Y6)*FtopOut)+dwy6*sigma_y)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=1))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=1))
    
    }else if (type == 'g'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)+dt*(Y1-x1)*a1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)+dt*(Y2-x2)*a2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)+dt*(Y3-x3)*a3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)+dt*(Y4-x4)*a4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)+dt*(Y5-x5)*a5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*(Y6-x6)*a6+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dwy6*sigma_y)


    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    
    }else if (type == 'g2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*(FbotIn*vbot+k1)+(x2-x1)*(FbotOut*f1+k1)-Fbot*Y1*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*Y2*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*Y3*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k5)+(Tmid-x4)*FmidIn*vmid-Fbot*Y4*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*Y5*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*(FtopIn*vtop+k6)+Ftop*Y6*a6)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((Tbot-x1)*(FbotIn*vbot+k1)+(x2-x1)*(FbotOut*f1+k1))+dwy1*sigma_x)
    model$addSystem(dY2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dwy2*sigma_x)
    model$addSystem(dY3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dwy3*sigma_x)
    model$addSystem(dY4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k5)+(Tmid-x4)*FmidIn*vmid)+dwy4*sigma_x)
    model$addSystem(dY5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dwy5*sigma_x)
    model$addSystem(dY6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*(FtopIn*vtop+k6))+dwy6*sigma_x)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=100))
    
    }else if (type == 'g3'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*(FbotIn*vbot+k1)+(x2-x1)*(FbotOut*f1+k1)-Fbot*Y1*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*Y2*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*Y3*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k5)+(Tmid-x4)*FmidIn*vmid-Fbot*Y4*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*Y5*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*(FtopIn*vtop+k6)+Ftop*Y6*a6)+dw6*sigma_x)

    model$addSystem(dY1~dt*(x1-Y1)+dwy1*sigma_x)
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*sigma_x)
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*sigma_x)
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*sigma_x)
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*sigma_x)
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*sigma_x)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=100))
    
    }else if (type == 'g4'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*(FbotIn*vbot+k1)+(x2-x1)*(FbotOut*f1+k1)-Fbot*Y1*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*Y2*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*Y3*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k5)+(Tmid-x4)*FmidIn*vmid-Fbot*Y4*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*Y5*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*(FtopIn*vtop+k6)+Ftop*Y6*a6)+dw6*sigma_x)

    model$addSystem(dY1~dt*(x1-Y1)*b1+dwy1*sigma_x)
    model$addSystem(dY2~dt*(x2-Y2)*b2+dwy2*sigma_x)
    model$addSystem(dY3~dt*(x3-Y3)*b3+dwy3*sigma_x)
    model$addSystem(dY4~dt*(x4-Y4)*b4+dwy4*sigma_x)
    model$addSystem(dY5~dt*(x5-Y5)*b5+dwy5*sigma_x)
    model$addSystem(dY6~dt*(x6-Y6)*b6+dwy6*sigma_x)

    model$setParameter(Y1 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y2 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y3 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y4 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y5 = c(init=0, lower=-1, upper=100))
    model$setParameter(Y6 = c(init=0, lower=-1, upper=100))
    
    }else if (type == 'glowpass'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1)*a)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2)*a)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3)*a)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*(Y4-x4)*a)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5)*a)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)*a+(ambientTemp-x6)*u)+dw6*sigma_x)

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

    }else if (type == 'glowpass2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1-Fbot*(Y1-x1)*a1)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2-Fbot*(Y2-x2)*a2)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3-Fbot*(Y3-x3)*a3)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid-Fbot*(Y4-x4)*a4)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+Ftop*(Y5-x5)*a5)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*FtopOut*f6+(Ttop-x6)*FtopIn*vtop+Ftop*(Y6-x6)*a6+(ambientTemp-x6)*u)+dw6*sigma_x)

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
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*f1+(Y1-x1)*u)+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*f2+(x3-x2)*FbotOut*f2+(Y2-x2)*u)+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*f3+(x4-x3)*FbotOut*f3+(Y3-x3)*u)+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*f4+(x5-x4)*FbotOut*f4+(Tmid-x4)*FmidIn*vmid+(Y4-x4)*u)+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*FtopOut*f5+(x6-x5)*FtopIn*f5+(Y5-x5)*u)+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k)+(Ttop-x6)*FtopIn*vtop+(Y6-x6)*u)+dw6*sigma_x)

    model$addSystem(dY1~dt/V1*((x2-Y1)*u+(Y2-Y1)*k1)+dwy1*sigma_y)
    model$addSystem(dY2~dt/V2*((x3-Y2)*u+(Y3-Y2)*k1+(Y1-Y2)*k1)+dwy2*sigma_y)
    model$addSystem(dY3~dt/V3*((x4-Y3)*u+(Y4-Y3)*k1+(Y2-Y3)*k1)+dwy3*sigma_y)
    model$addSystem(dY4~dt/V4*((x5-Y4)*u+(Y5-Y4)*k1+(Y3-Y4)*k1)+dwy4*sigma_y)
    model$addSystem(dY5~dt/V5*((x6-Y5)*u+(Y6-Y5)*k1+(Y4-Y5)*k1)+dwy5*sigma_y)
    model$addSystem(dY6~dt/V6*((x5-Y6)*u+(Y5-Y6)*k1+(ambientTemp-Y6)*k2)+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$ambientTemp[1])

    } else if (type == 'random'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*(Y1-x1)*a+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dt*(Y2-x2)*a+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dt*(Y3-x3)*a+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dt*(Y4-x4)*a+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dt*(Y5-x5)*a+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*(Y6-x6)*a+dw6*sigma_x)

    
    model$addSystem(dY1~dt*(x1-Y1)*b+dwy1*sigma_y)
    model$addSystem(dY2~dt*(x2-Y2)*b+dwy2*sigma_y)
    model$addSystem(dY3~dt*(x3-Y3)*b+dwy3*sigma_y)
    model$addSystem(dY4~dt*(x4-Y4)*b+dwy4*sigma_y)
    model$addSystem(dY5~dt*(x5-Y5)*b+dwy5*sigma_y)
    model$addSystem(dY6~dt*(x6-Y6)*b+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])


    model$setParameter(a = c(init=0, lower=-1, upper=1))
    model$setParameter(b = c(init=0, lower=-1, upper=1))
    model$setParameter(f1=0)
    model$setParameter(f2=0)
    model$setParameter(f3=0)
    model$setParameter(f4=0)
    model$setParameter(f5=0)
    model$setParameter(f6=0)
    model$setParameter(k1=0)
    model$setParameter(k2=0)
    model$setParameter(k3=0)
    model$setParameter(k4=0)
    model$setParameter(k5=0)
    model$setParameter(k6=0)
    model$setParameter(vbot=0)
    model$setParameter(vmid=0)
    model$setParameter(vtop=0)

    
    }else if (type == 'random2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*Y1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dt*Y2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dt*Y3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dt*Y4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dt*Y5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*Y6+dw6*sigma_x)

    model$addSystem(dY1~dwy1*sigma_y)
    model$addSystem(dY2~dwy2*sigma_y)
    model$addSystem(dY3~dwy3*sigma_y)
    model$addSystem(dY4~dwy4*sigma_y)
    model$addSystem(dY5~dwy5*sigma_y)
    model$addSystem(dY6~dwy6*sigma_y)

    model$setParameter(Y1 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y2 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y3 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y4 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y5 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y6 = c(init = 0, lb = -1, ub = 1))

    
    }else if (type == 'random3'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1))+dt*Y1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2)+(x3-x2)*(FbotOut*f2))+dt*Y2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3))+dt*Y3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4)+(Tmid-x4)*FmidIn*vmid)+dt*Y4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5))+dt*Y5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6)+(Ttop-x6)*FtopIn*vtop)+dt*Y6+dw6*sigma_x)

    model$addSystem(dY1~dwy1*sigma_y)
    model$addSystem(dY2~dwy2*sigma_y)
    model$addSystem(dY3~dwy3*sigma_y)
    model$addSystem(dY4~dwy4*sigma_y)
    model$addSystem(dY5~dwy5*sigma_y)
    model$addSystem(dY6~dwy6*sigma_y)

    model$setParameter(Y1 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y2 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y3 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y4 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y5 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y6 = c(init = 0, lb = -1, ub = 1))

    
    }else if (type == 'random4'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1))*Y1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2)+(x3-x2)*(FbotOut*f2))*Y2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3)+(x4-x3)*(FbotOut*f3))*Y3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4)+(x5-x4)*(FbotOut*f4)+(Tmid-x4)*FmidIn*vmid)*Y4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5)+(x6-x5)*(FtopIn*f5))*Y5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6)+(Ttop-x6)*FtopIn*vtop)*Y6+dw6*sigma_x)

    model$addSystem(dY1~dwy1*sigma_y)
    model$addSystem(dY2~dwy2*sigma_y)
    model$addSystem(dY3~dwy3*sigma_y)
    model$addSystem(dY4~dwy4*sigma_y)
    model$addSystem(dY5~dwy5*sigma_y)
    model$addSystem(dY6~dwy6*sigma_y)

    model$setParameter(Y1 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y2 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y3 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y4 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y5 = c(init = 0, lb = -1, ub = 1))
    model$setParameter(Y6 = c(init = 0, lb = -1, ub = 1))



    }else if (type == 'ar1'){

    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*(x1-Y1)*a+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dt*(x2-Y2)*a+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dt*(x3-Y3)*a+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dt*(x4-Y4)*a+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dt*(x5-Y5)*a+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*(x6-Y6)*a+dw6*sigma_x)

    
    model$addSystem(dY1~dt*(x1-Y1)+dwy1*0)
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*0)
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*0)
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*0)
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*0)
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*0)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])
    }else if (type == 'b3_ar1'){

    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*(x1-Y1)*a1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*v)+dt*(x2-Y2)*a2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*v)+dt*(x3-Y3)*a3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*v)+dt*(x4-Y4)*a4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)-Fbot*(x6-x5)*(x5-x4)*v)+dt*(x5-Y5)*a5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*(x6-Y6)*a+dw6*sigma_x)

    
    model$addSystem(dY1~dt*(x1-Y1)+dwy1*0)
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*0)
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*0)
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*0)
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*0)
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*0)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])

    }else if (type == 'b3_ar1_2'){

    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*Fbot*(x1-Y1)*a1/V1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*v)+dt*Fbot*(x2-Y2)*a2/V2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*v)+dt*Fbot*(x3-Y3)*a3/V3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*v)+dt*Fbot*(x4-Y4)*a4/V4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)-Fbot*(x6-x5)*(x5-x4)*v)+dt*Ftop*(x5-Y5)*a5/V5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*Ftop*(x6-Y6)*a6/V6+dw6*sigma_x)

    
    model$addSystem(dY1~dt*(x1-Y1)+dwy1*0)
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*0)
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*0)
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*0)
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*0)
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*0)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])
    }else if (type == 'test'){

    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*(x1-Y1)*h1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dt*(x2-Y2)*h2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dt*(x3-Y3)*h3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dt*(x4-Y4)*h4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dt*(x5-Y5)*h5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop+(ambientTemp-x6)*u)+dt*(x6-Y6)*h6+dw6*sigma_x)

    model$addSystem(dY1~dt*(x1-Y1)+dwy1*sigma_y)
    model$addSystem(dY2~dt*(x2-Y2)+dwy2*sigma_y)
    model$addSystem(dY3~dt*(x3-Y3)+dwy3*sigma_y)
    model$addSystem(dY4~dt*(x4-Y4)+dwy4*sigma_y)
    model$addSystem(dY5~dt*(x5-Y5)+dwy5*sigma_y)
    model$addSystem(dY6~dt*(x6-Y6)+dwy6*sigma_y)

    model$setParameter(Y1 = data$X1[1])
    model$setParameter(Y2 = data$X2[1])
    model$setParameter(Y3 = data$X3[1])
    model$setParameter(Y4 = data$X4[1])
    model$setParameter(Y5 = data$X5[1])
    model$setParameter(Y6 = data$X6[1])


    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)

    model$setParameter(k1 = 0)
    model$setParameter(k2 = 0)
    model$setParameter(k3 = 0)
    model$setParameter(k4 = 0)
    model$setParameter(k5 = 0)
    model$setParameter(k6 = 0)

    model$setParameter(vtop = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vbot = 0)
    
    model$setParameter(a1 = 0)
    model$setParameter(a2 = 0)
    model$setParameter(a3 = 0)
    model$setParameter(a4 = 0)
    model$setParameter(u = 0)


    model$setParameter(h1 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h2 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h3 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h4 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h5 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h6 = c(init = 0, lb = -10, ub = 10))
    

    }else if (type == 'test2'){

    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*dX1*h1+dw1*sigma_x)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dt*dX2*h2+dw2*sigma_x)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dt*dX3*h3+dw3*sigma_x)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dt*dX4*h4+dw4*sigma_x)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dt*dX5*h5+dw5*sigma_x)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*dX6*h6+dw6*sigma_x)

    # model$addSystem(dY1~dt*(x2-x1)*a1+dwy1*sigma_y)
    # model$addSystem(dY2~dt*(x3-2*x2+x1)*a2+dwy2*sigma_y)
    # model$addSystem(dY3~dt*(x4-2*x3+x2)*a3+dwy3*sigma_y)
    # model$addSystem(dY4~dt*(x5-2*x4+x3)*a4+dwy4*sigma_y)
    # model$addSystem(dY5~dt*(x6-2*x5+x4)*a5+dwy5*sigma_y)
    # model$addSystem(dY6~dt*(x6-x5)*a6+dwy6*sigma_y)

    # model$setParameter(Y1 = c(init = 0, lb = -1, ub = 1))
    # model$setParameter(Y2 = c(init = 0, lb = -1, ub = 1))
    # model$setParameter(Y3 = c(init = 0, lb = -1, ub = 1))
    # model$setParameter(Y4 = c(init = 0, lb = -1, ub = 1))
    # model$setParameter(Y5 = c(init = 0, lb = -1, ub = 1))
    # model$setParameter(Y6 = c(init = 0, lb = -1, ub = 1))

    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)

    model$setParameter(k1 = 0)
    model$setParameter(k2 = 0)
    model$setParameter(k3 = 0)
    model$setParameter(k4 = 0)
    model$setParameter(k5 = 0)
    model$setParameter(k6 = 0)

    model$setParameter(vtop = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vbot = 0)
    
    model$setParameter(a1 = 0)
    model$setParameter(a2 = 0)
    model$setParameter(a3 = 0)
    model$setParameter(a4 = 0)
    model$setParameter(u = 0)

    model$setParameter(h1 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h2 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h3 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h4 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h5 = c(init = 0, lb = -10, ub = 10))
    model$setParameter(h6 = c(init = 0, lb = -10, ub = 10))


    }else if (type == 'random5'){
    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dY2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dw2*sigma_x)
    model$addSystem(dY3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dw3*sigma_x)
    model$addSystem(dY4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dw4*sigma_x)
    model$addSystem(dY5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dw5*sigma_x)
    model$addSystem(dY6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)


    model$addSystem(dx1~dt*Y1+dwy1*sigma_y)
    model$addSystem(dx2~dt*Y2+dwy2*sigma_y)
    model$addSystem(dx3~dt*Y3+dwy3*sigma_y)
    model$addSystem(dx4~dt*Y4+dwy4*sigma_y)
    model$addSystem(dx5~dt*Y5+dwy5*sigma_y)
    model$addSystem(dx6~dt*Y6+dwy6*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)


    }else if (type == 'random5_1'){
    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dw1*sigma_x)
    model$addSystem(dY2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2))+dw2*sigma_x)
    model$addSystem(dY3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3))+dw3*sigma_x)
    model$addSystem(dY4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid)+dw4*sigma_x)
    model$addSystem(dY5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5))+dw5*sigma_x)
    model$addSystem(dY6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dw6*sigma_x)


    model$addSystem(dx1~dt*(Y1-x1)*a1+dwy1*sigma_y)
    model$addSystem(dx2~dt*(Y2-x2)*a2+dwy2*sigma_y)
    model$addSystem(dx3~dt*(Y3-x3)*a3+dwy3*sigma_y)
    model$addSystem(dx4~dt*(Y4-x4)*a4+dwy4*sigma_y)
    model$addSystem(dx5~dt*(Y5-x5)*a5+dwy5*sigma_y)
    model$addSystem(dx6~dt*(Y6-x6)*a6+dwy6*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    
    }else if (type == 'random6'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*(FbotOut*f1+k1))+dt*Y1+dw1*sigma_x1)
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*f2+k2)+(x3-x2)*(FbotOut*f2+k2)-Fbot*(x3-x2)*(x2-x1)*a1)+dt*Y2+dw2*sigma_x2)
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*f3+k3)+(x4-x3)*(FbotOut*f3+k3)-Fbot*(x4-x3)*(x3-x2)*a2)+dt*Y3+dw3*sigma_x3)
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*f4+k4)+(x5-x4)*(FbotOut*f4+k4)+(Tmid-x4)*FmidIn*vmid-Fbot*(x5-x4)*(x4-x3)*a3)+dt*Y4+dw4*sigma_x4)
    model$addSystem(dx5~dt/V5*((x4-x5)*(FtopOut*f5+k5)+(x6-x5)*(FtopIn*f5+k5)+Ftop*(x6-x5)*(x5-x4)*a4)+dt*Y5+dw5*sigma_x5)
    model$addSystem(dx6~dt/V6*((x5-x6)*(FtopOut*f6+k6)+(Ttop-x6)*FtopIn*vtop)+dt*Y6+dw6*sigma_x6)


    model$addSystem(dY1~dwy1*sigma_y)
    model$addSystem(dY2~dwy2*sigma_y)
    model$addSystem(dY3~dwy3*sigma_y)
    model$addSystem(dY4~dwy4*sigma_y)
    model$addSystem(dY5~dwy5*sigma_y)
    model$addSystem(dY6~dwy6*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)

    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)

    model$setParameter(k1 = 0)
    model$setParameter(k2 = 0)
    model$setParameter(k3 = 0)
    model$setParameter(k4 = 0)
    model$setParameter(k5 = 0)
    model$setParameter(k6 = 0)

    model$setParameter(vtop = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vbot = 0)
    
    model$setParameter(a1 = 0)
    model$setParameter(a2 = 0)
    model$setParameter(a3 = 0)
    model$setParameter(a4 = 0)
    }

    # model$setParameter(sigma_X = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_X = 2.083e-5)

    model$setParameter(sigma_x = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x1 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x2 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x3 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x4 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x5 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x6 = c(init=0.1, lower=0, upper=0.5))


    model$setParameter(sigma_y = c(init=1e-2, lower=0, upper=1))
    model$setParameter(sigma_y1 = c(init=1e-2, lower=0, upper=0.5))
    model$setParameter(sigma_y2 = c(init=1e-2, lower=0, upper=0.5))
    model$setParameter(sigma_y3 = c(init=1e-2, lower=0, upper=0.5))
    model$setParameter(sigma_y4 = c(init=1e-2, lower=0, upper=0.5))
    model$setParameter(sigma_y5 = c(init=1e-2, lower=0, upper=0.5))
    model$setParameter(sigma_y6 = c(init=1e-2, lower=0, upper=0.5))


    model$setParameter(sigma_x1 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x2 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x3 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x4 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x5 = c(init=0.1, lower=0, upper=0.5))
    model$setParameter(sigma_x6 = c(init=0.1, lower=0, upper=0.5))


    model$setParameter(g1 = c(init = 1e-3, lb = 0, ub =1))
    model$setParameter(g2 = c(init = 1e-3, lb = 0, ub = 1))
    model$setParameter(g3 = c(init= 1e-3, lb = 0, ub = 1))
    model$setParameter(g4 = c(init = 1e-3, lb = 0, ub = 1))
    model$setParameter(u = c(init = -0.1, lb = -1, ub = 1))
    model$setParameter(v = c(init = 1e-3, lb = 0, ub = 1))

    model$setParameter(f1 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f2 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f3 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f4 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f5 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f6 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f7 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f8 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f9 = c(init=0.1, lower=0, upper=100))
    model$setParameter(f10 = c(init=0.1, lower=0, upper=100))

    model$setParameter(k1 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k2 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k3 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k4 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k5 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k6 = c(init=0.1, lower=-100, upper=1000))
    model$setParameter(k = c(init=0.1, lower=0, upper=1000))

    model$setParameter(vtop = c(init=0.1, lower=0, upper=100))
    model$setParameter(vmid = c(init=0.1, lower=0, upper=100))
    model$setParameter(vbot = c(init=0.1, lower=0, upper=100))

    model$setParameter(a = c(init=0.1, lower=0, upper=10))
    model$setParameter(b = c(init=0.1, lower=0, upper=10))
    model$setParameter(u=c(init=0.1, lower=0, upper=100))

    model$setParameter(a1 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(a2 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(a3 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(a4 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(a5 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(a6 = c(init=1e-3, lower=-10, upper=10))

    model$setParameter(b1 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(b2 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(b3 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(b4 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(b5 = c(init=1e-3, lower=-10, upper=10))
    model$setParameter(b6 = c(init=1e-3, lower=-10, upper=10))




    model$addInput('Ttop','Tmid','Tbot','FtopIn','FtopOut','FmidIn','FmidOut','FbotIn','FbotOut','Ftop','Fmid','Fbot', 'ambientTemp')
    model$addInput('dX1','dX2','dX3','dX4','dX5','dX6')

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

data$dX1 <- c(0,diff(data$X1))
data$dX2 <- c(0,diff(data$X2))
data$dX3 <- c(0,diff(data$X3))
data$dX4 <- c(0,diff(data$X4))
data$dX5 <- c(0,diff(data$X5))
data$dX6 <- c(0,diff(data$X6))

basic3 <- load_fit('models/ctsm/basic2.RData')


D <- data[2000:8000,]

type <- 'random5_1'
model <- makemodel(type, data = D)
model <- setInitialState(model, D)
# model <- setParams(model, basic3)
model$ParameterValues

fit <- model$estimate(data = D, firstorder = TRUE)
model_dir <- 'models/ctsm/'
dir.create(model_dir, showWarnings = FALSE)

save(fit, file = paste0(model_dir, type ,'.RData'))



