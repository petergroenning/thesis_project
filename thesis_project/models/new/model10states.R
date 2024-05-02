library(ctsmr)

data <- read.csv('data/processed/data2m.csv')


makemodel <- function(type = 'basic', data = NULL){
    model <- ctsm$new()

    model$addObs(X1 ~ x1) # Input Layer (2m)
    model$addObs(X2 ~ x2) # (2m)
    model$addObs(X3 ~ x3) # (2m)
    model$addObs(X4 ~ x4) # (2m)
    model$addObs(X5 ~ x5) # (2m)
    model$addObs(X6 ~ x6) # Input Layer (2m)
    model$addObs(X7 ~ x7) # (2m)
    model$addObs(X8 ~ x8) # (1m)
    model$addObs(X9 ~ x9) # Input (1m)

    model$setVariance(X1 ~ sigma_X^2)
    model$setVariance(X2 ~ sigma_X^2)
    model$setVariance(X3 ~ sigma_X^2)
    model$setVariance(X4 ~ sigma_X^2)
    model$setVariance(X5 ~ sigma_X^2)
    model$setVariance(X6 ~ sigma_X^2)
    model$setVariance(X7 ~ sigma_X^2)
    model$setVariance(X8 ~ sigma_X^2)
    model$setVariance(X9 ~ sigma_X^2)

    model$setParameter(V1 = 1859)
    model$setParameter(V2 = 2959)
    model$setParameter(V3 = 4317)
    model$setParameter(V4 = 5929)
    model$setParameter(V5 = 7798)
    model$setParameter(V6 = 9923)
    model$setParameter(V7 = 12303)
    model$setParameter(V8 = 7125)
    model$setParameter(V9 = 7816)


    if (type == 'randomwalk'){
    model$addSystem(dx1~dw1*exp(sigma_x))
    model$addSystem(dx2~dw2*exp(sigma_x))
    model$addSystem(dx3~dw3*exp(sigma_x))
    model$addSystem(dx4~dw4*exp(sigma_x))
    model$addSystem(dx5~dw5*exp(sigma_x))
    model$addSystem(dx6~dw6*exp(sigma_x))
    model$addSystem(dx7~dw7*exp(sigma_x))
    model$addSystem(dx8~dw8*exp(sigma_x))
    model$addSystem(dx9~dw9*exp(sigma_x))


   }else if (type == 'basic1'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn+(x2-x1)*FbotOut*exp(f))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f)+(x3-x2)*FbotOut*exp(f))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f)+(x4-x3)*FbotOut*exp(f))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f)+(x5-x4)*FbotOut*exp(f))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f)+(x6-x5)*FbotOut*exp(f))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f)+(x7-x6)*FbotOut*exp(f)+(Tmid-x7)*FmidIn)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f)+(x8-x7)*FtopIn*exp(f))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f)+(x9-x8)*FtopIn*exp(f))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f)+(Ttop-x9)*FtopIn)+dw9*exp(sigma_x))


    } else if (type == 'basic2'){
    # Adding K to the model
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn+(x2-x1)*(FbotOut*exp(f)+exp(k)))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f)+exp(k))+(x3-x2)*(FbotOut*exp(f)+exp(k)))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f)+exp(k))+(x4-x3)*(FbotOut*exp(f)+exp(k)))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f)+exp(k))+(x5-x4)*(FbotOut*exp(f)+exp(k)))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FbotIn*exp(f)+exp(k))+(x6-x5)*(FbotOut*exp(f)+exp(k)))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FbotIn*exp(f)+exp(k))+(x7-x6)*(FbotOut*exp(f)+exp(k))+(Tmid-x7)*FmidIn)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*(FtopOut*exp(f)+exp(k))+(x8-x7)*(FtopIn*exp(f)+exp(k)))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*(FtopOut*exp(f)+exp(k))+(x9-x8)*(FtopIn*exp(f)+exp(k)))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*(FtopOut*exp(f)+exp(k))+(Ttop-x9)*FtopIn)+dw9*exp(sigma_x))
    } else if (type == 'basic3'){
    # Adding 4 k's to the model
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn+(x2-x1)*(FbotOut*exp(f)+exp(k1)))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn*exp(f)+exp(k1))+(x3-x2)*(FbotOut*exp(f)+exp(k2)))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn*exp(f)+exp(k2))+(x4-x3)*(FbotOut*exp(f)+exp(k3)))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn*exp(f)+exp(k3))+(x5-x4)*(FbotOut*exp(f)+exp(k4)))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FbotIn*exp(f)+exp(k4))+(x6-x5)*(FbotOut*exp(f)+exp(k5)))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FbotIn*exp(f)+exp(k5))+(x7-x6)*(FbotOut*exp(f)+exp(k6))+(Tmid-x7)*FmidIn)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*(FtopOut*exp(f)+exp(k6))+(x8-x7)*(FtopIn*exp(f)+exp(k7)))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*(FtopOut*exp(f)+exp(k7))+(x9-x8)*(FtopIn*exp(f)+exp(k8)))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*(FtopOut*exp(f)+exp(k8))+(Ttop-x9)*FtopIn)+dw9*exp(sigma_x))

    }else if (type == 'basic4'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f)+(x3-x2)*FbotOut*exp(f))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f)+(x4-x3)*FbotOut*exp(f))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f)+(x5-x4)*FbotOut*exp(f))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f)+(x6-x5)*FbotOut*exp(f))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f)+(x7-x6)*FbotOut*exp(f)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f)+(x8-x7)*FtopIn*exp(f))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f)+(x9-x8)*FtopIn*exp(f))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))

    }else if (type == 'basic5'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f1)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f2)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f3)+(x5-x4)*FbotOut*exp(f4))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f4)+(x6-x5)*FbotOut*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f5)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f6)+(x8-x7)*FtopIn*exp(f7))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f7)+(x9-x8)*FtopIn*exp(f8))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f8)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))



    }else if (type == 'basic5_1'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))

    }else if (type == 'basic5_2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*(FbotIn+FbotInk1)*vbot+(x2-x1)*(FbotOut+FbotOutk1)*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*(FbotIn+FbotInk1)*exp(f2)+(x3-x2)*(FbotOut+FbotOutk1)*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*(FbotIn+FbotInk1)*exp(f3)+(x4-x3)*(FbotOut+FbotOutk1)*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*(FbotIn+FbotInk1)*exp(f4)+(x5-x4)*(FbotOut+FbotOutk1)*exp(f4))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*(FbotIn+FbotInk1)*exp(f5)+(x6-x5)*(FbotOut+FbotOutk1)*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*(FbotIn+FbotInk1)*exp(f6)+(x7-x6)*(FbotOut+FbotOutk1)*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*(FtopOut+FtopOutk1)*exp(f7)+(x8-x7)*(FtopIn+FtopInk1)*exp(f7))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*(FtopOut+FtopOutk1)*exp(f8)+(x9-x8)*(FtopIn+FtopInk1)*exp(f8))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*(FtopOut+FtopOutk1)*exp(f9)+(Ttop-x9)*(FtopIn+FtopInk1)*vtop)+dw9*exp(sigma_x))
    
    }else if (type == 'nonlinear'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f1)+(x3-x2)*FbotOut*exp(f2)-exp(a)*Fbot*(x3-x2)*(x2-x1))+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f2)+(x4-x3)*FbotOut*exp(f3)-exp(a)*Fbot*(x4-x3)*(x3-x2))+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f3)+(x5-x4)*FbotOut*exp(f4)-exp(a)*Fbot*(x5-x4)*(x4-x3))+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f4)+(x6-x5)*FbotOut*exp(f5)-exp(a)*Fbot*(x6-x5)*(x5-x4))+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f5)+(x7-x6)*FbotOut*exp(f6)-exp(a)*Fbot*(x7-x6)*(x6-x5)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f6)+(x8-x7)*FtopIn*exp(f7)+exp(a)*Ftop*(x8-x7)*(x7-x6))+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f7)+(x9-x8)*FtopIn*exp(f8)+exp(a)*Ftop*(x9-x8)*(x8-x7))+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f8)+exp(a)*Ftop*(x10-x9)*(x9-x8)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))

    model$setParameter(a = log(c(init = 1e-3, lb = 1e-12, ub = 10)))

    }else if (type == 'model'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dt*Y1+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*Y2+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*Y3+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))+dt*Y4+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))+dt*Y5+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dt*Y6+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))+dt*Y7+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))+dt*Y8+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)+dt*Y9+dw9*exp(sigma_x))

    model$addSystem(dY1~dwy1*sigma_y)
    model$addSystem(dY2~dwy2*sigma_y)
    model$addSystem(dY3~dwy3*sigma_y)
    model$addSystem(dY4~dwy4*sigma_y)
    model$addSystem(dY5~dwy5*sigma_y)
    model$addSystem(dY6~dwy6*sigma_y)
    model$addSystem(dY7~dwy7*sigma_y)
    model$addSystem(dY8~dwy8*sigma_y)
    model$addSystem(dY9~dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)


    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)
    model$setParameter(f7 = 0)
    model$setParameter(f8 = 0)
    model$setParameter(f9 = 0)

    model$setParameter(vbot = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vtop = 0)


    }else if (type == 'model2'){
    model$addSystem(dx1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dt*Y1+dw1*exp(sigma_x))
    model$addSystem(dx2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dt*Y2+dw2*exp(sigma_x))
    model$addSystem(dx3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dt*Y3+dw3*exp(sigma_x))
    model$addSystem(dx4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))+dt*Y4+dw4*exp(sigma_x))
    model$addSystem(dx5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))+dt*Y5+dw5*exp(sigma_x))
    model$addSystem(dx6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dt*Y6+dw6*exp(sigma_x))
    model$addSystem(dx7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))+dt*Y7+dw7*exp(sigma_x))
    model$addSystem(dx8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))+dt*Y8+dw8*exp(sigma_x))
    model$addSystem(dx9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)+dt*Y9+dw9*exp(sigma_x))

    model$addSystem(dY1~dt*a*Y1+dwy1*sigma_y)
    model$addSystem(dY2~dt*a*Y2+dwy2*sigma_y)
    model$addSystem(dY3~dt*a*Y3+dwy3*sigma_y)
    model$addSystem(dY4~dt*a*Y4+dwy4*sigma_y)
    model$addSystem(dY5~dt*a*Y5+dwy5*sigma_y)
    model$addSystem(dY6~dt*a*Y6+dwy6*sigma_y)
    model$addSystem(dY7~dt*a*Y7+dwy7*sigma_y)
    model$addSystem(dY8~dt*a*Y8+dwy8*sigma_y)
    model$addSystem(dY9~dt*a*Y9+dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)


    model$setParameter(f1 = 0)
    model$setParameter(f2 = 0)
    model$setParameter(f3 = 0)
    model$setParameter(f4 = 0)
    model$setParameter(f5 = 0)
    model$setParameter(f6 = 0)
    model$setParameter(f7 = 0)
    model$setParameter(f8 = 0)
    model$setParameter(f9 = 0)

    model$setParameter(vbot = 0)
    model$setParameter(vmid = 0)
    model$setParameter(vtop = 0)
    model$setParameter(a = c(init = 0, lb = -1, ub = 1))

    }else if (type == 'model3'){
    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dY2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dY3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dY4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))+dw4*exp(sigma_x))
    model$addSystem(dY5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dY6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dY7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))+dw7*exp(sigma_x))
    model$addSystem(dY8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))+dw8*exp(sigma_x))
    model$addSystem(dY9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))

    model$addSystem(dx1~dt*(Y1)+dwy1*sigma_y)
    model$addSystem(dx2~dt*(Y2)+dwy2*sigma_y)
    model$addSystem(dx3~dt*(Y3)+dwy3*sigma_y)
    model$addSystem(dx4~dt*(Y4)+dwy4*sigma_y)
    model$addSystem(dx5~dt*(Y5)+dwy5*sigma_y)
    model$addSystem(dx6~dt*(Y6)+dwy6*sigma_y)
    model$addSystem(dx7~dt*(Y7)+dwy7*sigma_y)
    model$addSystem(dx8~dt*(Y8)+dwy8*sigma_y)
    model$addSystem(dx9~dt*(Y9)+dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)


    # model$setParameter(f1 = 0)
    # model$setParameter(f2 = 0)
    # model$setParameter(f3 = 0)
    # model$setParameter(f4 = 0)
    # model$setParameter(f5 = 0)
    # model$setParameter(f6 = 0)
    # model$setParameter(f7 = 0)
    # model$setParameter(f8 = 0)
    # model$setParameter(f9 = 0)

    # model$setParameter(vbot = 0)
    # model$setParameter(vmid = 0)
    # model$setParameter(vtop = 0)
    # model$setParameter(a = c(init = 0, lb = -1, ub = 1))

    }else if (type == 'model4'){
    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))+dw1*exp(sigma_x))
    model$addSystem(dY2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))+dw2*exp(sigma_x))
    model$addSystem(dY3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))+dw3*exp(sigma_x))
    model$addSystem(dY4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))+dw4*exp(sigma_x))
    model$addSystem(dY5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))+dw5*exp(sigma_x))
    model$addSystem(dY6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)+dw6*exp(sigma_x))
    model$addSystem(dY7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))+dw7*exp(sigma_x))
    model$addSystem(dY8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))+dw8*exp(sigma_x))
    model$addSystem(dY9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)+dw9*exp(sigma_x))

    model$addSystem(dx1~dt*(Y1-x1)*a1+dwy1*sigma_y)
    model$addSystem(dx2~dt*(Y2-x2)*a2+dwy2*sigma_y)
    model$addSystem(dx3~dt*(Y3-x3)*a3+dwy3*sigma_y)
    model$addSystem(dx4~dt*(Y4-x4)*a4+dwy4*sigma_y)
    model$addSystem(dx5~dt*(Y5-x5)*a5+dwy5*sigma_y)
    model$addSystem(dx6~dt*(Y6-x6)*a6+dwy6*sigma_y)
    model$addSystem(dx7~dt*(Y7-x7)*a7+dwy7*sigma_y)
    model$addSystem(dx8~dt*(Y8-x8)*a8+dwy8*sigma_y)
    model$addSystem(dx9~dt*(Y9-x9)*a9+dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)

    }else if (type == 'model5'){
    model$addSystem(dY1~dt/V1*((Tbot-x1)*FbotIn*vbot+(x2-x1)*FbotOut*exp(f1))-dt*Y1*b+dw1*exp(sigma_x))
    model$addSystem(dY2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))-dt*Y2*b+dw2*exp(sigma_x))
    model$addSystem(dY3~dt/V3*((x2-x3)*FbotIn*exp(f3)+(x4-x3)*FbotOut*exp(f3))-dt*Y3*b+dw3*exp(sigma_x))
    model$addSystem(dY4~dt/V4*((x3-x4)*FbotIn*exp(f4)+(x5-x4)*FbotOut*exp(f4))-dt*Y4*b+dw4*exp(sigma_x))
    model$addSystem(dY5~dt/V5*((x4-x5)*FbotIn*exp(f5)+(x6-x5)*FbotOut*exp(f5))-dt*Y5*b+dw5*exp(sigma_x))
    model$addSystem(dY6~dt/V6*((x5-x6)*FbotIn*exp(f6)+(x7-x6)*FbotOut*exp(f6)+(Tmid-x7)*FmidIn*vmid)-dt*Y6*b+dw6*exp(sigma_x))
    model$addSystem(dY7~dt/V7*((x6-x7)*FtopOut*exp(f7)+(x8-x7)*FtopIn*exp(f7))-dt*Y7*b+dw7*exp(sigma_x))
    model$addSystem(dY8~dt/V8*((x7-x8)*FtopOut*exp(f8)+(x9-x8)*FtopIn*exp(f8))-dt*Y8*b+dw8*exp(sigma_x))
    model$addSystem(dY9~dt/V9*((x8-x9)*FtopOut*exp(f9)+(Ttop-x9)*FtopIn*vtop)-dt*Y9*b+dw9*exp(sigma_x))

    model$addSystem(dx1~dt*(Y1)*a+dwy1*sigma_y)
    model$addSystem(dx2~dt*(Y2)*a+dwy2*sigma_y)
    model$addSystem(dx3~dt*(Y3)*a+dwy3*sigma_y)
    model$addSystem(dx4~dt*(Y4)*a+dwy4*sigma_y)
    model$addSystem(dx5~dt*(Y5)*a+dwy5*sigma_y)
    model$addSystem(dx6~dt*(Y6)*a+dwy6*sigma_y)
    model$addSystem(dx7~dt*(Y7)*a+dwy7*sigma_y)
    model$addSystem(dx8~dt*(Y8)*a+dwy8*sigma_y)
    model$addSystem(dx9~dt*(Y9)*a+dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)



    model$setParameter(a = log(c(init = 1e-2, lb = 1e-12, ub = 1)))
    model$setParameter(b = log(c(init = 1e-3, lb = 1e-12, ub = 1)))

    }else if (type == 'model6'){
    model$addSystem(dY1~dt/V1*((x2-x1)*FbotOut*exp(f1))-dt*Y1*b+dw1*exp(sigma_x))
    model$addSystem(dY2~dt/V2*((x1-x2)*FbotIn*exp(f2)+(x3-x2)*FbotOut*exp(f2))-dt*Y2*b+dw2*exp(sigma_x))
    model$addSystem(dY3~dt/V3*((x2-x3)*FbotIn*exp(f2)+(x4-x3)*FbotOut*exp(f2))-dt*Y3*b+dw3*exp(sigma_x))
    model$addSystem(dY4~dt/V4*((x3-x4)*FbotIn*exp(f3)+(x5-x4)*FbotOut*exp(f3))-dt*Y4*b+dw4*exp(sigma_x))
    model$addSystem(dY5~dt/V5*((x4-x5)*FbotIn*exp(f3)+(x6-x5)*FbotOut*exp(f3))-dt*Y5*b+dw5*exp(sigma_x))
    model$addSystem(dY6~dt/V6*((x5-x6)*FbotIn*exp(f4)+(x7-x6)*FbotOut*exp(f4))-dt*Y6*b+dw6*exp(sigma_x))
    model$addSystem(dY7~dt/V7*((x6-x7)*FtopOut*exp(f4)+(x8-x7)*FtopIn*exp(f4))-dt*Y7*b+dw7*exp(sigma_x))
    model$addSystem(dY8~dt/V8*((x7-x8)*FtopOut*exp(f5)+(x9-x8)*FtopIn*exp(f5))-dt*Y8*b+dw8*exp(sigma_x))
    model$addSystem(dY9~dt/V9*((x8-x9)*FtopOut*exp(f6))-dt*Y9*b+dw9*exp(sigma_x))

    model$addSystem(dx1~dt*(Y1)*exp(a)+dt/V1*((Tbot-x1)*FbotIn*vbot)+dwy1*sigma_y)
    model$addSystem(dx2~dt*(Y2)*exp(a)+dwy2*sigma_y)
    model$addSystem(dx3~dt*(Y3)*exp(a)+dwy3*sigma_y)
    model$addSystem(dx4~dt*(Y4)*exp(a)+dwy4*sigma_y)
    model$addSystem(dx5~dt*(Y5)*exp(a)+dwy5*sigma_y)
    model$addSystem(dx6~dt*(Y6)*exp(a)+dt/V6*((Tmid-x6)*FmidIn*vmid)+dwy6*sigma_y)
    model$addSystem(dx7~dt*(Y7)*exp(a)+dwy7*sigma_y)
    model$addSystem(dx8~dt*(Y8)*exp(a)+dwy8*sigma_y)
    model$addSystem(dx9~dt*(Y9)*exp(a)+dt/V9*((Ttop-x9)*FtopIn*vtop)+dwy9*sigma_y)

    model$setParameter(Y1 = 0)
    model$setParameter(Y2 = 0)
    model$setParameter(Y3 = 0)
    model$setParameter(Y4 = 0)
    model$setParameter(Y5 = 0)
    model$setParameter(Y6 = 0)
    model$setParameter(Y7 = 0)
    model$setParameter(Y8 = 0)
    model$setParameter(Y9 = 0)

    model$setParameter(a = log(c(init = 1e-2, lb = 1e-10, ub = 1)))
    model$setParameter(b = log(c(init = 1e-3, lb = 1e-10, ub = 1)))
    }
    
    model$addInput('Ttop','Tmid','Tbot','FtopIn','FtopOut','FmidIn','FmidOut','FbotIn','FbotOut','Ftop','Fmid','Fbot', 'ambientTemp')
    model$addInput('dX1','dX2','dX3','dX4','dX5','dX6')
    model$addInput('FbotInk1','FbotOutk1','FmidInk1', 'FtopInk1', 'FtopOutk1','Ttopk1','Tmidk1','Tbotk1')

    obs_noise <- 1/24*1.25e-4
    model$setParameter(sigma_X= sqrt(obs_noise))

    model$setParameter(sigma_x= log(c(init = 1e-1, lb = 1e-10, ub = 1)))
    model$setParameter(sigma_y= log(c(init = 1e-1, lb = 1e-10, ub = 1)))

    if ((type != 'model') | (type != 'model2') | (type != 'model5')){
    model$setParameter(f=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f1 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f2 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f3 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f4 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f5 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f6 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f7 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f8 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(f9 = log(c(init = 1, lb = 1e-12, ub = 10)))



    model$setParameter(k=log(c(init = 1, lb = 1e-12, ub = 10)))

    model$setParameter(k1=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k2=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k3=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k4=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k5=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k6=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k7=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k8=log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(k9=log(c(init = 1, lb = 1e-12, ub = 10)))

    model$setParameter(vbot = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vmid = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vtop = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vtop = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vtop1 = log(c(init = 1, lb = 1e-12, ub = 10)))
    model$setParameter(vtop2 = log(c(init = 1, lb = 1e-12, ub = 10)))


    model$setParameter(a1 = log(c(init = 1.1, lb = 1, ub = 10)))
    model$setParameter(a2 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a3 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a4 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a5 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a6 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a7 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a8 = log(c(init = 1.2, lb = 1, ub = 10)))
    model$setParameter(a9 = log(c(init = 1.2, lb = 1, ub = 10)))
    }
    


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
    model$setParameter(x8 = data$X8[1])
    model$setParameter(x9 = data$X9[1])

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

fit <- load_fit('models/ctsm9states/basic5_1.RData')
data$dX1 <- c(0,diff(data$X1))
data$dX2 <- c(0,diff(data$X2))
data$dX3 <- c(0,diff(data$X3))
data$dX4 <- c(0,diff(data$X4))
data$dX5 <- c(0,diff(data$X5))
data$dX6 <- c(0,diff(data$X6))

data$FbotInk1 <- c(0, data$FbotIn[1:(length(data$FbotIn)-1)])
data$FbotOutk1 <- c(0, data$FbotOut[1:(length(data$FbotOut)-1)])
data$FmidInk1 <- c(0, data$FmidIn[1:(length(data$FmidIn)-1)])
data$FtopInk1 <- c(0, data$FtopIn[1:(length(data$FtopIn)-1)])
data$FtopOutk1 <- c(0, data$FtopOut[1:(length(data$FtopOut)-1)])
data$Ttopk1 <- c(0, data$Ttop[1:(length(data$Ttop)-1)])
data$Tmidk1 <- c(0, data$Tmid[1:(length(data$Tmid)-1)])
data$Tbotk1 <- c(0, data$Tbot[1:(length(data$Tbot)-1)])


D <- data[3000:3500,]

type <- 'model6'
model <- makemodel(type, data = D)
model <- setInitialState(model, D)
# model <- setParams(model, fit)
fit <- model$estimate(data = D, firstorder = TRUE)
model_dir <- 'models/ctsm9states/'
dir.create(model_dir, showWarnings = FALSE)
save(fit, file = paste0(model_dir, type ,'.RData'))
