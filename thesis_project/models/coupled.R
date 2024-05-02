make_model<-function(model, type = 'model1'){
    # Add observation equations   
    model$addObs(X0~x0m)
    model$addObs(X1~x1m)
    model$addObs(X2~x2m)
    model$addObs(X3~x3m)
    model$addObs(X4~x4m)
    model$addObs(X5~x5m)
    model$addObs(X6~x6m)
    model$addObs(X7~x7m)
    model$addObs(X8~x8m)
    model$addObs(X9~x9m)
    model$addObs(X10~x10m)
    model$addObs(X11~x11m)
    model$addObs(X12~x12m)
    model$addObs(X13~x13m)
    model$addObs(X14~x14m)
    model$addObs(X15~x15m)

    # Set observation equation variances
    model$setVariance(X0 ~ sigma_X^2)
    model$setVariance(X1 ~ sigma_X^2)
    model$setVariance(X2 ~ sigma_X^2)
    model$setVariance(X3 ~ sigma_X^2)
    model$setVariance(X4 ~ sigma_X^2)
    model$setVariance(X5 ~ sigma_X^2)
    model$setVariance(X6 ~ sigma_X^2)
    model$setVariance(X7 ~ sigma_X^2)
    model$setVariance(X8 ~ sigma_X^2)
    model$setVariance(X9 ~ sigma_X^2)
    model$setVariance(X10 ~ sigma_X^2)
    model$setVariance(X11 ~ sigma_X^2)
    model$setVariance(X12 ~ sigma_X^2)
    model$setVariance(X13 ~ sigma_X^2)
    model$setVariance(X14 ~ sigma_X^2)
    model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-2,lb=0,ub=1))

    if (type == 'model1'){
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*v*FbotIn+(x1m-x0m)*(f0*FbotOut))/V0)+sigma_x*dw0)  
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((x0m-x1m)*f1*FbotIn+(x2m-x1m)*(f2*FbotOut))/V1)+sigma_x*dw1)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((x1m-x2m)*f3*FbotIn+(X3-x2m)*(f4*FbotOut))/V2)+sigma_x*dw2)

        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dw3)
        model$addSystem(dY1m~dt*(x1m-Y1m)*b+sigma_y*dw4)
        model$addSystem(dY2m~dt*(x2m-Y2m)*b+sigma_y*dw5)
    } else if (type == 'simple'){
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*v*FbotIn+(x1m-x0m)*(f2*FbotOut))/V0)+sigma_x*dwx0)  
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((x0m-x1m)*f1*FbotIn+(x2m-x1m)*(f2*FbotOut))/V1)+sigma_x*dwx1)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((x1m-x2m)*f1*FbotIn+(x3m-x2m)*(f2*FbotOut))/V2)+sigma_x*dwx2)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((x2m-x3m)*f3*FbotIn+(x4m-x3m)*(f4*FbotOut))/V3)+sigma_x*dwx3)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((x3m-x4m)*f3*FbotIn+(x5m-x4m)*(f4*FbotOut))/V4)+sigma_x*dwx4)
        model$addSystem(dx5m~dt*((Y5m-x5m)*a+((x4m-x5m)*f4*FbotIn+(X6-x5m)*(f4*FbotOut))/V5)+sigma_x*dwx5)


        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dwy0)
        model$addSystem(dY1m~dt*(x1m-Y1m)*b+sigma_y*dwy1)
        model$addSystem(dY2m~dt*(x2m-Y2m)*b+sigma_y*dwy2)
        model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwy3)
        model$addSystem(dY4m~dt*(x4m-Y4m)*b+sigma_y*dwy4)
        model$addSystem(dY5m~dt*(x5m-Y5m)*b+sigma_y*dwy5)


    } else if (type == 'nonlinear'){
        model$addSystem(dx0m~dt*(((Tbot-x0m)*v*FbotIn+(x1m-x0m)*(f2*FbotOut))/V0)+sigma_x*dwx0)  
        model$addSystem(dx1m~dt*(((x0m-x1m)*f1*FbotIn+(x2m-x1m)*(f2*FbotOut)-V*Fbot*(x2m-x1m)*(x1m-x0m))/V1)+sigma_x*dwx1)
        model$addSystem(dx2m~dt*(((x1m-x2m)*f1*FbotIn+(x3m-x2m)*(f2*FbotOut)-V*Fbot*(x2m-x1m)*(x1m-x0m))/V2)+sigma_x*dwx2)
        model$addSystem(dx3m~dt*(((x2m-x3m)*f3*FbotIn+(x4m-x3m)*(f4*FbotOut)-V*Fbot*(x2m-x1m)*(x1m-x0m))/V3)+sigma_x*dwx3)
        model$addSystem(dx4m~dt*(((x3m-x4m)*f3*FbotIn+(x5m-x4m)*(f4*FbotOut)-V*Fbot*(x2m-x1m)*(x1m-x0m))/V4)+sigma_x*dwx4)
        model$addSystem(dx5m~dt*(((x4m-x5m)*f4*FbotIn+(X6-x5m)*(f4*FbotOut)-V*Fbot*(x2m-x1m)*(x1m-x0m))/V5)+sigma_x*dwx5)



    }else if (type == 'model2'){
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*v*FbotIn+(x1m-x0m)*(f2*FbotOut))/V0)+sigma_x*dwx0)  
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((x0m-x1m)*f1*FbotIn+(x2m-x1m)*(f2*FbotOut))/V1)+sigma_x*dwx1)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((x1m-x2m)*f1*FbotIn+(x3m-x2m)*(f2*FbotOut))/V2)+sigma_x*dwx2)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((x2m-x3m)*f1*FbotIn+(x4m-x3m)*(f2*FbotOut))/V3)+sigma_x*dwx3)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((x3m-x4m)*f1*FbotIn+(x5m-x4m)*(f2*FbotOut))/V4)+sigma_x*dwx4)
        model$addSystem(dx5m~dt*((Y5m-x5m)*a+((x4m-x5m)*f1*FbotIn+(x6m-x5m)*(f2*FbotOut))/V5)+sigma_x*dwx5)
        model$addSystem(dx6m~dt*((Y6m-x6m)*a+((x5m-x6m)*f1*FbotIn+(x7m-x6m)*(f2*FbotOut))/V6)+sigma_x*dwx6)
        model$addSystem(dx7m~dt*((Y7m-x7m)*a+((x6m-x7m)*f1*FbotIn+(x8m-x7m)*(f2*FbotOut))/V7)+sigma_x*dwx7)
        model$addSystem(dx8m~dt*((Y8m-x8m)*a+((x7m-x8m)*f1*FbotIn+(x9m-x8m)*(f2*FbotOut))/V8)+sigma_x*dwx8)
        model$addSystem(dx9m~dt*((Y9m-x9m)*a+((x8m-x9m)*f1*FbotIn+(X10-x9m)*(f2*FbotOut))/V9)+sigma_x*dwx9)



        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dwy0)
        model$addSystem(dY1m~dt*(x1m-Y1m)*b+sigma_y*dwy1)
        model$addSystem(dY2m~dt*(x2m-Y2m)*b+sigma_y*dwy2)
        model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwy3)
        model$addSystem(dY4m~dt*(x4m-Y4m)*b+sigma_y*dwy4)
        model$addSystem(dY5m~dt*(x5m-Y5m)*b+sigma_y*dwy5)
        model$addSystem(dY6m~dt*(x6m-Y6m)*b+sigma_y*dwy6)
        model$addSystem(dY7m~dt*(x7m-Y7m)*b+sigma_y*dwy7)
        model$addSystem(dY8m~dt*(x8m-Y8m)*b+sigma_y*dwy8)
        model$addSystem(dY9m~dt*(x9m-Y9m)*b+sigma_y*dwy9)

    }else if (type == 'model3'){
        model$addSystem(dx0m~dt*(((Tbot-x0m)*v*FbotIn+(x1m-x0m)*(f2*FbotOut))/V0)+sigma_x*dwx0)  
        model$addSystem(dx1m~dt*(((x0m-x1m)*f1*FbotIn+(x2m-x1m)*(f2*FbotOut))/V1)+sigma_x*dwx1)
        model$addSystem(dx2m~dt*(((x1m-x2m)*f1*FbotIn+(x3m-x2m)*(f2*FbotOut))/V2)+sigma_x*dwx2)
        model$addSystem(dx3m~dt*(((x2m-x3m)*f1*FbotIn+(x4m-x3m)*(f2*FbotOut))/V3)+sigma_x*dwx3)
        model$addSystem(dx4m~dt*(((x3m-x4m)*f1*FbotIn+(x5m-x4m)*(f2*FbotOut))/V4)+sigma_x*dwx4)
        model$addSystem(dx5m~dt*(((x4m-x5m)*f1*FbotIn+(x6m-x5m)*(f2*FbotOut))/V5)+sigma_x*dwx5)
        model$addSystem(dx6m~dt*(((x5m-x6m)*f1*FbotIn+(x7m-x6m)*(f2*FbotOut))/V6)+sigma_x*dwx6)
        model$addSystem(dx7m~dt*(((x6m-x7m)*f1*FbotIn+(x8m-x7m)*(f2*FbotOut))/V7)+sigma_x*dwx7)
        model$addSystem(dx8m~dt*(((x7m-x8m)*f1*FbotIn+(x9m-x8m)*(f2*FbotOut))/V8)+sigma_x*dwx8)
        model$addSystem(dx9m~dt*(((x8m-x9m)*f1*FbotIn+(x10m-x9m)*(f2*FbotOut))/V9)+sigma_x*dwx9)
        model$addSystem(dx10m~dt*(((x9m-x10m)*f1*FbotIn+(x11m-x10m)*(f2*FbotOut)+(Tmid-x10m)*fmid*FmidIn)/V10)+sigma_x*dwx10)
        model$addSystem(dx11m~dt*(((x10m-x11m)*f3*FtopOut+(x12m-x11m)*(f4*FtopIn))/V11)+sigma_x*dwx11)
        model$addSystem(dx12m~dt*(((x11m-x12m)*f3*FtopOut+(x13m-x12m)*(f4*FtopIn))/V12)+sigma_x*dwx12)
        model$addSystem(dx13m~dt*(((x12m-x13m)*f3*FtopOut+(x14m-x13m)*(f4*FtopIn))/V13)+sigma_x*dwx13)
        model$addSystem(dx14m~dt*(((x13m-x14m)*f3*FtopOut+(x15m-x14m)*(f4*FtopIn))/V14)+sigma_x*dwx14)
        model$addSystem(dx15m~dt*(((Ttop-x15m)*(ftop)*(FtopIn)+(x14m-x15m)*(f5*FtopOut)+(ambientTemp-x15m)*u)/V15)+sigma_x*dwx15)
    }
  
    # Inputs
    model$addInput("X10",'FbotIn','FbotOut', 'Tbot', 'Fbot')
    #####SETPARAMETERS########
    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    # model$setParameter(Y0 = c(init=10,lb=0,ub=100))
    # model$setParameter(Y1 = c(init=10,lb=0,ub=100))
    # model$setParameter(Y2 = c(init=10,lb=0,ub=100))
    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-5,lb=0,ub=1))

    # Parameters
    model$setParameter(v = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f0 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f1 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f2 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f3 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f4 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(f5 = c(init=1e-3,lb=0,ub=2))
    model$setParameter(a = c(init=1e-3,lb=0,ub=2))
    model$setParameter(b = c(init=1e-3,lb=0,ub=2))
    model$setParameter(V= c(init=1e-3,lb=0,ub=2))
    model$setParameter(u = c(init=1e-3,lb=0,ub=2))

    model$setParameter(fmid = c(init=1e-3,lb=0,ub=2))
    model$setParameter(ftop = c(init=1e-3,lb=0,ub=2))
    return(model)}  