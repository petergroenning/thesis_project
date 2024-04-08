make_model<-function(model){
    # model$addObs(X0 ~ x0m)
    model$addObs(X1 ~ x1m)
    model$addObs(X2 ~ x2m)
    model$addObs(X3 ~ x3m)
    model$addObs(X4 ~ x4m)
    model$addObs(X5 ~ x5m)
    model$addObs(X6 ~ x6m)
    model$addObs(X7 ~ x7m)
    model$addObs(X8 ~ x8m)
    model$addObs(X9 ~ x9m)

    # model$addObs(X10 ~ x10m)
    # model$addObs(X15 ~ x15m)

    # Set observation equation variances
    # model$setVariance(X0 ~ sigma_X^2)
    model$setVariance(X1 ~ sigma_X^2)
    model$setVariance(X2 ~ sigma_X^2)
    model$setVariance(X3 ~ sigma_X^2)
    model$setVariance(X4 ~ sigma_X^2)
    model$setVariance(X5 ~ sigma_X^2)
    model$setVariance(X6 ~ sigma_X^2)
    model$setVariance(X7 ~ sigma_X^2)
    model$setVariance(X8 ~ sigma_X^2)
    model$setVariance(X9 ~ sigma_X^2)

    # model$setVariance(X10 ~ sigma_X^2)
    # model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    obs_std <- 0.15
    model$setParameter(sigma_X = obs_std)

    #Addsystemequations
    # model$addSystem(dx0m~dt*(Y)+dwX0*sigma_x)
    model$addSystem(dx1m~dt*((x2m*(e+f*FbotOut)-X0*(d+g*FbotIn))/2-x1m)/V1 +dwX1*sigma_x)
    model$addSystem(dx2m~dt*((x3m*(e+f*FbotOut)-x1m*(d+g*FbotIn))/2-x2m)/V2 +dwX2*sigma_x)
    model$addSystem(dx3m~dt*((x4m*(e+f*FbotOut)-x2m*(d+g*FbotIn))/2-x3m)/V3 +dwX3*sigma_x)
    model$addSystem(dx4m~dt*((x5m*(e+f*FbotOut)-x3m*(d+g*FbotIn))/2-x4m)/V4 +dwX4*sigma_x)
    model$addSystem(dx5m~dt*((x6m*(e+f*FbotOut)-x4m*(d+g*FbotIn))/2-x5m)/V5 +dwX5*sigma_x)
    model$addSystem(dx6m~dt*((x7m*(e+f*FbotOut)-x5m*(d+g*FbotIn))/2-x6m)/V6 +dwX6*sigma_x)
    model$addSystem(dx7m~dt*((x8m*(e+f*FbotOut)-x6m*(d+g*FbotIn))/2-x7m)/V7 +dwX7*sigma_x)
    model$addSystem(dx8m~dt*((x9m*(e+f*FbotOut)-x7m*(d+g*FbotIn))/2-x8m)/V8 +dwX8*sigma_x)
    model$addSystem(dx9m~dt*((x10m*(e+f*FbotOut)-x8m*(d+g*FbotIn))/2-x9m)/V9 +dwX9*sigma_x)

    # model$addSystem(dx0m~dt/V0*(FbotIn*(Tbot-x0m)*Vbot+(X1-x0m)*(kbot+fbot*FbotOut)+a*(20-x0m))+sigma_x*dw0)
    # model$addSystem(dx10m~dt/V10*(FmidIn*(Tmid-x10m)*Vmid+(X9-x10m)*(kmid1+fmid1*FbotIn)+(X11-x10m)*(kmid2+fmid2*FbotOut))+sigma_x*dw10)
    # model$addSystem(dx15m~dt/V15*(FtopIn*(Ttop-x15m)*Vtop+(X14-x15m)*(ktop+ftop*dFtop)+(ambientTemp-x15m)*Utop)+sigma_x*dw15)
    # model$addSystem(dY~dt/V0*(FbotIn*(Tbot-x0m)*Vbot+(X1-x0m)*(kbot+fbot*FbotOut))+dwY*sigma_y)

    # model$addSystem(dZ~dt*(x0m-Z)*a+dwZ*sigma_z)
    # model$setParameter(Z=c(init=8,lb=0,ub=100))
    model$addInput("X0","X10","X11","X12","X13","X14","dFtopIn","dFtopOut","dFtop","dFtopVol","dFmid","dFbot","Fbotk1","Fbotk2","Fbotk3","dFbotIn")
    #####SETPARAMETERS########
    #SystemNoise
    model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))
    model$setParameter(sigma_y=c(init=0.01,lb=0,ub=10))
    model$setParameter(sigma_z=c(init=0.01,lb=0,ub=10))
    #SystemParameters
    model$setParameter(Vmid=c(init=1,lb=0,ub=10))
    model$setParameter(Vbot=c(init=1,lb=0,ub=10))
    model$setParameter(Vtop=c(init=1,lb=0,ub=10))
    model$setParameter(Vsoil=c(init=1,lb=0,ub=10))

    model$setParameter(Ubot=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(Umid=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(Utop=c(init=0,lb=-1e3,ub=1e3))

    model$setParameter(kbot=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(kmid=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(kmid1=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(kmid2=c(init=1,lb=-1e3,ub=1e3))

    model$setParameter(ktop=c(init=1e-6,lb=-1e3,ub=1e3))


    model$setParameter(ftop=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(fbot=c(init=1,lb=0,ub=1e3))
    model$setParameter(fmid1=c(init=0,lb=-1e3,ub=1e3))
    model$setParameter(fmid2=c(init=0,lb=-1e3,ub=1e3))
    model$setParameter(V=c(init=1e3,lb=0,ub=1e6))

    model$setParameter(a=c(init=1,lb=-1e3,ub=1e1))
    model$setParameter(b=c(init=1,lb=-1e3,ub=1e1))

    model$setParameter(c=c(init=1e-2,lb=-1e3,ub=1e3))
    model$setParameter(d=c(init=1,lb=0,ub=1e3))
    model$setParameter(e=c(init=1,lb=0,ub=1e3))
    model$setParameter(f=c(init=1e-3,lb=0,ub=1e3))
    model$setParameter(g=c(init=1e-3,lb=0,ub=1e3))

    return(model)}