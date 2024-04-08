make_model<-function(model){

    model$addObs(X15 ~ x15m)

    # Set observation equation variances

    model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    obs_std <- 0.15
    model$setParameter(sigma_X = obs_std)

    #Addsystemequations
    model$addSystem(dx15m~dt/V15*(FtopIn*(Ttop-x15m)*Vtop+(X14-x15m)*(ktop+ftop*dFtop)+(ambientTemp-x15m)*Utop)+sigma_x*dw15)

    

    model$addInput("X1","X2","X3","X4","X5","X6","X7","X8","X9","X11","X12","X13","X14","dFtopIn","dFtopOut","dFtop","dFtopVol")
    #####SETPARAMETERS########
    #SystemNoise
    model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))
    model$setParameter(sigma_y=c(init=0.001,lb=0,ub=10))
    #SystemParameters
    model$setParameter(Vmid=c(init=1,lb=0,ub=10))
    model$setParameter(Vbot=c(init=1,lb=0,ub=10))
    model$setParameter(Vtop=c(init=1,lb=0,ub=10))
    model$setParameter(Vsoil=c(init=1,lb=0,ub=10))

    model$setParameter(Ubot=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(Umid=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(Utop=c(init=0,lb=-1e2,ub=1e2))

    model$setParameter(kbot=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(kmid=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(ktop=c(init=1e-6,lb=-1e2,ub=1e2))


    model$setParameter(ftop=c(init=1,lb=-1e3,ub=1e3))
    model$setParameter(fmid=c(init=0,lb=-1e3,ub=1e3))
    model$setParameter(v=c(init=1e3,lb=0,ub=1e6))

    model$setParameter(a=c(init=0,lb=-1e3,ub=1e1))
    model$setParameter(b=c(init=0,lb=-1e3,ub=1e1))

    model$setParameter(c=c(init=1e-2,lb=-1e2,ub=1e2))
    model$setParameter(d=c(init=1e-2,lb=-1e2,ub=1e2))

    return(model)}