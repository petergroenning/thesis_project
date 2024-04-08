make_model<-function(model){
    # model$addObs(X0 ~ x0m)
    # model$addObs(X10 ~ x10m)
    model$addObs(X15 ~ x15m)

    # Set observation equation variances
    # model$setVariance(X0 ~ sigma_X^2)
    # model$setVariance(X10 ~ sigma_X^2)
    model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    obs_std <- 0.15
    model$setParameter(sigma_X = obs_std)

    #Addsystemequations
    # model$addSystem(dx0m~dt/V0*(FbotIn*(Tbot-x0m)*Vbot)+sigma_x*dw0)
    # model$addSystem(dx10m~dt/V10*(FmidIn*(Tmid-x10m)*Vmid)+sigma_x*dw10)
    model$addSystem(dx15m~dt/V15*(FtopIn*(Ttop-x15m)*Vtop)+sigma_x*dw15)


    model$addInput("X1","X2","X3","X4","X5","X6","X7","X8","X9","X11","X12","X13","X14")
    #####SETPARAMETERS########
    #SystemNoise
    model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

    #SystemParameters
    model$setParameter(Vmid=c(init=1,lb=0,ub=10))
    model$setParameter(Vbot=c(init=1,lb=0,ub=10))
    model$setParameter(Vtop=c(init=1,lb=0,ub=10))
    model$setParameter(Vsoil=c(init=1,lb=0,ub=10))

    model$setParameter(Ubot=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(Umid=c(init=1,lb=-1e2,ub=1e2))
    model$setParameter(Utop=c(init=1,lb=-1e2,ub=1e2))




    return(model)}