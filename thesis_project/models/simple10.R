make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X10~x10m)
    # Set observation equation variances
    model$setVariance(X10 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx10m~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx10m~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx10m~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx10m~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y-x10m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y+((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)/V10-x10m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y-x10m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y+((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))+(Tmid-x10m)*FmidIn*fmid)/V10-x10m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y-x10m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y+((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k2+f2*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)/V10-x10m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y-x10m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V10*((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY)
        model$addSystem(dx10m~dt*(Y+((x10m-X9)*(k1+f1*(FbotIn+FmidOut))+(x10m-X11)*(k1+f1*(FbotOut+FmidIn))-(X11-x10m)*(x10m-X9)*Fbot*v+(Tmid-x10m)*FmidIn*fmid)/V10-x10m)+sigma_x*dw1)

    }

    # Add Inputs
    model$addInput("X9","X11",'FmidIn','FmidOut','Fmid','Tmid')

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-3,lb=-1e-33,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters
    model$setParameter(k1 = c(init=-20,lb=-30,ub=30))
    model$setParameter(k2 = c(init=-20,lb=-30,ub=30))
    model$setParameter(f1 = c(init=0,lb=-20,ub=1))
    model$setParameter(f2 = c(init=5e-1,lb=-2,ub=2))
    model$setParameter(fmid=c(init=1,lb=0,ub=10))
    model$setParameter(v = c(init=1e-7,lb=-1,ub=2))
    return(model)}