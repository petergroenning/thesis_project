make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X14~x14m)
    # Set observation equation variances
    model$setVariance(X14 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx14m~dt/4*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx14m~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx14m~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx14m~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y-x14m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y+((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut))/V14-x14m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y-x14m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y+((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut))/V14-x14m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y-x14m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y+((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k2+f2*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)/V14-x14m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1+f1*FtopIn)+(x14m-X15)*(k1+f1*FtopOut)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y-x14m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V14*((x14m-X13)*(k1-f1*FtopOut)+(x14m-X15)*(k1-f1*FtopIn)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(Y+((x14m-X13)*(k1-f1*FtopOut)+(x14m-X15)*(k1-f1*FtopIn)+(X15-x14m)*(x14m-X13)*Ftop*v)/V14-x14m)+sigma_x*dw1)

    }else if (type == 'model'){
        model$addSystem(dY14m~dt/V14*((x14m-X13)*(k-f*FtopOut)+(x14m-X15)*(k-f*FtopIn)+(X15-x14m)*(x14m-X13)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx14m~dt*((Y14m-x14m)*a+((x14m-X13)*(k-f*FtopOut)+(x14m-X15)*(k-f*FtopIn)+(X15-x14m)*(x14m-X13)*Ftop*v)/V14)+sigma_x*dw1)
    }else if (type == 'model1'){
        model$addSystem(dY14m~dt*(((x14m-X13)*(k3-f3*(FtopOut))+(x14m-X15)*(k4-f4*(FtopIn)))/V14)+sigma_y*dwY)
        model$addSystem(dx14m~dt*(((x14m-X13)*(k1-f1*(FtopOut))+(x14m-X15)*(k2-f2*(FtopIn)))/V14+(Y14m-x14m)*a)+(sigma_x)*dw1)

    }else if (type == 'model2'){
        model$addSystem(dY14m~dt*(x14m-Y14m)+sigma_y*dwY)
        model$addSystem(dx14m~dt*((Y14m-x14m)*a+((x14m-X13)*(k1-f1*(FtopOut))+(x14m-X15)*(k2-f2*(FtopIn)))/V14)+(sigma_x)*dw1)

    }else if (type == 'model3'){
        model$addSystem(dY14m~dt*(x14m-Y14m)*b+sigma_y*dwY)
        model$addSystem(dx14m~dt*((Y14m-x14m)*a+((X13-x14m)*(k1+f1*(FtopOut))+(X15-x14m)*(k2+f2*FtopIn))/V14)+(sigma_x)*dw1)

    }else if (type == 'model4'){
        model$addSystem(dY14m~dt*((x14m-Y14m)*b)+sigma_y*dwY)
        model$addSystem(dx14m~dt*((Y14m-x14m)*a+((X13-x14m)*(f1*(FtopOut))+(X15-x14m)*(f2*FtopIn))/V14)+(sigma_x)*dw1)

    }

    # Add Inputs
    model$addInput("X13","X15")
    model$addInput('FtopInk1','FtopInk2')
    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=1e-30,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(u = c(init=1e-6,lb=0,ub=40))
    # Parameters
    model$setParameter(k = c(init=10,lb=0,ub=50))
    model$setParameter(f = c(init=1,lb=0,ub=10))

    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))
    return(model)}