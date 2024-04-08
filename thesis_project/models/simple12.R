make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X12~x12m)
    # Set observation equation variances
    model$setVariance(X12 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx12m~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx12m~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx12m~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx12m~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y-x12m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y+((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut))/V12-x12m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y-x12m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y+((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut))/V12-x12m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y-x12m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y+((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k2+f2*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)/V12-x12m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y-x12m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V12*((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx12m~dt*(Y+((x12m-X11)*(k1+f1*FtopIn)+(x12m-X13)*(k1+f1*FtopOut)+(X13-x12m)*(x12m-X11)*Ftop*v)/V12-x12m)+sigma_x*dw1)

    }

    # Add Inputs
    model$addInput("X11","X13")

     #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=1e-30,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters
    model$setParameter(k1 = c(init=-20,lb=-30,ub=30))
    model$setParameter(k2 = c(init=-20,lb=-30,ub=30))
    model$setParameter(f1 = c(init=0,lb=-20,ub=1))
    model$setParameter(f2 = c(init=5e-1,lb=-2,ub=2))
    model$setParameter(v = c(init=1e-2,lb=0,ub=2))
    return(model)}