make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X11~x11m)
    # Set observation equation variances
    model$setVariance(X11 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y-x11m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y+((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut))/V11-x11m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y-x11m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y+((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut))/V11-x11m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y-x11m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y11m+((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k2+f2*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)/V11-x11m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY11m~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y11m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V11*((x11m-X10)*(k1+f1*FtopIn)+(x11m-X12)*(k1+f1*FtopOut)+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y-x11m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY11m~dt/V11*((x11m-x10m)*(k1-f1*(FtopOut))+(x11m-X12)*(k1-f1*(FtopIn))+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(Y11m+((x11m-x10m)*(k1-f1*(FtopOut))+(x11m-X12)*(k1-f1*(FtopIn))+(X12-x11m)*(x11m-X10)*Ftop*v)/V11-x11m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY11m~dt/V11*((x11m-X10)*(k-f*(FtopOut))+(x11m-X12)*(k-f*(FtopIn))+(X12-x11m)*(x11m-X10)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx11m~dt*((Y11m-x11m)*a+((x11m-X10)*(k-f*(FtopOut))+(x11m-X12)*(k-f*(FtopIn))+(X12-x11m)*(x11m-X10)*Ftop*v)/V11)+sigma_x*dw1)
        
    }else if (type == 'model1'){
        model$addSystem(dY11m~dt*(((x11m-X10)*(k3-f3*(FtopOut))+(x11m-X12)*(k4-f4*(FtopIn)))/V11)+sigma_y*dwY)
        model$addSystem(dx11m~dt*(((x11m-X10)*(k1-f1*(FtopOut))+(x11m-X12)*(k2-f2*(FtopIn)))/V11+(Y11m-x11m)*a)+(sigma_x)*dw1)

    }else if (type == 'model2'){
        model$addSystem(dY11m~dt*(x11m-Y11m)+sigma_y*dwY)
        model$addSystem(dx11m~dt*((Y11m-x11m)*a+((x11m-X10)*(k1-f1*(FtopOut))+(x11m-X12)*(k2-f2*(FtopIn)))/V11)+(sigma_x)*dw1)

    }else if (type == 'model3'){
        model$addSystem(dY11m~dt*(x11m-Y11m)*b+sigma_y*dwY)
        model$addSystem(dx11m~dt*((Y11m-x11m)*a+((X10-x11m)*(k1+f1*(FtopOut))+(X12-x11m)*(k2+f2*(FtopIn)))/V11)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY11m~dt*(x11m-Y11m)*b+sigma_y*dwY)
        model$addSystem(dx11m~dt*((Y11m-x11m)*a+((X10-x11m)*(f1*(FtopOut))+(X12-x11m)*(f2*(FtopIn)))/V11)+(sigma_x)*dw1)
    }

    # Add Inputs
    model$addInput("X10","X12")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=1e-30,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    # Parameters
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))
    return(model)}