make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X13~x13m)
    # Set observation equation variances
    model$setVariance(X13 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx13m~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx13m~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx13m~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx13m~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y-x13m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y+((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut))/V13-x13m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y-x13m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut))+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y+((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut))/V13-x13m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y-x13m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y+((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k2+f2*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)/V13-x13m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1+f1*FtopIn)+(x13m-X14)*(k1+f1*FtopOut)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y-x13m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V13*((x13m-X12)*(k1-f1*FtopOut)+(x13m-X14)*(k1-f1*FtopIn)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(Y+((x13m-X12)*(k1-f1*FtopOut)+(x13m-X14)*(k1-f1*FtopIn)+(X14-x13m)*(x13m-X12)*Ftop*v)/V13-x13m)+sigma_x*dw1)

    }else if (type == 'model'){
        # model$addSystem(dY13m~dt/V13*((X12-x13m)*(k+f*(FtopOut))+(X14-x13m)*(k+f*(FtopIn))+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        # model$addSystem(dx13m~dt*((Y13m-x13m)*a+((X14-x13m)*(k+f*(FtopOut))+(X14-x13m)*(k+f*(FtopIn))+(X14-x13m)*(x13m-X12)*Ftop*v)/V13)+sigma_x*dw1)
        model$addSystem(dY13m~dt/V13*((x13m-X12)*(k-f*FtopOut)+(x13m-X14)*(k-f*FtopIn)+(X14-x13m)*(x13m-X12)*Ftop*v)+sigma_y*dwY)
        model$addSystem(dx13m~dt*((Y13m-x13m)*a+((x13m-X12)*(k-f*FtopOut)+(x13m-X14)*(k-f*FtopIn)+(X14-x13m)*(x13m-X12)*Ftop*v)/V13)+sigma_x*dw1)
    }else if (type == 'model1'){
        model$addSystem(dY13m~dt*(((x13m-X12)*(k3-f3*(FtopOut))+(x13m-X14)*(k4-f4*(FtopIn)))/V13)+sigma_y*dwY)
        model$addSystem(dx13m~dt*(((x13m-X12)*(k1-f1*(FtopOut))+(x13m-X14)*(k2-f2*(FtopIn)))/V13+(Y13m-x13m)*a)+(sigma_x)*dw1)

    }else if (type == 'model2'){
        model$addSystem(dY13m~dt*(x13m-Y13m)+sigma_y*dwY)
        model$addSystem(dx13m~dt*((Y13m-x13m)*a+((x13m-X12)*(k1-f1*(FtopOut))+(x13m-X14)*(k2-f2*(FtopIn)))/V13)+(sigma_x)*dw1)

    }else if (type == 'model3'){
        model$addSystem(dY13m~dt*(x13m-Y13m)*b+sigma_y*dwY)
        model$addSystem(dx13m~dt*((Y13m-x13m)*a+((X12-x13m)*(k1+f1*(FtopOut))+(X14-x13m)*(k2+f2*(FtopIn)))/V13)+(sigma_x)*dw1)

    }else if (type == 'model4'){
        model$addSystem(dY13m~dt*(x13m-Y13m)*b+sigma_y*dwY)
        model$addSystem(dx13m~dt*((Y13m-x13m)*a+((X12-x13m)*(f1*(FtopOut))+(X14-x13m)*(f2*(FtopIn)))/V13)+(sigma_x)*dw1)

    }

    # Add Inputs
    model$addInput("X12","X14")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=1e-30,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    # Parameters
    # Parameters
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))

    return(model)}