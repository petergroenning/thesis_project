make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X9~x9m)
    # Set observation equation variances
    model$setVariance(X9 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx9m~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx9m~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx9m~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx9m~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y-x9m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y+((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut))/V9-x9m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y-x9m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y+((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut))/V9-x9m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y-x9m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y+((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k2+f2*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)/V9-x9m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*FbotIn)+(x9m-X10)*(k1-f1*FbotOut)-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y-x9m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V9*((x9m-X8)*(k1-f1*(FbotIn))+(x9m-X10)*(k1-f1*(FbotOut))-(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx9m~dt*(Y+((x9m-X8)*(k1-f1*(FbotIn))+(x9m-X10)*(k1-f1*(FbotOut))-(X10-x9m)*(x9m-X8)*Fbot*v)/V9-x9m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY9m~dt/V9*((x9m-X8)*(k-f*FbotIn)+(x9m-X10)*(k-f*(FbotOut))-(X10-x9m)*(x9m-X8)*Fbot*(v))+sigma_y*dwY9)
        model$addSystem(dx9m~dt*((Y9m-x9m)*a+((x9m-X8)*(k-f*FbotIn)+(x9m-X10)*(k-f*(FbotOut))-(X10-x9m)*(x9m-X8)*Fbot*(v))/V9)+sigma_x*dw9)
    }

    # Add Inputs
    model$addInput("X8","X10","FmidVol")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-3,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters
    model$setParameter(k = c(init=10,lb=0,ub=50))
    model$setParameter(f = c(init=1,lb=0,ub=10))
    model$setParameter(v = c(init=1e-4,lb=0,ub=1))
    model$setParameter(a = c(init=1,lb=0,ub=2))
    return(model)}