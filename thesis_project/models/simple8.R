make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X8~x8m)
    # Set observation equation variances
    model$setVariance(X8 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx8m~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx8m~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx8m~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx8m~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y-x8m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y+((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut))/V8-x8m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y-x8m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y+((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut))/V8-x8m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y-x8m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y+((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k2+f2*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)/V8-x8m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y-x8m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V8*((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*(Y+((x8m-X7)*(k1-f1*FbotIn)+(x8m-X9)*(k1-f1*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)/V8-x8m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY8m~dt/V8*((x8m-X7)*(k-f*FbotIn)+(x8m-X9)*(k-f*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx8m~dt*((Y8m-x8m)*a+((x8m-X7)*(k-f*FbotIn)+(x8m-X9)*(k-f*FbotOut)-(X9-x8m)*(x8m-X7)*Fbot*v)/V8)+sigma_x*dw1) }

    # Add Inputs
    model$addInput("X7","X9")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))

    # Parameters
    model$setParameter(k = c(init=1,lb=0,ub=50))
    model$setParameter(f = c(init=5e-1,lb=0,ub=10))
    model$setParameter(v = c(init=1e-4,lb=0,ub=1))
    model$setParameter(a = c(init=0.5,lb=0,ub=2))
    return(model)}