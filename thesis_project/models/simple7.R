make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X7~x7m)
    # Set observation equation variances
    model$setVariance(X7 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))/V7-x7m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut))/V7-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)/V7-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1+f1*FbotIn)+(x7m-X8)*(k1+f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)/V7-x7m)+sigma_x*dw1)

    }

    # Add Inputs
    model$addInput("X6","X8")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))

    # Parameters
    model$setParameter(k1 = c(init=-20,lb=-30,ub=30))
    model$setParameter(k2 = c(init=2,lb=-30,ub=30))
    model$setParameter(f1 = c(init=5e-1,lb=-4,ub=4))
    model$setParameter(f2 = c(init=-5e-1,lb=-2,ub=2))
    model$setParameter(v = c(init=1e-1,lb=-2,ub=2))
    return(model)}