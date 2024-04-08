make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X6~x6m)
    # Set observation equation variances
    model$setVariance(X6 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))/V6-x6m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut))/V6-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)/V6-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1+f1*FbotIn)+(x6m-X7)*(k1+f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)/V6-x6m)+sigma_x*dw1)

    }

    # Add Inputs
    model$addInput("X5","X7")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=50))

    # Parameters
    model$setParameter(k1 = c(init=-20,lb=-30,ub=30))
    model$setParameter(k2 = c(init=2,lb=-30,ub=30))
    model$setParameter(f1 = c(init=5e-1,lb=-4,ub=4))
    model$setParameter(f2 = c(init=-5e-1,lb=-2,ub=2))
    model$setParameter(v = c(init=1e-1,lb=-2,ub=2))
    return(model)}