make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X2~x2m)
    # Set observation equation variances
    model$setVariance(X2 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx2m~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx2m~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx2m~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx2m~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y-x2m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y+((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut))/V2-x2m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y-x2m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y+((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut))/V2-x2m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y-x2m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y+((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k2+f2*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)/V2-x2m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y-x2m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y+((x2m-X1)*(k1+f1*FbotIn)+(x2m-X3)*(k1+f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)/V2-x2m)+sigma_x*dw1)
 
    }

    # Add Inputs
    model$addInput("X1","X3")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-3,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))

    # Parameters
    model$setParameter(k1 = c(init=-20,lb=-50,ub=50))
    model$setParameter(k2 = c(init=-20,lb=-50,ub=50))
    model$setParameter(f1 = c(init=-5e-1,lb=-2,ub=2))
    model$setParameter(f2 = c(init=5e-1,lb=-2,ub=2))
    model$setParameter(v = c(init=1e-1,lb=-1,ub=1))
    return(model)}