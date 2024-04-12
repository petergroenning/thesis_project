make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X7~x7m)
    # Set observation equation variances
    model$setVariance(X7 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx7m~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut))/V7-x7m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut))/V7-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k2+f2*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)/V7-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y-x7m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V7*((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(Y+((x7m-X6)*(k1-f1*FbotIn)+(x7m-X8)*(k1-f1*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)/V7-x7m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY7m~dt/V7*((x7m-X6)*(k-f*FbotIn)+(x7m-X8)*(k-f*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx7m~dt*((Y7m-x7m)*a+((x7m-X6)*(k-f*FbotIn)+(x7m-X8)*(k-f*FbotOut)-(X8-x7m)*(x7m-X6)*Fbot*v)/V7)+sigma_x*dw1)

    }else if (type == 'model1'){
        model$addSystem(dY7m~dt*(((x7m-X6)*(k3-f3*(FbotIn))+(x7m-X8)*(k4-f4*(FbotOut)))/V7)+sigma_y*dwY)
        model$addSystem(dx7m~dt*(((x7m-X6)*(k1-f1*(FbotIn))+(x7m-X8)*(k2-f2*(FbotOut)))/V7+(Y7m-x7m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY7m~dt*(x7m-Y7m)+sigma_y*dwY)
        model$addSystem(dx7m~dt*((Y7m-x7m)*a+((x7m-X6)*(k1-f1*(FbotIn))+(x7m-X8)*(k2-f2*(FbotOut)))/V7)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY7m~dt*(x7m-Y7m)*b+sigma_y*dwY)
        model$addSystem(dx7m~dt*((Y7m-x7m)*a+((X6-x7m)*(k1+f1*(FbotIn))+(X8-x7m)*(k2+f2*(FbotOut)))/V7)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY7m~dt*(x7m-Y7m)*b+sigma_y*dwY)
        model$addSystem(dx7m~dt*((Y7m-x7m)*a+((X6-x7m)*(f1*(FbotIn))+(X8-x7m)*(f2*(FbotOut)))/V7)+(sigma_x)*dw1)
    }
    
    
    # Parameters
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))
    # Add Inputs
    model$addInput("X6","X8")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))

    # Parameters

    return(model)}