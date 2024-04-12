make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X4~x4m)
    # Set observation equation variances
    model$setVariance(X4 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx4m~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx4m~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx4m~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx4m~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y-x4m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y+((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut))/V4-x4m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y-x4m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y+((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut))/V4-x4m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y-x4m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y+((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k2+f2*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)/V4-x4m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y-x4m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V4*((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(Y+((x4m-X3)*(k1-f1*FbotIn)+(x4m-X5)*(k1-f1*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)/V4-x4m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY4m~dt/V4*((x4m-X3)*(k-f*FbotIn)+(x4m-X5)*(k-f*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((x4m-X3)*(k-f*FbotIn)+(x4m-X5)*(k-f*FbotOut)-(X5-x4m)*(x4m-X3)*Fbot*v)/V4)+sigma_x*dw1)
    }else if (type == 'model1'){
        model$addSystem(dY4m~dt*(((x4m-X3)*(k3-f3*(FbotIn))+(x4m-X5)*(k4-f4*(FbotOut)))/V4)+sigma_y*dwY)
        model$addSystem(dx4m~dt*(((x4m-X3)*(k1-f1*(FbotIn))+(x4m-X5)*(k2-f2*(FbotOut)))/V4+(Y4m-x4m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY4m~dt*(x4m-Y4m)+sigma_y*dwY)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((x4m-X3)*(k1-f1*(FbotIn))+(x4m-X5)*(k2-f2*(FbotOut)))/V4)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY4m~dt*(x4m-Y4m)*b+sigma_y*dwY)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((X3-x4m)*(k1+f1*(FbotIn))+(X5-x4m)*(k2+f2*(FbotOut)))/V4)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY4m~dt*(x4m-Y4m)*b+sigma_y*dwY)
        model$addSystem(dx4m~dt*((Y4m-x4m)*a+((X3-x4m)*(f1*(FbotIn))+(X5-x4m)*(f2*(FbotOut)))/V4)+(sigma_x)*dw1)
    }
    
    # Parameters
    # model$setParameter(Y4m = c(init=8,lb=0,ub=40))
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))

    # Add Inputs
    model$addInput("X3","X5")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))


    return(model)}