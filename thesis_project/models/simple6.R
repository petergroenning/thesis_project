make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X6~x6m)
    # Set observation equation variances
    model$setVariance(X6 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx6m~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut))/V6-x6m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut))/V6-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k2+f2*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)/V6-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y-x6m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V6*((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(Y+((x6m-X5)*(k1-f1*FbotIn)+(x6m-X7)*(k1-f1*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)/V6-x6m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY6m~dt/V6*((x6m-X5)*(k-f*FbotIn)+(x6m-X7)*(k-f*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx6m~dt*((Y6m-x6m)*a+((x6m-X5)*(k-f*FbotIn)+(x6m-X7)*(k-f*FbotOut)-(X7-x6m)*(x6m-X5)*Fbot*v)/V6)+sigma_x*dw1)
    }else if (type == 'model1'){
        model$addSystem(dY6m~dt*(((x6m-X5)*(k3-f3*(FbotIn))+(x6m-X7)*(k4-f4*(FbotOut)))/V6)+sigma_y*dwY)
        model$addSystem(dx6m~dt*(((x6m-X5)*(k1-f1*(FbotIn))+(x6m-X7)*(k2-f2*(FbotOut)))/V6+(Y6m-x6m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY6m~dt*(x6m-Y6m)+sigma_y*dwY)
        model$addSystem(dx6m~dt*((Y6m-x6m)*a+((x6m-X5)*(k1-f1*(FbotIn))+(x6m-X7)*(k2-f2*(FbotOut)))/V6)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY6m~dt*(x6m-Y6m)*b+sigma_y*dwY)
        model$addSystem(dx6m~dt*((Y6m-x6m)*a+((X5-x6m)*(k1+f1*(FbotIn))+(X7-x6m)*(k2+f2*(FbotOut)))/V6)+(sigma_x)*dw1)
    } else if (type == 'model4'){
        model$addSystem(dY6m~dt*(x6m-Y6m)*b+sigma_y*dwY)
        model$addSystem(dx6m~dt*((Y6m-x6m)*a+((X5-x6m)*(f1*(FbotIn))+(X7-x6m)*(f2*(FbotOut)))/V6)+(sigma_x)*dw1)
    } 
    # Parameters
    # model$setParameter(Y6m = c(init=8,lb=0,ub=40))
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))
    # Add Inputs
    model$addInput("X5","X7")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))

    return(model)}