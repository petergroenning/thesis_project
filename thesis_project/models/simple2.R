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
        model$addSystem(dY~dt/V2*((x2m-X1)*(k1-f1*FbotIn)+(x2m-X3)*(k1-f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(Y+((x2m-X1)*(k1-f1*FbotIn)+(x2m-X3)*(k1-f1*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)/V2-x2m)+sigma_x*dw1)
 
    } else if (type == 'model'){
        model$addSystem(dY2m~dt/V2*((x2m-X1)*(k-f*FbotIn)+(x2m-X3)*(k-f*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((x2m-X1)*(k-f*FbotIn)+(x2m-X3)*(k-f*FbotOut)-(X3-x2m)*(x2m-X1)*Fbot*v)/V2)+sigma_x*dw1)

    }else if (type == 'model1'){
        model$addSystem(dY2m~dt*(((x2m-X1)*(k3-f3*(FbotIn))+(x2m-X3)*(k4-f4*(FbotOut)))/V2)+sigma_y*dwY)
        model$addSystem(dx2m~dt*(((x2m-X1)*(k1-f1*(FbotIn))+(x2m-X3)*(k2-f2*(FbotOut)))/V2+(Y2m-x2m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY2m~dt*(x2m-Y2m)+sigma_y*dwY)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((x2m-X1)*(k1-f1*(FbotIn))+(x2m-X3)*(k2-f2*(FbotOut)))/V2)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY2m~dt*(x2m-Y2m)*b+sigma_y*dwY)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((X1-x2m)*(k1+f1*(FbotIn))+(X3-x2m)*(k2+f2*(FbotOut)))/V2)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY2m~dt*(x2m-Y2m)*b+sigma_y*dwY)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((X1-x2m)*(f1*(FbotIn))+(X3-x2m)*(f2*(FbotOut)))/V2)+(sigma_x)*dw1)
    }else if (type == 'model5'){
        model$addSystem(dY2m~dt*((x2m-Y2m)*b+(X1-Y2m)*k1+(X3-Y2m)*k2)+sigma_y*dwY)
        model$addSystem(dx2m~dt*((Y2m-x2m)*a+((X1-x2m)*(f1*(FbotIn))+(X3-x2m)*(f2*(FbotOut)))/V2)+(sigma_x)*dw1)
    }
    
    
    # Parameters
    # Parameters
    model$setParameter(k1 = c(init=10,lb=0,ub=50))
    model$setParameter(k2 = c(init=10,lb=0,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))

    # Add Inputs
    model$addInput("X1","X3")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-5,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    




    return(model)}