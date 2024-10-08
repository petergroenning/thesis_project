make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X3~x3m)
    # Set observation equation variances
    model$setVariance(X3 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx3m~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'simple'){
        model$addSystem(dx3m~dt/V3*((X2-x3m)*(f1*FbotIn)+(X4-x3m)*(f2*FbotOut))+sigma_x*dw1)
    }
    else if (type == 'linear1'){
        model$addSystem(dx3m~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx3m~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx3m~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y-x3m)+sigma_x*dw1)


    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y+((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut))/V3-x3m)+sigma_x*dw1)


    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y-x3m)+sigma_x*dw1)


    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y+((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut))/V3-x3m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y-x3m)+sigma_x*dw1)


    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y+((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)/V3-x3m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y-x3m)+sigma_x*dw1)


    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V3*((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(Y+((x3m-X2)*(k1-f1*FbotIn)+(x3m-X4)*(k1-f1*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)/V3-x3m)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY3m~dt/V3*((x3m-X2)*(k-f*FbotIn)+(x3m-X4)*(k-f*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+(sigma_y)*dwY)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((x3m-X2)*(k-f*FbotIn)+(x3m-X4)*(k-f*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)/V3)+(sigma_x)*dw1)
    }else if (type == 'model1'){
        model$addSystem(dY3m~dt*(((x3m-X2)*(k3-f3*(FbotIn))+(x3m-X4)*(k4-f4*(FbotOut)))/V3)+sigma_y*dwY)
        model$addSystem(dx3m~dt*(((x3m-X2)*(k1-f1*(FbotIn))+(x3m-X4)*(k2-f2*(FbotOut)))/V3+(Y3m-x3m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY3m~dt*(x3m-Y3m)+sigma_y*dwY)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((x3m-X2)*(k1-f1*(FbotIn))+(x3m-X4)*(k2-f2*(FbotOut)))/V3)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwY)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((X2-x3m)*(k1+f1*(FbotIn))+(X4-x3m)*(k2+f2*(FbotOut)))/V3)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwY)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((X2-x3m)*(f1*(FbotIn))+(X4-x3m)*(f2*(FbotOut)))/V3)+(sigma_x)*dw1)
    }else if (type == 'model5'){
        model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwY)
        model$addSystem(dx3m~dt*((Y3m-x3m)*a+((X2-x3m)*(f1*(FbotIn))+(X4-x3m)*(f2*(FbotOut)))/V3)+(sigma_x)*dw1)
    }else if (type == 'nl') {
         model$addSystem(dx3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nl1') {
         model$addSystem(dx3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+dt*(Y3m-x3m)*a+sigma_x*dw1)
         model$addSystem(dY3m~dt*(x3m-Y3m)*b+sigma_y*dwY)
    } else if (type == 'nl2') {
         model$addSystem(dx3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut)-(X4-x3m)*(x3m-X2)*Fbot*v)+dt*(Y3m-x3m)*a+sigma_x*dw1)
         model$addSystem(dY3m~dt*(x3m-Y3m)+0*dwY)
    } else if (type == 'test1'){
        model$addSystem(dY3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut))+sigma_x*dw1)
         model$addSystem(dx3m~dt*Y3m+sigma_y*dwY)
    }else if (type == 'test2'){
        model$addSystem(dY3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut))+sigma_x*dw1)
         model$addSystem(dx3m~dt*Y3m*a+sigma_y*dwY)
    }else if (type == 'test3'){
        model$addSystem(dY3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut))+sigma_x*dw1)
         model$addSystem(dx3m~dt*(Y3m-x3m)*a+sigma_y*dwY)
    }else if (type == 'test4'){
        model$addSystem(dY3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut))-dt*Y3m*b+sigma_x*dw1)
         model$addSystem(dx3m~dt*(Y3m)*a+sigma_y*dwY)
    }
    
    
    # Parameters
    # model$setParameter(Y3m = c(init=8,lb=0,ub=40))
    # Parameters
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=-10,ub=10))
    model$setParameter(b = c(init=0.1,lb=0,ub=1))
    model$setParameter(v = c(init=1e-4,lb=0,ub=10))


    # Add Inputs
    model$addInput("X2","X4")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-4,lb=0,ub=1))



    return(model)}