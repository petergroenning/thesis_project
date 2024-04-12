make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X1~x1m)
    # Set observation equation variances
    model$setVariance(X1 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    # MODEL NO DELAY 
    if (type == 'linear'){
        model$addSystem(dx1m~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut))+sigma_x*dw1)

    } else if (type == 'linear1'){
        model$addSystem(dx1m~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut))+sigma_x*dw1)

    } else if (type == 'nonlinear') {
         model$addSystem(dx1m~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'nonlinear1') {
        model$addSystem(dx1m~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_x*dw1)

    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y)+sigma_x*dw1)

    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y-x1m)+sigma_x*dw1)

    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut))/V1-x1m)+sigma_x*dw1)

    } else if (type == 'linear1_lag1'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y)+sigma_x*dw1)

    } else if (type == 'linear1_lag2'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y-x1m)+sigma_x*dw1)

    } else if (type == 'linear1_lag3'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut))+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut))/V1-x1m)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y-x1m)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag3'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1-x1m)+sigma_x*dw1)

    } else if (type == 'nonlinear1_lag1'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y)+sigma_x*dw1)

    } else if (type == 'nonlinear1_lag2'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y-x1m)+sigma_x*dw1)

    } else if (type == 'nonlinear1_lag3'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1-f1*FbotIn)+(x1m-X2)*(k1-f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1-f1*FbotIn)+(x1m-X2)*(k1-f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1-x1m)+sigma_x*dw1)

    } else if (type == 'nonlinear1_lag4'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1)+sigma_x*dw1)

    } else if (type == 'nonlinear1_lag5'){
        model$addSystem(dY~dt/V1*((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(Y+((x1m-X0)*(k1+f1*FbotIn)+(x1m-X2)*(k1+f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1-u*x1m)+sigma_x*dw1)

    # } 
    # else if (type == 'model2'){
    #     model$addSystem(dY1m~dt/V1*((x1m-X0)*(k1-f1*FbotIn)+(x1m-X2)*(k1-f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
    #     model$addSystem(dx1m~dt*(Y1m-x1m+((x1m-X0)*(k1-f1*FbotIn)+(x1m-X2)*(k1-f1*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1)+sigma_x*dw1)

    } else if (type == 'model'){
        model$addSystem(dY1m~dt/V1*((x1m-X0)*(k-f*FbotIn)+(x1m-X2)*(k-f*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((x1m-X0)*(k-f*FbotIn)+(x1m-X2)*(k-f*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)/V1)+sigma_x*dw1)

    }else if (type == 'model1'){
        model$addSystem(dY1m~dt*(((x1m-X0)*(k3-f3*(FbotIn))+(x1m-X2)*(k4-f4*(FbotOut)))/V1)+sigma_y*dwY)
        model$addSystem(dx1m~dt*(((x1m-X0)*(k1-f1*(FbotIn))+(x1m-X2)*(k2-f2*(FbotOut)))/V1+(Y1m-x1m)*a)+(sigma_x)*dw1)
    }else if (type == 'model2'){
        model$addSystem(dY1m~dt*(x1m-Y1m)+sigma_y*dwY)
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((x1m-X0)*(k1-f1*(FbotIn))+(x1m-X2)*(k2-f2*(FbotOut)))/V1)+(sigma_x)*dw1)
    }else if (type == 'model3'){
        model$addSystem(dY1m~dt*(x1m-Y1m)*b+sigma_y*dwY)
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((X0-x1m)*(k1+f1*(FbotIn))+(X2-x1m)*(k2+f2*(FbotOut)))/V1)+(sigma_x)*dw1)
    }else if (type == 'model4'){
        model$addSystem(dY1m~dt*(x1m-Y1m)*b+sigma_y*dwY)
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((X0-x1m)*(f1*(FbotIn))+(X2-x1m)*(f2*(FbotOut)))/V1)+(sigma_x)*dw1)
    }else if (type == 'model5'){
        model$addSystem(dY1m~dt*(x1m-Y1m)*b-dt*u*Y1m+sigma_y*dwY)
        model$addSystem(dx1m~dt*((Y1m-x1m)*a+((X0-x1m)*(f1*(FbotIn))+(X2-x1m)*(f2*(FbotOut)))/V1)+(sigma_x)*dw1)
    }



    # Add Inputs
    model$addInput("X0","X2")

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-5,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    # Parameters
    model$setParameter(k1 = c(init=1,lb=-50,ub=50))
    model$setParameter(k2 = c(init=1,lb=-50,ub=50))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=50))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=50))

    model$setParameter(a = c(init=1,lb=0,ub=10))
    model$setParameter(b = c(init=1e-4,lb=0,ub=10))
    model$setParameter(u = c(init=1e-6,lb=0,ub=40))

    return(model)}