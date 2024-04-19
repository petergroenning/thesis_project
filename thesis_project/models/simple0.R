make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X0~x0m)
    # Set observation equation variances
    model$setVariance(X0 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    if (type == 'linear'){
        model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(v1)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))+sigma_x*dw0)

    } else if (type == 'nonlinear'){
        model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_x*dw0)
        
    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))+sigma_y*dw0)
        model$addSystem(dx0m~dt/V0*(Y)+sigma_x*dw1)

    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))+sigma_y*dw0)
        model$addSystem(dx0m~dt/V0*(Y-x0m)+sigma_x*dw1)

    } else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))+sigma_y*dw0)
        model$addSystem(dx0m~dt/V0*(Y-x0m+((Tbot-x0m)*(v1)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))/V0)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag1'){
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
        model$addSystem(dx0m~dt/V0*(Y)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag2'){
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
        model$addSystem(dx0m~dt/V0*(Y-x0m)+sigma_x*dw1)

    } else if (type == 'nonlinear_lag3'){
        # model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
        # model$addSystem(dx0m~dt/V0*(Y-x0m+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))/V0)+sigma_x*dw1)
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(x0m-X1)*(k1-f1*(FbotOut)))+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y-x0m)+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(x0m-X1)*(k1-f1*(FbotOut)))/V0)+sigma_x*dw1)
    } else if (type == 'model'){
        model$addSystem(dY0m~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(x0m-X1)*(k-f*(FbotOut)))+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(x0m-X1)*(k-f*(FbotOut)))/V0)+sigma_x*dw1)  
    } else if (type == 'model2'){
        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*v*(FbotIn)+(x0m-X1)*(k1-f1*(FbotOut)))/V0)+sigma_x*dw1)  
    }else if (type == 'model3'){
        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*(v)*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)))/V0)+sigma_x*dw1)  
    }else if (type == 'model4'){
        model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*(v)*(FbotIn)+(X1-x0m)*(f1*(FbotOut)))/V0)+sigma_x*dw1)  
    }else if (type == 'model5'){
         model$addSystem(dY0m~dt*(x0m-Y0m)*b+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*(v)*(FbotIn)+(X1-x0m)*(f1*(FbotOut)+f2*FbotIn))/V0)+sigma_x*dw1)   
    }else if (type == 'model6'){
        model$addSystem(dY0m~dt*(x0m-Y0m)*b-dt*u*Y0m+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y0m-x0m)*a+((Tbot-x0m)*(v)*(FbotIn)+(X1-x0m)*(f1*(FbotOut)))/V0)+sigma_x*dw1)  
    }
  
    
  
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=20))

    model$addInput("X1",'FbotIn','FbotOut', 'Tbot', 'Fbot')
    #####SETPARAMETERS########
    #SystemNoise
    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    # model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters

    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-5,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))

    model$setParameter(k1 = c(init=500,lb=0,ub=10000))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=100))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=100))
    model$setParameter(v = c(init=1,lb=0,ub=100))
    model$setParameter(a = c(init=1,lb=0,ub=4))
    model$setParameter(b= c(init=1,lb=0,ub=4))
    model$setParameter(u = c(init=1e-6,lb=0,ub=40))
    # Parameters    
    # model$setParameter(k = c(init=1,lb=0,ub=50))
    # model$setParameter(f = c(init=5e-1,lb=0,ub=10))
    # model$setParameter(v = c(init=1e-4,lb=0,ub=10))
    # model$setParameter(a = c(init=0.5,lb=0,ub=1))
    # # model$setParameter(f1=c(init=1e-6,lb=-3,ub=10))
    # # model$setParameter(f2=c(init=1e-6,lb=-3,ub=10))
    # model$setParameter(v1 = c(init=1,lb=-2,ub=10))
    # model$setParameter(v2 = c(init=1,lb=-2,ub=10))
    # model$setParameter(k1=c(init=-1.2,lb=-40,ub=40))
    # model$setParameter(k2=c(init=-1.2,lb=-40,ub=40))



    # model$setParameter(sigma_x=c(init=2e-2,lb=0,ub=10))
    # model$setParameter(f1 = c(init=5e-1,lb=0,ub=2))
    # model$setParameter(f2 = c(init=5e-1,lb=0,ub=2))
    # model$setParameter(v1 = c(init=1e-1,lb=1e-2,ub=2))
    # model$setParameter(v2 = c(init=1e-1,lb=1e-2,ub=2))

    # model$setParameter(f1=c(init=1e-1,lb=-2,ub=1))
    # model$setParameter(v=c(init=1e-3,lb=-1,ub=1))

    # model$setParameter(k1=c(init=-1.2,lb=-10,ub=1))
    # model$setParameter(k2=c(init=-1.2,lb=-10,ub=1))

    return(model)}