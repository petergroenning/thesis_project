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
        model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
        model$addSystem(dx0m~dt*((Y-x0m)+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))/V0)+sigma_x*dw1)
    } 
    # else if (type == 'nonlinear1_lag1'){
    #     model$addSystem(dY~dt/V0*((x0m-X1)*(k1+f1*(FbotOut)+f2*(x0m)*FbotOut)+(Tbot-x0m)*FbotIn*(f1+f2*(x0m-X1)))+sigma_y*dwY)
    #     model$addSystem(dx0m~dt/V0*(Y)+sigma_x*dw1)
    # } else if (type == 'nonlinear1_lag2'){
    #     model$addSystem(dY~dt/V0*((x0m-X1)*(k1+f1*(FbotOut)+f2*(x0m)*FbotOut)+(Tbot-x0m)*FbotIn*(f1+f2*(x0m-X1)))+sigma_y*dwY)
    #     model$addSystem(dx0m~dt/V0*(Y-x0m)+sigma_x*dw1)
    # } else if (type == 'nonlinear1_lag3'){

    #     model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+k2*(x0m-X1)+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
    #      model$addSystem(dx0m~dt*((Y))+sigma_x*dw1)
    # } else if (type == 'nonlinear1_lag4'){
    #     model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+k2*(x0m-X1)+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
    #     model$addSystem(dx0m~dt*((Y-x0m)+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+k2*(x0m-X1)+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))/V0)+sigma_x*dw1)
    # }else if (type == 'nonlinear1_lag5'){
    #     model$addSystem(dY~dt/V0*((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
    #     model$addSystem(dx0m~dt*((Y-x0m)+((Tbot-x0m)*(v1+v2*(X1-x0m))*(FbotIn)+(X1-x0m)*(k1+f1*(FbotOut)+f2*(Fbot)*(X1-x0m)))/V0)+sigma_x*dw1)
    # }
    
    
  
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=20))

    model$addInput("X1",'FbotIn','FbotOut', 'Tbot', 'Fbot')
    #####SETPARAMETERS########
    #SystemNoise
    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters

    # model$setParameter(f1=c(init=1e-6,lb=-3,ub=10))
    # model$setParameter(f2=c(init=1e-6,lb=-3,ub=10))
    # model$setParameter(v1 = c(init=1e-1,lb=-2,ub=10))
    # model$setParameter(v2 = c(init=1e-1,lb=-2,ub=10))
    # model$setParameter(k1=c(init=-1.2,lb=-40,ub=40))
    # model$setParameter(k2=c(init=-1.2,lb=-40,ub=40))



    # model$setParameter(sigma_x=c(init=2e-2,lb=0,ub=10))
    model$setParameter(f1 = c(init=5e-1,lb=0,ub=2))
    model$setParameter(f2 = c(init=5e-1,lb=0,ub=2))
    model$setParameter(v1 = c(init=1e-1,lb=1e-2,ub=2))
    model$setParameter(v2 = c(init=1e-1,lb=1e-2,ub=2))

    model$setParameter(f1=c(init=1e-1,lb=-2,ub=1))
    model$setParameter(f2=c(init=1e-3,lb=-1,ub=1))

    model$setParameter(k1=c(init=-1.2,lb=-10,ub=1))
    model$setParameter(k2=c(init=-1.2,lb=-10,ub=1))

    return(model)}