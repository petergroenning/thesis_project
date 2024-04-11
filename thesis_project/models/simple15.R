make_model<-function(model, type = 'linear'){
    # Add observation equations   
    model$addObs(X15~x15m)
    # Set observation equation variances
    model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    if (type == 'linear'){
        model$addSystem(dx15m~dt/V15*((Ttop-x15m)*(v1)*FtopIn+(X14-x15m)*(k1+FtopOut*f1)+(ambientTemp-x15m)*u)+sigma_x*dw15)

    } else if (type == 'nonlinear1'){
        model$addSystem(dx15m~dt/V15*((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1)+(ambientTemp-x15m)*u)+sigma_x*dw15)

    } else if (type == 'nonlinear2'){
        model$addSystem(dx15m~dt/V15*((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1+f2*Ftop*(x15m-X14))+(ambientTemp-x15m)*u)+sigma_x*dw15)


    } else if (type == 'linear_lag1'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1)*FtopIn+(X14-x15m)*(k1+FtopOut*f1)+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*Y+sigma_x*dw15)
    } else if (type == 'linear_lag2'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1)*FtopIn+(X14-x15m)*(k1+FtopOut*f1))+sigma_y*dwY)
        model$addSystem(dx15m~dt*(Y-x15m)+sigma_x*dw15)
    }else if (type == 'linear_lag3'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1)*FtopIn+(X14-x15m)*(k1+FtopOut*f1)+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*(Y-x15m+((Ttop-x15m)*(v1)*FtopIn+(X14-x15m)*(k1+FtopOut*f1)+(ambientTemp-x15m)*u)/V15)+sigma_x*dw15)

    }else if (type == 'nonlinear2_lag1'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1+f2*Ftop*(x15m-X14))+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*Y+sigma_x*dw15)
    } else if (type == 'nonlinear2_lag2'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1+f2*Ftop*(x15m-X14))+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*(Y-x15m)+sigma_x*dw15)
    }else if (type == 'nonlinear2_lag3'){
        model$addSystem(dY~dt/V15*((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1+f2*Ftop*(x15m-X14))+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*(Y-x15m+((Ttop-x15m)*(v1+v2*(X14-x15m))*FtopIn+(X14-x15m)*(k1+FtopOut*f1+f2*Ftop*(x15m-X14))+(ambientTemp-x15m)*u)/V15)+sigma_x*dw15)
    } else if (type == 'model'){
        model$addSystem(dY15m~dt/V15*((Ttop-x15m)*(v1+v2*(x15m-X14))*FtopIn+(X14-x15m)*(k+f*FtopOut)+(ambientTemp-x15m)*u)+sigma_y*dwY)
        model$addSystem(dx15m~dt*((Y15m-x15m)*a+((Ttop-x15m)*(v1+v2*(x15m-X14))*FtopIn+(X14-x15m)*(k+f*FtopOut)+(ambientTemp-x15m)*u)/V15)+sigma_x*dw15)

    }
    
    
  
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=20))

    model$addInput("X14",'FtopIn','FtopOut', 'Tbot', 'Ftop')
    #####SETPARAMETERS########
    #SystemNoise
    #SystemNoise
    model$setParameter(sigma_x = c(init=2e-3,lb=0,ub=1))

    # Hidden State
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=40))
    # Parameters
    model$setParameter(k = c(init=10,lb=0,ub=50))
    model$setParameter(f = c(init=1,lb=0,ub=10))
    model$setParameter(v = c(init=1e-4,lb=0,ub=1))
    model$setParameter(a = c(init=1,lb=0,ub=2))

    model$setParameter(v1 = c(init=0,lb=-4,ub=5))
    model$setParameter(v2 = c(init=0,lb=-4,ub=5))
    model$setParameter(u = c(init=1,lb=0,ub=5))
  

    return(model)}