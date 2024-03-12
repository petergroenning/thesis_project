library(ctsmr)

make_model<-function(model){
    # Add observation equations and variances
    model$addSystem(dx1m~dt/V1*((x2m-x1m)*A2*k/Cp+(Uin-Uout)*f*(x2m-x1m)+Uout*(Tdischarge-x1m)*vbot)+sigma_x*dw1)
    model$addSystem(dx2m~dt/V2*((x1m-x2m)*A2*k/Cp+(x3m-x2m)*A3*k/Cp+(Uin-Uout)*f*(x3m-x1m)/2)+sigma_x*dw2)
    model$addSystem(dx3m~dt/V3*((x2m-x3m)*A3*k/Cp+(x4m-x3m)*A4*k/Cp+(Uin-Uout)*f*(x4m-x2m)/2)+sigma_x*dw3)
    model$addSystem(dx4m~dt/V4*((x3m-x4m)*A4*k/Cp+(x5m-x4m)*A5*k/Cp+(Uin-Uout)*f*(x5m-x3m)/2)+sigma_x*dw4)
    model$addSystem(dx5m~dt/V5*((x4m-x5m)*A5*k/Cp+(x6m-x5m)*A6*k/Cp+(Uin-Uout)*f*(x6m-x4m)/2)+sigma_x*dw5)
    model$addSystem(dx6m~dt/V6*((x5m-x6m)*A6*k/Cp+(x7m-x6m)*A7*k/Cp+(Uin-Uout)*f*(x7m-x5m)/2)+sigma_x*dw6)
    model$addSystem(dx7m~dt/V7*((x6m-x7m)*A7*k/Cp+(x8m-x7m)*A8*k/Cp+(Uin-Uout)*f*(x8m-x6m)/2)+sigma_x*dw7)
    model$addSystem(dx8m~dt/V8*((x7m-x8m)*A8*k/Cp+(x9m-x8m)*A9*k/Cp+(Uin-Uout)*f*(x9m-x7m)/2)+sigma_x*dw8)
    model$addSystem(dx9m~dt/V9*((x8m-x9m)*A9*k/Cp+(x10m-x9m)*A10*k/Cp+(Uin-Uout)*f*(x10m-x8m)/2)+sigma_x*dw9)
    model$addSystem(dx10m~dt/V10*((x9m-x10m)*A10*k/Cp+(x11m-x10m)*A11*k/Cp+(Uin-Uout)*f*(x11m-x9m)/2)+sigma_x*dw10)
    model$addSystem(dx11m~dt/V11*((x10m-x11m)*A11*k/Cp+(x12m-x11m)*A12*k/Cp+(Uin-Uout)*f*(x12m-x10m)/2)+sigma_x*dw11)
    model$addSystem(dx12m~dt/V12*((x11m-x12m)*A12*k/Cp+(x13m-x12m)*A13*k/Cp+(Uin-Uout)*f*(x13m-x11m)/2)+sigma_x*dw12)
    model$addSystem(dx13m~dt/V13*((x12m-x13m)*A13*k/Cp+(x14m-x13m)*A14*k/Cp+(Uin-Uout)*f*(x14m-x12m)/2)+sigma_x*dw13)
    model$addSystem(dx14m~dt/V14*((x13m-x14m)*A14*k/Cp+(Uin-Uout)*f*(x14m-x13m)+Uin*(Tcharge-x14m)*vtop)+sigma_x*dw14)

    #SystemNoise
    model$setParameter(sigma_x=c(init=1,lb=0,ub=10))

    #SystemParameters
    model$setParameter(vbot=c(init=1.4323e-1,lb=0,ub=1e6))
    model$setParameter(vtop=c(init=1.7673,lb=0,ub=1e6))
    model$setParameter(k = c(init = 6.5176, lb = 0, ub = 1e6))
    model$setParameter(f = c(init = 4e-6 * 3751, lb = 0, ub = 1e6))

    return(model)
}