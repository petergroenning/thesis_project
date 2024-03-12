make_model <- function(model){
     # Add system equations

    model$addSystem(dx0m ~ dt/V0*((x1m-x0m)*(A1*k1/Cp+FbotOut*f1)+(Tbot-x0m)*FbotIn*vbot)+sigma_x*dw0)
    model$addSystem(dx1m ~ dt/V1*((x0m-x1m)*(A1*k1/Cp+FbotIn*f1)+(x2m-x1m)*(A2*k1/Cp+FbotOut*f1)+(Tbot-x1m)*(FbotIn)*vbot1)+sigma_x*dw1)
    model$addSystem(dx2m ~ dt/V2*((x1m-x2m)*(A2*k1/Cp+FbotIn*f1)+(x3m-x2m)*(A3*k1/Cp+FbotOut*f1)+(Tbot-x2m)*(FbotIn)*vbot1)+sigma_x*dw2)
    model$addSystem(dx3m ~ dt/V3*((x2m-x3m)*(A3*k1/Cp+FbotIn*f1)+(x4m-x3m)*(A4*k1/Cp+FbotOut*f1)+(Tbot-x3m)*(FbotIn)*vbot1)+sigma_x*dw3)
    model$addSystem(dx4m ~ dt/V4*((x3m-x4m)*(A4*k1/Cp+FbotIn*f1)+(x5m-x4m)*(A5*k1/Cp+FbotOut*f1)+(Tbot-x4m)*(FbotIn)*vbot1)+sigma_x*dw4)
    model$addSystem(dx5m ~ dt/V5*((x4m-x5m)*(A5*k1/Cp+FbotIn*f1)+(x6m-x5m)*(A6*k1/Cp+FbotOut*f1)+(Tbot-x5m)*(FbotIn)*vbot1)+sigma_x*dw5)
    model$addSystem(dx6m ~ dt/V6*((x5m-x6m)*(A6*k1/Cp+FbotIn*f1)+(x7m-x6m)*(A7*k1/Cp+FbotOut*f1))+sigma_x*dw6)
    model$addSystem(dx7m ~ dt/V7*((x6m-x7m)*(A7*k1/Cp+FbotIn*f1)+(x8m-x7m)*(A8*k1/Cp+FbotOut*f1))+sigma_x*dw7)
    model$addSystem(dx8m ~ dt/V8*((x7m-x8m)*(A8*k1/Cp+FbotIn*f1)+(x9m-x8m)*(A9*k1/Cp+FbotOut*f1))+sigma_x*dw8)
    model$addSystem(dx9m ~ dt/V9*((x8m-x9m)*(A9*k1/Cp+FbotIn*f1)+(x10m-x9m)*(A10*k1/Cp+FbotOut*f1)+(Tmid-x9m)*FmidIn*vmid9)+sigma_x*dw9)
    model$addSystem(dx10m ~ dt/V10*((x9m-x10m)*(A10*k1/Cp+FtopOut*f2)+(x11m-x10m)*(A11*k1/Cp+FtopIn*f2)+(Tmid-x10m)*FmidIn*vmid10)+sigma_x*dw10)
    model$addSystem(dx11m ~ dt/V11*((x10m-x11m)*(A11*k1/Cp+FtopOut*f2)+(x12m-x11m)*(A12*k1/Cp+FtopIn*f2)+(Tmid-x11m)*FmidIn*vmid11)+sigma_x*dw11)
    model$addSystem(dx12m ~ dt/V12*((x11m-x12m)*(A12*k1/Cp+FtopOut*f2)+(x13m-x12m)*(A13*k1/Cp+FtopIn*f2)+(Tmid-x12m)*FmidIn*vmid12)+sigma_x*dw12)
    model$addSystem(dx13m ~ dt/V13*((x12m-x13m)*(A13*k1/Cp+FtopOut*f2)+(x14m-x13m)*(A14*k1/Cp+FtopIn*f2))+sigma_x*dw13)
    model$addSystem(dx14m ~ dt/V14*((x13m-x14m)*(A14*k1/Cp+FtopOut*f2)+(x15m-x14m)*(A15*k1/Cp+FtopIn*f2))+sigma_x*dw14)
    model$addSystem(dx15m ~ dt/V15*((x14m-x15m)*(A15*k1/Cp+FtopOut*f3)+(ambientTemp-x15m)*Atop*k2/Cp+(Ttop-x15m)*FtopIn*vtop) + sigma_x*dw15)
   
   
       
    ###### SET PARAMETERS ########
    # System Noise
    model$setParameter(sigma_x = c(init = 1, lb = 0, ub = 10))
    model$setParameter(k1 = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(k2 = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(k3 = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(vbot = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(vbot1 = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(vmid = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(vtop = c(init = 1, lb = 0, ub = 1000))
    model$setParameter(f1 = c(init = 1.8483e-04, lb =0, ub = 1000))
    model$setParameter(f2 = c(init = 1.8483e-04, lb =0, ub = 1000))
    model$setParameter(f3 = c(init = 1.8483e-04, lb =0, ub = 1000))

    model$setParameter(vmid8 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid9 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid10 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid11= c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid12 = c(init=1, lb = 0, ub = 1000))

    return(model)
}
