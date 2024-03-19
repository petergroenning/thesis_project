make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*A1*k/Cp+FbotOut*(x1m-x0m)+(Tbot-x0m)*FbotIn*v0)+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*((x0m-x1m)*A1*k/Cp+(x2m-x1m)*A2*k/Cp-FbotIn*(x1m-x0m)+FbotOut*(x2m-x1m)+(Tbot-x1m)*FbotIn*v1)+sigma_x*dw1)
model$addSystem(dx2m~dt/V2*((x1m-x2m)*A2*k/Cp+(x3m-x2m)*A3*k/Cp-FbotIn*(x2m-x1m)+FbotOut*(x3m-x2m)+(Tbot-x2m)*FbotIn*v2)+sigma_x*dw2)
model$addSystem(dx3m~dt/V3*((x2m-x3m)*A3*k/Cp+(x4m-x3m)*A4*k/Cp-FbotIn*(x3m-x2m)+FbotOut*(x4m-x3m)+(Tbot-x3m)*FbotIn*v3)+sigma_x*dw3)
model$addSystem(dx4m~dt/V4*((x3m-x4m)*A4*k/Cp+(x5m-x4m)*A5*k/Cp-FbotIn*(x4m-x3m)+FbotOut*(x5m-x4m)+(Tbot-x4m)*FbotIn*v4)+sigma_x*dw4)
model$addSystem(dx5m~dt/V5*((x4m-x5m)*A5*k/Cp+(x6m-x5m)*A6*k/Cp-FbotIn*(x5m-x4m)+FbotOut*(x6m-x5m)+(Tbot-x5m)*FbotIn*v5)+sigma_x*dw5)
model$addSystem(dx6m~dt/V6*((x5m-x6m)*A6*k/Cp+(x7m-x6m)*A7*k/Cp-FbotIn*(x6m-x5m)+FbotOut*(x7m-x6m))+sigma_x*dw6)
model$addSystem(dx7m~dt/V7*((x6m-x7m)*A7*k/Cp+(x8m-x7m)*A8*k/Cp-FbotIn*(x7m-x6m)+FbotOut*(x8m-x7m))+sigma_x*dw7)
model$addSystem(dx8m~dt/V8*((x7m-x8m)*A8*k/Cp+(x9m-x8m)*A9*k/Cp-FbotIn*(x8m-x7m)+FbotOut*(x9m-x8m))+sigma_x*dw8)
model$addSystem(dx9m~dt/V9*((x8m-x9m)*A9*k/Cp+(x10m-x9m)*A10*k/Cp-FbotIn*(x9m-x8m)+FbotOut*(x10m-x9m))+sigma_x*dw9)
model$addSystem(dx10m~dt/V10*((x9m-x10m)*A10*k/Cp+(x11m-x10m)*A11*k/Cp-FtopOut*(x10m-x9m)+FtopIn*(x11m-x10m)+(Tmid-x10m)*FmidIn*v10)+sigma_x*dw10)
model$addSystem(dx11m~dt/V11*((x10m-x11m)*A11*k/Cp+(x12m-x11m)*A12*k/Cp-FtopOut*(x11m-x10m)+FtopIn*(x12m-x11m)+(Tmid-x11m)*FmidIn*v11)+sigma_x*dw11)
model$addSystem(dx12m~dt/V12*((x11m-x12m)*A12*k/Cp+(x13m-x12m)*A13*k/Cp-FtopOut*(x12m-x11m)+FtopIn*(x13m-x12m)+(Tmid-x12m)*FmidIn*v12)+sigma_x*dw12)
model$addSystem(dx13m~dt/V13*((x12m-x13m)*A13*k/Cp+(x14m-x13m)*A14*k/Cp-FtopOut*(x13m-x12m)+FtopIn*(x14m-x13m))+sigma_x*dw13)
model$addSystem(dx14m~dt/V14*((x13m-x14m)*A14*k/Cp+(x15m-x14m)*A15*k/Cp-FtopOut*(x14m-x13m)+FtopIn*(x15m-x14m))+sigma_x*dw14)
model$addSystem(dx15m~dt/V15*((x14m-x15m)*A15*k/Cp+(ambientTemp-x15m)*Atop*ktop/Cp-FtopOut*(x15m-x14m)+(Ttop-x15m)*FtopIn)+sigma_x*dw15)




######SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=1,lb=0,ub=10))

#SystemParameters
model$setParameter(vbot=c(init=1.4323e-1,lb=0,ub=1e6))
model$setParameter(vmid=c(init=3.8134e-1,lb=0,ub=1e6))
model$setParameter(vtop=c(init=1.7673,lb=0,ub=1e6))
model$setParameter(k = c(init = 6.5176, lb = 0, ub = 1e6))
model$setParameter(ktop = c(init = 2.3303, lb = 0, ub = 1e6))
model$setParameter(f = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(f1 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(f2 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(fbot = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(ftop = c(init = 1, lb = 0, ub = 1e6))

model$setParameter(v0 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v1 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v2 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v3 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v4 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v5 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v10 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v11 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v12 = c(init = 1, lb = 0, ub = 1e6))

model$setParameter(vmid11 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(vmid12 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(vmid13 = c(init = 1, lb = 0, ub = 1e6))




return(model)
}
