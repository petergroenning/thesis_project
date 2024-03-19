make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*(A1*k/Cp)+f1*FbotOut*(x1m-x0m)+(Tbot-x0m)*FbotIn)+dt/c*(-log(1+exp(c*(x0m-x1m))))+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*((x0m-x1m)*(A1*k/Cp)+(x2m-x1m)*A2*k/Cp-Fbot*(x2m-x0m)/2)+dt/c*(log(1+exp(c*(x0m-x1m)))-log(1+exp(c*(x1m-x2m))))+sigma_x*dw1)
model$addSystem(dx2m~dt/V2*((x1m-x2m)*A2*k/Cp+(x3m-x2m)*A3*k/Cp-Fbot*(x3m-x1m)/2)+dt/c*(log(1+exp(c*(x1m-x2m)))-log(1+exp(c*(x2m-x3m))))+sigma_x*dw2)
model$addSystem(dx3m~dt/V3*((x2m-x3m)*A3*k/Cp+(x4m-x3m)*A4*k/Cp-Fbot*(x4m-x2m)/2)+dt/c*(log(1+exp(c*(x2m-x3m)))-log(1+exp(c*(x3m-x4m))))+sigma_x*dw3)
model$addSystem(dx4m~dt/V4*((x3m-x4m)*A4*k/Cp+(x5m-x4m)*A5*k/Cp-Fbot*(x5m-x3m)/2)+dt/c*(log(1+exp(c*(x3m-x4m)))-log(1+exp(c*(x4m-x5m))))+sigma_x*dw4)
model$addSystem(dx5m~dt/V5*((x4m-x5m)*A5*k/Cp+(x6m-x5m)*A6*k/Cp-Fbot*(x6m-x4m)/2)+dt/c*(log(1+exp(c*(x4m-x5m)))-log(1+exp(c*(x5m-x6m))))+sigma_x*dw5)
model$addSystem(dx6m~dt/V6*((x5m-x6m)*A6*k/Cp+(x7m-x6m)*A7*k/Cp-Fbot*(x7m-x5m)/2)+dt/c*(log(1+exp(c*(x5m-x6m)))-log(1+exp(c*(x6m-x7m))))+sigma_x*dw6)
model$addSystem(dx7m~dt/V7*((x6m-x7m)*A7*k/Cp+(x8m-x7m)*A8*k/Cp-Fbot*(x8m-x6m)/2)+dt/c*(log(1+exp(c*(x6m-x7m)))-log(1+exp(c*(x7m-x8m))))+sigma_x*dw7)
model$addSystem(dx8m~dt/V8*((x7m-x8m)*A8*k/Cp+(x9m-x8m)*A9*k/Cp-Fbot*(x9m-x7m)/2)+dt/c*(log(1+exp(c*(x7m-x8m)))-log(1+exp(c*(x8m-x9m))))+sigma_x*dw8)
model$addSystem(dx9m~dt/V9*((x8m-x9m)*A9*k/Cp+(x10m-x9m)*(A10*k/Cp+FmidVol*f2)-Fbot*(x10m-x8m)/2)+dt/c*(log(1+exp(c*(x8m-x9m)))-log(1+exp(c*(x9m-x10m))))+sigma_x*dw9)
model$addSystem(dx10m~dt/V10*((x9m-x10m)*(A10*k/Cp+FmidVol*f2)+(x11m-x10m)*(A11*k/Cp+FmidVol*f2)+Ftop*(x11m-x9m)/2+(Tmid-x10m)*FmidIn)+dt/c*(log(1+exp(c*(x9m-x10m)))-log(1+exp(c*(x10m-x11m))))+sigma_x*dw10)
model$addSystem(dx11m~dt/V11*((x10m-x11m)*(A11*k/Cp+FmidVol*f2)+(x12m-x11m)*A12*k/Cp+Ftop*(x12m-x10m)/2)+dt/c*(log(1+exp(c*(x10m-x11m)))-log(1+exp(c*(x11m-x12m))))+sigma_x*dw11)
model$addSystem(dx12m~dt/V12*((x11m-x12m)*A12*k/Cp+(x13m-x12m)*A13*k/Cp+Ftop*(x13m-x11m)/2)+dt/c*(log(1+exp(c*(x11m-x12m)))-log(1+exp(c*(x12m-x13m))))+sigma_x*dw12)
model$addSystem(dx13m~dt/V13*((x12m-x13m)*A13*k/Cp+(x14m-x13m)*A14*k/Cp+Ftop*(x14m-x12m)/2)+dt/c*(log(1+exp(c*(x12m-x13m)))-log(1+exp(c*(x13m-x14m))))+sigma_x*dw13)
model$addSystem(dx14m~dt/V14*((x13m-x14m)*A14*k/Cp+(x15m-x14m)*(A15*k/Cp)+Ftop*(x15m-x13m)/2)+dt/c*(log(1+exp(c*(x13m-x14m)))-log(1+exp(c*(x14m-x15m))))+sigma_x*dw14)
model$addSystem(dx15m~dt/V15*((x14m-x15m)*(A15*k/Cp)+(ambientTemp-x15m)*Atop*ktop/Cp-f3*FtopOut*(x15m-x14m)+(Ttop-x15m)*FtopIn)+dt/c*(log(1+exp(c*(x14m-x15m))))+sigma_x*dw15)




######SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=1,lb=0,ub=10))

#SystemParameters
model$setParameter(k = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(ktop = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(c = 1e2)
model$setParameter(f1 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(f2 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(f3 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(d = c(init = 0.5, lb = 0, ub = 1))




return(model)
}
