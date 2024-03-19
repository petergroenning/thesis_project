make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/C*(K*(x1m-x0m)+FbotOut*(x1m-x0m)+(Tbot-x0m)*FbotIn*v0)+sigma_x*dw0)
model$addSystem(dx1m~dt/C*(K*(x0m-2*x1m+x2m)-Fbot*(x2m-x0m)/2+(Tbot-x1m)*FbotIn*v1)+sigma_x*dw1)
model$addSystem(dx2m~dt/C*(K*(x1m-2*x2m+x3m)-Fbot*(x3m-x1m)/2+(Tbot-x2m)*FbotIn*v2)+sigma_x*dw2)
model$addSystem(dx3m~dt/C*(K*(x2m-2*x3m+x4m)-Fbot*(x4m-x2m)/2+(Tbot-x3m)*FbotIn*v3)+sigma_x*dw3)
model$addSystem(dx4m~dt/C*(K*(x3m-2*x4m+x5m)-Fbot*(x5m-x3m)/2+(Tbot-x4m)*FbotIn*v4)+sigma_x*dw4)
model$addSystem(dx5m~dt/C*(K*(x4m-2*x5m+x6m)-Fbot*(x6m-x4m)/2+(Tbot-x5m)*FbotIn*v5)+sigma_x*dw5)
model$addSystem(dx6m~dt/C*(K*(x5m-2*x6m+x7m)-Fbot*(x7m-x5m)/2)+sigma_x*dw6)
model$addSystem(dx7m~dt/C*(K*(x6m-2*x7m+x8m)-Fbot*(x8m-x6m)/2)+sigma_x*dw7)
model$addSystem(dx8m~dt/C*(K*(x7m-2*x8m+x9m)-Fbot*(x9m-x7m)/2)+sigma_x*dw8)
model$addSystem(dx9m~dt/C*(K*(x8m-2*x9m+x10m)-Fbot*(x10m-x8m)/2)+sigma_x*dw9)
model$addSystem(dx10m~dt/C*(K*(x9m-2*x10m+x11m)+Ftop*(x11m-x9m)/2+(Tmid-x10m)*FmidIn*v10)+sigma_x*dw10)
model$addSystem(dx11m~dt/C*(K*(x10m-2*x11m+x12m)+Ftop*(x12m-x10m)/2+(Tmid-x11m)*FmidIn*v11)+sigma_x*dw11)
model$addSystem(dx12m~dt/C*(K*(x11m-2*x12m+x13m)+Ftop*(x13m-x11m)/2+(Tmid-x12m)*FmidIn*v12)+sigma_x*dw12)
model$addSystem(dx13m~dt/C*(K*(x12m-2*x13m+x14m)+Ftop*(x14m-x12m)/2+(Tmid-x13m)*FmidIn*v13)+sigma_x*dw13)
model$addSystem(dx14m~dt/C*(K*(x13m-2*x14m+x15m)+Ftop*(x15m-x13m)/2+(Ttop-x14m)*FtopIn*v14)+sigma_x*dw14)
model$addSystem(dx15m~dt/C*(K*(x14m-x15m)+(ambientTemp-x15m)*Ktop-FtopOut*(x15m-x14m)+(Ttop-x15m)*FtopIn*v15)+sigma_x*dw15)




######SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=1,lb=0,ub=10))

#SystemParameters
model$setParameter(K= c(init = 1, lb = 0, ub = 1e6))
model$setParameter(Ktop= c(init = 1, lb = 0, ub = 1e6))
model$setParameter(C= c(init=100,lb=0,ub=1e6))
model$setParameter(v0 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v1 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v2 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v3 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v4 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v5 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v10 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v11 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v12 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v13 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v14 = c(init = 1, lb = 0, ub = 1e6))
model$setParameter(v15 = c(init = 1, lb = 0, ub = 1e6))




return(model)
}
