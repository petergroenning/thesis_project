make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*bot1*FbotOut+(Tbot-x0m)*bot2*FbotIn)+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*((x2m-x1m)^3*(FbotOut*a)+(x0m-x1m)*FbotIn*b)+sigma_x*dw1)
model$addSystem(dx2m~dt/V2*((x3m-x2m)^3*(FbotOut*a)+(x1m-x2m)*FbotIn*b)+sigma_x*dw2)
model$addSystem(dx3m~dt/V3*((x4m-x3m)^3*(FbotOut*a)+(x2m-x3m)*FbotIn*b)+sigma_x*dw3)
model$addSystem(dx4m~dt/V4*((x5m-x4m)^3*(FbotOut*a)+(x3m-x4m)*FbotIn*b)+sigma_x*dw4)
model$addSystem(dx5m~dt/V5*((x6m-x5m)^3*(FbotOut*a)+(x4m-x5m)*FbotIn*b)+sigma_x*dw5)
model$addSystem(dx6m~dt/V6*((x7m-x6m)^3*(FbotOut*a)+(x5m-x6m)*FbotIn*b)+sigma_x*dw6)
model$addSystem(dx7m~dt/V7*((x8m-x7m)^3*(FbotOut*a)+(x6m-x7m)*FbotIn*b)+sigma_x*dw7)
model$addSystem(dx8m~dt/V8*((x9m-x8m)^3*(FbotOut*a)+(x7m-x8m)*FbotIn*b)+sigma_x*dw8) 
model$addSystem(dx9m~dt/V9*((x10m-x9m)^3*(FbotOut*a)+(x8m-x9m)*FbotIn*b)+sigma_x*dw9)

model$addSystem(dx10m~dt/V10*((x9m-x10m)*(FtopOut+FmidIn)*c+(x11m-x10m)*(FtopIn+FmidIn)*d+(Tmid-x10m)*FmidIn*vmid)+sigma_x*dw10)
model$addSystem(dx11m~dt/V11*((x10m-x11m)*(FtopOut*c)+(x12m-x11m)*(FtopIn*d))+sigma_x*dw11)
model$addSystem(dx12m~dt/V12*((x11m-x12m)*(FtopOut*c)+(x13m-x12m)*(FtopIn*d))+sigma_x*dw12)
model$addSystem(dx13m~dt/V13*((x12m-x13m)*(FtopOut*c)+(x14m-x13m)*(FtopIn*d))+sigma_x*dw13)
model$addSystem(dx14m~dt/V14*((x13m-x14m)*(FtopOut*c)+(x15m-x14m)*(FtopIn*d))+sigma_x*dw14)
model$addSystem(dx15m~dt/V15*((x14m-x15m)*top1*FtopOut+(Ttop-x15m)*top2*FtopIn)+sigma_x*dw15)



#####SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

#SystemParameters
model$setParameter(a=c(init=2,lb=0,ub=1e6))
model$setParameter(b=c(init=2,lb=0,ub=1e6))
model$setParameter(c=c(init=2,lb=0,ub=1e6))
model$setParameter(d=c(init=2,lb=0,ub=1e6))
model$setParameter(vmid=c(init=1,lb=0,ub=1e6))
model$setParameter(bot1=c(init=1,lb=0,ub=1e3))
model$setParameter(bot2=c(init=1,lb=0,ub=1e3))
model$setParameter(top1=c(init=1,lb=0,ub=1e3))
model$setParameter(top2=c(init=1,lb=0,ub=1e3))





return(model)
}