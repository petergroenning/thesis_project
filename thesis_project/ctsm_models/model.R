make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*Vbot*FbotOut+(Tbot-x0m)*Rbot*FbotIn+k_bot*(x1m-x0m))+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*((x2m-x1m)*(x1m-x0m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x0m-2*x1m+x2m)+((x2m-x1m)*FbotOut*Rbot3-(x1m-x0m)*FbotIn*Rbot4))+sigma_x*dw1+dt/s*(log(1+exp(s*(x0m-x1m)))))
model$addSystem(dx2m~dt/V2*((x3m-x2m)*(x2m-x1m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x1m-2*x2m+x3m)+((x3m-x2m)*FbotOut*Rbot3-(x2m-x1m)*FbotIn*Rbot4))+sigma_x*dw2+dt/s*(log(1+exp(s*(x1m-x2m)))))
model$addSystem(dx3m~dt/V3*((x4m-x3m)*(x3m-x2m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x2m-2*x3m+x4m)+((x4m-x3m)*FbotOut*Rbot3-(x3m-x2m)*FbotIn*Rbot4))+sigma_x*dw3+dt/s*(log(1+exp(s*(x2m-x3m)))))
model$addSystem(dx4m~dt/V4*((x5m-x4m)*(x4m-x3m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x3m-2*x4m+x5m)+((x5m-x4m)*FbotOut*Rbot3-(x4m-x3m)*FbotIn*Rbot4))+sigma_x*dw4+dt/s*(log(1+exp(s*(x3m-x4m)))))
model$addSystem(dx5m~dt/V5*((x6m-x5m)*(x5m-x4m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x4m-2*x5m+x6m)+((x6m-x5m)*FbotOut*Rbot3-(x5m-x4m)*FbotIn*Rbot4))+sigma_x*dw5+dt/s*(log(1+exp(s*(x4m-x5m)))))
model$addSystem(dx6m~dt/V6*((x7m-x6m)*(x6m-x5m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x5m-2*x6m+x7m)+((x7m-x6m)*FbotOut*Rbot3-(x6m-x5m)*FbotIn*Rbot4))+sigma_x*dw6+dt/s*(log(1+exp(s*(x5m-x6m)))))
model$addSystem(dx7m~dt/V7*((x8m-x7m)*(x7m-x6m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x6m-2*x7m+x8m)+((x8m-x7m)*FbotOut*Rbot3-(x7m-x6m)*FbotIn*Rbot4))+sigma_x*dw7+dt/s*(log(1+exp(s*(x6m-x7m)))))
model$addSystem(dx8m~dt/V8*((x9m-x8m)*(x8m-x7m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x7m-2*x8m+x9m)+((x9m-x8m)*FbotOut*Rbot3-(x8m-x7m)*FbotIn*Rbot4))+sigma_x*dw8+dt/s*(log(1+exp(s*(x7m-x8m)))))
model$addSystem(dx9m~dt/V9*((x10m-x9m)*(x9m-x8m)*(FbotOut*Rbot1-FbotIn*Rbot2)+k*(x8m-2*x9m+x10m)+((x10m-x9m)*FbotOut*Rbot3-(x9m-x8m)*FbotIn*Rbot4))+sigma_x*dw9+dt/s*(log(1+exp(s*(x8m-x9m)))))

model$addSystem(dx10m~dt/V10*((x11m-x10m)*(x10m-x9m)*((FtopIn+FmidIn)*Rtop2-(FtopOut+FmidIn)*Rtop1)+(Tmid-x10m)*FmidIn*Vmid+k*(x9m-2*x10m+x11m)+((x11m-x10m)*(FtopIn+FmidIn)*Rtop3-(x10m-x9m)*(FtopOut+FmidIn)*Rtop4))+sigma_x*dw10+dt/s*(log(1+exp(s*(x9m-x10m)))))
model$addSystem(dx11m~dt/V11*((x12m-x11m)*(x11m-x10m)*(FtopIn*Rtop2-FtopOut*Rtop1)+k*(x10m-2*x11m+x12m)+((x12m-x11m)*FtopIn*Rtop3-(x11m-x10m)*FtopOut*Rtop4))+sigma_x*dw11+dt/s*(log(1+exp(s*(x10m-x11m)))))
model$addSystem(dx12m~dt/V12*((x13m-x12m)*(x12m-x11m)*(FtopIn*Rtop2-FtopOut*Rtop1)+k*(x11m-2*x12m+x13m)+((x13m-x12m)*FtopIn*Rtop3-(x12m-x11m)*FtopOut*Rtop4))+sigma_x*dw12+dt/s*(log(1+exp(s*(x11m-x12m)))))
model$addSystem(dx13m~dt/V13*((x14m-x13m)*(x13m-x12m)*(FtopIn*Rtop2-FtopOut*Rtop1)+k*(x12m-2*x13m+x14m)+((x14m-x13m)*FtopIn*Rtop3-(x13m-x12m)*FtopOut*Rtop4))+sigma_x*dw13+dt/s*(log(1+exp(s*(x12m-x13m)))))
model$addSystem(dx14m~dt/V14*((x15m-x14m)*(x14m-x13m)*(FtopIn*Rtop2-FtopOut*Rtop1)+k*(x13m-2*x14m+x15m)+((x15m-x14m)*FtopIn*Rtop3-(x14m-x13m)*FtopOut*Rtop4))+sigma_x*dw14+dt/s*(log(1+exp(s*(x13m-x14m)))))
model$addSystem(dx15m~dt/V15*((Ttop-x15m)*Vtop*FtopIn-(x15m-x14m)*Rtop*FtopOut+k_top*(x14m-x15m))+sigma_x*dw15+dt/s*(log(1+exp(s*(x14m-x15m)))))



#####SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

#SystemParameters
model$setParameter(Rbot1=c(init=2,lb=0,ub=1e6))
model$setParameter(Rbot2=c(init=2,lb=0,ub=1e6))
model$setParameter(Rbot3=c(init=2,lb=0,ub=1e6))
model$setParameter(Rbot4=c(init=2,lb=0,ub=1e6))

model$setParameter(Rtop1=c(init=1,lb=0,ub=1e6))
model$setParameter(Rtop2=c(init=1,lb=0,ub=1e6))
model$setParameter(Rtop3=c(init=1,lb=0,ub=1e6))
model$setParameter(Rtop4=c(init=1,lb=0,ub=1e6))

model$setParameter(Vmid=c(init=1,lb=0,ub=1e6))
model$setParameter(Vbot=c(init=1,lb=0,ub=1e3))
model$setParameter(Vtop=c(init=1,lb=0,ub=1e3))


model$setParameter(k=c(init=1,lb=0,ub=1e6))
model$setParameter(k_bot=c(init=1,lb=0,ub=1e6))
model$setParameter(k_top=c(init=1,lb=0,ub=1e6))

model$setParameter(Rtop=c(init=1,lb=0,ub=1e6))
model$setParameter(Rbot=c(init=1,lb=0,ub=1e6))

model$setParameter(s = 50)



return(model)
}