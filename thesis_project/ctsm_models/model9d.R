make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(1/b11+(x-x0m)/b21)*FbotIn+(x1m-x0m)*(1/b12+(x1m-x0m)/b22)*(FbotOut))+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*((x2m-x1m)*(1/k11+(x1m-x0m)/k21)*(FbotOut)+(x0m-x1m)*(1/k12+(x2m-x1m)/k22)*(FbotIn))+sigma_x*dw1)
model$addSystem(dx2m~dt/V2*((x3m-x2m)*(1/k11+(x2m-x1m)/k21)*(FbotOut)+(x1m-x2m)*(1/k12+(x3m-x2m)/k22)*(FbotIn))+sigma_x*dw2)
model$addSystem(dx3m~dt/V3*((x4m-x3m)*(1/k11+(x3m-x2m)/k21)*(FbotOut)+(x2m-x3m)*(1/k12+(x4m-x3m)/k22)*(FbotIn))+sigma_x*dw3)
model$addSystem(dx4m~dt/V4*((x5m-x4m)*(1/k11+(x4m-x3m)/k21)*(FbotOut)+(x3m-x4m)*(1/k12+(x5m-x4m)/k22)*(FbotIn))+sigma_x*dw4)
model$addSystem(dx5m~dt/V5*((x6m-x5m)*(1/k11+(x5m-x4m)/k21)*(FbotOut)+(x4m-x5m)*(1/k12+(x6m-x5m)/k22)*(FbotIn))+sigma_x*dw5)
model$addSystem(dx6m~dt/V6*((x7m-x6m)*(1/k11+(x6m-x5m)/k21)*(FbotOut)+(x5m-x6m)*(1/k12+(x7m-x6m)/k22)*(FbotIn))+sigma_x*dw6)
model$addSystem(dx7m~dt/V7*((x8m-x7m)*(1/k11+(x7m-x6m)/k21)*(FbotOut)+(x6m-x7m)*(1/k12+(x8m-x7m)/k22)*(FbotIn))+sigma_x*dw7)
model$addSystem(dx8m~dt/V8*((x9m-x8m)*(1/k11+(x8m-x7m)/k21)*(FbotOut)+(x7m-x8m)*(1/k12+(x9m-x8m)/k22)*(FbotIn))+sigma_x*dw8) 
model$addSystem(dx9m~dt/V9*((x10m-x9m)*(1/k11+(x9m-x8m)/k21)*(FbotOut)+(x8m-x9m)*(1/k12+(x10m-x9m)/k22)*(FbotIn))+sigma_x*dw9)

model$addSystem(dx10m~dt/V10*((x9m-x10m)*(1/h11+(x11m-x10m)/h21)*(FtopOut)+(x11m-x10m)*(1/h12+(x10m-x9m)/h22)*(FtopIn)+(Tmid-x10m)*FmidIn)+sigma_x*dw10)
model$addSystem(dx11m~dt/V11*((x10m-x11m)*(1/h11+(x12m-x11m)/h21)*(FtopOut)+(x12m-x11m)*(1/h12+(x11m-x10m)/h22)*(FtopIn))+sigma_x*dw11)
model$addSystem(dx12m~dt/V12*((x11m-x12m)*(1/h11+(x13m-x12m)/h21)*(FtopOut)+(x13m-x12m)*(1/h12+(x12m-x11m)/h22)*(FtopIn))+sigma_x*dw12)
model$addSystem(dx13m~dt/V13*((x12m-x13m)*(1/h11+(x14m-x13m)/h21)*(FtopOut)+(x14m-x13m)*(1/h12+(x13m-x12m)/h22)*(FtopIn))+sigma_x*dw13)
model$addSystem(dx14m~dt/V14*((x13m-x14m)*(1/h11+(x15m-x14m)/h21)*(FtopOut)+(x15m-x14m)*(1/h12+(x14m-x13m)/h22)*(FtopIn))+sigma_x*dw14)
model$addSystem(dx15m~dt/V15*((Ttop-x15m)*(1/t11+(x15m-x14m)/t21)*(FtopIn)+(x14m-x15m)*(1/t12+(x14m-x15m)/t22)*FtopOut)+sigma_x*dw15)




######SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

#SystemParameters
model$setParameter(k11=c(init=4.3685e-01,lb=0,ub=2000))
model$setParameter(k21=c(init=2.6533e+00,lb=0,ub=2000))
model$setParameter(k12=c(init=2.1770e-01,lb=0,ub=2000))
model$setParameter(k22=c(init=1.9720e+02,lb=0,ub=2000))

model$setParameter(h11=c(init=3.0,lb=0,ub=2000))
model$setParameter(h21=c(init=6.638,lb=0,ub=2000))
model$setParameter(h12=c(init=1.40,lb=0,ub=2000))
model$setParameter(h22=c(init=7.8672,lb=0,ub=2000))

model$setParameter(b11=c(init= 4.0655e-01,lb=0,ub=2000))
model$setParameter(b21=c(init= 6.2723e-02,lb=0,ub=2000))
model$setParameter(b12=c(init=5.4108e-02,lb=0,ub=2000))
model$setParameter(b22=c(init=1.9971e+02,lb=0,ub=2000))

model$setParameter(t11=c(init=4.9876e-01,lb=0,ub=2000))
model$setParameter(t21=c(init= 4.4134e+01,lb=0,ub=2000))
model$setParameter(t12=c(init= 1.1060,lb=0,ub=2000))
model$setParameter(t22=c(init=1.7180e+02,lb=0,ub=2000))


model$setParameter(V=c(init=1.2483e+04,lb=1000,ub=20000))


return(model)
}
