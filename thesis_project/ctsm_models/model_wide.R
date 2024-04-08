make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*bot1*FbotOut+(Tbot-x0m)*Vbot*FbotIn+(x0m-x1m)*pbot)+sigma_x*dw0)
model$addSystem(dx1m~dt/V1*(k1*(x2m-x1m)*FbotOut+h1*(x0m-x1m)*FbotIn+(x1m-x2m)*p1+(x1m-x0m)*pbot)+sigma_x*dw1)
model$addSystem(dx2m~dt/V2*(k1*(x3m-x2m)*FbotOut+h1*(x1m-x2m)*FbotIn+(x2m-x3m)*p1+(x2m-x1m)*p1)+sigma_x*dw2)
model$addSystem(dx3m~dt/V3*(k2*(x4m-x3m)*FbotOut+h2*(x2m-x3m)*FbotIn+(x3m-x4m)*p1+(x3m-x2m)*p1)+sigma_x*dw3)
model$addSystem(dx4m~dt/V4*(k2*(x5m-x4m)*FbotOut+h2*(x3m-x4m)*FbotIn+(x4m-x5m)*p2+(x4m-x3m)*p1)+sigma_x*dw4)
model$addSystem(dx5m~dt/V5*(k3*(x6m-x5m)*FbotOut+h3*(x4m-x5m)*FbotIn+(x5m-x6m)*p2+(x5m-x4m)*p2)+sigma_x*dw5)
model$addSystem(dx6m~dt/V6*(k3*(x7m-x6m)*FbotOut+h3*(x5m-x6m)*FbotIn+(x6m-x7m)*p2+(x6m-x5m)*p2)+sigma_x*dw6)
model$addSystem(dx7m~dt/V7*(k4*(x8m-x7m)*FbotOut+h4*(x6m-x7m)*FbotIn+(x7m-x8m)*p3+(x7m-x6m)*p2)+sigma_x*dw7)
model$addSystem(dx8m~dt/V8*(k4*(x9m-x8m)*FbotOut+h4*(x7m-x8m)*FbotIn+(x8m-x9m)*p3+(x8m-x7m)*p3)+sigma_x*dw8)
model$addSystem(dx9m~dt/V9*(k5*(x10m-x9m)*FbotOut+h5*(x8m-x9m)*FbotIn+(x9m-x10m)*p3+(x9m-x8m)*q3)+sigma_x*dw9)

model$addSystem(dx10m~dt/V10*(k5*(x11m-x10m)*FtopIn+h5*(x9m-x10m)*FtopOut+(Tmid-x10m)*FmidIn*Vmid+(x10m-x11m)*p4+(x10m-x9m)*p3)+sigma_x*dw10)
model$addSystem(dx11m~dt/V11*(k6*(x12m-x11m)*FtopIn+h6*(x10m-x11m)*FtopOut+(x11m-x12m)*p4+(x11m-x10m)*p4)+sigma_x*dw11)
model$addSystem(dx12m~dt/V12*(k6*(x13m-x12m)*FtopIn+h6*(x11m-x12m)*FtopOut+(x12m-x13m)*p5+(x12m-x11m)*q4)+sigma_x*dw12)
model$addSystem(dx13m~dt/V13*(k7*(x14m-x13m)*FtopIn+h7*(x12m-x13m)*FtopOut+(x13m-x14m)*p5+(x13m-x12m)*p5)+sigma_x*dw13)
model$addSystem(dx14m~dt/V14*(k7*(x15m-x14m)*FtopIn+h7*(x13m-x14m)*FtopOut+(x14m-x15m)*ptop+(x14m-x13m)*p5)+sigma_x*dw14)
model$addSystem(dx15m~dt/V15*((Ttop-x15m)*Vtop*FtopIn-(x15m-x14m)*top1*FtopOut+(x15m-ambientTemp)*ptop+(x15m-x14m)*qtop)+sigma_x*dw15)



#####SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

#SystemParameters

model$setParameter(Vmid=c(init=1,lb=0,ub=10))
model$setParameter(Vbot=c(init=1,lb=0,ub=10))
model$setParameter(Vtop=c(init=1,lb=0,ub=10))
model$setParameter(top1=c(init=1,lb=-10,ub=10))
model$setParameter(bot1=c(init=1,lb=0,ub=10))

model$setParameter(p1=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p2=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p3=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p4=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p5=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p6=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(p7=c(init=1e-2,lb=-1e2,ub=1e2))

model$setParameter(q1=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q2=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q3=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q4=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q5=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q6=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(q7=c(init=1e-2,lb=-1e2,ub=1e2))


model$setParameter(ptop=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(qtop=c(init=1e-2,lb=-1e2,ub=1e2))
model$setParameter(pbot=c(init=1e-2,lb=-1e2,ub=1e2))


model$setParameter(k1=c(init=1,lb=0,ub=10))
model$setParameter(k2=c(init=1,lb=0,ub=10))
model$setParameter(k3=c(init=1,lb=0,ub=10))
model$setParameter(k4=c(init=1,lb=0,ub=10))
model$setParameter(k5=c(init=1,lb=0,ub=10))
model$setParameter(k6=c(init=1,lb=0,ub=10))
model$setParameter(k7=c(init=1,lb=0,ub=10))

model$setParameter(h1=c(init=1,lb=0,ub=10))
model$setParameter(h2=c(init=1,lb=0,ub=10))
model$setParameter(h3=c(init=1,lb=0,ub=10))
model$setParameter(h4=c(init=1,lb=0,ub=10))
model$setParameter(h5=c(init=1,lb=0,ub=10))
model$setParameter(h6=c(init=1,lb=0,ub=10))
model$setParameter(h7=c(init=1,lb=0,ub=10))







model$setParameter(ktop=c(init=1,lb=0,ub=10))

model$setParameter(s = 50)



return(model)}