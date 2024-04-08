make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/V0*((x1m-x0m)*bot1*FbotOut+(Tbot-x0m)*Vbot*FbotIn)+sigma_x*dw0-l1*dt/s*V0/(V0+V1)*(log(1+exp(s*(x0m-x1m)))))
model$addSystem(dx1m~dt/V1*(k1*(x2m-x1m)*FbotOut+h1*(x0m-x1m)*FbotIn)+sigma_x*dw1+l1*dt/s*(V0/(V0+V1)*log(1+exp(s*(x0m-x1m)))-V1/(V1+V2)*log(1+exp(s*(x1m-x2m)))))
model$addSystem(dx2m~dt/V2*(k1*(x3m-x2m)*FbotOut+h1*(x1m-x2m)*FbotIn)+sigma_x*dw2+l2*dt/s*(V1/(V1+V2)*log(1+exp(s*(x1m-x2m)))-V2/(V2+V3)*log(1+exp(s*(x2m-x3m)))))
model$addSystem(dx3m~dt/V3*(k2*(x4m-x3m)*FbotOut+h2*(x2m-x3m)*FbotIn)+sigma_x*dw3+l2*dt/s*(V2/(V2+V3)*log(1+exp(s*(x2m-x3m)))-V3/(V3+V4)*log(1+exp(s*(x3m-x4m)))))
model$addSystem(dx4m~dt/V4*(k2*(x5m-x4m)*FbotOut+h2*(x3m-x4m)*FbotIn)+sigma_x*dw4+l3*dt/s*(V3/(V3+V4)*log(1+exp(s*(x3m-x4m)))-V4/(V4+V5)*log(1+exp(s*(x4m-x5m)))))
model$addSystem(dx5m~dt/V5*(k3*(x6m-x5m)*FbotOut+h3*(x4m-x5m)*FbotIn)+sigma_x*dw5+l3*dt/s*(V4/(V4+V5)*log(1+exp(s*(x4m-x5m)))-V5/(V5+V6)*log(1+exp(s*(x5m-x6m)))))
model$addSystem(dx6m~dt/V6*(k3*(x7m-x6m)*FbotOut+h3*(x5m-x6m)*FbotIn)+sigma_x*dw6+l4*dt/s*(V5/(V5+V6)*log(1+exp(s*(x5m-x6m)))-V6/(V6+V7)*log(1+exp(s*(x6m-x7m)))))
model$addSystem(dx7m~dt/V7*(k4*(x8m-x7m)*FbotOut+h4*(x6m-x7m)*FbotIn)+sigma_x*dw7+l4*dt/s*(V6/(V6+V7)*log(1+exp(s*(x6m-x7m)))-V7/(V7+V8)*log(1+exp(s*(x7m-x8m)))))
model$addSystem(dx8m~dt/V8*(k4*(x9m-x8m)*FbotOut+h4*(x7m-x8m)*FbotIn)+sigma_x*dw8+l5*dt/s*(V7/(V7+V8)*log(1+exp(s*(x7m-x8m)))-V8/(V8+V9)*log(1+exp(s*(x8m-x9m)))))
model$addSystem(dx9m~dt/V9*(k5*(x10m-x9m)*FbotOut+h5*(x8m-x9m)*FbotIn)+sigma_x*dw9+l5*dt/s*(V8/(V8+V9)*log(1+exp(s*(x8m-x9m)))-V9/(V9+V10)*log(1+exp(s*(x9m-x10m)))))

model$addSystem(dx10m~dt/V10*(k5*(x11m-x10m)*FtopIn+h5*(x9m-x10m)*FtopOut+(Tmid-x10m)*FmidIn*Vmid)+sigma_x*dw10+l6*dt/s*(V9/(V9+V10)*log(1+exp(s*(x9m-x10m)))-V10/(V10+V11)*log(1+exp(s*(x10m-x11m)))))
model$addSystem(dx11m~dt/V11*(k6*(x12m-x11m)*FtopIn+h6*(x10m-x11m)*FtopOut)+sigma_x*dw11+l6*dt/s*(V10/(V10+V11)*log(1+exp(s*(x10m-x11m)))-V11/(V11+V2)*log(1+exp(s*(x11m-x12m)))))
model$addSystem(dx12m~dt/V12*(k6*(x13m-x12m)*FtopIn+h6*(x11m-x12m)*FtopOut)+sigma_x*dw12+l6*dt/s*(V11/(V11+V12)*log(1+exp(s*(x11m-x12m)))-V12/(V12+V13)*log(1+exp(s*(x12m-x13m)))))
model$addSystem(dx13m~dt/V13*(k7*(x14m-x13m)*FtopIn+h7*(x12m-x13m)*FtopOut)+sigma_x*dw13+l7*dt/s*(V12/(V12+V13)*log(1+exp(s*(x12m-x13m)))-V13/(V13+V14)*log(1+exp(s*(x13m-x14m)))))
model$addSystem(dx14m~dt/V14*(k7*(x15m-x14m)*FtopIn+h7*(x13m-x14m)*FtopOut)+sigma_x*dw14+l7*dt/s*(V13/(V13+V14)*log(1+exp(s*(x13m-x14m)))-V14/(V14+V15)*log(1+exp(s*(x14m-x15m)))))
model$addSystem(dx15m~dt/V15*((Ttop-x15m)*Vtop*FtopIn-(x15m-x14m)*top1*FtopOut+ktop*(ambientTemp-x15m))+sigma_x*dw15+l8*dt/s*(V14/(V14+V15)*log(1+exp(s*(x14m-x15m)))))



#####SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=0.15,lb=0,ub=10))

#SystemParameters

model$setParameter(Vmid=c(init=1,lb=0,ub=10))
model$setParameter(Vbot=c(init=1,lb=0,ub=10))
model$setParameter(Vtop=c(init=1,lb=0,ub=10))
model$setParameter(top1=c(init=1,lb=0,ub=10))
model$setParameter(bot1=c(init=1,lb=0,ub=10))

model$setParameter(f1=c(init=1,lb=0,ub=10))
model$setParameter(f2=c(init=1,lb=0,ub=10))
model$setParameter(f3=c(init=1,lb=0,ub=10))
model$setParameter(f4=c(init=1,lb=0,ub=10))
model$setParameter(f5=c(init=1,lb=0,ub=10))
model$setParameter(f6=c(init=1,lb=0,ub=10))
model$setParameter(f7=c(init=1,lb=0,ub=10))

model$setParameter(h1=c(init=1,lb=0,ub=10))
model$setParameter(h2=c(init=1,lb=0,ub=10))
model$setParameter(h3=c(init=1,lb=0,ub=10))
model$setParameter(h4=c(init=1,lb=0,ub=10))
model$setParameter(h5=c(init=1,lb=0,ub=10))
model$setParameter(h6=c(init=1,lb=0,ub=10))
model$setParameter(h7=c(init=1,lb=0,ub=10))


model$setParameter(k1=c(init=1,lb=0,ub=10))
model$setParameter(k2=c(init=1,lb=0,ub=10))
model$setParameter(k3=c(init=1,lb=0,ub=10))
model$setParameter(k4=c(init=1,lb=0,ub=10))
model$setParameter(k5=c(init=1,lb=0,ub=10))
model$setParameter(k6=c(init=1,lb=0,ub=10))
model$setParameter(k7=c(init=1,lb=0,ub=10))

model$setParameter(l1=c(init=1,lb=0,ub=10))
model$setParameter(l2=c(init=1,lb=0,ub=10))
model$setParameter(l3=c(init=1,lb=0,ub=10))
model$setParameter(l4=c(init=1,lb=0,ub=10))
model$setParameter(l5=c(init=1,lb=0,ub=10))
model$setParameter(l6=c(init=1,lb=0,ub=10))
model$setParameter(l7=c(init=1,lb=0,ub=10))
model$setParameter(l8=c(init=1,lb=0,ub=10))




model$setParameter(ktop=c(init=1,lb=0,ub=10))

model$setParameter(s = 100)



return(model)
}