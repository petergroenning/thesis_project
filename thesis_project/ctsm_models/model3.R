make_model<-function(model){
#Addsystemequations

model$addSystem(dx0m~dt/(V0*Cp)*((x1m-x0m)*A1*k1+(Tsoil-x0m)*(A0+S0)*k2)+dt*((Tbot-x0m)*FbotIn/V0*vbot+(x1m-x0m)*(FbotOut)*A1/A0*f3)+sigma_x*dw0)
model$addSystem(dx1m~dt/(V1*Cp)*((x0m-x1m)*A1*k1+(x2m-x1m)*A2*k1+(Tsoil-x1m)*S1*k2)+dt*f3*A2/A1*((x1m-Tbot)*(-FbotIn)+(x2m-x0m)*FbotOut)/2+sigma_x*dw1)
model$addSystem(dx2m~dt/(V2*Cp)*((x1m-x2m)*A2*k1+(x3m-x2m)*A3*k1+(Tsoil-x2m)*S2*k2)+dt*f3*A3/A2*((x2m-Tbot)*(-FbotIn)+(x3m-x1m)*FbotOut)/2+sigma_x*dw2)
model$addSystem(dx3m~dt/(V3*Cp)*((x2m-x3m)*A3*k1+(x4m-x3m)*A4*k1+(Tsoil-x3m)*S3*k2)+dt*f3*A4/A3*((x3m-Tbot)*(-FbotIn)+(x4m-x2m)*FbotOut)/2+sigma_x*dw3)
model$addSystem(dx4m~dt/(V4*Cp)*((x3m-x4m)*A4*k1+(x5m-x4m)*A5*k1+(Tsoil-x4m)*S4*k2)+dt*f3*A5/A4*((x4m-Tbot)*(-FbotIn)+(x5m-x3m)*FbotOut)/2+sigma_x*dw4)
model$addSystem(dx5m~dt/(V5*Cp)*((x4m-x5m)*A5*k1+(x6m-x5m)*A6*k1+(Tsoil-x5m)*S5*k2)+dt*f3*A6/A5*((x5m-Tbot)*(-FbotIn)+(x6m-x4m)*FbotOut)/2+sigma_x*dw5)
model$addSystem(dx6m~dt/(V6*Cp)*((x5m-x6m)*A6*k1+(x7m-x6m)*A7*k1+(Tsoil-x6m)*S6*k2)+dt*(x7m-x5m)/2*(-Fbot)*f1*A7/A6+sigma_x*dw6)
model$addSystem(dx7m~dt/(V7*Cp)*((x6m-x7m)*A7*k1+(x8m-x7m)*A8*k1+(Tsoil-x7m)*S7*k2)+dt*(x8m-x6m)/2*(-Fbot)*f1*A8/A7+sigma_x*dw7)
model$addSystem(dx8m~dt/(V8*Cp)*((x7m-x8m)*A8*k1+(x9m-x8m)*A9*k1+(Tsoil-x8m)*S8*k2)+dt*(x9m-x7m)/2*(-Fbot)*f1*A9/A8+sigma_x*dw8)
model$addSystem(dx9m~dt/(V9*Cp)*((x8m-x9m)*A9*k1+(x10m-x9m)*A10*k1+(Tsoil-x9m)*S9*k2)+dt*(x10m-x8m)/2*(-Fbot)*f1*A10/A9+dt*(Tmid-x9m)*FmidIn/V9*vmid9+sigma_x*dw9)
model$addSystem(dx10m~dt/(V10*Cp)*((x9m-x10m)*A10*k1+(x11m-x10m)*A11*k1+(Tsoil-x10m)*S10*k2)+dt*(x11m-x9m)/2*Ftop*f2*A11/A10+dt*(Tmid-x10m)*FmidIn/V10*vmid10+sigma_x*dw10)
model$addSystem(dx11m~dt/(V11*Cp)*((x10m-x11m)*A11*k1+(x12m-x11m)*A12*k1+(Tsoil-x11m)*S11*k2)+dt*(x12m-x10m)/2*Ftop*f2*A12/A11+dt*(Tmid-x11m)*FmidIn/V11*vmid11+sigma_x*dw11)
model$addSystem(dx12m~dt/(V12*Cp)*((x11m-x12m)*A12*k1+(x13m-x12m)*A13*k1+(Tsoil-x12m)*S12*k2)+dt*(x13m-x11m)/2*Ftop*f2*A13/A12+dt*(Tmid-x12m)*FmidIn/V12*vmid12+sigma_x*dw12)
model$addSystem(dx13m~dt/(V13*Cp)*((x12m-x13m)*A13*k1+(x14m-x13m)*A14*k1+(Tsoil-x13m)*S13*k2)+dt*(x14m-x12m)/2*Ftop*f2*A14/A13+sigma_x*dw13)
model$addSystem(dx14m~dt/(V14*Cp)*((x13m-x14m)*A14*k1+(x15m-x14m)*A15*k1+(Tsoil-x14m)*S14*k2)+dt*(x15m-x13m)/2*Ftop*f2*A15/A14+sigma_x*dw14)
model$addSystem(dx15m~dt/(V15*Cp)*((x14m-x15m)*A15*k1+(ambientTemp-x15m)*Atop*k3)+dt*(Ttop-x15m)*FtopIn/V15*vtop+dt*(x15m-x14m)*(-FtopOut)*f4*Atop/A15+sigma_x*dw15)




######SETPARAMETERS########
#SystemNoise
model$setParameter(sigma_x=c(init=1,lb=0,ub=10))

#SystemParameters
model$setParameter(k1=c(init=1,lb=0,ub=1000))
model$setParameter(k2=c(init=1,lb=0,ub=1000))
model$setParameter(k3=c(init=1,lb=0,ub=1000))
model$setParameter(vbot=c(init=1,lb=0,ub=1000))
model$setParameter(vmid=c(init=1,lb=0,ub=1000))
model$setParameter(vtop=c(init=1,lb=0,ub=1000))
model$setParameter(f1=c(init=1.8483e-04,lb=0,ub=1000))
model$setParameter(f2=c(init=1.8483e-04,lb=0,ub=1000))
model$setParameter(f3=c(init=1.8483e-04,lb=0,ub=1000))
model$setParameter(f4=c(init=1.8483e-04,lb=0,ub=1000))

model$setParameter(vmid8=c(init=1,lb=0,ub=1000))
model$setParameter(vmid9=c(init=1,lb=0,ub=1000))
model$setParameter(vmid10=c(init=1,lb=0,ub=1000))
model$setParameter(vmid11=c(init=1,lb=0,ub=1000))
model$setParameter(vmid12=c(init=1,lb=0,ub=1000))



return(model)
}
