make_model<-function(model){
  model$addObs(X0~x0m)
  model$addObs(X1~x1m)
  model$addObs(X2~x2m)
  model$addObs(X3~x3m)
  model$addObs(X4~x4m)
  model$addObs(X5~x5m)
  model$addObs(X6~x6m)
  model$addObs(X7~x7m)
  model$addObs(X8~x8m)
  model$addObs(X9~x9m)
  model$addObs(X10~x10m)
  model$addObs(X11~x11m)
  model$addObs(X12~x12m)
  model$addObs(X13~x13m)
  model$addObs(X14~x14m)
  model$addObs(X15~x15m)

  # Set observation equation variances
  model$setVariance(X0 ~ sigma_X^2)
  model$setVariance(X1 ~ sigma_X^2)
  model$setVariance(X2 ~ sigma_X^2)
  model$setVariance(X3 ~ sigma_X^2)
  model$setVariance(X4 ~ sigma_X^2)
  model$setVariance(X5 ~ sigma_X^2)
  model$setVariance(X6 ~ sigma_X^2)
  model$setVariance(X7 ~ sigma_X^2)
  model$setVariance(X8 ~ sigma_X^2)
  model$setVariance(X9 ~ sigma_X^2)
  model$setVariance(X10 ~ sigma_X^2)
  model$setVariance(X11 ~ sigma_X^2)
  model$setVariance(X12 ~ sigma_X^2)
  model$setVariance(X13 ~ sigma_X^2)
  model$setVariance(X14 ~ sigma_X^2)
  model$setVariance(X15 ~ sigma_X^2)

  # Observation Noise
  model$setParameter(sigma_X = c(init=2e-6,lb=0,ub=1))

  # System Noise
  model$setParameter(sigma_x = c(init=2e-2,lb=0,ub=1))

  # ADD SYSTEM EQUATIONS
  # Bottom Input/Output Layer
  model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(vbot+vbot*x0m)*(FbotIn)+(x1m-x0m)*(kbot+fbot1*FbotOut+fbot2*x0m*FbotOut)+ubot*x0m)+sigma_x*dw0)
  
  # Bottom layers
  model$addSystem(dx1m~dt/V1*((x0m-x1m)*(k1+f1*FbotIn)+(x2m-x1m)*(k2+f2*FbotOut)-(x2m-x1m)*(x1m-x0m)*Fbot*v)+sigma_x*dw1)
  model$addSystem(dx2m~dt/V2*((x1m-x2m)*(k1+f1*FbotIn)+(x3m-x2m)*(k2+f2*FbotOut)-(x3m-x2m)*(x2m-x1m)*Fbot*v)+sigma_x*dw2)
  model$addSystem(dx3m~dt/V3*((x2m-x3m)*(k1+f1*FbotIn)+(x4m-x3m)*(k2+f2*FbotOut)-(x4m-x3m)*(x3m-x2m)*Fbot*v)+sigma_x*dw3)
  model$addSystem(dx4m~dt/V4*((x3m-x4m)*(k1+f1*FbotIn)+(x5m-x4m)*(k2+f2*FbotOut)-(x5m-x4m)*(x4m-x3m)*Fbot*v)+sigma_x*dw4)
  model$addSystem(dx5m~dt/V5*((x4m-x5m)*(k1+f1*FbotIn)+(x6m-x5m)*(k2+f2*FbotOut)-(x6m-x5m)*(x5m-x4m)*Fbot*v)+sigma_x*dw5)
  model$addSystem(dx6m~dt/V6*((x5m-x6m)*(k1+f1*FbotIn)+(x7m-x6m)*(k2+f2*FbotOut)-(x7m-x6m)*(x6m-x5m)*Fbot*v)+sigma_x*dw6)
  model$addSystem(dx7m~dt/V7*((x6m-x7m)*(k1+f1*FbotIn)+(x8m-x7m)*(k2+f2*FbotOut)-(x8m-x7m)*(x7m-x6m)*Fbot*v)+sigma_x*dw7)
  model$addSystem(dx8m~dt/V8*((x7m-x8m)*(k1+f1*FbotIn)+(x9m-x8m)*(k2+f2*FbotOut)-(x9m-x8m)*(x8m-x7m)*Fbot*v)+sigma_x*dw8)
  model$addSystem(dx9m~dt/V9*((x8m-x9m)*(k1+f1*FbotIn)+(x10m-x9m)*(k2+f2*FbotOut)-(x10m-x9m)*(x9m-x8m)*Fbot*v)+sigma_x*dw9)


  # Middle Input/Output Layer
  model$addSystem(dx10m~dt/V10*((Tmid-x10m)*(vmid1+vmid2*(x11m-x10m))*FmidIn+(x9m-x10m)*(k1+f1*FbotIn)+(x11m-x10m)*(k2+f2*FbotOut)+(x11m-x10m)*(x10m-x9m)*v*Fbot)+sigma_x*dw10)

  # Middle layers
  model$addSystem(dx11m~dt/V11*((x10m-x11m)*(k1+f1*FbotIn)+(x12m-x11m)*(k1+f1*FbotOut)+(x12m-x11m)*(x11m-x10m)*Ftop*v)+sigma_x*dw11)
  model$addSystem(dx12m~dt/V12*((x11m-x12m)*(k1+f1*FbotIn)+(x13m-x12m)*(k1+f1*FbotOut)+(x13m-x12m)*(x12m-x11m)*Ftop*v)+sigma_x*dw12)
  model$addSystem(dx13m~dt/V13*((x12m-x13m)*(k1+f1*FbotIn)+(x14m-x13m)*(k1+f1*FbotOut)+(x14m-x13m)*(x13m-x12m)*Ftop*v)+sigma_x*dw13)
  model$addSystem(dx14m~dt/V14*((x13m-x14m)*(k1+f1*FbotIn)+(x15m-x14m)*(k1+f1*FbotOut)+(x15m-x14m)*(x14m-x13m)*Ftop*v)+sigma_x*dw14)

  # Top Input/Output Layer  
  model$addSystem(dx15m~dt/V15*(((Ttop+a*dTtop)-x15m)*(vtop1+(x15m-x14m)*vtop2)*(FtopIn+b*dFtopIn)+(x14m-x15m)*(ktop+(FtopOut)*ftop1+FtopOut*ftop2*x15m)+(ambientTemp-x15m)*utop1+utop2*x15m)+sigma_x*dw15)

  # Add Inputs
  model$addInput('dTtop','dFtopIn')

  # Set Parameters
  model$setParameter(k1 = c(init=-10,lb=-20,ub=0))
  model$setParameter(k2 = c(init=-15,lb=-20,ub=0))
  model$setParameter(f1 = c(init=5e-1,lb=0,ub=1))
  model$setParameter(f2 = c(init=5e-1,lb=0,ub=1))
  model$setParameter(v = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vbot = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vmid1 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vmid2 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vtop1 = c(init=1,lb=0,ub=2))
  model$setParameter(vtop2 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(a=c(init=2.4,lb=0,ub=3))
  model$setParameter(b=c(init=3.7e-2,lb=0,ub=5e-1))
  model$setParameter(ubot=c(init=-1e-1,lb=-1,ub=0))
  model$setParameter(utop1=c(init=1e-4,lb=0,ub=1e-2))
  model$setParameter(utop2=c(init=-1,lb=-10,ub=0))
  model$setParameter(fbot1=c(init=1e-1,lb=0,ub=1))
  model$setParameter(fbot2=c(init=1e-3,lb=0,ub=5e-1))
  model$setParameter(ftop1=c(init=2,lb=0,ub=10))
  model$setParameter(ftop2=c(init=-2e-2,lb=-1,ub=0))
  model$setParameter(kbot=c(init=-1.2,lb=-10,ub=0))
  model$setParameter(ktop=c(init=-5e-4,lb=-1,ub=0))


  return(model)}