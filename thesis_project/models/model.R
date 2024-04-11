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
  # model$addObs(X10~x10m)
  # model$addObs(X11~x11m)
  # model$addObs(X12~x12m)
  # model$addObs(X13~x13m)
  # model$addObs(X14~x14m)
  # model$addObs(X15~x15m)


  # model$addObs(X0~Y0m)
  # model$addObs(X1~Y1m)
  # model$addObs(X2~Y2m)
  # model$addObs(X3~Y3m)
  # model$addObs(X4~Y4m)
  # model$addObs(X5~Y5m)
  # model$addObs(X6~Y6m)
  # model$addObs(X7~Y7m)
  # model$addObs(X8~Y8m)
  # model$addObs(X9~Y9m)
  # model$addObs(X10~Y10m)
  # model$addObs(X11~Y11m)
  # model$addObs(X12~Y12m)
  # model$addObs(X13~Y13m)
  # model$addObs(X14~Y14m)
  # model$addObs(X15~Y15m)


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
  # model$setVariance(X10 ~ sigma_X^2)
  # model$setVariance(X11 ~ sigma_X^2)
  # model$setVariance(X12 ~ sigma_X^2)
  # model$setVariance(X13 ~ sigma_X^2)
  # model$setVariance(X14 ~ sigma_X^2)
  # model$setVariance(X15 ~ sigma_X^2)




  model$addSystem(dY0m~dt/V0*((Tbot-x0m)*(vbot1+vbot2*(x1m-x0m))*(FbotIn)+(x0m-x1m)*(k*1+f*(FbotOut)+fbot*(Fbot)*(x1m-x0m)))+sigma_y*dwY0)
  model$addSystem(dx0m~dt*(Y0m-x0m+((Tbot-x0m)*(vbot1+vbot2*(x1m-x0m))*(FbotIn)+(x0m-x1m)*(k*1+f*(FbotOut)+fbot*(Fbot)*(x1m-x0m)))/V0)+sigma_x*dw0)

  model$addSystem(dY1m~dt/V1*((x1m-x0m)*(k*1-f*FbotIn)+(x1m-x2m)*(k*2-f*FbotOut)-(x2m-x1m)*(x1m-x0m)*Fbot*v*14)+sigma_y*dwY)
  model$addSystem(dx1m~dt*(Y1m+((x1m-x0m)*(k*1-f*FbotIn)+(x1m-x2m)*(k*2-f*FbotOut)-(x2m-x1m)*(x1m-x0m)*Fbot*v*14)/V1-x1m)+sigma_x*dw1)

  model$addSystem(dY2m~dt/V2*((x2m-x1m)*(k*2-f*FbotIn)+(x2m-x3m)*(k*3-f*FbotOut)-(x3m-x2m)*(x2m-x1m)*Fbot*(v*13))+sigma_y*dwY2)
  model$addSystem(dx2m~dt*(Y2m-x2m+((x2m-x1m)*(k*2-f*FbotIn)+(x2m-x3m)*(k*3-f*FbotOut)-(x3m-x2m)*(x2m-x1m)*Fbot*(v*13))/V2)+sigma_x*dw2)

  model$addSystem(dY3m~dt/V3*((x3m-x2m)*(k*4-f*FbotIn)+(x3m-x4m)*(k*5-f*FbotOut)-(x4m-x3m)*(x3m-x2m)*Fbot*(v*12))+sigma_y*dwY3)
  model$addSystem(dx3m~dt*(Y3m-x3m+((x3m-x2m)*(k*4-f*FbotIn)+(x3m-x4m)*(k*5-f*FbotOut)-(x4m-x3m)*(x3m-x2m)*Fbot*(v*12))/V3)+sigma_x*dw3)

  model$addSystem(dY4m~dt/V4*((x4m-x3m)*(k*5-f*FbotIn)+(x4m-x5m)*(k*6-f*FbotOut)-(x5m-x4m)*(x4m-x3m)*Fbot*(v*11))+sigma_y*dwY4)
  model$addSystem(dx4m~dt*(Y4m-x4m+((x4m-x3m)*(k*5-f*FbotIn)+(x4m-x5m)*(k*6-f*FbotOut)-(x5m-x4m)*(x4m-x3m)*Fbot*(v*11))/V4)+sigma_x*dw4)

  model$addSystem(dY5m~dt/V5*((x5m-x4m)*(k*6-f*FbotIn)+(x5m-x6m)*(k*7-f*FbotOut)-(x6m-x5m)*(x5m-x4m)*Fbot*(v*10))+sigma_y*dwY5)
  model$addSystem(dx5m~dt*(Y5m-x5m+((x5m-x4m)*(k*6-f*FbotIn)+(x5m-x6m)*(k*7-f*FbotOut)-(x6m-x5m)*(x5m-x4m)*Fbot*(v*10))/V5)+sigma_x*dw5)

  model$addSystem(dY6m~dt/V6*((x6m-x5m)*(k*7-f*FbotIn)+(x6m-x7m)*(k*8-f*FbotOut)-(x7m-x6m)*(x6m-x5m)*Fbot*(v*9))+sigma_y*dwY6)
  model$addSystem(dx6m~dt*(Y6m-x6m+((x6m-x5m)*(k*7-f*FbotIn)+(x6m-x7m)*(k*8-f*FbotOut)-(x7m-x6m)*(x6m-x5m)*Fbot*(v*9))/V6)+sigma_x*dw6)

  model$addSystem(dY7m~dt/V7*((x7m-x6m)*(k*8-f*FbotIn)+(x7m-x8m)*(k*9-f*FbotOut)-(x8m-x7m)*(x7m-x6m)*Fbot*(v*8))+sigma_y*dwY7)
  model$addSystem(dx7m~dt*(Y7m-x7m+((x7m-x6m)*(k*8-f*FbotIn)+(x7m-x8m)*(k*9-f*FbotOut)-(x8m-x7m)*(x7m-x6m)*Fbot*(v*8))/V7)+sigma_x*dw7)

  model$addSystem(dY8m~dt/V8*((x8m-x7m)*(k*9-f*FbotIn)+(x8m-x9m)*(k*10-f*FbotOut)-(x9m-x8m)*(x8m-x7m)*Fbot*(v*7))+sigma_y*dwY8)
  model$addSystem(dx8m~dt*(Y8m-x8m+((x8m-x7m)*(k*9-f*FbotIn)+(x8m-x9m)*(k*10-f*FbotOut)-(x9m-x8m)*(x8m-x7m)*Fbot*(v*7))/V8)+sigma_x*dw8)

  model$addSystem(dY9m~dt/V9*((x9m-x8m)*(k*10-f*FbotIn)+(x9m-X10)*(k*11-f*(FbotOut+FmidIn))-(X10-x9m)*(x9m-x8m)*Fbot*(v*6))+sigma_y*dwY9)
  model$addSystem(dx9m~dt*(Y9m-x9m+((x9m-x8m)*(k*10-f*FbotIn)+(x9m-X10)*(k*11-f*(FbotOut+FmidIn))-(X10-x9m)*(x9m-x8m)*Fbot*(v*6))/V9)+sigma_x*dw9)

  # model$addSystem(dY10m~dt/V10*((x10m-x9m)*(k*A10-f*(FbotIn+FmidOut))+(x10m-x11m)*(k*A11-f*(FbotOut+FmidIn))-(x11m-x10m)*(x10m-x9m)*Fbot*(v*5)+(Tmid-x10m)*FmidIn*fmid)+sigma_y*dwY10)
  # model$addSystem(dx10m~dt*(Y10m-x10m+((x10m-x9m)*(k*A10-f*(FbotIn+FmidOut))+(x10m-x11m)*(k*A11-f*(FbotOut+FmidIn))-(x11m-x10m)*(x10m-x9m)*Fbot*(v*5)+(Tmid-x10m)*FmidIn*fmid)/V10)+sigma_x*dw10)

  # model$addSystem(dY11m~dt/V11*((x10m-x11m)*(k*A11+f*(FtopIn+FmidIn))+(x12m-x11m)*(k*A12+f*FtopOut)+(x12m-x11m)*(x11m-x10m)*Ftop*(v*4))+sigma_y*dwY11)
  # model$addSystem(dx11m~dt*(Y11m-x11m+((x10m-x11m)*(k*A11+f*(FtopIn+FmidIn))+(x12m-x11m)*(k*A12+f*FtopOut)+(x12m-x11m)*(x11m-x10m)*Ftop*(v*4))/V11)+sigma_x*dw11)

  # model$addSystem(dY12m~dt/V12*((x11m-x12m)*(k*A12+f*FtopIn)+(x11m-x12m)*(k*A13+f*FtopOut)+(x13m-x12m)*(x12m-x11m)*Ftop*(v*3))+sigma_y*dwY12)
  # model$addSystem(dx12m~dt*(Y12m-x12m+((x11m-x12m)*(k*A12+f*FtopIn)+(x11m-x12m)*(k*A13+f*FtopOut)+(x13m-x12m)*(x12m-x11m)*Ftop*(v*3))/V12)+sigma_x*dw12)

  # model$addSystem(dY13m~dt/V13*((x12m-x13m)*(k*A13+f*FtopIn)+(x12m-x13m)*(k*A14+f*FtopOut)+(x14m-x13m)*(x13m-x12m)*Ftop*(v*2))+sigma_y*dwY13)
  # model$addSystem(dx13m~dt*(Y13m-x13m+((x12m-x13m)*(k*A13+f*FtopIn)+(x12m-x13m)*(k*A14+f*FtopOut)+(x14m-x13m)*(x13m-x12m)*Ftop*(v*2))/V13)+sigma_x*dw13)

  # model$addSystem(dY14m~dt/V14*((x13m-x14m)*(k*A14+f*FtopIn)+(X15-x14m)*(k*A15+f*FtopOut)+(X15-x14m)*(x14m-x13m)*Ftop*(v*1))+sigma_y*dwY14)
  # model$addSystem(dx14m~dt*(Y14m-x14m+((x13m-x14m)*(k*A14+f*FtopIn)+(X15-x14m)*(k*A15+f*FtopOut)+(X15-x14m)*(x14m-x13m)*Ftop*(v*1))/V14)+sigma_x*dw14)

  # model$addSystem(dY15m~dt/V15*((Ttop-x15m)*(vtop1+vtop2*(x15m-x14m))*FtopIn+(x14m-x15m)*(ktop+FtopOut*f+ftop*Ftop*(x14m-x15m)))+sigma_y*dwY15)
  # model$addSystem(dx15m~dt*(Y15m-x15m+((Ttop-x15m)*(vtop1+vtop2*(x15m-x14m))*FtopIn+(x14m-x15m)*(ktop+FtopOut*f+ftop*Ftop*(x14m-x15m)))/V15)+sigma_x*dw15)

  model$addInput('X15','X10')
  model$addInput('Ttop','Ftop','FtopIn','FtopOut','Tmid','Fmid','FmidIn','FmidOut','Tbot','Fbot','FbotIn','FbotOut','FtopVol','FmidVol','FbotVol','ambientTemp')

  # Observation Noise
  model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=1e-2))

  # System Noise
  model$setParameter(sigma_x = c(init=2e-3,lb=0,ub=1))
  model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))


  model$setParameter(k = c(init=1,lb=0,ub=5))
  model$setParameter(f = c(init=1e-1,lb=0,ub=3))
  model$setParameter(v = c(init=1e-5,lb=0,ub=1e-3))

  # model$setParameter(k=c(init=5,lb=0,ub=30))
  # model$setParameter(f=c(init=0.5,lb=0,ub=2))
  # # model$setParameter(v=c(init=7e-3,lb=0,ub=5)) 
  
  model$setParameter(vbot1=c(init=1.35e-1,lb=0,ub=2))
  model$setParameter(vbot2=c(init=1.97e-1,lb=0,ub=2))
  model$setParameter(fbot=c(init=1.94e-2,lb=0,ub=2))

  # model$setParameter(vtop1=c(init=1.35,lb=0,ub=2))
  # model$setParameter(vtop2=c(init=2.19e-2,lb=0,ub=1))
  # model$setParameter(ftop=c(init=-1.94e-2,lb=-1,ub=1))
  # model$setParameter(ktop=c(init=12,lb=1,ub=20))

  # model$setParameter(vmid=c(init=3.17e-1,lb=0,ub=2))
  # model$setParameter(fmid=c(init=3.1745e-01,lb=0,ub=2))
  
  

  return(model)}