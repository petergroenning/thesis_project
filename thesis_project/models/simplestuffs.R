make_model<-function(model){
 model$addObs(X0~x0m)
#  model$addObs(X1~x1m)
#  model$addObs(X2~x2m)
#  model$addObs(X3~x3m)
#  model$addObs(X4~x4m)
#  model$addObs(X5~x5m)
#  model$addObs(X6~x6m)
#  model$addObs(X7~x7m)
#  model$addObs(X8~x8m)
#  model$addObs(X9~x9m)
    # model$addObs(X10 ~ x10m)
    # model$addObs(X11 ~ x11m)
    # model$addObs(X12 ~ x12m)
    # model$addObs(X13 ~ x13m)
    # model$addObs(X14 ~ x14m)
    # model$addObs(X15 ~ x15m)

    # Set observation equation variances
    model$setVariance(X0 ~ sigma_X^2)
    # model$setVariance(X1 ~ sigma_X^2)
    # model$setVariance(X2 ~ sigma_X^2)
    # model$setVariance(X3 ~ sigma_X^2)
    # model$setVariance(X4 ~ sigma_X^2)
      # model$setVariance(X5 ~ sigma_X^2)
    #  model$setVariance(X6 ~ sigma_X^2)
    # model$setVariance(X7 ~ sigma_X^2)
    # model$setVariance(X8 ~ sigma_X^2)
    # model$setVariance(X9 ~ sigma_X^2)
    # model$setVariance(X10 ~ sigma_X^2)
    #  model$setVariance(X11 ~ sigma_X^2)
    # model$setVariance(X12 ~ sigma_X^2)
    # model$setVariance(X13 ~ sigma_X^2)
    # model$setVariance(X14 ~ sigma_X^2)
    # model$setVariance(X15 ~ sigma_X^2)

    # Observation Noise
    # model$setParameter(sigma_X = obs_std)
    model$setParameter(sigma_X = c(init=2e-4,lb=0,ub=10))

    #Addsystemequations
    
    
    # model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(vbot1+vbot2*(X1-x0m))*FbotIn+(X1-x0m)*(kbot+fbot1*FbotOut+fbot2*x0m*FbotOut)+u*x0m)+sigma_x*dw0)
   #  model$addSystem(dx0m~dt/V0*(FbotIn*(Tbot-x0m)*Vbot*(X1-x0m)+(X1-x0m)*(kbot+fbot*FbotOut))+dwX*sigma_x)
  # model$addSystem(dY~dt/V1*((X0-x1m)*(k1+f1*FbotIn)+(X2-x1m)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v)+sigma_y*dwY)
  # model$addSystem(dx1m~dt*(Y-x1m)+sigma_x*dw1)
  # model$addSystem(dx1m~dt/V1*(((X0-x1m)*(k1+f1*FbotIn)+(X2-x1m)*(k2+f2*FbotOut)-(X2-x1m)*(x1m-X0)*Fbot*v))+sigma_x*dw1)



  # model$addSystem(dx2m~dt/V2*((X1-x2m)*(k1+f1*FbotIn)+(X3-x2m)*(k2+f2*FbotOut)+(X3-x2m)*(x2m-X1)*Fbot*v)+sigma_x*dw2)
  # model$addSystem(dx3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k2+f2*FbotOut)+(X4-x3m)*(x3m-X2)*Fbot*v)+sigma_x*dw3)
  # model$addSystem(dx4m~dt/V4*((X3-x4m)*(k1+f1*FbotIn)+(X5-x4m)*(k2+f2*FbotOut)+(X5-x4m)*(x4m-X3)*Fbot*v)+sigma_x*dw4)
  # model$addSystem(dx5m~dt/V5*((X4-x5m)*(k1+f1*FbotIn)+(X6-x5m)*(k2+f2*FbotOut)+(X6-x5m)*(x5m-X4)*Fbot*v)+sigma_x*dw5)
  # model$addSystem(dx6m~dt/V6*((X5-x6m)*(k1+f1*FbotIn)+(X7-x6m)*(k2+f2*FbotOut)+(X7-x6m)*(x6m-X5)*Fbot*v)+sigma_x*dw6)
  # model$addSystem(dx7m~dt/V7*((X6-x7m)*(k1+f1*FbotIn)+(X8-x7m)*(k2+f2*FbotOut)+(X8-x7m)*(x7m-X6)*Fbot*v)+sigma_x*dw7)
  # model$addSystem(dx8m~dt/V8*((X7-x8m)*(k1+f1*FbotIn)+(X9-x8m)*(k2+f2*FbotOut)+(X9-x8m)*(x8m-X7)*Fbot*v)+sigma_x*dw8)
  # model$addSystem(dx9m~dt/V9*((X8-x9m)*(k1+f1*FbotIn)+(X10-x9m)*(k2+f2*FbotOut)+(X10-x9m)*(x9m-X8)*Fbot*v)+sigma_x*dw9)
   #  model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(Vbot1*(X1-x0m)+Vbot2)*(FbotIn)+(X1-x0m)*(kbot+fbot1*(FbotOut)+(FbotOut)*x0m*fbot2)+a*x0m)+sigma_x*dw0)



  # model$addSystem(dx11m~dt/V11*((X10-x11m)*(k1+f1*FtopOut)+(X12-x11m)*(k2+f2*FtopIn))+sigma_x*dw11)
  # model$addSystem(dx12m~dt/V12*((X11-x12m)*(k1+f1*FtopOut)+(X13-x12m)*(k2+f2*FtopIn))+sigma_x*dw12)
  # model$addSystem(dx13m~dt/V13*((X12-x13m)*(k1+f1*FtopOut)+(X14-x13m)*(k2+f2*FtopIn))+sigma_x*dw13)
  # model$addSystem(dx14m~dt/V14*((X13-x14m)*(k1+f1*FtopOut)+(X15-x14m)*(k2+f2*FtopIn))+sigma_x*dw14)

    # model$addSystem(dx1m~dt/V1*((x1m-x0m)*(FbotIn*Vbot1-kbot)+(X2-x1m)*(FbotOut*Vbot2+kbot)+(X2-x1m)*(x1m-x0m)*Fbot*fbot1+a*(Tbot-x1m)*FbotIn)+sigma_x*dw1)
    # model$addSystem(dx10m~dt/V10*((Tmid-x10m)*(vmid1+(X11-x10m)*vmid3)*(FmidIn)+(X9-x10m)*(k1+f1*FbotIn)+(X11-x10m)*(k2+f2*FbotOut)+(X11-x10m)*(x10m-X9)*v*Fbot)+sigma_x*dw10)

  #  model$addSystem(dx5m~dt/V5*((X4-x5m)*(k1+f1*(FbotIn+a*dFbotIn))+(X6-x5m)*(k1+f1*(FbotOut+b*dFbotOut))+(X6-x5m)*(x5m-X4)*f2*(Fbot+c*dFbot))+sigma_x*dw5)
   # model$addSystem(dx3m~dt/V3*((X2-x3m)*(k1+f1*FbotIn)+(X4-x3m)*(k1+f1*FbotOut)+(X4-x3m)*(x3m-X2)*f2*Fbot)+sigma_x*dw3)

  # model$addSystem(dY~dt/V15*(((Ttop)-x15m)*(vtop1+(x15m-X14)*vtop2)*(FtopIn)+(X14-x15m)*(ktop+(FtopOut)*ftop1+FtopOut*ftop2*x15m)+(ambientTemp-x15m)*utop1+utop2*x15m)+sigma_y*dwY)
  # model$addSystem(dx15m~dt*(Y-x15m+((Ttop)-x15m)*vtop1*FtopIn/V15)+sigma_x*dw15)

  # model$addSystem(dx15m~dt/V15*(((Ttop)-x15m)*(vtop1+(x15m-X14)*vtop2)*(FtopIn)+(X14-x15m)*(ktop+(FtopOut)*ftop1+FtopOut*ftop2*x15m)+(ambientTemp-x15m)*utop1+utop2*x15m)+sigma_x*dw15)


   model$addSystem(dY~dt/V0*((Tbot-x0m)*(vbot1+vbot2*(X1-x0m))*(FbotIn)+(X1-x0m)*(kbot1+kbot2*(x0m-X1)+fbot1*(FbotOut)+fbot2*(Fbot)*(X1-x0m)))+sigma_y*dw0)
   model$addSystem(dx0m~dt*((Y-x0m))+sigma_x*dw1)


  #  model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(vbot1+vbot2*(X1-x0m))*(FbotIn)+(X1-x0m)*(kbot1+kbot2*(X1-x0m)+fbot1*(FbotOut)+fbot2*(Fbot)*(X1-x0m)))+sigma_x*dw0)


   # model$addSystem(dx0m~dt/V0*((Tbot-x0m)*(vbot1+vbot2*(X1-x0m))*(FbotIn)+(X1-x0m)*(kbot+fbot1*(FbotOut)+FbotOut*(x0m-X1)*fbot2)+ubot*x0m)+sigma_x*dw0)

    # model$addSystem(dZ~dt*(x0m+Z)*K+sigma_z*dwz)
    # model$setParameter(Z = c(init=8,lb=0,ub=10))
    # model$setParameter(sigma_z = c(init=2e-4,lb=0,ub=1))
    # model$setParameter(K = c(init=1e-3,lb=0,ub=1e-1))

    # model$addSystem(dY~dwy*sigma_y)
    model$setParameter(sigma_y = c(init=1e-2,lb=0,ub=1))
    model$setParameter(Y = c(init=10,lb=5,ub=20))

    model$addInput('dFtopIn','dTtop','dTbot','dFbotIn', 'dFbotOut', 'FbotInk1', 'Tbotk1')
    model$addInput("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14")
    #####SETPARAMETERS########
    #SystemNoise
  model$setParameter(sigma_x=c(init=2e-2,lb=0,ub=10))
  model$setParameter(k1 = c(init=-10,lb=-20,ub=0))
  model$setParameter(k2 = c(init=-15,lb=-20,ub=0))
  model$setParameter(f1 = c(init=5e-1,lb=0,ub=2))
  model$setParameter(f2 = c(init=5e-1,lb=0,ub=2))
  model$setParameter(v = c(init=1e-1,lb=0,ub=2))
  model$setParameter(vbot1 = c(init=1e-1,lb=1e-2,ub=1))
  model$setParameter(vbot2 = c(init=1e-1,lb=1e-2,ub=1))
  model$setParameter(vmid1 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vmid2 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(vtop1 = c(init=1,lb=0,ub=2))
  model$setParameter(vtop2 = c(init=1e-1,lb=0,ub=1))
  model$setParameter(a=c(init=-1e-3,lb=-1,ub=1))
  model$setParameter(b=c(init=3.7e-2,lb=-10,ub=5e-1))
  model$setParameter(c=c(init=3e-2,lb=-1,ub=1))
  model$setParameter(ubot=c(init=-1e-1,lb=-1,ub=0))
  model$setParameter(utop1=c(init=1e-4,lb=0,ub=1e-2))
  model$setParameter(utop2=c(init=-1,lb=-10,ub=0))
  model$setParameter(fbot1=c(init=1e-1,lb=0,ub=1))
  model$setParameter(fbot2=c(init=1e-3,lb=0,ub=5e-1))
  model$setParameter(ftop1=c(init=2,lb=0,ub=10))
  model$setParameter(ftop2=c(init=-2e-2,lb=-1,ub=0))
  model$setParameter(kbot1=c(init=-1.2,lb=-10,ub=0))
    model$setParameter(kbot2=c(init=-1.2,lb=-10,ub=0))

  model$setParameter(ktop=c(init=-5e-4,lb=-1,ub=1))
    return(model)}