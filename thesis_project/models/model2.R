make_model<-function(model){
    # Add observation equations
    model$addObs(X0 ~ x0m)
    model$addObs(X1 ~ x1m)
    model$addObs(X2 ~ x2m)
    model$addObs(X3 ~ x3m)
    model$addObs(X4 ~ x4m)
    model$addObs(X5 ~ x5m)
    model$addObs(X6 ~ x6m)
    model$addObs(X7 ~ x7m)
    model$addObs(X8 ~ x8m)
    model$addObs(X9 ~ x9m)
    model$addObs(X10 ~ x10m)
    model$addObs(X11 ~ x11m)
    model$addObs(X12 ~ x12m)
    model$addObs(X13 ~ x13m)
    model$addObs(X14 ~ x14m)
    model$addObs(X15 ~ x15m)

    # Set observation equation variances
    model$setVariance(X0 ~ sigma_X0^2)
    model$setVariance(X1 ~ sigma_X1^2)
    model$setVariance(X2 ~ sigma_X2^2)
    model$setVariance(X3 ~ sigma_X3^2)
    model$setVariance(X4 ~ sigma_X4^2)
    model$setVariance(X5 ~ sigma_X5^2)
    model$setVariance(X6 ~ sigma_X6^2)
    model$setVariance(X7 ~ sigma_X7^2)
    model$setVariance(X8 ~ sigma_X8^2)
    model$setVariance(X9 ~ sigma_X9^2)
    model$setVariance(X10 ~ sigma_X10^2)
    model$setVariance(X11 ~ sigma_X11^2)
    model$setVariance(X12 ~ sigma_X12^2)
    model$setVariance(X13 ~ sigma_X13^2)
    model$setVariance(X14 ~ sigma_X14^2)
    model$setVariance(X15 ~ sigma_X15^2)

    model$setParameter(sigma_X0 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X1 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X2 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X3 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X4 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X5 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X6 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X7 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X8 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X9 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X10 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X11 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X12 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X13 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X14 = c(init=2e-4,lb=0,ub=10))
    model$setParameter(sigma_X15 = c(init=2e-4,lb=0,ub=10))

    model$addInput('FbotIn','FbotOut','Fbot','Tbot','FmidIn','FmidOut','Fmid','Tmid','FtopIn','FtopOut','Ftop','Ttop')
    #Addsystemequations
    model$addSystem(dY0m~dt/V0*((Tbot-x0m)*(v1bot+v2bot*(x1m-x0m))*(FbotIn)+(x0m-x1m)*(k0-f0*(FbotOut)))+sigma_y0*dwy0)
    model$addSystem(dx0m~dt*((Y0m-x0m)*a0+((Tbot-x0m)*(v1bot+v2bot*(x1m-x0m))*(FbotIn)+(x0m-x1m)*(k0-f0*(FbotOut)))/V0)+sigma_x0*dwx0)  

    model$addSystem(dY1m~dt/V1*((x1m-x0m)*(k1-f1*FbotIn)+(x1m-x2m)*(k1-f1*FbotOut)-(x2m-x1m)*(x1m-x0m)*Fbot*v1)+sigma_y1*dwy1)
    model$addSystem(dx1m~dt*((Y1m-x1m)*a1+((x1m-x0m)*(k1-f1*FbotIn)+(x1m-x2m)*(k1-f1*FbotOut)-(x2m-x1m)*(x1m-x0m)*Fbot*v1)/V1)+sigma_x1*dwx1)

    model$addSystem(dY2m~dt/V2*((x2m-x1m)*(k2-f2*FbotIn)+(x2m-x3m)*(k2-f2*FbotOut)-(x3m-x2m)*(x2m-x1m)*Fbot*v2)+sigma_y2*dwy2)
    model$addSystem(dx2m~dt*((Y2m-x2m)*a2+((x2m-x1m)*(k2-f2*FbotIn)+(x2m-x3m)*(k2-f2*FbotOut)-(x3m-x2m)*(x2m-x1m)*Fbot*v2)/V2)+sigma_x2*dwx2)

    model$addSystem(dY3m~dt/V3*((x3m-x2m)*(k3-f3*FbotIn)+(x3m-x4m)*(k3-f3*FbotOut)-(x4m-x3m)*(x3m-x2m)*Fbot*v3)+sigma_y3*dwy3)
    model$addSystem(dx3m~dt*((Y3m-x3m)*a3+((x3m-x2m)*(k3-f3*FbotIn)+(x3m-x4m)*(k3-f3*FbotOut)-(x4m-x3m)*(x3m-x2m)*Fbot*v3)/V3)+sigma_x3*dwx3)

    model$addSystem(dY4m~dt/V4*((x4m-x3m)*(k4-f4*FbotIn)+(x4m-x5m)*(k4-f4*FbotOut)-(x5m-x4m)*(x4m-x3m)*Fbot*v4)+sigma_y4*dwy4)
    model$addSystem(dx4m~dt*((Y4m-x4m)*a4+((x4m-x3m)*(k4-f4*FbotIn)+(x4m-x5m)*(k4-f4*FbotOut)-(x5m-x4m)*(x4m-x3m)*Fbot*v4)/V4)+sigma_x4*dwx4)

    model$addSystem(dY5m~dt/V5*((x5m-x4m)*(k5-f5*FbotIn)+(x5m-x6m)*(k5-f5*FbotOut)-(x6m-x5m)*(x5m-x4m)*Fbot*v5)+sigma_y5*dwy5)
    model$addSystem(dx5m~dt*((Y5m-x5m)*a5+((x5m-x4m)*(k5-f5*FbotIn)+(x5m-x6m)*(k5-f5*FbotOut)-(x6m-x5m)*(x5m-x4m)*Fbot*v5)/V5)+sigma_x5*dwx5)

    model$addSystem(dY6m~dt/V6*((x6m-x5m)*(k6-f6*FbotIn)+(x6m-x7m)*(k6-f6*FbotOut)-(x7m-x6m)*(x6m-x5m)*Fbot*v6)+sigma_y6*dwy6)
    model$addSystem(dx6m~dt*((Y6m-x6m)*a6+((x6m-x5m)*(k6-f6*FbotIn)+(x6m-x7m)*(k6-f6*FbotOut)-(x7m-x6m)*(x6m-x5m)*Fbot*v6)/V6)+sigma_x6*dwx6)

    model$addSystem(dY7m~dt/V7*((x7m-x6m)*(k7-f7*FbotIn)+(x7m-x8m)*(k7-f7*FbotOut)-(x8m-x7m)*(x7m-x6m)*Fbot*v7)+sigma_y7*dwy7)
    model$addSystem(dx7m~dt*((Y7m-x7m)*a7+((x7m-x6m)*(k7-f7*FbotIn)+(x7m-x8m)*(k7-f7*FbotOut)-(x8m-x7m)*(x7m-x6m)*Fbot*v7)/V7)+sigma_x7*dwx7)

    model$addSystem(dY8m~dt/V8*((x8m-x7m)*(k8-f8*FbotIn)+(x8m-x9m)*(k8-f8*FbotOut)-(x9m-x8m)*(x8m-x7m)*Fbot*v8)+sigma_y8*dwy8)
    model$addSystem(dx8m~dt*((Y8m-x8m)*a8+((x8m-x7m)*(k8-f8*FbotIn)+(x8m-x9m)*(k8-f8*FbotOut)-(x9m-x8m)*(x8m-x7m)*Fbot*v8)/V8)+sigma_x8*dwx8)

    model$addSystem(dY9m~dt/V9*((x9m-x8m)*(k9-f9*FbotIn)+(x9m-x10m)*(k9-f9*FbotOut)-(x10m-x9m)*(x9m-x8m)*Fbot*v9)+sigma_y9*dwy9)
    model$addSystem(dx9m~dt*((Y9m-x9m)*a9+((x9m-x8m)*(k9-f9*FbotIn)+(x9m-x10m)*(k9-f9*FbotOut)-(x10m-x9m)*(x9m-x8m)*Fbot*v9)/V9)+sigma_x9*dwx9)

    model$addSystem(dY10m~dt/V10*((x10m-x9m)*(k10-f10*(FbotIn+FmidOut))+(x10m-x11m)*(k10-f10*(FbotOut+FmidIn))-(x11m-x10m)*(x10m-x9m)*Fbot*v10+(Tmid-x10m)*FmidIn*fmid)+sigma_y10*dwy10)
    model$addSystem(dx10m~dt*((Y10m-x10m)*a10+((x10m-x9m)*(k10-f10*(FbotIn+FmidOut))+(x10m-x11m)*(k10-f10*(FbotOut+FmidIn))-(x11m-x10m)*(x10m-x9m)*Fbot*v10++(Tmid-x10m)*FmidIn*fmid)/V10)+sigma_x10*dw10)
    
    model$addSystem(dY11m~dt/V11*((x11m-x10m)*(k11-f11*(FtopOut))+(x11m-x12m)*(k11-f11*(FtopIn))+(x12m-x11m)*(x11m-x10m)*Ftop*v11)+sigma_y11*dwy11)
    model$addSystem(dx11m~dt*((Y11m-x11m)*a11+((x11m-x10m)*(k11-f11*(FtopOut))+(x11m-x12m)*(k11-f11*(FtopIn))+(x12m-x11m)*(x11m-x10m)*Ftop*v11)/V11)+sigma_x11*dw11)

    model$addSystem(dY12m~dt/V12*((x12m-x11m)*(k12-f12*(FtopOut))+(x12m-x13m)*(k12-f12*(FtopIn))+(x13m-x12m)*(x12m-x11m)*Ftop*v12)+sigma_y12*dwy12)
    model$addSystem(dx12m~dt*((Y12m-x12m)*a12+((x12m-x11m)*(k12-f12*(FtopOut))+(x12m-x13m)*(k12-f12*(FtopIn))+(x13m-x12m)*(x12m-x11m)*Ftop*v12)/V12)+sigma_x12*dw12)

    model$addSystem(dY13m~dt/V13*((x13m-x12m)*(k13-f13*FtopOut)+(x13m-x14m)*(k13-f13*FtopIn)+(x14m-x13m)*(x13m-x12m)*Ftop*v13)+sigma_y13*dwy13)
    model$addSystem(dx13m~dt*((Y13m-x13m)*a13+((x13m-x12m)*(k13-f13*FtopOut)+(x13m-x14m)*(k13-f13*FtopIn)+(x14m-x13m)*(x13m-x12m)*Ftop*v13)/V13)+sigma_x13*dw13)

    model$addSystem(dY14m~dt/V14*((x13m-x14m)*(k14-f14*FtopOut)+(x14m-x15m)*(k14-f14*FtopIn)+(x15m-x14m)*(x14m-x13m)*Ftop*v14)+sigma_y14*dwy14)
    model$addSystem(dx14m~dt*((Y14m-x14m)*a14+((x13m-x14m)*(k14-f14*FtopOut)+(x14m-x15m)*(k14-f14*FtopIn)+(x15m-x14m)*(x14m-x13m)*Ftop*v14)/V14)+sigma_x14*dw14)

    model$addSystem(dY15m~dt/V15*((Ttop-x15m)*(v1top+v2top*(x15m-x14m))*FtopIn+(x14m-x15m)*(k15+f15*FtopOut)+(ambientTemp-x15m)*u)+sigma_y15*dwy15)
    model$addSystem(dx15m~dt*((Y15m-x15m)*a15+((Ttop-x15m)*(v1top+v2top*(x15m-x14m))*FtopIn+(x14m-x15m)*(k15+f15*FtopOut)+(ambientTemp-x15m)*u)/V15)+sigma_x15*dw15)


    # Set parameters
    model$setParameter(sigma_y0 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y1 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y2 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y3 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y4 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y5 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y6 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y7 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y8 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y9 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y10 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y11 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y12 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y13 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y14 = c(init=1e-2,lb=0,ub=1))
    model$setParameter(sigma_y15 = c(init=1e-2,lb=0,ub=1))

    model$setParameter(sigma_x0 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x1 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x2 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x3 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x4 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x5 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x6 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x7 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x8 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x9 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x10 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x11 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x12 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x13 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x14 = c(init=1e-3,lb=0,ub=1))
    model$setParameter(sigma_x15 = c(init=1e-3,lb=0,ub=1))

    model$setParameter(a0 = c(init=1,lb=0,ub=2))
    model$setParameter(a1 = c(init=1,lb=0,ub=2))
    model$setParameter(a2 = c(init=1,lb=0,ub=2))
    model$setParameter(a3 = c(init=1,lb=0,ub=2))
    model$setParameter(a4 = c(init=1,lb=0,ub=2))
    model$setParameter(a5 = c(init=1,lb=0,ub=2))
    model$setParameter(a6 = c(init=1,lb=0,ub=2))
    model$setParameter(a7 = c(init=1,lb=0,ub=2))
    model$setParameter(a8 = c(init=1,lb=0,ub=2))
    model$setParameter(a9 = c(init=1,lb=0,ub=2))
    model$setParameter(a10 = c(init=1,lb=0,ub=2))
    model$setParameter(a11 = c(init=1,lb=0,ub=2))
    model$setParameter(a12 = c(init=1,lb=0,ub=2))
    model$setParameter(a13 = c(init=1,lb=0,ub=2))
    model$setParameter(a14 = c(init=1,lb=0,ub=2))
    model$setParameter(a15 = c(init=1,lb=0,ub=2))

    model$setParameter(v1bot = c(init=0,lb=-4,ub=5))
    model$setParameter(v2bot = c(init=0,lb=-4,ub=5))
    model$setParameter(v1top = c(init=0,lb=-4,ub=5))
    model$setParameter(v2top = c(init=0,lb=-4,ub=5))
    model$setParameter(fmid = c(init=1,lb=0,ub=10))

    model$setParameter(k0 = c(init=10,lb=0,ub=50))
    model$setParameter(f0 = c(init=1,lb=0,ub=10))
    
    model$setParameter(k1 = c(init=10,lb=0,ub=50))
    model$setParameter(f1 = c(init=1,lb=0,ub=10))
    model$setParameter(v1 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k2 = c(init=10,lb=0,ub=50))
    model$setParameter(f2 = c(init=1,lb=0,ub=10))
    model$setParameter(v2 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k3 = c(init=10,lb=0,ub=50))
    model$setParameter(f3 = c(init=1,lb=0,ub=10))
    model$setParameter(v3 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k4 = c(init=10,lb=0,ub=50))
    model$setParameter(f4 = c(init=1,lb=0,ub=10))
    model$setParameter(v4 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k5 = c(init=10,lb=0,ub=50))
    model$setParameter(f5 = c(init=1,lb=0,ub=10))
    model$setParameter(v5 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k6 = c(init=10,lb=0,ub=50))
    model$setParameter(f6 = c(init=1,lb=0,ub=10))
    model$setParameter(v6 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k7 = c(init=10,lb=0,ub=50))
    model$setParameter(f7 = c(init=1,lb=0,ub=10))
    model$setParameter(v7 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k8 = c(init=10,lb=0,ub=50))
    model$setParameter(f8 = c(init=1,lb=0,ub=10))
    model$setParameter(v8 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k9 = c(init=10,lb=0,ub=50))
    model$setParameter(f9 = c(init=1,lb=0,ub=10))
    model$setParameter(v9 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k10 = c(init=10,lb=0,ub=50))
    model$setParameter(f10 = c(init=1,lb=0,ub=10))
    model$setParameter(v10  = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k11 = c(init=10,lb=0,ub=50))
    model$setParameter(f11 = c(init=1,lb=0,ub=10))
    model$setParameter(v11 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k12 = c(init=10,lb=0,ub=50))
    model$setParameter(f12 = c(init=1,lb=0,ub=10))
    model$setParameter(v12 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k13 = c(init=10,lb=0,ub=50))
    model$setParameter(f13 = c(init=1,lb=0,ub=10))
    model$setParameter(v13 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k14 = c(init=10,lb=0,ub=50))
    model$setParameter(f14 = c(init=1,lb=0,ub=10))
    model$setParameter(v14 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(k15 = c(init=10,lb=0,ub=50))
    model$setParameter(f15 = c(init=1,lb=0,ub=10))
    model$setParameter(v15 = c(init=1e-4,lb=0,ub=1))

    model$setParameter(u = c(init=1,lb=0,ub=5))
    return(model)}