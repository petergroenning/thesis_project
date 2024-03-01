make_model <- function(model){
     # Add system equations

    # model$addSystem(dx0m ~ dt/(V0 * Cp) * ((x1m - x0m) * A1 * k1) + dt * ((Tbot - x0m) * FbotIn/V0 * vbot + (x1m - x0m) * (FbotOut) * 1/V0 * f3) + sigma_x * dw0)
    model$addSystem(dx0m ~ 1/V0*((x1m - x0m)*(Cp*A1*k1+FbotOut*f3) + (Tbot - x0m)*FbotIn*vbot)*dt + sigma_x * dw0)
    # model$addSystem(dx1m ~ dt/(V1 * Cp) * ((x0m - x1m) * A1 * k1 + (x2m - x1m) * A2 * k1) + dt * f3 * 1/V1 * ((x1m - Tbot) * (-FbotIn) + (x2m - x0m) * FbotOut)/2 + sigma_x * dw1)
    model$addS


    model$addSystem(dx2m ~ dt/(V2 * Cp) * ((x1m - x2m) * A2 * k1 + (x3m - x2m) * A3 * k1) + dt * f3 * 1/V2 * ((x2m - Tbot) * (-FbotIn) + (x3m - x1m) * FbotOut)/2 + sigma_x * dw2)
    model$addSystem(dx3m ~ dt/(V3 * Cp) * ((x2m - x3m) * A3 * k1 + (x4m - x3m) * A4 * k1) + dt * f3 * 1/V3 * ((x3m - Tbot) * (-FbotIn) + (x4m - x2m) * FbotOut)/2 + sigma_x * dw3)
    model$addSystem(dx4m ~ dt/(V4 * Cp) * ((x3m - x4m) * A4 * k1 + (x5m - x4m) * A5 * k1) + dt * f3 * 1/V4 * ((x4m - Tbot) * (-FbotIn) + (x5m - x3m) * FbotOut)/2 + sigma_x * dw4)
    model$addSystem(dx5m ~ dt/(V5 * Cp) * ((x4m - x5m) * A5 * k1 + (x6m - x5m) * A6 * k1) + dt * f3 * 1/V5 * ((x5m - Tbot) * (-FbotIn) + (x6m - x4m) * FbotOut)/2 + sigma_x * dw5)
    model$addSystem(dx6m ~ dt/(V6 * Cp) * ((x5m - x6m) * A6 * k1 + (x7m - x6m) * A7 * k1) + dt * (x7m - x5m) /2* (-Fbot) * f1 * 1/V6 + sigma_x * dw6)
    model$addSystem(dx7m ~ dt/(V7 * Cp) * ((x6m - x7m) * A7 * k1 + (x8m - x7m) * A8 * k1) + dt * (x8m - x6m) /2* (-Fbot) * f1 * 1/V7 + sigma_x * dw7)
    model$addSystem(dx8m ~ dt/(V8 * Cp) * ((x7m - x8m) * A8 * k1 + (x9m - x8m) * A9 * k1) + dt * (x9m - x7m) /2* (-Fbot) * f1 * 1/V8 + sigma_x * dw8)
    model$addSystem(dx9m ~ dt/(V9 * Cp) * ((x8m - x9m) * A9 * k1 + (x10m - x9m) * A10 * k1 ) + dt * (x10m - x8m)/2 * (-Fbot) * f1 * 1/V9 + dt * (Tmid - x9m) * FmidIn/V9 * vmid9 + sigma_x * dw9)
    model$addSystem(dx10m ~ dt/(V10 * Cp) * ((x9m - x10m) * A10 * k1 + (x11m - x10m) * A11 * k1 ) + dt *  (x11m - x9m)/2 * Ftop * f2 * 1/V10 + dt * (Tmid - x10m) * FmidIn/V10 * vmid10 + sigma_x * dw10)
    model$addSystem(dx11m ~ dt/(V11 * Cp) * ((x10m - x11m) * A11 * k1 + (x12m - x11m) * A12 * k1) + dt * (x12m - x10m)/2 * Ftop * f2 * 1/V11 + dt * (Tmid - x11m) * FmidIn/V11 * vmid11+ sigma_x * dw11)
    model$addSystem(dx12m ~ dt/(V12 * Cp) * ((x11m - x12m) * A12 * k1 + (x13m - x12m) * A13 * k1) + dt * (x13m - x11m)/2 * Ftop * f2 * 1/V12 + dt * (Tmid - x12m) * FmidIn/V12 * vmid12 +sigma_x * dw12)
    model$addSystem(dx13m ~ dt/(V13 * Cp) * ((x12m - x13m) * A13 * k1 + (x14m - x13m) * A14 * k1) + dt * (x14m - x12m)/2 * Ftop * f2 * 1/V13 + sigma_x * dw13)
    model$addSystem(dx14m ~ dt/(V14 * Cp) * ((x13m - x14m) * A14 * k1 + (x15m - x14m) * A15 * k1) + dt * (x15m - x13m)/2 * Ftop * f2 * 1/V14 + sigma_x * dw14)
    model$addSystem(dx15m ~ dt/(V15 * Cp) * ((x14m - x15m) * A15 * k1 + (ambientTemp - x15m) * Atop * k3) + dt * (Ttop - x15m) * FtopIn/V15 * vtop + dt * (x15m - x14m) * (-FtopOut) * f4 * 1 / V15 + sigma_x * dw15)



     
    ###### SET PARAMETERS ########
    # System Noise
    model$setParameter(sigma_x = c(init = 1, lb = 0, ub = 10))

    # System Parameters
    model$setParameter(k1 = c(init = 1, lb = 0, ub = 1e6))
    model$setParameter(k2 = c(init = 1e-7, lb = 0, ub = 1e6))
    model$setParameter(k3 = c(init = 1, lb = 0, ub = 1e6))
    model$setParameter(vbot = c(init = 1, lb = 0, ub = 1e6))
    model$setParameter(vmid = c(init = 1, lb = 0, ub = 1e6))
    model$setParameter(vtop = c(init = 1, lb = 0, ub = 1e6))
    model$setParameter(f1 = c(init = 1, lb =0, ub = 1e6))
    model$setParameter(f2 = c(init = 1, lb =0, ub = 1e6))
    model$setParameter(f3 = c(init = 1, lb =0, ub = 1e6))
    model$setParameter(f4 = c(init = 1, lb =0, ub = 1e6))

    model$setParameter(vmid8 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid9 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid10 = c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid11= c(init=1, lb = 0, ub = 1000))
    model$setParameter(vmid12 = c(init=1, lb = 0, ub = 1000))



    return(model)
}
