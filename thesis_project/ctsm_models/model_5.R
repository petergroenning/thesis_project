make_model <- function(model){
     # Add system equations
    model$addSystem(dx0m ~ 1/(c1*b1x0+c2*b2x0+c3*b3x0+c4*b4x0+c5*b5x0) *  ((x1m - x0m) * (R + FbotVol * Cvol0)) * dt + (Tbot - x0m) * FbotIn * Vbot * dt  + sigma_x0 * dw0)
    model$addSystem(dx1m ~ 1/(c1*b1x1+c2*b2x1+c3*b3x1+c4*b4x1+c5*b5x1) *  ((x0m - x1m) * (R + FbotVol * Cvol0) + (x2m - x1m) * (R + FbotVol * Cvol0)) * dt + sigma_x1 * dw1)
    model$addSystem(dx2m ~ 1/(c1*b1x2+c2*b2x2+c3*b3x2+c4*b4x2+c5*b5x2) *  ((x1m - x2m) * (R + FbotVol * Cvol0) + (x3m - x2m) * (R + FbotVol * Cvol0)) * dt + sigma_x2 * dw2)
    model$addSystem(dx3m ~ 1/(c1*b1x3+c2*b2x3+c3*b3x3+c4*b4x3+c5*b5x3) *  ((x2m - x3m) * (R + FbotVol * Cvol0) + (x4m - x3m) * (R + FbotVol * Cvol0)) * dt + sigma_x3 * dw3)
    model$addSystem(dx4m ~ 1/(c1*b1x4+c2*b2x4+c3*b3x4+c4*b4x4+c5*b5x4) *  ((x3m - x4m) * (R + FbotVol * Cvol0) + (x5m - x4m) * (R + FbotVol * Cvol1)) * dt + sigma_x4 * dw4)
    model$addSystem(dx5m ~ 1/(c1*b1x5+c2*b2x5+c3*b3x5+c4*b4x5+c5*b5x5) *  ((x4m - x5m) * (R + FbotVol * Cvol1) + (x6m - x5m) * (R + FbotVol * Cvol1)) * dt + sigma_x5 * dw5)
    model$addSystem(dx6m ~ 1/(c1*b1x6+c2*b2x6+c3*b3x6+c4*b4x6+c5*b5x6) *  ((x5m - x6m) * (R + FbotVol * Cvol1) + (x7m - x6m) * (R + FbotVol * Cvol1)) * dt + sigma_x6 * dw6)
    model$addSystem(dx7m ~ 1/(c1*b1x7+c2*b2x7+c3*b3x7+c4*b4x7+c5*b5x7) *  ((x6m - x7m) * (R + FbotVol * Cvol1) + (x8m - x7m) * (R + FbotVol * Cvol1)) * dt + sigma_x7 * dw7)
    model$addSystem(dx8m ~ 1/(c1*b1x8+c2*b2x8+c3*b3x8+c4*b4x8+c5*b5x8) *  ((x7m - x8m) * (R + FbotVol * Cvol1) + (x9m - x8m) * (R + FbotVol * Cvol2)) * dt + sigma_x8 * dw8)
    model$addSystem(dx9m ~ 1/(c1*b1x9+c2*b2x9+c3*b3x9+c4*b4x9+c5*b5x9) *  ((x8m - x9m) * (R + FbotVol * Cvol2) + (x10m - x9m) * (R + FbotVol * Cvol2)) * dt + sigma_x9 * dw9)
    model$addSystem(dx10m ~ 1/(c1*b1x10+c2*b2x10+c3*b3x10+c4*b4x10+c5*b5x10) * ((x9m - x10m) * (R + FbotVol * Cvol2) + (x11m - x10m) *  (R + FtopVol * Cvol2)) * dt + sigma_x10 * dw10)
    model$addSystem(dx11m ~ 1/(c1*b1x11+c2*b2x11+c3*b3x11+c4*b4x11+c5*b5x11) * ((x10m - x11m) * (R + FtopVol * Cvol2) + (x12m - x11m) * (R + FtopVol * Cvol3)) * dt + sigma_x11 * dw11)
    model$addSystem(dx12m ~ 1/(c1*b1x12+c2*b2x12+c3*b3x12+c4*b4x12+c5*b5x12) * ((x11m - x12m) * (R + FtopVol * Cvol3) + (x13m - x12m) * (R + FtopVol * Cvol3)) * dt + sigma_x12 * dw12)
    model$addSystem(dx13m ~ 1/(c1*b1x13+c2*b2x13+c3*b3x13+c4*b4x13+c5*b5x13) * ((x12m - x13m) * (R + FtopVol * Cvol3) + (x14m - x13m) * (R + FtopVol * Cvol3)) * dt + sigma_x13 * dw13)
    model$addSystem(dx14m ~ 1/(c1*b1x14+c2*b2x14+c3*b3x14+c4*b4x14+c5*b5x14) * ((x13m - x14m) * (R + FtopVol * Cvol3) + (x15m - x14m) * (R + FtopVol * Cvol3)) * dt + sigma_x14 * dw14)
    model$addSystem(dx15m ~ 1/(c1*b1x15+c2*b2x15+c3*b3x15+c4*b4x15+c5*b5x15) * ((x14m - x15m) * (R + FtopVol * Cvol3) + (ambientTemp-x15m)*Utop) * dt + (Ttop - x15m) * FtopIn * Vtop * dt + sigma_x15 * dw15)

    ###### SET PARAMETERS ########
    # System Noise
    model$setParameter(sigma_x = c(init = 1, lb = 0, ub = 2))

    # System Parameters
    lb <- 1e-8
    ub <- 1000
    init <- 1
    model$setParameter(R = c(init=init, lb = lb, ub = ub))
    model$setParameter(C = c(init=init, lb = lb, ub = ub))
    model$setParameter(Vbot = c(init=1e-4, lb = lb, ub = ub))
    model$setParameter(Vbot1 = c(init=1e-4, lb = lb, ub = ub))
    model$setParameter(Vtop =c(init=1e-4, lb = lb, ub = ub))
    model$setParameter(Vtop1 =c(init=1e-4, lb = lb, ub = ub))
    model$setParameter(Utop = c(init=init, lb = lb, ub = ub))
    model$setParameter(Cbot = c(init=init, lb = lb, ub = ub))
    model$setParameter(Ctop = c(init=init, lb = lb, ub = ub))
    model$setParameter(Cvol0 = c(init=1e-3, lb = 1e-10, ub = 1e2))
    model$setParameter(Cvol1 = c(init=1e-3, lb = 1e-10, ub = 1e2))
    model$setParameter(Cvol2 = c(init=1e-3, lb = 1e-10, ub = 1e2))
    model$setParameter(Cvol3 = c(init=1e-3, lb = 1e-10, ub = 1e2))

    # Spline params
    init <- 10
    lb <- -1000
    ub <- 1000
    model$setParameter(c1 = c(init=init, lb = lb, ub = ub))
    model$setParameter(c2 = c(init=init, lb = lb, ub = ub))
    model$setParameter(c3 = c(init=init, lb = lb, ub = ub))
    model$setParameter(c4 = c(init=init, lb = lb, ub = ub))
    model$setParameter(c5 = c(init=init, lb = lb, ub = ub))

    model$setParameter(sigma_x0 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x1 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x2 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x3 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x4 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x5 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x6 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x7 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x8 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x9 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x10 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x11 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x12 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x13 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x14 = c(init = 1e-1, lb = 0, ub = 2))
    model$setParameter(sigma_x15 = c(init = 1e-1, lb = 0, ub = 2))
    


    return(model)
}
