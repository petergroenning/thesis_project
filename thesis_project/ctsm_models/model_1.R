make_model <- function(model){
     # Add system equations



    model$addSystem(dx0m ~ 1/C *  ((x1m - x0m) * R) * dt + (Tbot - x0m) * Vbot * FbotIn * dt +  sigma_x * dw0)
    model$addSystem(dx1m ~ 1/C *  ((x0m - x1m) * R + (x2m - x1m) * R ) * dt + sigma_x * dw1)
    model$addSystem(dx2m ~ 1/C *  ((x1m - x2m) * R + (x3m - x2m) * R) * dt + sigma_x * dw2)
    model$addSystem(dx3m ~ 1/C *  ((x2m - x3m) * R + (x4m - x3m) * R) * dt + sigma_x * dw3)
    model$addSystem(dx4m ~ 1/C *  ((x3m - x4m) * R + (x5m - x4m) * R) * dt + sigma_x * dw4)
    model$addSystem(dx5m ~ 1/C *  ((x4m - x5m) * R + (x6m - x5m) * R) * dt + sigma_x * dw5)
    model$addSystem(dx6m ~ 1/C *  ((x5m - x6m) * R + (x7m - x6m) * R) * dt + sigma_x * dw6)
    model$addSystem(dx7m ~ 1/C *  ((x6m - x7m) * R + (x8m - x7m) * R) * dt + sigma_x * dw7)
    model$addSystem(dx8m ~ 1/C *  ((x7m - x8m) * R + (x9m - x8m) * R) * dt + sigma_x * dw8)
    model$addSystem(dx9m ~ 1/C *  ((x8m - x9m) * R + (x10m - x9m) * R) * dt + sigma_x * dw9)
    model$addSystem(dx10m ~ 1/C * ((x9m - x10m) * R + (x11m - x10m) *  R) * dt + sigma_x * dw10)
    model$addSystem(dx11m ~ 1/C * ((x10m - x11m) * R + (x12m - x11m) * R) * dt + sigma_x * dw11)
    model$addSystem(dx12m ~ 1/C * ((x11m - x12m) * R + (x13m - x12m) * R) * dt + sigma_x * dw12)
    model$addSystem(dx13m ~ 1/C * ((x12m - x13m) * R + (x14m - x13m) * R) * dt + sigma_x * dw13)
    model$addSystem(dx14m ~ 1/C * ((x13m - x14m) * R + (x15m - x14m) * R) * dt + sigma_x * dw14)
    model$addSystem(dx15m ~ 1/C* ((x14m - x15m) * R+(ambientTemp-x15m)*Utop) * dt + (Ttop - x15m) * Vtop * FtopIn * dt + sigma_x * dw15)

     
    ###### SET PARAMETERS ########
    # System Noise
    model$setParameter(sigma_x = c(init = 1, lb = 0, ub = 10))

    # System Parameters
    lb <- 1e-5
    ub <- 1000
    init <- 1
    model$setParameter(R = c(init=init, lb = lb, ub = ub))
    model$setParameter(C = c(init=init, lb = lb, ub = ub))
    model$setParameter(Vbot = c(init=init, lb = lb, ub = ub))
    model$setParameter(Vtop =c(init=init, lb = lb, ub = ub))
    model$setParameter(Utop = c(init=init, lb = lb, ub = ub))
    model$setParameter(Cbot = c(init=init, lb = lb, ub = ub))
    model$setParameter(Ctop = c(init=init, lb = lb, ub = ub))

    return(model)
}
