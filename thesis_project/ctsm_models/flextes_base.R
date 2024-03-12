library(ctsmr)

model <- ctsm$new()
    # Add observation equations and variances
    model$addObs(T14 ~ x1m)
    model$addObs(T13 ~ x2m)
    model$addObs(T12 ~ x3m)
    model$addObs(T11 ~ x4m)
    model$addObs(T10 ~ x5m)
    model$addObs(T09 ~ x6m)
    model$addObs(T08 ~ x7m)
    model$addObs(T07 ~ x8m)
    model$addObs(T06 ~ x9m)
    model$addObs(T05 ~ x10m)
    model$addObs(T04 ~ x11m)
    model$addObs(T03 ~ x12m)
    model$addObs(T02 ~ x13m)
    model$addObs(T01 ~ x14m)


    # Set observation equation variances
    model$setVariance(T01 ~ sigma_T^2)
    model$setVariance(T02 ~ sigma_T^2)
    model$setVariance(T03 ~ sigma_T^2)
    model$setVariance(T04 ~ sigma_T^2)
    model$setVariance(T05 ~ sigma_T^2)
    model$setVariance(T06 ~ sigma_T^2)
    model$setVariance(T07 ~ sigma_T^2)
    model$setVariance(T08 ~ sigma_T^2)
    model$setVariance(T09 ~ sigma_T^2)
    model$setVariance(T10 ~ sigma_T^2)
    model$setVariance(T11 ~ sigma_T^2)
    model$setVariance(T12 ~ sigma_T^2)
    model$setVariance(T13 ~ sigma_T^2)
    model$setVariance(T14 ~ sigma_T^2)


    # Observation Noise
    obs_std <- 0.15
    model$setParameter(sigma_T = obs_std)

    model$addInput("Uin", "Uout")


    #### Physical Parameters #####
    # VOLUMES
    model$setParameter(V1 = 744)
    model$setParameter(V2 = 1294)
    model$setParameter(V3 = 1872)
    model$setParameter(V4 = 2481)
    model$setParameter(V5 = 3122)
    model$setParameter(V6 = 3795)
    model$setParameter(V7 = 4499)
    model$setParameter(V8 = 5235)
    model$setParameter(V9 = 6004)
    model$setParameter(V10 = 6804)
    model$setParameter(V11 = 7636)
    model$setParameter(V12 = 8500)
    model$setParameter(V13 = 9396)
    model$setParameter(V14 = 10324)


    # AREAS
    model$setParameter(A1 = 496)
    model$setParameter(A2 = 1024)
    model$setParameter(A3 = 1584)
    model$setParameter(A4 = 2176)
    model$setParameter(A5 = 2800)
    model$setParameter(A6 = 3456)
    model$setParameter(A7 = 4144)
    model$setParameter(A8 = 4864)
    model$setParameter(A9 = 5616)
    model$setParameter(A10 = 6400)
    model$setParameter(A11 = 7216)
    model$setParameter(A12 = 8064)
    model$setParameter(A13 = 8944)
    model$setParameter(A14 = 9856)
    
    # SPECIFIC HEATS
    model$setParameter(Cp = 4148) # kJ / (m^3 * K)
    model$setParameter(Tcharge = 90)
    model$setParameter(Tdischarge = 50)


 

setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)

    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]

    model$ParameterValues$initial <- newparams
    return(model)
}

setInitialState <- function(model, data){
    model$setParameter(x1m = data$T14[1])
    model$setParameter(x2m = data$T13[1])
    model$setParameter(x3m = data$T12[1])
    model$setParameter(x4m = data$T11[1])
    model$setParameter(x5m = data$T10[1])
    model$setParameter(x6m = data$T09[1])
    model$setParameter(x7m = data$T08[1])
    model$setParameter(x8m = data$T07[1])
    model$setParameter(x9m = data$T06[1])
    model$setParameter(x10m = data$T05[1])
    model$setParameter(x11m = data$T04[1])
    model$setParameter(x12m = data$T03[1])
    model$setParameter(x13m = data$T02[1])
    model$setParameter(x14m = data$T01[1])
    return(model)
}