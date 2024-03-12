library(ctsmr)

model <- ctsm$new()
    # Add observation equations and variances
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

    # Observation Noise
    obs_std <- 0.15
    model$setParameter(sigma_X0 = obs_std)
    model$setParameter(sigma_X1 = obs_std)
    model$setParameter(sigma_X2 = obs_std)
    model$setParameter(sigma_X3 = obs_std)
    model$setParameter(sigma_X4 = obs_std)
    model$setParameter(sigma_X5 = obs_std)
    model$setParameter(sigma_X6 = obs_std)
    model$setParameter(sigma_X7 = obs_std)
    model$setParameter(sigma_X8 = obs_std)
    model$setParameter(sigma_X9 = obs_std)
    model$setParameter(sigma_X10 = obs_std)
    model$setParameter(sigma_X11 = obs_std)
    model$setParameter(sigma_X12 = obs_std)
    model$setParameter(sigma_X13 = obs_std)
    model$setParameter(sigma_X14 = obs_std)
    model$setParameter(sigma_X15 = obs_std)


    #### Physical Parameters #####
    # VOLUMES
    model$setParameter(V0 = 808)
    model$setParameter(V1 = 1051)
    model$setParameter(V2 = 1326)
    model$setParameter(V3 = 1633)
    model$setParameter(V4 = 1973)
    model$setParameter(V5 = 2344)
    model$setParameter(V6 = 2747)
    model$setParameter(V7 = 3182)
    model$setParameter(V8 = 3649)
    model$setParameter(V9 = 4149)
    model$setParameter(V10 = 4680)
    model$setParameter(V11 = 5243)
    model$setParameter(V12 = 5838)
    model$setParameter(V13 = 6465)
    model$setParameter(V14 = 7125)
    model$setParameter(V15 = 7816)

    # AREAS
    model$setParameter(A0 = 697)
    model$setParameter(A1 = 924)
    model$setParameter(A2 = 1183)
    model$setParameter(A3 = 1475)
    model$setParameter(A4 = 1798)
    model$setParameter(A5 = 2153)
    model$setParameter(A6 = 2540)
    model$setParameter(A7 = 2959)
    model$setParameter(A8 = 3410)
    model$setParameter(A9 = 3894)
    model$setParameter(A10 = 4409)
    model$setParameter(A11 = 4956)
    model$setParameter(A12 = 5535)
    model$setParameter(A13 = 6147)
    model$setParameter(A14 = 6790)
    model$setParameter(A15 = 7465)
    model$setParameter(Atop = 8172)

    # SIDE AREAS
    model$setParameter(S0 = 114)
    model$setParameter(S1 = 130)
    model$setParameter(S2 = 146)
    model$setParameter(S3 = 162)
    model$setParameter(S4 = 178)
    model$setParameter(S5 = 194)
    model$setParameter(S6 = 210)
    model$setParameter(S7 = 226)
    model$setParameter(S8 = 242)
    model$setParameter(S9 = 258)
    model$setParameter(S10 = 274)
    model$setParameter(S11 = 290)
    model$setParameter(S12 = 306)
    model$setParameter(S13 = 322)
    model$setParameter(S14 = 338)
    model$setParameter(S15 = 354)

    # SPECIFIC HEATS
    model$setParameter(Cp = 4148) # kJ / (m^3 * K)

    # TEMPERATURES
    model$setParameter(Tsoil = 15)

    model$addInput("Ttop", "Tbot", "Tmid", "FtopIn", "FmidIn", "FbotIn", "FtopOut", "FbotOut","FmidOut", "FtopVol", "FmidVol", "FbotVol", "ambientTemp", "Ftop", "Fbot")
    model$addInput("b1x0","b1x1","b1x1","b1x2","b1x3","b1x4","b1x5","b1x6","b1x7","b1x8","b1x9","b1x10","b1x11","b1x12","b1x13","b1x14","b1x15")
    model$addInput("b2x0","b2x1","b2x1","b2x2","b2x3","b2x4","b2x5","b2x6","b2x7","b2x8","b2x9","b2x10","b2x11","b2x12","b2x13","b2x14","b2x15")
    model$addInput("b3x0","b3x1","b3x1","b3x2","b3x3","b3x4","b3x5","b3x6","b3x7","b3x8","b3x9","b3x10","b3x11","b3x12","b3x13","b3x14","b3x15")
    model$addInput("b4x0","b4x1","b4x1","b4x2","b4x3","b4x4","b4x5","b4x6","b4x7","b4x8","b4x9","b4x10","b4x11","b4x12","b4x13","b4x14","b4x15")
    model$addInput("b5x0","b5x1","b5x1","b5x2","b5x3","b5x4","b5x5","b5x6","b5x7","b5x8","b5x9","b5x10","b5x11","b5x12","b5x13","b5x14","b5x15")
    model$addInput("b6x0","b6x1","b6x1","b6x2","b6x3","b6x4","b6x5","b6x6","b6x7","b6x8","b6x9","b6x10","b6x11","b6x12","b6x13","b6x14","b6x15")
    model$addInput("b7x0","b7x1","b7x1","b7x2","b7x3","b7x4","b7x5","b7x6","b7x7","b7x8","b7x9","b7x10","b7x11","b7x12","b7x13","b7x14","b7x15")


 

setParams <- function(model, fit){
    oldparams <- fit$xm[row.names(model$ParameterValues)]
    newparams <- model$ParameterValues$initial
    names(newparams) <- row.names(model$ParameterValues)

    newparams[!is.na(oldparams)] <- oldparams[!is.na(oldparams)]

    model$ParameterValues$initial <- newparams
    return(model)
}

setInitialState <- function(model, data){
    model$setParameter(x0m = data$X0[1])
    model$setParameter(x1m = data$X1[1])
    model$setParameter(x2m = data$X2[1])
    model$setParameter(x3m = data$X3[1])
    model$setParameter(x4m = data$X4[1])
    model$setParameter(x5m = data$X5[1])
    model$setParameter(x6m = data$X6[1])
    model$setParameter(x7m = data$X7[1])
    model$setParameter(x8m = data$X8[1])
    model$setParameter(x9m = data$X9[1])
    model$setParameter(x10m = data$X10[1])
    model$setParameter(x11m = data$X11[1])
    model$setParameter(x12m = data$X12[1])
    model$setParameter(x13m = data$X13[1])
    model$setParameter(x14m = data$X14[1])
    model$setParameter(x15m = data$X15[1])
    return(model)
}