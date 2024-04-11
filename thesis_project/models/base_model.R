library(ctsmr)

model <- ctsm$new()


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

    model$addInput("Ttop", "Tbot", "Tmid", "FtopIn", "FmidIn", "FbotIn", "FtopOut", "FbotOut","FmidOut", "FtopVol", "FmidVol", "FbotVol", "ambientTemp", "Ftop", "Fbot", "Fmid")


 

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

    # model$setParameter(Y0m = data$X0[1])
    # model$setParameter(Y1m = data$X1[1])
    # model$setParameter(Y2m = data$X2[1])
    # model$setParameter(Y3m = data$X3[1])
    # model$setParameter(Y4m = data$X4[1])
    # model$setParameter(Y5m = data$X5[1])
    # model$setParameter(Y6m = data$X6[1])
    # model$setParameter(Y7m = data$X7[1])
    # model$setParameter(Y8m = data$X8[1])
    # model$setParameter(Y9m = data$X9[1])
    # model$setParameter(Y10m = data$X10[1])
    # model$setParameter(Y11m = data$X11[1])
    # model$setParameter(Y12m = data$X12[1])
    # model$setParameter(Y13m = data$X13[1])
    # model$setParameter(Y14m = data$X14[1])
    # model$setParameter(Y15m = data$X15[1])
    return(model)
}