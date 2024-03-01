setprm <- function(model){
    ## Init state
    model$setParameter(  Tret1 = c(init=30 ,lb=0     ,ub=60 ) ) # C
    model$setParameter(  Tret2 = c(init=30 ,lb=0     ,ub=60 ) ) # C
    model$setParameter(  Tret3 = c(init=30 ,lb=0     ,ub=60 ) ) # C
    model$setParameter(  Ti = c(init=20 ,lb=0     ,ub=25 ) ) # C
    ## Set the initial value for the optimization
    model$setParameter(  Ch = c(init=0.8 ,lb=1E-5  ,ub=500 ) ) # kWh/C
    model$setParameter(  Ci = c(init=1.36,lb=1E-5  ,ub=500 ) ) # kWh/C
    model$setParameter( Rih = c(init=4  ,lb=0.1   ,ub=1E4) )   # C/kW (was 0.6 in paper)
    model$setParameter( Rie = c(init=5  ,lb=1     ,ub=1E4) )   # C/kW
    model$setParameter(  Aw = c(init=8   ,lb=1E-5  ,ub=500 ) ) # m^2
    model$setParameter( p11 = c(init=1   ,lb=-30   ,ub=10 ) )
    model$setParameter( p22 = c(init=1   ,lb=-30   ,ub=10 ) )
    model$setParameter( p33 = c(init=1   ,lb=-30   ,ub=10 ) )
    model$setParameter( p44 = c(init=1   ,lb=-30   ,ub=10 ) )
    model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
    model$setParameter( e22 = c(init=-1  ,lb=-50   ,ub=10 ) )
    model$setParameter( e33 = c(init=-1  ,lb=-50   ,ub=10 ) )
    ## From data report p. 23
    ##60.95+1862.09+713.56+732.93+424.31+451.85+425.92+573.94+778.77
    ## PhMax = c(init=6.024)) # The max heating power
    model$setParameter(  cw = c(init=0.001118)  ) # kWh/(l C) The specific heat capacity of water
    model$setParameter(flowMax = c(init=180)  ) # l/h. What could it be!? Should for this house (cooling 30 C): 6.024 / (0.001118 * 30) = 180
    model$setParameter(slope = c(init=1)) # The slope parameter of the Sigmoid function tur
    invisible(model)
}
