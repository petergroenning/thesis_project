sdeTi <- function(X, rerun = FALSE){
    ## Try if the fit was cached
    cache_val <- cache_load(!rerun & !pg$rerun, X, pg$cache_dir)
    if(class(cache_val) == "ctsmr"){ return(cache_val) }

    ## Generate a new object of class ctsm
    model <- ctsm$new()
    ## Add system equations and thereby also states
    model$addSystem(dTi ~ ( 1/(Ci*Rie)*(Te-Ti) + gA/Ci*Gv + 1/Ci*Qi )*dt + exp(p11)*dw1)
    ## Set the names of the inputs
    model$addInput(Te,Gv,Qi)

    ## Set the observation equation: Ti is the state, yTi is the measured output
    model$addObs(yTi ~ Ti)
    ## Set the variance of the measurement error
    model$setVariance(yTi ~ exp(e11))

    ## Set the initial value (for the optimization) of the value of the state at the starting time point
    model$setParameter(  Ti = c(init=35  ,lb=10    ,ub=45) )
    ## Set the initial value of the parameters for the optimization
    model$setParameter(  Ci = c(init=1E5 ,lb=1E4   ,ub=1E7) )
    model$setParameter( Rie = c(init=0.1 ,lb=1E-5  ,ub=10) )
    model$setParameter(  gA = c(init=1   ,lb=0.001 ,ub=10))
    model$setParameter( p11 = c(init=1   ,lb=-50   ,ub=10) )
    model$setParameter( p22 = c(init=1   ,lb=-50   ,ub=10) )
    model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10) )

    ## Run the parameter optimization
    fit <- model$estimate(data = X, firstorder=pg$firstorder, threads = pg$threads)
    fit$data[[1]] <- X
    fit$Rnames <- c("Rie")

    ## Save the fit, filename is generated above and kept in pg$cache_dir folder
    cache_save(fit, pg$cache_dir, cache_val)

    ## Return
    return(fit)
}
