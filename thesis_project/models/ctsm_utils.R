save_ctsm_fit <- function(fit, name){
    logLik <- round(fit$loglik)
    save(fit, file = paste0("models/", name, "/fit_ll_", logLik, ".RData"))
}

save_TMB_fit <- function(fit, name){
    logLik <- round(-fit$nll)
    save(fit, file = paste0("models/", name, "/fit_ll_", logLik, ".RData"))
}


load_best_fit <- function(name){
    load(paste0("models/", name, "/", max(list.files(paste0("models/", name)), na.rm = TRUE)))
    return(fit)
}
get_CTSM_parameter_value <- function(fit, param, value){
    # make param to string
    param <- as.character(param)
      if (param %in% names(fit$xm)){
            return(fit$xm[[param]][1])
        } else {
            return(value)
        }
    return (value)
}
    

get_TMB_parameter_value <- function(fit, param, value){
    # make param to string
    if (class(fit) == 'ctsmrTMB.fit'){
        param <- as.character(param)
        
        if (param %in% names(fit$par.fixed)){
            value <- fit$par.fixed[[param]][1]
            if (grepl("log", param, fixed = TRUE)){
                
                value <- exp(value)
            }

            return(value)
        } else {
            return(value)
        }
    } else {
        return(value)
    }
}

getParam <- function(param, value, fit){
    if (class(fit) == 'ctsmrTMB.fit'){
        return(get_TMB_parameter_value(fit, param, value))
    } else {
        return(get_CTSM_parameter_value(fit, param, value))
    }
}