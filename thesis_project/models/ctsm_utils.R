save_ctsm_fit <- function(fit, name){
    logLik <- round(fit$loglik)
    save(fit, file = paste0("models/", name, "/fit_ll_", logLik, ".RData"))
}

load_best_fit <- function(name){
    load(paste0("models/", name, "/", max(list.files(paste0("models/", name)), na.rm = TRUE)))
    return(fit)
}
get_fit_parameter_value <- function(fit, param, value){
    # make param to string
    param <- as.character(param)
    if (param %in% names(fit$xm)){
        return(fit$xm[[param]][1])
    } else {
        return(value)
    }
}