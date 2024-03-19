makefit <- function(model){
    ## Needs to be runned for simulate to work


    model$AnalyseModel()
    ## Make a fit modelect to pass to simulate.ctsmr
    fit <- list()
    fit$model <- model
    fit$xm <- model$ParameterValues$initial
    names(fit$xm) <- row.names(model$ParameterValues)
    ## Can be nessary such that xm has the right length
    fit$xm <- fit$xm[model$pars]
    ##
    class(fit) <- "ctsmr"
    ##
    return(fit)
}


