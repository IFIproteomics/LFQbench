### Common functions for formatting Software inputs 

## Depends on: fswe.variables.R 

guessExperiment <- function(exp, injections){
    # exp: one of the experiments: i.e. experiments[[1]]
    # injections: list of injections (preferably a unique list!) of the experiment
    
    #remove extensions
    injections <- as.vector(sapply(injections, file_path_sans_ext))
    all(sapply(injections, is.element, exp ))
}

guessInjection <- function(varname, exp){
    injections <- experiments[[exp]]
    as.character(names(which(sapply(injections,  grepl, varname, ignore.case = T))))
}

guessExperiment_wide <- function(exp, varnames){
    # exp: one of the experiments: i.e. experiments[[1]]
    # varnames: column names of the data frame
    any(sapply(exp, grepl, varnames, ignore.case = T))
}




