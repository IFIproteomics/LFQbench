### Common functions for formatting Software inputs 

## Depends on: fswe.variables.R 

substrRight <- function(x, n) { substr(x, nchar(x)-n+1, nchar(x)) }
rmlastchars <- function(x, n) { substr(x, 1, nchar(x) - n) }

stopHere = function() {
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt))
    stop()
}


guessSoftwareSource <- function(filename){
    softsource <- NULL
    grepsw <- function(a, filen) { return(substring(filen, 1, nchar(a)) == a) }
    softsource <- software_sources[unlist(lapply(software_sources, grepsw, filename))]
    
    if(length(softsource) == 0) {
        stop("Software source can not be guessed by filename! Review file names: they should start by the software source.")
    }
    return(softsource)    
}


take1stentry <- function(entries){
    first_entry <- strsplit(as.character(entries), "\\||\\/", fixed=F, perl=T)
    first_entry <- first_entry[[1]][4]
    if(substrRight(first_entry, 3)=="/sp") first_entry <- rmlastchars(first_entry, 3)
    return(first_entry)
}

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

guessSpecie <- function(proteinid){
    sp <- names(which(sapply(species, grepl, proteinid)))
    if(length(sp) == 0) sp <- "NA"
    if(length(sp) > 1) sp <- "multiple"
    sp
}

sum_top_n <- function(values, n, minimum = 1){
    # This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
    # choose the top N peptides.
    if(length(which(!is.na(values))) < minimum) {return (NA)}
    if (n > length(values)) n <- length(values)
    sum(sort(values, decreasing=T)[1:n], na.rm=T)
}

avg_top_n <- function(values, n, minimum = 1){
    # This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
    # choose the top N peptides.
    if(length(which(!is.na(values))) < minimum) {return (NA)}
    if (n > length(values)) n <- length(values)
    mean(sort(values, decreasing=T)[1:n], na.rm=T)
}

sumNA <- function(values){ 
    sumv <- sum(values, na.rm=T)
    if(sumv == 0) sumv <- NA
    sumv
}

avgNA <- function(values){ 
    avgv <- mean(values, na.rm=T)
    if(avgv == 0) avgv <- NA
    avgv
}

single_hits <- function(values){
    # choose single hit proteins.
    if(length(which(!is.na(values))) > 1) {return (NA)}
    sum(values)
}

guessSep <- function(filename){
    filename <- gsub(".gz$", "", filename)
    extension <- file_ext(filename)
    if(extension == "tsv") return ("\t")
    else if(extension == "csv") return (",")
    else if(extension == "txt") return ("\t")
    else return(NA)
}

