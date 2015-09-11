
#' guessSep
#' This function assigns a separator for separated values files readers
#' depending on file extension
#' @param filename 
guessSep <- function(filename){
  filename <- gsub(".gz$", "", filename)
  extension <- file_ext(filename)
  if(extension == "tsv") return ("\t")
  else if(extension == "csv") return (",")
  else if(extension == "txt") return ("\t")
  else return(NA)
}

#' guessSoftwareSource
#' @param filename The file name should start by one of the software sources
#' @param software_sources a list defined at the variables file.
guessSoftwareSource <- function(filename, software_sources, allowNA = FALSE){
  softsource <- NA
  grepsw <- function(a, filen) { return(substring(filen, 1, nchar(a)) == a) }
  softsource <- software_sources[unlist(lapply(software_sources, grepsw, filename))]
  
  if(length(softsource) == 0 & !allowNA) {
    stop("Software source can not be guessed by filename! Review file names: they should start by the software source.")
  }
  return(softsource)    
}

#' stopHere
#' A stop with no error messages
stopHere = function() {
  opt <- options(show.error.messages=FALSE) 
  on.exit(options(opt))
  stop()
}


#' take1stentry
#' 
#' this function takes the first entry of a swissprot database header
#' @param entries char variable with the swissprot header
take1stentry <- function(entries){
  
  substrRight <- function(x, n) { substr(x, nchar(x)-n+1, nchar(x)) }
  rmlastchars <- function(x, n) { substr(x, 1, nchar(x) - n) }
  
  first_entry <- strsplit(as.character(entries), "\\||\\/", fixed=F, perl=T)
  first_entry <- first_entry[[1]][4]
  if(substrRight(first_entry, 3)=="/sp") first_entry <- rmlastchars(first_entry, 3)
  return(first_entry)
}

#' guessOrganism
#' 
#' this function uses Swissprot annotations i.e. QXXXXX_HUMAN
#' to guess what's the species of the entry. The header (proteinid)
#' may have several entries, and thus it can belong to several
#' organisms. In that case, the function returns 'multiple'.
#' 
#' @param proteinid char variable containing the Swissprot header.
#' @param species vector defining the organisms that are part of the experiment. It is defined at fswe.datasets.R
#' 
#' An example of species vector:
#' species <- vector(mode="list", length=3)
#' names(species) <- c("HUMAN", "YEAST", "ECOLI")
#' species[[1]] <- "_HUMAN"
#' species[[2]] <- "_YEAS"
#' species[[3]] <- "_ECOLI"
#' 
guessOrganism <- function(proteinid, species){
  sp <- names(which(sapply(species, grepl, proteinid)))
  if(length(sp) == 0) sp <- "NA"
  if(length(sp) > 1) sp <- "multiple"
  sp
}

#' sum_top_n
#' 
#' This function sums the values of the top N values of a vector 'values'. It also allows 
#' to have a minimum number of values to be summed. If the minimum of values is not met, it
#' returns NA. 
#' 
#' This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
#' choose the top N peptides.
#' 
#' @param values the vector with the values. 
#' @param n number of top N values to be summed. 
#' @param minimum minimum number of valid (not NA) values allowed in order to return a value.
#' 
sum_top_n <- function(values, n, minimum = 1){
  if(length(which(!is.na(values))) < minimum) {return (NA)}
  if (n > length(values)) n <- length(values)
  sum(sort(values, decreasing=T)[1:n], na.rm=T)
}

#' sum_top_n
#' 
#' This function averages the values of the top N values of a vector 'values'. It also allows 
#' to have a minimum number of values to be averaged. If the minimum of values is not met, it
#' returns NA. 
#' 
#' This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
#' choose the top N peptides.
#' 
#' @param values the vector with the values. 
#' @param n number of top N values to be averaged. 
#' @param minimum minimum number of valid (not NA) values allowed in order to return a value.
#' 
avg_top_n <- function(values, n, minimum = 1){
  # This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
  # choose the top N peptides.
  if(length(which(!is.na(values))) < minimum) {return (NA)}
  if (n > length(values)) n <- length(values)
  mean(sort(values, decreasing=T)[1:n], na.rm=T)
}


#' sumNA
#' 
#' sums values substituting NAs by zeroes
#' This is useful for summary functions in dplyr, for example
#' 
#' @param values
sumNA <- function(values){ 
  sumv <- sum(values, na.rm=T)
  if(sumv == 0) sumv <- NA
  sumv
}

#' avgNA
#' 
#' averages values substituting NAs by zeroes
#' This is useful for summary functions in dplyr, for example
#' 
#' @param values
avgNA <- function(values){ 
  avgv <- mean(values, na.rm=T)
  if(avgv == 0) avgv <- NA
  avgv
}

#' single_hits
#' 
#' This function filters vectors by returning NA if the vector has more than one valid element.
#' This is useful for summary functions in dplyr, for example
#' 
#' @param values
single_hits <- function(values){
  # choose single hit proteins.
  if(length(which(!is.na(values))) > 1) {return (NA)}
  sum(values)
}

### Common functions for formatting Software inputs 

## Depends on: fswe.variables.R 

guessExperiment <- function(exp, injections){
    # exp: one of the experiments: i.e. experiments[[1]]
    # injections: list of injections (preferably a unique list!) of the experiment
    
    #remove extensions
    injections <- as.vector(sapply(injections, file_path_sans_ext))
    all(sapply(injections, is.element, exp ))
}

guessInjection <- function(varname, injections){
    as.character(names(which(sapply(injections,  grepl, varname, ignore.case = T))))
}

guessExperiment_wide <- function(exp, varnames){
    # exp: one of the experiments: i.e. experiments[[1]]
    # varnames: column names of the data frame
    any(sapply(exp, grepl, varnames, ignore.case = T))
}

