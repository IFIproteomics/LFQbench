#' loadLibrary
#' This function load an specific library for the package.
#' @param x The library to be load
#' @export

loadLibrary <- function(x) 
{
  if(!require(x, character.only=T, quietly = T)) 
  { 
    install.packages(x); require(x, character.only=T) 
  } 
}

# +
# This define + for String concatenation
# @export

#"+"=function(...) UseMethod("+")
#"+.default"=.Primitive("+")
#"+.character"=function(...) paste(...,sep="")

#' mkdir
#' 
#'  This recursively create a directory if it does not exist
#'  @param dirName the directory to create the directory
#'  @export
 
mkdir = function(dirName){
	if(!file.exists(dirName))
	{
    	cat(paste("folder ", dirName, " does not exist, creating ... ", sep = ""))
	    dir.create(dirName, recursive=T, showWarnings=F)
    	cat("done!\n")
  	}
}


#' createNumericPairs
#' 
#' This function creates pairs for a range of integer numbers (e.g. for from=1, to=4 the result is rbind(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
#' 
#' @param from
#' @param to
#' @export
 
createNumericPairs = function(from, to) {
  pairs = rbind(from:(from+1))
  for(num1 in from:(to-1)) {
    for(num2 in (num1+1):to) {
      if(num2<=(from+1)) next
      pairs = as.matrix(rbind(pairs, c(num1, num2)))
    }
  }
  rownames(pairs) = 1:nrow(pairs)
  return(pairs)
}
################################################################################

#' rms
#' 
#' This root mean square
#' @param x 
#' @export
#' 
rms = function(x) sqrt( sum(x^2) / length(x) )

#'qboxplot
#' This boxplots for matrix columns or list of vectors using controllable quantile based whiskers
# 
#' @param vals a dataset as list or matrix of values
#' @param labs a vector of names for each dataset
#' @param lims overrides yLim parameter of boxplot
#' @param whiskerQuantile whiskers range from (this value) to (1 - this value)
#' @export 
#' 
qboxplot = function( vals, labs=NULL, lims=NULL, whiskerQuantile=0.05, horizontal=F, ... ){
  if( is.matrix(vals) ){
    cnames = colnames(vals)
    vals = lapply( 1:ncol(vals), function(i) vals[,i] )
    names(vals) = cnames
  }
  
  vals = lapply(vals, function(x) x[!is.na(x)] )
  
  d = as.list( vals )
  nCols = length( d )
  
  bp = boxplot( d, plot=F, na.rm=T )
  
  # fill labels
  if( is.null(labs) ) {
    if( !is.null(colnames(vals)) )
    {
      labs = colnames(vals)
    }
    else
      if( !is.null(names(d)) )
      {
        labs = names(d)
      }
    else
    {
      labs = 1:nCols
    }            
  }
  
  # fill limits
  if( is.null(lims) )
  {
    low = min(unlist(vals), na.rm=T)
    hi = max(unlist(vals), na.rm=T)
    rng = (hi-low)*0.1
    lims = c(low-rng, hi+rng)
  }
  
  for(i in 1:nCols)
  {
    y = unlist( d[i] )
    
    qnt=quantile( y, probs=c( whiskerQuantile, 0.25, 0.50, 0.75, 1-whiskerQuantile ), na.rm=T )
    o = y[ y<qnt[1] | y>qnt[5] ]
    g = rep.int(i, length(o))
    if(i==1)
    {
      bp$stats=qnt
      bp$out=o
      bp$group=g
    }
    else
    {
      bp$stats=cbind( bp$stats, qnt )
      bp$out=c( bp$out, o )
      bp$group=c( bp$group, g )
    }
  }
  
  bxp(bp,
      h = 0.15
      , ylim=lims
      , names=labs
      , horizontal=horizontal
      , ...
  )
}


#' spiner
#' 
#' This function compute the spiner 
#' @param samples
#' @param sampleNames 
#' @param partNames 
#' @param plotTitle 
#' @param bgCols
#' @param fgCols 
#' @export 
#' 
spiner = function(samples, sampleNames=colnames(samples), partNames=rownames(samples), plotTitle="", bgCols=NULL, fgCols="black"){
  # calculate relative values of a vector
  rel = function(v) v/sum(v)
  # additive middle positions for vector values
  mipo = function(v)
  { 
    c = 1:length(v); 
    for(i in 1:length(v)) c[i] = v[i]/2.0 + sum(v[0:(i-1)]); 
    return(c);
  }
  ###########################
  
  nSamples = ncol(samples)
  nParts = nrow(samples)
  
  sampleBGCols = bgCols
  if(is.null(bgCols)) sampleBGCols = gray.colors( nParts );
  
  sampleFGCols = fgCols
  
  # labels rotation
  par(las = 1)
  # margins: bottom, left, top, right
  # par( mar=c(2,2.5,2,6)+0.5 )
  # plot data
  spineplot( 
    t(samples), 
    xlab="", 
    ylab="", 
    main=plotTitle, 
    axes=F, 
    col=sampleBGCols
  )
  
  sampleX = 1:nSamples; sampleX[] = 1; sampleX = mipo( rel(sampleX) );
  
  # write sample names below their data
  axis(1, labels=sampleNames, at=sampleX, tick=F)
  # write left ruler
  axis(2, lwd = 1, lwd.ticks = 1, pos = -0.03)
  
  sy = 1:nParts
  for(si in 1:nSamples)
  {
    sx = 1:nParts
    sx[] = sampleX[si]
    sr = rel( samples[,si] )
    sy = mipo( sr )
    text( sx, sy, paste("", round(sr*100, 1)," %",sep = ""), col=sampleFGCols )
  }
  
  axis(4, labels=partNames, at=sy, tick=F)
}

#' evalCommandLineArguments
#' This function evaluates the commandLine Arguments 
#' @export
evalCommandLineArguments = function(){
  args=commandArgs(trailingOnly = T)
  sapply(args[grep("=", args)], 
         function(arg) eval(expr = parse(text=arg), envir = globalenv())
  ) 
}

#' guessSep
#' This function assigns a separator for separated values files readers
#' depending on file extension
#' @param filename 
#' @export
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
#' @export
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
#' @export
stopHere = function() {
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt))
    stop()
}


#' take1stentry
#' 
#' this function takes the first entry of a swissprot database header
#' @param entries char variable with the swissprot header
#' @export 
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
single_hits <- function(values){
    # choose single hit proteins.
    if(length(which(!is.na(values))) > 1) {return (NA)}
    sum(values)
}

