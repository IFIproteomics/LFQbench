#' loadLibrary
#' install and load a library
#' @param x the name of the library to load
#' @export
loadLibrary <- function(x) 
{
  if(!require(x, character.only=T, quietly = T)) 
  { 
    install.packages(x)
    require(x, character.only=T) 
  } 
}

# Define + for String concatenation
# "+"=function(...) UseMethod("+")
# "+.default"=.Primitive("+")
# "+.character"=function(...) paste(...,sep="")

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

#' moveFiles
#' 
#'  move all files (with any file extension) to a new directory
#'  @param srcDir the origin folder
#'  @param tarDir the target folder
#'  @param namePattern the file name pattern for files to be moved
#'  @param rmSrc if the origin folder should be removed, it can only be removed if it is empty
#'  @export
moveFiles = function(srcDir, tarDir, namePattern="\\..+", rmSrc=F)
{
    mkdir(tarDir)
    files = list.files( srcDir, pattern=namePattern, include.dirs=F, full.names=F, recursive=F )
    nix = sapply(files, function(f) file.rename( file.path(srcDir, f), file.path(tarDir, f) ) )
    if(rmSrc) file.remove(srcDir)
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
         function(arg) 
             eval(expr = parse(text=arg), envir = globalenv())
  )
}

