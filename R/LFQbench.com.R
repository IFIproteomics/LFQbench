#' Load and install function
#' 
#' Load function upload different libraries and dependencies in the package
#' 
#' @param x A name of a library to be used or uploaded
#' @return This funtion only install the dependecies libraries
#' 
#' @export
#' 

loadLibrary <- function(x) 
{
    if(!require(x, character.only=T, quietly = T)) 
    { 
        install.packages(x); require(x, character.only=T) 
    } 
}

################################################################################

#' + function 
#' 
#' Define + for String concatenation
#' @param x set of String values to be concatenated
#' @return a resulted String 
#' @export
#'

"+"=function(...) UseMethod("+")
"+.default"=.Primitive("+")
"+.character"=function(...) paste(...,sep="")

################################################################################

#' mkdir to create folders
#' 
#'  Create a folder recursivelly if it does not exist
#'  
#'  @param dirName the directory to be created 
#'  @return do not return any information, only create the folder. 
#'  @export
#'  

mkdir <- function(dirName)
{
    if(!file.exists(dirName))
    {
        if(DEBUG) cat("folder '"+dirName+"' does not exist, creating ... ")
        
        dir.create(dirName, recursive=T, showWarnings=F)
        
        if(DEBUG) cat("done!\n")
    }
}

################################################################################

#' createNumericPairs
#' 
#' Create pairs for a range of integer numbers
#' @param from variable from the start number of the range
#' @param to variable end of the number range
#'  
#' @export
#' 

createNumericPairs <- function(from, to) 
{
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

#' root mean square
#' 
#' return the mean square for a number
#' @param x number 
#' @return return the mena square for a number
#' @export
#' 
rms <- function(x) sqrt( sum(x^2) / length(x) )

################################################################################

#' qboxplots for matrix columns or list of vectors
#'
#' qboxplots for matrix columns or list of vectorsusing controllable quantile based whiskers
# 
#' @param vals a dataset as list or matrix of values
#' @param labs a vector of names for each dataset
#' @param lims overrides yLim parameter of boxplot
#' @param whiskerQuantile whiskers range from (this value) to (1 - this value)
#' @return a boxplot
#' @export 

qboxplot <- function( vals, labs=NULL, lims=NULL, whiskerQuantile=0.05, horizontal=F, ... )
{
    if( is.matrix(vals) ){
        cnames = colnames(vals)
        vals = lapply( 1:ncol(vals), function(i) vals[,i] )
        names(vals) = cnames
    }
    
    # remove na
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

################################################################################

#' spiner
#' 
#' generate an spineplot using different paramters. Spineplot are a special cases of mosaic plots, and can be seen as a generalization of stacked (or highlighted) bar plots. Analogously, spinograms are an extension of histograms.
#' 
#' @param samples 
#' @param sampleNames samples names 
#' @param partNames Row names 
#' @param plotTitle plot Title
#' @return spineplot  
#' @export
#' 

spiner <- function(samples, sampleNames=colnames(samples), partNames=rownames(samples), plotTitle="", bgCols=NULL, fgCols="black")
{
    # calculate relative values of a vector
    rel = function(v) v/sum(v)
    # additive middle positions for vector values
    mipo = function(v)
    { 
        c = 1:length(v); 
        for(i in 1:length(v)) c[i] = v[i]/2.0 + sum(v[0:(i-1)]); 
        return(c);
    }
    
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
        text( sx, sy, ""+ round(sr*100, 1) + " %", col=sampleFGCols )
    }
    
    axis(4, labels=partNames, at=sy, tick=F)
}

################################################################################

