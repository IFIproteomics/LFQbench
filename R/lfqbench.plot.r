################################################################################
#' LFQbench.plotToFile
#' 
#'  Plot user content to file by using LFQbench canvas configuration.
#'  
#' @param file the pdfFile
#' @param plotFunc the plot function
#' @export
LFQbench.plotToFile = function( file, plotFunc, ... )
{
  pdf(file=file, onefile=T, width=LFQbench.Config$PlotWidth, height=LFQbench.Config$PlotHeight, family="Helvetica", pointsize=9)
  par(LFQbench.Config$par)
  if(is.function(plotFunc)) res = plotFunc(...)
  par(LFQbench.Config$parBackup)
  dev.off()
  if(exists("res")) return(res)
}
################################################################################

################################################################################
# calculate qc function
getQCFunction = function( ratios, ensureValueRange=c(0, 1) )
{
  # make absolute values and sort them in ascending order
  absRatios = sort( abs(ratios) )
  # ensure value range
  if( absRatios[length( absRatios )]<ensureValueRange[2] ) absRatios = c(absRatios, ensureValueRange[2])
  if( absRatios[1]>ensureValueRange[1] ) absRatios = c(ensureValueRange[1], absRatios)
  return( approxfun( absRatios, (1:length( absRatios ))/length( absRatios ), method="linear" ) )
}
################################################################################

################################################################################
#' expandRange
#' expand a range by a relative factor.
#' @param r the range as a vector of two numbers in increasing order
#' @param f the ralative expansion factor, e.g. 
#'          f=0 will not affect the original range,
#'          f=0.2 will expand the range by 20%, t.m. 10% to each side
#' @export
expandRange = function(r=c(0,1), f=0)
{
    margin = abs(r[2] - r[1]) * f / 2
    res = r + c(-margin, margin)
    return(res)
}
################################################################################

################################################################################
#' emptyPlot
#' show an empty plot with the plotting area in given ranges
#' x/y-Axes are plotted on the bottom/left side with given line width
#' grid can be switched off
#' @param xRange=c(0,1)
#' @param yRange=c(0,1)
#' @param lwd=1
#' @param grid=T
#' @param showXlab=T
#' @param showYlab=T
#' @param axes=T
#' @param cex.axis=1
#' @param plotAreaExpansionFactor=0.02 # this will expand the diplayed plot area by 1% to each side
#' @export
emptyPlot = function(xRange=c(0,1), yRange=c(0,1), lwd=1, grid=T, showXlab=T, showYlab=T, axes=T, cex.axis=1, at.x = NULL, at.y = NULL, plotAreaExpansionFactor=0.02)
{
    xRange = expandRange(xRange, plotAreaExpansionFactor)
    yRange = expandRange(yRange, plotAreaExpansionFactor)
	plot( 0, 0, xlim=xRange, ylim=yRange, type="n", axes=F, xlab="", ylab="", main="" )
	if(axes) addAxes(lwd, showXlab, showYlab, cex.axis=cex.axis, at.x = at.x, at.y = at.y)
	if(grid) grid()
}
################################################################################

################################################################################
#' addAxes
#' 
#' add normal axes having full lengths of plotting area
#' 
#' @param lwd=1
#' @param showXlab=T
#' @param showYlab=T
#' @param showXAxis=T
#' @param showYAxis=T
#' @param cex.axis=1
#' @param at.x=NULL 
#' @param at.y=NULL
#' @export
addAxes = function(lwd=1, showXlab=T, showYlab=T, showXAxis=T, showYAxis=T, cex.axis=1, at.x = NULL, at.y = NULL)
{
  usr = par()$usr
  if(showXAxis) lines(x=usr[c(1,2)], y=usr[c(3,3)], lwd=lwd, xpd=T)
  if(showYAxis) lines(x=usr[c(1,1)], y=usr[c(3,4)], lwd=lwd, xpd=T)
  
  if(showXlab) {
      at = NULL
      if(!is.null(LFQbench.Config$AxisXLabelNumDiv)) at = seq(usr[1], usr[2], by = LFQbench.Config$AxisXLabelNumDiv)
      if(!is.null(at.x)) at = at.x
      axis(1, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis, at=at)  
  } 
  if(showYlab){
      at=NULL
      if(!is.null(LFQbench.Config$AxisYLabelNumDiv)) at = seq(usr[3], usr[4], by = LFQbench.Config$AxisYLabelNumDiv)
      if(!is.null(at.y)) at = at.y
      axis(2, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis, at=at)
  } 
      
}
################################################################################

################################################################################
#' addXLab
#' add an x-axis label shifted by marginShift from par()$mgp[1]
#' @param xlab=""
#' @param marginShift=-0.2
#' @param cex.lab=LFQbench.Config$AxisLabelSize
#' @param ... other parameter passed to title(...)
#' @export
addXLab = function(xlab="", marginShift=-0.2, cex.lab = LFQbench.Config$AxisLabelSize, ...)
{
  m = par()$mgp
  m[1] = m[1]+marginShift
  title(xlab=xlab, mgp=m, cex.lab=cex.lab, ...)
}

#' addYLab
#' add an y-axis label shifted by marginShift from par()$mgp[1]
#' @param ylab=""
#' @param marginShift=-0.2
#' @param cex.lab=LFQbench.Config$AxisLabelSize
#' @param ... other parameter passed to title(...)
#' @export
addYLab = function(ylab="", marginShift=0.3, cex.lab = LFQbench.Config$AxisLabelSize, ...)
{
  m = par()$mgp
  m[1] = m[1]+marginShift
  title(ylab=ylab, mgp=m, cex.lab=cex.lab, ...)
}
################################################################################

################################################################################
# add a horizontal line to plot area
hLine = function(y, ...)
{
  lines(x=par()$usr[c(1,2)], y=c(y,y), ...)
}
################################################################################

################################################################################
# add a vertical line to plot area
vLine = function(x, ...)
{
  lines(x=c(x,x), y=par()$usr[c(3,4)], ...)
}
################################################################################

################################################################################
# make a color darker or lighter by scaling
scaleColor=function(col, scale=1.5, alpha=1)
{
  rgbVec = t( col2rgb(col)/255*scale )
  rgbVec[rgbVec>255] = 255
  return( rgb(rgbVec, alpha=alpha) )
}
################################################################################

################################################################################
#' plotSpeciesLegend
#' 
#'  This function add the species legend to a plot
#' @param pos the legend position on plot area
plotSpeciesLegend = function( pos="top", ... )
{
  legend(x=pos, legend=LFQbench.Config$AllSpeciesNames, col=LFQbench.Config$SpeciesColors, lwd=LFQbench.Config$PlotLegendLineWidth, inset=.02, bg="white", ...)	
}
################################################################################

################################################################################
#' LFQbench.plotSpeciesLegends
#' 
#'  plot vertical and horizontal species legends
#' @param plotFolder the folder, where plots should go to
#' @export
LFQbench.plotSpeciesLegends = function(plotFolder = LFQbench.Config$PlotFilesLocation, ...)
{
  pdf(file = paste(plotFolder,"/species_legend_vertical.pdf", sep = ""), width = 1.05, height = 0.66, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=F, ... )
  dev.off()
  pdf(file = paste(plotFolder,"/species_legend_horizontal.pdf",sep = ""), width = 2.9, height = 0.35, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=T, ... )
  dev.off()
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
addScatterPointsForSpecies = function(species, samplePairResult, minAlpha=.2, showExpectationLine=T, showRegressionLine=T, useCfgColor=T, rampColors=T, ...)
{
  ds = samplePairResult$data[[species]]
  theCol = ifelse(useCfgColor, LFQbench.Config$SpeciesColors[which(LFQbench.Config$AllSpeciesNames==ds$species)], ds$col)
  cols=theCol
  if(rampColors)
  {
    # get density for a log-ratio value
    lr2d = approxfun(ds$density$x, ds$density$y, method = "linear")
    # densities for all log-ratios of the species
    densities = lr2d(ds$y)
    # scale densities between 0 and 1
    densities = (densities - min(densities))
    densities = densities / max(densities)
    densities[densities<minAlpha] = minAlpha
    # set alpha for each log-ratio by its scaled density
    cols = alpha(theCol, densities)
  }
  # draw points
  points(ds, pch=20, col=cols, ... )
  # draw expectation line
  if(showExpectationLine) 
  {
    abline(h = ds$expectation, col=scaleColor(theCol, .8), lty="dashed", lwd=LFQbench.Config$PlotCurveLineWidth)
  }
  # draw lowess regression line
  if(showRegressionLine)
  {
    regLine = lowess( ds$x, ds$y )
    lines( regLine, col=scaleColor(theCol,.5), lty=5, lwd=LFQbench.Config$PlotCurveLineWidth )
  }
  return(ds$y)
}
################################################################################

################################################################################
# scale a vector to values ranging from 0 to 1
scaleTo01 = function(v)
{
  v = v - min(v)
  v = v / max(v)
  return( v )
}
################################################################################

################################################################################
addScatterPointsForSpecies2 = function(species, samplePairResult, minAlpha=.2, showExpectationLine=T, showRegressionLine=T, useCfgColor=T, rampColors=T, ...)
{
  ds = samplePairResult$data[[species]]
  theCol = ifelse(useCfgColor, LFQbench.Config$SpeciesColors[which(LFQbench.Config$AllSpeciesNames==ds$species)], ds$col)
  cols=theCol
  if(rampColors)
  {
    # get density for log2(B)
    xDens = density( ds$x )
    x2d = approxfun(xDens$x, xDens$y, method = "linear")
    # normalized densities for log2-intensities
    d4x = scaleTo01( x2d( ds$x ) )
      
    # get density for a log-ratio value
    y2d = approxfun(ds$density$x, ds$density$y, method = "linear")
    # normalized densities for log-ratios
    d4y = scaleTo01( y2d(ds$y) )

    # calculate 2d density
    densities = d4x * d4y 
  
    # limit by min threshold
    densities[ densities<minAlpha ] = minAlpha
    
    # set alpha for each log-ratio by its scaled density
    cols = alpha(theCol, densities)
  }
  # draw points
  points(ds, pch=20, col=cols, ... )
  # draw expectation line
  if(showExpectationLine) 
  {
    abline(h = ds$expectation, col=scaleColor(theCol, .8), lty="dashed", lwd=LFQbench.Config$PlotCurveLineWidth)
  }
  # draw lowess regression line
  if(showRegressionLine)
  {
    regLine = lowess( ds$x, ds$y, f=LFQbench.Config$ScatterPlotLowessBandWidth )
    lines( regLine, col=scaleColor(theCol,.5), lty=5, lwd=LFQbench.Config$PlotCurveLineWidth )
  }
  return(ds$y)
}
################################################################################

################################################################################
#' makeScatter
#' 
#' make a scatter plot with smoothed colors
#' @param samplePair the result set data of a sample pair
#' @param showLegend display the legend
#' @param showRegLines display regression curves
#' @param showExpLines display expectation lines
#' @param useCfgColor ignore result set colors and use colors from current configuration settings
makeScatter = function( samplePair, showLegend=F, showRegLines=F, showExpLines=T, useCfgColor=T, at.x = NULL, at.y = NULL) {
        
  emptyPlot(samplePair$xlim, samplePair$ylim, grid=F, 
            lwd = LFQbench.Config$AxisLineThickness, 
            cex.axis = LFQbench.Config$AxisAnnotationSize, 
            at.x = at.x, at.y = at.y)

  if( !any(names(LFQbench.Config)=="DisplayAxisLabels") ) LFQbench.Config[["DisplayAxisLabels"]] <<- TRUE
  if(LFQbench.Config$DisplayAxisLabels) {
      addXLab( as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ) )
      addYLab( as.expression( bquote( Log[2]~"("~.(paste(samplePair$name1, ":", samplePair$name2, sep=""))~")" ) ) )
  }
  logRatios = sapply(LFQbench.Config$AllSpeciesNames, 
                     addScatterPointsForSpecies2, samplePair, 
                      minAlpha=LFQbench.Config$ScatterPlotPointMinAlpha, showExpectationLine=showExpLines, showRegressionLine=showRegLines, 
                      useCfgColor=T, rampColors=ifelse(LFQbench.Config$ScatterPlotPointMinAlpha==1, F, T), cex=LFQbench.Config$PlotPointSize
                        )
  if(showLegend) plotSpeciesLegend( horiz=T )
  return(logRatios)
}
################################################################################

################################################################################
#' LFQbench.showScatterPlot
#'
#' draw scatterplot
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @param showLegend display the legend
#' @param showRegLines display regression curves
#' @export
LFQbench.showScatterPlot = function( samplePair, showLegend=F, showRegLines=F, at.x = NULL, at.y = NULL )
{
  par(LFQbench.Config$par)
  logRatios = makeScatter(samplePair, showLegend, showRegLines, at.x = at.x, at.y = at.y )
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.showLogRatioBoxPlot
#'
#' draw log ratio quiantile based boxplot
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @export
LFQbench.showLogRatioBoxPlot = function( samplePair )
{
  par(LFQbench.Config$par)
  logRatios = sapply( samplePair$data, function(d) { return(d$y) } )
  qboxplot( logRatios, pch=20, ylab="",
            whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F, lims=LFQbench.Config$LogRatioPlotRange,
      cex.lab = LFQbench.Config$AxisLabelSize
    )
  addAxes(showXlab=F, lwd = LFQbench.Config$AxisLineThickness, cex.axis = LFQbench.Config$AxisAnnotationSize)
  axis(1, labels=LFQbench.Config$AllSpeciesNames, at=(1:LFQbench.Config$NumberOfSpecies), 
       lwd.ticks=LFQbench.Config$AxisLineThickness, lwd=0, cex.axis = LFQbench.Config$AxisAnnotationSize )
  addYLab( as.expression( bquote( Log[2]~"("~.(paste(samplePair$name1, ":", samplePair$name2, sep=""))~")" ) ))
  grid()
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.showQuantBarPlot
#'
#' draw barplot showing quantification bars
#' @param samplePair the result set data of a hybrid proteome sample pair
LFQbench.showQuantBarPlot = function( samplePair )
{
  par(LFQbench.Config$par)
  values = unlist(sapply( samplePair$data, function(d) { return(d$y) } ))
  colors = unlist(sapply( samplePair$data, function(d) { return(rep(d$col,length(d$y)) ) } ))
  cat(paste("",length(values)," values ...\n", sep=""))
  sortedIndices = order(values, decreasing=T)
  colors = colors[sortedIndices]
  values = values[sortedIndices]
  barplot(values, col=colors, border=NA,
        ylim=samplePair$ylim,
        main="", names.arg=NA, xlab="", ylab="",
        cex.lab = LFQbench.Config$AxisLabelSize, cex.axis = LFQbench.Config$AxisAnnotationSize
  )
  addXLab( as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ) )
  addYLab( as.expression( bquote( Log[2]~"("~.(paste(samplePair$name1, ":", samplePair$name2, sep=""))~")" ) ) )
  legendLabels = sapply( samplePair$data, function(d) d$species )
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="topright", legend=legendLabels, col=legendColors, lwd=LFQbench.Config$PlotLegendLineWidth, inset=.02, bg="white")
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.showScatterAndDensityPlot
#'
#' draw scatter plot with an attached density plot
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @param showLegend display the legend
#' @param showRegLines display regression curves
#' @param scatterPlotWidth the window portion used for scatter plot
#' @export
LFQbench.showScatterAndDensityPlot = function(samplePair, showLegend=F, showRegLines=F, scatterPlotWidth=0.8, at.x = NULL, at.y = NULL)
{
	par(LFQbench.Config$par)
	par(fig=c(0,scatterPlotWidth,0,1), new=F)
	logRatios = makeScatter(samplePair, showLegend, showRegLines, at.x = at.x, at.y = at.y)
	# density plot
	par(fig=c(scatterPlotWidth,1,0,1), new=T)
	pm = par()$mar
	pm[c(2,4)] = 0
	par(mar=pm)
	xLim=samplePair$ylim
	yLim=range( lapply(samplePair$data, function(d) range(d$density$y) ) )
	emptyPlot(yLim, xLim, axes=F, grid=F)
  mkLines = function(d)
  {
    idx = d$density$x > min(d$y) & d$density$x < max(d$y)
    theCol = LFQbench.Config$SpeciesColors[which(LFQbench.Config$AllSpeciesNames==d$species)]
    lines(d$density$y[idx], d$density$x[idx] , type="l", col=theCol, lwd=LFQbench.Config$PlotCurveLineWidth ) 
    hLine(y=d$expectation, lty=2, lwd=LFQbench.Config$PlotCurveLineWidth, col=scaleColor(theCol,.8) )
  }
	nix = sapply( samplePair$data, mkLines )
	par(LFQbench.Config$parBackup)
	par(fig=c(0,1,0,1), new=F)
}
################################################################################

################################################################################
#' LFQbench.showScatterAndBoxPlot
#'
#' draw scatter plot with an attached box plot
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @param showLegend display the legend
#' @param showRegLines display regression curves
#' @param scatterPlotWidth the window portion used for scatter plot
#' @export
LFQbench.showScatterAndBoxPlot = function(samplePair, showLegend=F, showRegLines=F, scatterPlotWidth=0.8, at.x = NULL, at.y = NULL)
{
  par(LFQbench.Config$par)
  par(fig=c(0,scatterPlotWidth,0,1), new=F)
  logRatios = makeScatter(samplePair, showLegend, showRegLines, at.x = at.x, at.y = at.y)
  # box plot
  par(fig=c(scatterPlotWidth,1,0,1), new=T)
  pm = par()$mar
  pm[c(2,4)] = 0
  par(mar=pm)
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"("~.(paste(samplePair$name1, ":", samplePair$name2, sep=""))~")" ) ),
            whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F, lims=samplePair$ylim, border=LFQbench.Config$SpeciesColors
  )
  # add horizontal lines
  makeExpLines = function(d) 
  {
    theCol = LFQbench.Config$SpeciesColors[which(LFQbench.Config$AllSpeciesNames==d$species)]
    hLine(y=d$expectation, lty=2, lwd=LFQbench.Config$PlotCurveLineWidth, col=scaleColor(theCol,.8) )
    return(d$expectation) 
  }  
  expRatios = sapply( samplePair$data, makeExpLines )
  par(LFQbench.Config$parBackup)
  par(fig=c(0,1,0,1), new=F)
}
################################################################################

################################################################################
#' LFQbench.showDistributionDensityPlot
#'
#' draw log-ratio distribution kernel density plot
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @param showLegend display the legend
#' @export
LFQbench.showDistributionDensityPlot = function(samplePair, showLegend=F)
{
  xLim=samplePair$ylim
  yLim=range( lapply(samplePair$data, function(d) range(d$density$y) ) )
  par(LFQbench.Config$par)
  emptyPlot(xLim, yLim, lwd = LFQbench.Config$AxisLineThickness, cex.axis = LFQbench.Config$AxisAnnotationSize )
  addXLab( as.expression( bquote( Log[2]~"("~.(paste(samplePair$name1, ":", samplePair$name2, sep=""))~")" ) ) )
  addYLab( "Density", .6 )
  sapply(samplePair$data, function(d) lines(d$density, type="l", col=d$col, lwd=LFQbench.Config$PlotCurveLineWidth ) )
  if(showLegend) plotSpeciesLegend( )
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.showAllSpeciesQC
#'
#' draw quantification curves
#' @param samplePair the result set data of a hybrid proteome sample pair
#' @export
LFQbench.showAllSpeciesQC = function(samplePair)
{
  par(LFQbench.Config$par)
  emptyPlot(xRange=samplePair$qcrange)
  x = (1:500)/500*samplePair$qcrange[2]
  title(
    main="",
    xlab=as.expression( bquote( "Absolute"~log[2]~"("~.(samplePair$name1)~":"~.(samplePair$name2)~")" ) ),
    ylab="Parts"
  )
  sapply(samplePair$data, function(d) lines(x, d$qcfunction(x), type="l", col=d$col, lwd=LFQbench.Config$PlotCurveLineWidth ) )
  legendLabels = sapply( samplePair$data, function(d) paste(d$species , ": " , d$auqc, sep="" ))
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=LFQbench.Config$PlotLegendLineWidth, inset=.02, bg="white")
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.showSingleSpeciesQC
#'
#' draw quantification curves
#' @param resultSet the result set data
#' @param speciesName the name of a species
#' @export
LFQbench.showSingleSpeciesQC = function(resultSet, speciesName)
{
  speciesIndex = which( LFQbench.Config$AllSpeciesNames == speciesName )
  par( LFQbench.Config$par )
  emptyPlot( xRange=LFQbench.Config$AUQCRatioRange )
  title(
    main="",
    xlab=as.expression( bquote( Absolute~log[2]-ratio~"for"~.(speciesName) ) ),
    ylab="Parts"
  )
  x = (1:500)/500*LFQbench.Config$AUQCRatioRange[2]
  sapply(rs$result, 
         function(sp){ 
             lines(x, sp$data[[speciesIndex]]$qcfunction(x), type="l", col=sp$col, lwd=LFQbench.Config$PlotCurveLineWidth ) 
        } )
  legendLabels = sapply( rs$result, function(sp) paste("AUQC(",sp$name1, ":", sp$name2, "): ", sp$data[[speciesIndex]]$auqc, sep = "" ))
  legendColors = sapply( rs$result, function(sp) sp$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=LFQbench.Config$PlotLegendLineWidth, inset=.02, bg="white")
  par( LFQbench.Config$parBackup )
}
################################################################################

################################################################################
#' LFQbench.plotProteinDispersionBySpecies
#'
#' draw protein dispersion as cv between replicates by species
#' @param d the data of a result set (resultSet$data)
#' @export
LFQbench.plotProteinDispersionBySpecies = function( d )
{
  d4s = sapply( LFQbench.Config$AllSpeciesNames, function(sp) as.vector( na.exclude( d$cv[ d$species==sp ] ) ) )
  par(LFQbench.Config$par)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=LFQbench.Config$AllSpeciesNames, at=(1:LFQbench.Config$NumberOfSpecies), lwd.ticks=1, lwd=0)
  grid()
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
#' LFQbench.plotProteinDispersionBySample
#'
#' draw protein dispersion as cv between replicates by sample
#' @param d the data of a result set (resultSet$data)
#' @export
LFQbench.plotProteinDispersionBySample = function( d )
{
  d4s = sapply( 1:LFQbench.Config$NumberOfSamples, function(si) as.vector( na.exclude( d$cv[,si] ) ) )
  par(LFQbench.Config$par)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=LFQbench.Config$AllSampleNames, at=(1:LFQbench.Config$NumberOfSamples), lwd.ticks=1, lwd=0)
  grid()
  par(LFQbench.Config$parBackup)
}
################################################################################

################################################################################
plotLogRatios = function( logRatios, File )
{
  cat("creating log ratio box plot ...")
  pdf(file=File, onefile=T, width=LFQbench.Config$PlotWidth*1.5, height=LFQbench.Config$PlotHeight*2, family="Helvetica", pointsize=9)
  
  par(LFQbench.Config$par)
  par( las=2 )
  mars = LFQbench.Config$par$mar
  mars[1] = 14
  par( mar = mars )
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"(A:B)" ) ),
            whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F, lims=LFQbench.Config$LogRatioPlotRange, lwd=1, cex=LFQbench.Config$PlotPointSize
  )
  addAxes(showXlab=F)
  par()
  axis(1, labels=names(logRatios), at=(1:length(logRatios)), lwd.ticks=1, lwd=0)
  # grid()
  par(LFQbench.Config$parBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
plotCVs = function( CVs, File )
{
  cat("creating CV boxplot ...")
  pdf(file=File, onefile=T, width=LFQbench.Config$PlotWidth*1.5, height=LFQbench.Config$PlotHeight*1.5, family="Helvetica", pointsize=9)
  
  par(LFQbench.Config$par)
  par( las=2 )
  mars = LFQbench.Config$par$mar
  mars[1] = 10
  par( mar = mars )
  qboxplot( CVs, pch=20,
            ylab="Protein quantification dispersion (CV)",
            whiskerQuantile=LFQbench.Config$BoxPlotWhiskerQuantile, axes=F, lims=range(CVs, na.rm = T), lwd=1, cex=LFQbench.Config$PlotPointSize
  )
  addAxes( showXlab=F )
  par()
  axis(1, labels=names(CVs), at=(1:length(CVs)), lwd.ticks=1, lwd=0)
  # grid()
  par(LFQbench.Config$parBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
#' LFQbench.plotSampleComposition
#' 
#' draw the sample composition as a spine plot
#' @param sampleAmounts
#' @param bgColors
#' @param fgColors
#' @param labMainSize
#' @param pdfFile
#' @param pdfWidth
#' @param pdfHeight
#' @param pdfFontSize
#' @param pdfFontFamily
#' @param pdfMar
#' @export
LFQbench.plotSampleComposition = function(
  sampleAmounts = data.frame( LFQbench.Config$SampleComposition, row.names = 1  ),
  bgColors=NULL, fgColors=NULL, 
  labMainSize=1.5,
  pdfFile=paste(LFQbench.Config$PlotFilesLocation, "/SampleComposition.pdf", sep = ""),
  pdfWidth=2, pdfHeight=3,
  pdfFontSize=9, pdfFontFamily="Helvetica",
  pdfMar=c(2.3,2.8,.5,4.6)
  )
{
  n = dim(sampleAmounts)[2]
  m = dim(sampleAmounts)[1]
  
  if(!is.null(pdfFile)) pdf(file=pdfFile, width=pdfWidth, height=pdfHeight, pointsize=pdfFontSize, family=pdfFontFamily)
  
  par(mfrow=c(1,1), cex.main=labMainSize, family=pdfFontFamily, mar=pdfMar)

  if(is.null(bgColors)) bgColors = gray.colors(m, 0.9, 0.2)
  if(is.null(fgColors)) fgColors = c("black", rep("white", m-1) )
  
  spiner( sampleAmounts, 
          sampleNames = colnames(sampleAmounts),
          bgCols=bgColors,
          fgCols=fgColors, 
          plotTitle= ""
  )
  
  if(!is.null(pdfFile)) dev.off()
}
################################################################################