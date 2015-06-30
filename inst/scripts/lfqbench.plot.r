################################################################################
# user defined R graphics parameters
if( !exists("cfg") ) cfg=list()
if( is.null(cfg$initialized) )
{
  cfg$par = list(
    # plot area margins: c(bottom, left, top, right)
    mar=c(3,3.2,.5,.5),
    # plot axis: c(title, label, line)
    mgp=c(2,0.6,0),
    # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
    las=1
  )
  cfg$PlotWidth  = 6
  cfg$PlotHeight = 4
  cfg$PlotCurveLineWidth = 2
  cfg$PlotLegendLineWidth = 4
  cfg$PlotPointSize = 2
  cfg$ScatterPlotPointType = 20
  cfg$PlotPointMinAlpha=.2
}
# backup original R graphics parameters
cfg$parBackup = par()[ names(cfg$par) ]
################################################################################

################################################################################
plotToFile = function( file, plotFunc, ... )
{
  pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
  par(cfg$par)
  if(is.function(plotFunc)) res = plotFunc(...)
  par(cfg$parBackup)
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
# show empty plot with coordinate system in given ranges
# x/y-Axes are plotted on the bottom/left side with given line width
# grid can be switched off
emptyPlot = function(xRange=c(0,1), yRange=c(0,1), lwd=1, grid=T, showXlab=T, showYlab=T, axes=T, cex.axis=1)
{
	plot( 0, 0, xlim=xRange, ylim=yRange, type="n", axes=F, xlab="", ylab="", main="" )
	if(axes) addAxes(lwd, showXlab, showYlab, cex.axis=cex.axis)
	if(grid) grid()
}
################################################################################

################################################################################
# add normal axes having full lengths of plotting area
addAxes = function(lwd=1, showXlab=T, showYlab=T, showXAxis=T, showYAxis=T, cex.axis=1)
{
  usr = par()$usr
  if(showXAxis) lines(x=usr[c(1,2)], y=usr[c(3,3)], lwd=lwd, xpd=T)
  if(showYAxis) lines(x=usr[c(1,1)], y=usr[c(3,4)], lwd=lwd, xpd=T)
  if(showXlab) axis(1, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis)
  if(showYlab) axis(2, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis)
}
################################################################################

################################################################################
# add an x-axis label shifted by marginShift from par()$mgp[1]
addXLab = function(xlab="", marginShift=-0.2, cex.lab = cfg$AxisLabelSize, ...)
{
  m = par()$mgp
  m[1] = m[1]+marginShift
  title(xlab=xlab, mgp=m, cex.lab=cex.lab, ...)
}
# add an y-axis label shifted by marginShift from par()$mgp[1]
addYLab = function(ylab="", marginShift=0.3, cex.lab = cfg$AxisLabelSize, ...)
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
plotSpeciesLegend = function( pos="top", ... )
{
  legend(x=pos, legend=cfg$AllSpeciesNames, col=cfg$SpeciesColors, lwd=cfg$PlotLegendLineWidth, inset=.02, bg="white", ...)	
}
################################################################################

################################################################################
plotSpeciesLegends = function()
{
  pdf(file = cfg$PlotFilesLocation + "/species_legend_vertical.pdf", width = 1.05, height = 0.66, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=F )
  dev.off()
  pdf(file = cfg$PlotFilesLocation + "/species_legend_horizontal.pdf", width = 2.9, height = 0.35, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=T )
  dev.off()
  par(cfg$parBackup)
}
################################################################################

################################################################################
addScatterPointsForSpecies = function(species, samplePairResult, minAlpha=.2, showExpectationLine=T, showRegressionLine=T, useCfgColor=T, rampColors=T, ...)
{
  ds = samplePairResult$data[[species]]
  theCol = ifelse(useCfgColor, cfg$SpeciesColors[which(cfg$AllSpeciesNames==ds$species)], ds$col)
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
    abline(h = ds$expectation, col=scaleColor(theCol, .8), lty="dashed", lwd=cfg$PlotCurveLineWidth)
  }
  # draw lowess regression line
  if(showRegressionLine)
  {
    regLine = lowess( ds$x, ds$y )
    lines( regLine, col=scaleColor(theCol,.5), lty=5, lwd=cfg$PlotCurveLineWidth )
  }
  return(ds$y)
}
################################################################################

################################################################################
# make a scatter plot with smoothed colors
makeScatter = function(samplePair, showLegend=F, showRegLines=F, showExpLines=T, useCfgColor=T )
{
  emptyPlot(samplePair$xlim, samplePair$ylim, grid=F, lwd = cfg$AxisLineThickness, cex.axis = cfg$AxisAnnotationSize)
  addXLab( as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ) )
  addYLab( as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ) )
  logRatios = sapply(cfg$AllSpeciesNames, 
                     addScatterPointsForSpecies, samplePair, 
                      minAlpha=cfg$PlotPointMinAlpha, showExpectationLine=showExpLines, showRegressionLine=showRegLines, 
                      useCfgColor=T, rampColors=ifelse(cfg$PlotPointMinAlpha==1, F, T), cex=cfg$PlotPointSize
  )
  if(showLegend) plotSpeciesLegend( horiz=T )
  return(logRatios)
}
################################################################################

################################################################################
# draw scatterplot
showScatterPlot = function( samplePair, showLegend=F, showRegLines=F )
{
  par(cfg$par)
  logRatios = makeScatter(samplePair, showLegend, showRegLines )
  par(cfg$parBackup)
}
################################################################################

################################################################################
# draw log ratio quiantile based boxplot
showLogRatioBoxPlot = function(samplePair)
{
  par(cfg$par)
  logRatios = sapply( samplePair$data, function(d) { return(d$y) } )
  qboxplot( logRatios, pch=20, ylab="",
            whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F, lims=cfg$LogRatioPlotRange,
      cex.lab = cfg$AxisLabelSize
    )
  addAxes(showXlab=F, lwd = cfg$AxisLineThickness, cex.axis = cfg$AxisAnnotationSize)
  axis(1, labels=cfg$AllSpeciesNames, at=(1:cfg$NumberOfSpecies), 
       lwd.ticks=cfg$AxisLineThickness, lwd=0, cex.axis = cfg$AxisAnnotationSize )
  addYLab( as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ))
  grid()
  par(cfg$parBackup)
}
################################################################################

################################################################################
# draw barplot showing quantification bars
showQuantBarPlot = function(samplePair)
{
  par(cfg$par)
  values = unlist(sapply( samplePair$data, function(d) { return(d$y) } ))
  colors = unlist(sapply( samplePair$data, function(d) { return(rep(d$col,length(d$y)) ) } ))
  cat(""+length(values)+" values ...\n")
  sortedIndices = order(values, decreasing=T)
  colors = colors[sortedIndices]
  values = values[sortedIndices]
  barplot(values, col=colors, border=NA,
        ylim=samplePair$ylim,
        main="", names.arg=NA, xlab="", ylab="",
        cex.lab = cfg$AxisLabelSize, cex.axis = cfg$AxisAnnotationSize
  )
  addXLab( as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ) )
  addYLab( as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ) )
  legendLabels = sapply( samplePair$data, function(d) d$species )
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="topright", legend=legendLabels, col=legendColors, lwd=cfg$PlotLegendLineWidth, inset=.02, bg="white")
  par(cfg$parBackup)
}
################################################################################

################################################################################
# draw scatter- and density plot
showScatterAndDensityPlot = function(samplePair, showLegend=F, showRegLines=F, scatterPlotWidth=0.8)
{
	par(cfg$par)
	par(fig=c(0,scatterPlotWidth,0,1), new=F)
	logRatios = makeScatter(samplePair, showLegend, showRegLines)
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
    theCol = cfg$SpeciesColors[which(cfg$AllSpeciesNames==d$species)]
    lines(d$density$y[idx], d$density$x[idx] , type="l", col=theCol, lwd=cfg$PlotCurveLineWidth ) 
    hLine(y=d$expectation, lty=2, lwd=cfg$PlotCurveLineWidth, col=scaleColor(theCol,.8) )
  }
	nix = sapply( samplePair$data, mkLines )
	par(cfg$parBackup)
	par(fig=c(0,1,0,1), new=F)
}
################################################################################

loadLibrary("scales")

################################################################################
# draw scatter- and boxplot
showScatterAndBoxPlot = function(samplePair, showLegend=F, showRegLines=F, scatterPlotWidth=0.8)
{
  par(cfg$par)
  par(fig=c(0,scatterPlotWidth,0,1), new=F)
  logRatios = makeScatter(samplePair, showLegend, showRegLines)
  # box plot
  par(fig=c(scatterPlotWidth,1,0,1), new=T)
  pm = par()$mar
  pm[c(2,4)] = 0
  par(mar=pm)
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ),
            whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F, lims=samplePair$ylim, border=cfg$SpeciesColors
  )
  # add horizontal lines
  makeExpLines = function(d) 
  {
    theCol = cfg$SpeciesColors[which(cfg$AllSpeciesNames==d$species)]
    hLine(y=d$expectation, lty=2, lwd=cfg$PlotCurveLineWidth, col=scaleColor(theCol,.8) )
    return(d$expectation) 
  }  
  expRatios = sapply( samplePair$data, makeExpLines )
  par(cfg$parBackup)
  par(fig=c(0,1,0,1), new=F)
}
################################################################################

################################################################################
# draw log-ratio distribution kernel density plot
showDistributionDensityPlot = function(samplePair, showLegend=F)
{
  xLim=samplePair$ylim
  yLim=range( lapply(samplePair$data, function(d) range(d$density$y) ) )
  par(cfg$par)
  emptyPlot(xLim, yLim, lwd = cfg$AxisLineThickness, cex.axis = cfg$AxisAnnotationSize )
  addXLab( as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ) )
  addYLab( "Density", .6 )
  sapply(samplePair$data, function(d) lines(d$density, type="l", col=d$col, lwd=cfg$PlotCurveLineWidth ) )
  if(showLegend) plotSpeciesLegend( )
  par(cfg$parBackup)
}
################################################################################

################################################################################
# draw quantification curve
showAllSpeciesQC = function(samplePair)
{
  par(cfg$par)
  emptyPlot(xRange=samplePair$qcrange)
  x = (1:500)/500*samplePair$qcrange[2]
  title(
    main="",
    xlab=as.expression( bquote( "Absolute"~log[2]~"("~.(samplePair$name1)~":"~.(samplePair$name2)~")" ) ),
    ylab="Parts"
  )
  sapply(samplePair$data, function(d) lines(x, d$qcfunction(x), type="l", col=d$col, lwd=cfg$PlotCurveLineWidth ) )
  legendLabels = sapply( samplePair$data, function(d) d$species + ": " + d$auqc )
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=cfg$PlotLegendLineWidth, inset=.02, bg="white")
  par(cfg$parBackup)
}
################################################################################

################################################################################
showSingleSpeciesQC = function(samplePairsData, speciesName="MOUSE")
{
  speciesIndex = which( cfg$AllSpeciesNames == speciesName )
  par( cfg$par )
  emptyPlot( xRange=cfg$AUQCRatioRange )
  title(
    main="",
    xlab=as.expression( bquote( Absolute~log[2]-ratio~"for"~.(speciesName) ) ),
    ylab="Parts"
  )
  x = (1:500)/500*cfg$AUQCRatioRange[2]
  sapply(samplePairsData, function(sp){ lines(x, sp$data[[speciesIndex]]$qcfunction(x), type="l", col=sp$col, lwd=cfg$PlotCurveLineWidth ) } )
  legendLabels = sapply( samplePairsData, function(sp) "auqc("+sp$name1+":"+sp$name2+"): "+sp$data[[speciesIndex]]$auqc )
  legendColors = sapply( samplePairsData, function(sp) sp$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=cfg$PlotLegendLineWidth, inset=.02, bg="white")
  par( cfg$parBackup )
}
################################################################################

################################################################################
# draw protein dispersion between replicates by species
plotProteinDispersionBySpecies = function( d )
{
  d4s = sapply( cfg$AllSpeciesNames, function(sp) as.vector( na.exclude( d$cv[ d$species==sp ] ) ) )
  par(cfg$par)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=cfg$AllSpeciesNames, at=(1:cfg$NumberOfSpecies), lwd.ticks=1, lwd=0)
  grid()
  par(cfg$parBackup)
}
################################################################################

################################################################################
# draw protein dispersion between replicates by sample
plotProteinDispersionBySample = function( d )
{
  d4s = sapply( 1:cfg$NumberOfSamples, function(si) as.vector( na.exclude( d$cv[,si] ) ) )
  par(cfg$par)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=cfg$AllSampleNames, at=(1:cfg$NumberOfSamples), lwd.ticks=1, lwd=0)
  grid()
  par(cfg$parBackup)
}
################################################################################

################################################################################
plotLogRatios = function( logRatios, File )
{
  cat("creating log ratio box plot ...")
  pdf(file=File, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*2, family="Helvetica", pointsize=9)
  
  par(cfg$par)
  par( las=2 )
  mars = cfg$par$mar
  mars[1] = 14
  par( mar = mars )
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"(A:B)" ) ),
            whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F, lims=cfg$LogRatioPlotRange, lwd=1, cex=cfg$PlotPointSize
  )
  addAxes(showXlab=F)
  par()
  axis(1, labels=names(logRatios), at=(1:length(logRatios)), lwd.ticks=1, lwd=0)
  # grid()
  par(cfg$parBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
plotCVs = function( CVs, File )
{
  cat("creating CV boxplot ...")
  pdf(file=File, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*1.5, family="Helvetica", pointsize=9)
  
  par(cfg$par)
  par( las=2 )
  mars = cfg$par$mar
  mars[1] = 10
  par( mar = mars )
  qboxplot( CVs, pch=20,
            ylab="Protein quantification dispersion (CV)",
            whiskerQuantile=cfg$BoxPlotWhiskerQuantile, axes=F, lims=range(CVs, na.rm = T), lwd=1, cex=cfg$PlotPointSize
  )
  addAxes( showXlab=F )
  par()
  axis(1, labels=names(CVs), at=(1:length(CVs)), lwd.ticks=1, lwd=0)
  # grid()
  par(cfg$parBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
plotSampleComposition = function(
  sampleAmounts = data.frame( cfg$SampleComposition, row.names = 1  ),
  bgColors=NULL, fgColors=NULL, 
  labMainSize=1.5,
  pdfFile=cfg$PlotFilesLocation+"/SampleComposition.pdf",
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