################################################################################
# user defined R graphics parameters
if(!exists("GraphicsParameters")) GraphicsParameters = list(
    # plot area margins: c(bottom, left, top, right)
    mar=c(3,3.2,.5,.5),
    # plot axis: c(title, label, line)
    mgp=c(2,0.6,0),
    # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
    las=1
  )
# backup original R graphics parameters
GraphicsParametersBackup = par()[ names(GraphicsParameters) ]
if(!exists("PlotCurveLineWidth")) PlotCurveLineWidth = 2
if(!exists("PlotLegendLineWidth")) PlotLegendLineWidth = 8
if(!exists("PlotPointSize")) PlotPointSize = 2
if(!exists("ScatterPlotPointType")) ScatterPlotPointType = "."
################################################################################

################################################################################
plotToFile = function( file, plotFunc, ... )
{
  pdf(file=file, onefile=T, width=PlotWidth, height=PlotHeight, family="Helvetica", pointsize=9)
  par(GraphicsParameters)
  if(is.function(plotFunc)) res = plotFunc(...)
  par(GraphicsParametersBackup)
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
# add normal 
addAxes = function(lwd=1, showXlab=T, showYlab=T, showXAxis=T, showYAxis=T, cex.axis=1)
{
  usr = par()$usr
  if(showXlab) axis(1, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis)
  if(showYlab) axis(2, lwd=0, lwd.ticks=lwd, cex.axis=cex.axis)
  if(showXAxis) lines(x=usr[c(1,2)], y=usr[c(3,3)], lwd=lwd*2)
  if(showYAxis) lines(x=usr[c(1,1)], y=usr[c(3,4)], lwd=lwd*2)
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
  legend(x=pos, legend=AllSpeciesNames, col=SpeciesColors, lwd=PlotLegendLineWidth, inset=.02, bg="white", ...)	
}
################################################################################

################################################################################
plotSpeciesLegends = function()
{
  pdf(file = PlotFilesLocation + "/species_legend_vertical.pdf", width = 1.05, height = 0.66, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=F )
  dev.off()
  pdf(file = PlotFilesLocation + "/species_legend_horizontal.pdf", width = 2.9, height = 0.35, family = "Helvetica", pointsize = 9)
  par(mar=c(0,0,0.1,0))
  emptyPlot( 0:1, 0:1, lwd=0, grid=F, showXlab=F, showYlab=F, axes=F )
  plotSpeciesLegend(pos="top", horiz=T )
  dev.off()
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw scatterplot
showScatterPlot = function( samplePair, showLegend=F, pointType=ScatterPlotPointType, pointSize=PlotPointSize, cex.lab=1.2, cex.axis=1 )
{
  par(GraphicsParameters)
  emptyPlot(samplePair$xlim, samplePair$ylim, grid=F, cex.axis=cex.axis)
  title(main="",
        xlab=as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ),
        ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ),
        cex.lab=cex.lab
  )
  logRatios = sapply( samplePair$data, 
    function(d) { 
      y = d$y
      y[y<samplePair$ylim[1]] = samplePair$ylim[1]
      y[y>samplePair$ylim[2]] = samplePair$ylim[2]
      points( d$x, y, pch=pointType, col= d$col, cex=pointSize )
      return(d$y) 
  } )
  expRatios = sapply( samplePair$data, 
                      function(d) 
                      {
                        theCol = scaleColor(d$col,.8)
                        hLine(y=d$expectation, lty=2, lwd=PlotCurveLineWidth, col=theCol )
                        return(d$expectation) 
                      }
              )
  if(showLegend) plotSpeciesLegend()
  
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw log ratio quiantile based boxplot
showLogRatioBoxPlot = function(samplePair)
{
  par(GraphicsParameters)
  logRatios = sapply( samplePair$data, function(d) { return(d$y) } )
  qboxplot( logRatios, pch=20,
      ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ),
      whiskerQuantile=BoxPlotWhiskerQuantile, axes=F, lims=LogRatioPlotRange,
    )
  addAxes(showXlab=F)
  axis(1, labels=AllSpeciesNames, at=(1:NumberOfSpecies), lwd.ticks=1, lwd=0)
  grid()
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw barplot showing quantification bars
showQuantBarPlot = function(samplePair)
{
  par(GraphicsParameters)
  # emptyPlot(samplePair$xlim, samplePair$ylim)
  
  values = unlist(sapply( samplePair$data, function(d) { return(d$y) } ))
  colors = unlist(sapply( samplePair$data, function(d) { return(rep(d$col,length(d$y)) ) } ))
  
  cat(""+length(values)+" values ...\n")
  
  sortedIndices = order(values, decreasing=T)
  colors = colors[sortedIndices]
  values = values[sortedIndices]
  
  barplot(values, col=colors, border=NA,
        ylim=samplePair$ylim,
        main="", names.arg=NA,
        xlab=as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ),
        ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) )
  )
  
  legendLabels = sapply( samplePair$data, function(d) d$species )
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="topright", legend=legendLabels, col=legendColors, lwd=PlotLegendLineWidth, inset=.02, bg="white")
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw scatter- and density plot
showScatterAndDensityPlot = function(samplePair, showLegend=F)
{
	par(GraphicsParameters)
	par(fig=c(0,0.7,0,1), new=F)
	
	# scatter plot
	emptyPlot(samplePair$xlim, samplePair$ylim, grid=F)
	title(main="",
		xlab=as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ),
		ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) )
	)
	logRatios = sapply( samplePair$data, function(d) { points( d$x, d$y, pch=ScatterPlotPointType, col=d$col, cex=PlotPointSize ); return(d$y) } )
	expRatios = sapply( samplePair$data, function(d) { hLine(y=d$expectation, lty=2, lwd=PlotCurveLineWidth, col=scaleColor(d$col,.8) ); return(d$expectation) } )
	
	if(showLegend) plotSpeciesLegend( horiz=T )
	
	# density plot
	par(fig=c(0.7,1,0,1), new=T)
	pm = par()$mar
	pm[c(2,4)] = 0
	par(mar=pm)
	
	xLim=samplePair$ylim
	yLim=range( lapply(samplePair$data, function(d) range(d$density$y) ) )
	emptyPlot(yLim, xLim, axes=F, grid=F)
	nix = sapply(
		samplePair$data,
		function(d) 
		{
			idx = d$density$x > min(d$y) & d$density$x < max(d$y)
			lines(d$density$y[idx], d$density$x[idx] , type="l", col=d$col, lwd=PlotCurveLineWidth ) 
			hLine(y=d$expectation, lty=2, lwd=PlotCurveLineWidth, col=scaleColor(d$col,.8) )
		}
	)
	
	par(GraphicsParametersBackup)
	par(fig=c(0,1,0,1), new=F)
}
################################################################################

################################################################################
# draw scatter- and boxplot
showScatterAndBoxPlot = function(samplePair, showLegend=F)
{
  par(GraphicsParameters)
  par(fig=c(0,0.7,0,1), new=F)
  
  # scatter plot
  emptyPlot(samplePair$xlim, samplePair$ylim, grid=F)
  title(main="",
        xlab=as.expression( bquote( Log[2]~"("~.(samplePair$name2)~")" ) ),
        ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) )
  )
  nSets = length(samplePair$data)
  logRatios = sapply( samplePair$data, function(d) { points( d$x, d$y, pch=ScatterPlotPointType, col=d$col, cex=PlotPointSize ); return(d$y) } )
  expRatios = sapply( samplePair$data, function(d) { hLine(y=d$expectation, lty=2, lwd=PlotCurveLineWidth, col=scaleColor(d$col,.8) ); return(d$expectation) } )
  
  if(showLegend) plotSpeciesLegend( horiz=T )
  
  # box plot
  par(fig=c(0.7,1,0,1), new=T)
  pm = par()$mar
  pm[c(2,4)] = 0
  par(mar=pm)
  
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ),
            whiskerQuantile=BoxPlotWhiskerQuantile, axes=F, lims=samplePair$ylim, border=SpeciesColors
  )
  
  # add horizontal lines
  expRatios = sapply( samplePair$data, function(d) { hLine(y=d$expectation, lty=2, lwd=PlotCurveLineWidth, col=scaleColor(d$col,.8) ); return(d$expectation) } )
  
  par(GraphicsParametersBackup)
  par(fig=c(0,1,0,1), new=F)
}
################################################################################

################################################################################
# draw log-ratio distribution kernel density plot
showDistributionDensityPlot = function(samplePair, showLegend=F)
{
  xLim=samplePair$ylim
  yLim=range( lapply(samplePair$data, function(d) range(d$density$y) ) )
  par(GraphicsParameters)
  emptyPlot(xLim, yLim)
  title(main="",
        xlab=as.expression( bquote( Log[2]~"("~.(samplePair$name1+":"+samplePair$name2)~")" ) ), 
        ylab="Density"
  )
  sapply(samplePair$data, function(d) lines(d$density, type="l", col=d$col, lwd=PlotCurveLineWidth ) )

  if(showLegend) plotSpeciesLegend( )
  
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw quantification curve
showAllSpeciesQC = function(samplePair)
{
  par(GraphicsParameters)
  emptyPlot(xRange=samplePair$qcrange)
  x = (1:500)/500*samplePair$qcrange[2]
  title(
    main="",
    xlab=as.expression( bquote( "Absolute"~log[2]~"("~.(samplePair$name1)~":"~.(samplePair$name2)~")" ) ),
    ylab="Parts"
  )
  sapply(samplePair$data, function(d) lines(x, d$qcfunction(x), type="l", col=d$col, lwd=PlotCurveLineWidth ) )
  legendLabels = sapply( samplePair$data, function(d) d$species + ": " + d$auqc )
  legendColors = sapply( samplePair$data, function(d) d$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=PlotLegendLineWidth, inset=.02, bg="white")
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
showSingleSpeciesQC = function(samplePairsData, speciesName="MOUSE")
{
  speciesIndex = which( AllSpeciesNames == speciesName )
  par( GraphicsParameters )
  emptyPlot( xRange=AUQCRatioRange )
  title(
    main="",
    xlab=as.expression( bquote( Absolute~log[2]-ratio~"for"~.(speciesName) ) ),
    ylab="Parts"
  )
  x = (1:500)/500*AUQCRatioRange[2]
  sapply(samplePairsData, function(sp){ lines(x, sp$data[[speciesIndex]]$qcfunction(x), type="l", col=sp$col, lwd=PlotCurveLineWidth ) } )
  legendLabels = sapply( samplePairsData, function(sp) "auqc("+sp$name1+":"+sp$name2+"): "+sp$data[[speciesIndex]]$auqc )
  legendColors = sapply( samplePairsData, function(sp) sp$col )
  legend(x="bottomright", legend=legendLabels, col=legendColors, lwd=PlotLegendLineWidth, inset=.02, bg="white")
  par( GraphicsParametersBackup )
}
################################################################################

################################################################################
# draw protein dispersion between replicates by species
plotProteinDispersionBySpecies = function( d )
{
  d4s = sapply( AllSpeciesNames, function(sp) as.vector( na.exclude( d$cv[ d$species==sp ] ) ) )
  par(GraphicsParameters)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=AllSpeciesNames, at=(1:NumberOfSpecies), lwd.ticks=1, lwd=0)
  grid()
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
# draw protein dispersion between replicates by sample
plotProteinDispersionBySample = function( d )
{
  d4s = sapply( 1:NumberOfSamples, function(si) as.vector( na.exclude( d$cv[,si] ) ) )
  par(GraphicsParameters)
  qboxplot(d4s, pch=20, ylab="Protein quantification dispersion (CV)", whiskerQuantile=BoxPlotWhiskerQuantile, axes=F )
  addAxes(showXlab=F)
  axis(1, labels=AllSampleNames, at=(1:NumberOfSamples), lwd.ticks=1, lwd=0)
  grid()
  par(GraphicsParametersBackup)
}
################################################################################

################################################################################
plotLogRatios = function( logRatios, File )
{
  cat("creating log ratio box plot ...")
  pdf(file=File, onefile=T, width=PlotWidth*1.5, height=PlotHeight*2, family="Helvetica", pointsize=9)
  
  par(GraphicsParameters)
  par( las=2 )
  mars = GraphicsParameters$mar
  mars[1] = 14
  par( mar = mars )
  qboxplot( logRatios, pch=20,
            ylab=as.expression( bquote( Log[2]~"(A:B)" ) ),
            whiskerQuantile=BoxPlotWhiskerQuantile, axes=F, lims=LogRatioPlotRange, lwd=1, cex=PlotPointSize
  )
  addAxes(showXlab=F)
  par()
  axis(1, labels=names(logRatios), at=(1:length(logRatios)), lwd.ticks=1, lwd=0)
  # grid()
  par(GraphicsParametersBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
plotCVs = function( CVs, File )
{
  cat("creating CV boxplot ...")
  pdf(file=File, onefile=T, width=PlotWidth*1.5, height=PlotHeight*1.5, family="Helvetica", pointsize=9)
  
  par(GraphicsParameters)
  par( las=2 )
  mars = GraphicsParameters$mar
  mars[1] = 10
  par( mar = mars )
  qboxplot( CVs, pch=20,
            ylab="Protein quantification dispersion (CV)",
            whiskerQuantile=BoxPlotWhiskerQuantile, axes=F, lims=range(CVs, na.rm = T), lwd=1, cex=PlotPointSize
  )
  addAxes( showXlab=F )
  par()
  axis(1, labels=names(CVs), at=(1:length(CVs)), lwd.ticks=1, lwd=0)
  # grid()
  par(GraphicsParametersBackup)
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
plotSampleComposition = function(
  sampleAmounts = data.frame( SampleComposition, row.names = 1  ),
  bgColors=NULL, fgColors=NULL, 
  labMainSize=1.5,
  pdfFile=PlotFilesLocation+"/SampleComposition.pdf",
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