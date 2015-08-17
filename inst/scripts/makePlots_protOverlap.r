rm( list=ls() )


DEBUG=T

library(LFQbench)
#source('lfqbench.com.r')
source('lfqbench.plot.r')
source('lfqbench.config.r')
source('lfqbench.defs.r')
#source('lfqbench.access.r')

####################################################################################################
plotWithCfg = function(func, ...)
{
    par(cfg$par)
    func(...)
    par(cfg$parBackup)
}
####################################################################################################
getFirstSamplePair = function(rsName) 
{
    sp = ResultSets[[rsName]]$result[[1]]
    sp$xlim=c(9, 26)
    return(sp)
}
####################################################################################################
plotRS = function(rsName, scatter=T, box=F, kde=F, showRegLines=F)
{
    sp = getFirstSamplePair(rsName)
    if(scatter) showScatterPlot( sp, showRegLines=showRegLines )
    if(box) showScatterAndBoxPlot( sp, showRegLines=showRegLines )
    if(kde) showScatterAndDensityPlot(sp, showRegLines=showRegLines)
}
####################################################################################################
plotRsToFile = function(rsName, ...)
{
    cat(paste("plotting dataset", rsName, " ... "))
    pdf(file.path(outDir, paste0(rsName,".pdf")) , family = "Helvetica", pointsize = 9, width = cfg$PlotWidth, height = cfg$PlotHeight)
    #   tiff(outDir+"/"+rsName+".tiff", 
    #        compression = "zip", family = "Helvetica", pointsize = 9, width = cfg$PlotWidth, height = cfg$PlotHeight, 
    #        units="in", res=600)
    #   png(outDir+"/"+rsName+".png", family = "Helvetica", pointsize = 9, width = cfg$PlotWidth, height = cfg$PlotHeight, 
    #      units="in", res=600)
    #   jpeg(outDir+"/"+rsName+".jpg", quality=75, family = "Helvetica", pointsize = 9, width = cfg$PlotWidth, height = cfg$PlotHeight, 
    #     units="in", res=600)
    plotRS(rsName, ...)
    dev.off()  
    cat("done!\n")
}
####################################################################################################
plotHComp <- function(file=NULL)
{
    if(!is.null(file))
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight*.8, family="Helvetica", pointsize=9)
        par(cfg$par)
    }
    # par(mar=c(4,3,1,2.5))
    par(mar=c(4,5.5,1,2.5))
    par(lwd=2)
    sampleAmounts = data.frame( cfg$SampleComposition, row.names = 1  )
    barData = as.matrix(sampleAmounts[,2:1])
    barplot(
        barData, 
        col=cfg$SpeciesColors, 
        names.arg = colnames(barData), 
        horiz=T, 
        axes=F,
        cex.names=cfg$AxisLabelSize,
        ylab = "Sample",
        cex.lab = cfg$AxisLabelSize      
    )
    # text( x = c(32.5, 72.5, 90 ), y = rep(2.7, 3), labels = cfg$AllSpeciesNames, xpd=T, cex = cfg$AxisLabelSize )
    axis(1, at = seq(0,100,20), labels = paste0(seq(0,100,20),"%"),line = 0.6, cex.axis = cfg$AxisAnnotationSize, lwd=cfg$AxisLineThickness)
    # legend(50, 3.8, xjust=0.5, cex=cfg$AxisLabelSize, xpd=T, horiz=T, legend=cfg$AllSpeciesNames, col=cfg$SpeciesColors, lwd=cfg$PlotLegendLineWidth, bg="white", bty="n")
    if(!is.null(file))
    {
        par(cfg$parBackup)
        dev.off()
    }
}
####################################################################################################
speciesLegendH <- function(file=NULL)
{
    if(!is.null(file))
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight*.2, family="Helvetica", pointsize=9)
        par(cfg$par)
    }
    par(mar=c(0,0,0,0))
    plot(1, type = "n", xlab = "", ylab = "", main = "", axes=F)
    legend("top", xjust=0.5, cex=cfg$AxisLabelSize, xpd=T, horiz=T, legend=cfg$AllSpeciesNames, col=cfg$SpeciesColors, lwd=cfg$PlotLegendLineWidth, bg="white", bty="n")
    if(!is.null(file))
    {
        par(cfg$parBackup)
        dev.off()
    }
}
####################################################################################################
speciesLegendV <- function(file=NULL)
{
    if(!is.null(file))
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
        par(cfg$par)
    }
    par(mar=c(0,0,0,0))
    plot(1, type = "n", xlab = "", ylab = "", main = "", axes=F)
    legend("top", xjust=0.5, cex=cfg$AxisLabelSize, xpd=T, horiz=F, legend=cfg$AllSpeciesNames, col=cfg$SpeciesColors, lwd=cfg$PlotLegendLineWidth, bg="white", bty="n")
    if(!is.null(file))
    {
        par(cfg$parBackup)
        dev.off()
    }
}
####################################################################################################
intensityRangeLabels = c("[0,2]","(2,4]","(4,6]","(6,8]","(8,10]", bquote("(10,"*infinity*")"))

showBarPlot = function(vals, cols, labs, showLegend=F, valLab="Value", showGrid=T, catLabs=as.expression(intensityRangeLabels), file=NULL)
{
    if(!is.null(file)) pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
    par(cfg$par)
    par(  mar = c( 6.3, 5.4, .7, 0 ),
          # plot axis: c(title, label, line)
          mgp = c( 4, 0.3, 0 )
    )
    
    # ylim = c(0, ceiling( max(vals, na.rm=T)*2 )/2 )
    ylim = c(0, max(vals, na.rm=T) )
    barx = barplot(vals, col=cols, names.arg=catLabs, 
                   ylim = ylim, cex.names = cfg$AxisAnnotationSize,
                   border = NA, beside = T, plot=T, axes=F, las=2, space=c(0.3, 3) 
    )
    
    # suppress 0.05 y-axis ticks
    yTicks = unique(floor((axTicks(2)+.05)*10)/10)
    
    if(showGrid) 
    {
        abline(h=yTicks, col="lightgray", lty="dashed")
    }
    
    if(showLegend)
    {
        legend("topright", legend = labs, lwd = 6, col=cols, bty = "n", ncol=2 ) #, cex = cfg$AxisAnnotationSize)
    }
    
    apply(barx, 2, function(v) lines(range(v), c(0,0), lwd=cfg$AxisLineThickness, xpd=T) )
    
    axis(2, pos=1, at = yTicks, labels = format(yTicks, 1, trim=T ),
         lwd=0, lwd.ticks=cfg$AxisLineThickness, cex.axis=cfg$AxisAnnotationSize, mgp=c( 6, 1, 0 )
    )
    lines(x=c(1,1), y=ylim, lwd=cfg$AxisLineThickness)
    
    title(xlab=as.expression( bquote( ""~Log[2]*"(B)" ) ), cex.lab=cfg$AxisLabelSize, mgp=c( 5.3, 1, 0 ) )
    title(ylab=valLab, cex.lab=cfg$AxisLabelSize, mgp=c( 2.8, 1, 0 ) )
    
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
    return(barx)
}
####################################################################################################
plotBarLegend = function(labs, cols, nCol=2, file=NULL)
{
    if(!is.null(file)) pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
    par(  mar = c( 0, 0.5, 0, 0 ))
    plot(0,0,type = "n", axes=F, xlab="", ylab="" )
    legend("top", legend = labs, lwd = 8, col=cols, bty = "n", ncol=nCol, cex = cfg$AxisAnnotationSize, xpd=T)
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}
####################################################################################################
c4s = function(m, spc) m[,colnames(m)==spc]
# accuracy: the median deviation of log-ratios to the expected value
getAcc = function(rs, spc) c4s(rs$result[[1]]$rangedAccuracy, spc)
# precision: the standard deviation of log-ratios
getPrc = function(rs, spc) c4s(rs$result[[1]]$rangedPrecision, spc)
####################################################################################################
plotRound = function(roundNames, species="ECOLI", filePrefix="", cols=NULL)
{
    labs = gsub("_r1", "", roundNames)
    labs = gsub("_r2", " (it.2)", labs)
    labs = gsub("_builtin", " (built-in)", labs)
    if( is.null(cols) ) cols = brewer.pal(length(roundNames), "Dark2")
    acc = t(abs(sapply(ResultSets[roundNames], getAcc, species)))
    prc = t(abs(sapply(ResultSets[roundNames], getPrc, species)))
    showBarPlot(
        acc, cols, labs, 
        file= file.path(outDir, paste0(filePrefix, "_accuracy_(", species, "_logratio_delta).pdf")), 
        valLab=bquote( group("|",Delta*"("~log[2]*"(A:B)"~")", "|") )
    )
    showBarPlot(
        prc, cols, labs, 
        file= file.path(outDir, paste0(filePrefix, "_precision_(", species, "_logratio_sd).pdf")),
        valLab=bquote( "SD("~log[2]*"(A:B) )" )
    )
    return(list( labs=labs, cols=cols ))
}
####################################################################################################

software.names = c(
    "DIAumpire",
    "OpenSWATH", 
    "PeakView", 
    "Skyline", 
    "Spectronaut"
)

protein.experiments.r1 = c(
    "DIAumpire_builtin_r1", "DIAumpire_r1",
    "OpenSWATH_r1", 
    "PeakView_builtin_r1", "PeakView_r1",
    "Skyline_r1", 
    "Spectronaut_r1"
)

protein.experiments.r2 = c(
    "OpenSWATH_r2",
    "PeakView_r2",
    "Skyline_r2",
    "Spectronaut_r2"
)

peptide.experiments.r1 = c(
    "DIAumpire_peptides_r1", 
    "OpenSWATH_peptides_r1", 
    "PeakView_peptides_r1", 
    "Skyline_peptides_r1", 
    "Spectronaut_peptides_r1"
)

peptide.experiments.r2 = c(
    "OpenSWATH_peptides_r2", 
    "PeakView_peptides_r2", 
    "Skyline_peptides_r2",
    "Spectronaut_peptides_r2"
)

experimentColors = brewer.pal(10, "Paired")
names(experimentColors) = c(
    "DIAumpire_r1",
    "OpenSWATH_r1", "OpenSWATH_r2",
    "PeakView_r1",  "PeakView_r2", 
    "Skyline_r1", "Skyline_r2",
    "Spectronaut_r1", "Spectronaut_r2"
)


#if(DEBUG) cat( "R working directory: " + getwd() + "\n")

working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w"
evalCommandLineArguments()
datafile <- file.path(working_dir, "log","ResultSets.rda")
outDir= file.path(working_dir, "output_figures")
load(file = datafile)
mkdir(outDir)


experimentNames = names(ResultSets)

# r1.names = experimentNames[grep("r1", experimentNames)]
# r2.names = experimentNames[grep("(OpenSWATH|PeakView|Spectronaut)", experimentNames)]
# r1.cols = experimentColors[r1.names]
# r2.cols = experimentColors[r2.names]

# plotRound(r1.names, species = "ECOLI", filePrefix = "round1", cols=r1.cols)
# plotRound(r1.names, species = "YEAST", filePrefix = "round1", cols=r1.cols)
# r1.labsAndCols = plotRound(r1.names, species = "HUMAN", filePrefix = "round1", cols=r1.cols)
# plotBarLegend(r1.labsAndCols$labs, r1.labsAndCols$cols, nCol=2, outDir+"/round1_legend_h.pdf")
# plotBarLegend(r1.labsAndCols$labs, r1.labsAndCols$cols, nCol=1, outDir+"/round1_legend_v.pdf")

# plotRound(r2.names, species = "ECOLI", filePrefix = "round2", cols=r2.cols)
# plotRound(r2.names, species = "YEAST", filePrefix = "round2", cols=r2.cols)
# r2.labsAndCols = plotRound(r2.names, species = "HUMAN", filePrefix = "round2", cols=r2.cols)
# plotBarLegend(r2.labsAndCols$labs, r2.labsAndCols$cols, nCol=2, outDir+"/round2_legend_h.pdf")
# plotBarLegend(r2.labsAndCols$labs, r2.labsAndCols$cols, nCol=1, outDir+"/round2_legend_v.pdf")

mvAnteil = function(spc="HUMAN", mv)
{
    spcs = names(mv)
    v = mv[ , spcs==spc ]
    nAll = sum(v)
    nMV = sum(v[-1])
    return( nMV/nAll*100 )
}
mv4rs = function(rs) sapply(cfg$AllSpeciesNames, mvAnteil, rs$data$missingvalues )

nams2labs = function(nams)
{
    labs = gsub("_r1", "", nams)
    labs = gsub("_r2", "\n(improved)", labs) # It was: labs = gsub("_r2", " (improved)", labs)
    labs = gsub("_peptides", "", labs)
    labs = gsub("_builtin", "\n(built-in)", labs)
    labs = gsub("DIAumpire", "DIA-Umpire", labs)
    labs = gsub("OpenSWATH", "OpenSWATH", labs)
    return(labs)
}

plotMissingValues = function( missingValues, file=NULL, labXSize=cfg$AxisAnnotationSize, labYSize=cfg$AxisAnnotationSize, catNameShift=0 )
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
    }
    
    par(cfg$par)
    par(  mar = c( 8.3, 5, .8, 0 ),
          # plot axis: c(title, label, line)
          # mgp = c( 4, 0.3, 0 )
          mgp=c( 3.4, 1, 0 )
    )
    labs = nams2labs(colnames(missingValues))
    barx = barplot(missingValues, names.arg = rep("", ncol(missingValues)), 
                   beside=T, ylim=c(0,100), col=cfg$SpeciesColors, cex.axis=cfg$AxisAnnotationSize, axes=F
    )
    labsx = colMeans(barx)
    text(labsx+.5, rep(-5, length(labs)) - catNameShift, labels = labs, xpd=T, cex=labXSize, srt=90, pos=2 )
    addYLab("Incomplete cases (%)", cex.lab=labYSize)
    addAxes(showYAxis = T, showXAxis = F, cex.axis = labYSize, lwd = cfg$AxisLineThickness, showXlab = F)
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

# figure 2: missing values in round 1 peptides
mvs = sapply(ResultSets[peptide.experiments.r1], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_peptides_r1.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

mvs = sapply(ResultSets[protein.experiments.r1], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_proteins_r1.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

# Supp. to figure 2: missing values in round 2 
mvs = sapply(ResultSets[peptide.experiments.r2], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_peptides_r2.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

mvs = sapply(ResultSets[protein.experiments.r2], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_proteins_r2.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)


# make scatterplot files for all datasets
sapply( experimentNames,  plotRsToFile,  scatter=F, box=T, kde=F, showRegLines=T )
speciesLegendH(file.path(outDir, "species_legend_h.pdf"))

# figure 5
# ECOLI log-ratios as boxplots combined for each software Built-in, R1, R2
# having empty space for missing datasets

getLogRatios = function( rsIndex=1, species="ECOLI" ) 
{
    if(length(rsIndex)<1)
        return(NULL)
    else
        return(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$y)
}

getProteinLogRatios = function( software="DIAumpire", species="ECOLI" )
{
    # built-in
    biLR = getLogRatios( grep( paste0(software, "_builtin_r1"), experimentNames), species)
    # r1
    r1LR = getLogRatios( grep( paste0(software, "_r1"), experimentNames), species)
    # r2
    r2LR = getLogRatios( grep( paste0(software, "_r2"), experimentNames), species)
    return( list("Built-in"=biLR, "Iteration 1"=r1LR, "Iteration 2"=r2LR))
}

getPartialLogRatios = function( rsIndex=1, species="ECOLI", fromPart=0.0, toPart=1.0 )
{
    if(length(rsIndex)<1)
    {
        return(NULL)
    }
    else
    {
        # increasing order of log intensities
        xorder = order(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$x)
        # sort log ratios in the order of increasing log intensities
        lrs = getLogRatios(rsIndex = rsIndex, species = species)
        lrs = lrs[xorder]
        n = length(lrs)
        fromIdx = floor(fromPart * n)
        toIdx = floor(toPart * n)
        return(lrs[fromIdx:toIdx])
    }		
}

getPartialProteinLogRatios = function( software="DIAumpire", species="ECOLI", fromPart=0, toPart=1 )
{
    cat(paste0("selecting log ratios for software: ", software, " ...\n"))	
    # built-in
    biLR = getPartialLogRatios( grep(paste0(software, "_builtin_r1"), experimentNames), species=species, fromPart, toPart)
    # r1
    r1LR = getPartialLogRatios( grep(paste0(software, "_r1"), experimentNames), species=species, fromPart, toPart)
    # r2
    r2LR = getPartialLogRatios( grep(paste0(software, "_r2"), experimentNames), species=species, fromPart, toPart)
    return( list("Built-in"=biLR, "Iteration 1"=r1LR, "Iteration 2"=r2LR) )
}

expectedLogRatios = log2( cfg$AllExpectedAmounts[,1] / cfg$AllExpectedAmounts[,2] )

# boxplot taking all log ratios into account
showBoxPlot4Softwares = function( species="ECOLI", file=NULL  )
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*2, height=cfg$PlotHeight*1.5, family="Helvetica", pointsize=9)
    }
    par(cfg$par)
    par(  mar = c( 10.5, 5.5, .5, 0.5 ),
          # plot axis: c(title, label, line)
          mgp=c( 3.6, 1.0, 0 )
    )
    lrs = unlist( lapply(software.names, getProteinLogRatios ), recursive = F )
    
    boxPos = c(1:3, 5:7, 9:11, 13:15, 17:19)
    sepPos = c(4,8,12,16)
    
    qboxplot( lrs, whiskerQuantile = 0.01, 
              ylab=as.expression( bquote( Log[2]~"(A:B)" ) ),
              horizontal = F, pch=20, cex=cfg$PlotPointSize, axes=F,
              cex.lab=cfg$AxisLabelSize,
              at=boxPos
    )
    
    plotArea = par()$usr
    
    # plot expectation line
    abline(h=expectedLogRatios[species], lty=2, lwd=cfg$PlotCurveLineWidth )
    abline(v=sepPos, lty=2, lwd=cfg$PlotCurveLineWidth, col="gray" )
    
    addAxes(lwd = cfg$AxisLineThickness, showXlab = F, showYlab = T, showXAxis = T, showYAxis = T, cex.axis = cfg$AxisAnnotationSize )
    
    yAxRange = par()$usr[3:4]
    yAxSize = max(yAxRange) - min(yAxRange)
    
    axis(1, labels=names(lrs), at=boxPos, lwd.ticks = cfg$AxisLineThickness, lwd=0, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=min(yAxRange)-0.36*yAxSize , labels=nams2labs(software.names), at=boxPos[(1:5)*3-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
    return(plotArea)
}

showBoxPlot4Softwares( species="ECOLI", file = file.path(outDir, "software_boxplots_ecoli.pdf" ))



# boxplots for low/mid/high intensity ranges
showPartialBoxPlot4Softwares = function( species="ECOLI", fromPart=0, toPart=1, file=NULL  )
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*2, height=cfg$PlotHeight*1.5, family="Helvetica", pointsize=9)
    }
    mar = c( 0.5, 5.5, .5, 0.5 )
    mgp=c( 3.6, 1.0, 0 )
    par( cfg$par )
    par( mar=mar, mgp=mgp )
    lrs = unlist( lapply( software.names, getPartialProteinLogRatios, species, fromPart, toPart ), recursive = F )
    
    boxPos = c( 1:3, 5:7, 9:11, 13:15, 17:19 )
    sepPos = c( 4, 8, 12, 16 )
    
    qboxplot( lrs, whiskerQuantile = cfg$BoxPlotWhiskerQuantile, lims=c(expectedLogRatios[species] - 3, expectedLogRatios[species] + 3.2), 
              ylab=as.expression( bquote( Log[2]~"(A:B)" ) ),
              horizontal = F, pch=20, cex=cfg$PlotPointSize, axes=F,
              cex.lab=cfg$AxisLabelSize, at=boxPos
    )
    
    plotArea = par()$usr
    
    # plot expectation line
    abline(h=expectedLogRatios[species], lty=2, lwd=cfg$PlotCurveLineWidth )
    abline(v=sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    
    addAxes(lwd = cfg$AxisLineThickness, showXlab = F, showYlab = T, showXAxis = T, showYAxis = T, cex.axis = cfg$AxisAnnotationSize )
    
    if(!is.null(file)) 
    {
        yAxRange = par()$usr[3:4]
        yAxSize = max(yAxRange) - min(yAxRange)
        
        axis(1, labels=names(lrs), at=boxPos, lwd.ticks = cfg$AxisLineThickness, lwd=0, cex.axis=cfg$AxisAnnotationSize, las=2)
        axis(1, pos=min(yAxRange)-0.36*yAxSize , labels=nams2labs(software.names), at=boxPos[(1:5)*3-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)	
    }
    
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
    return( list(
        usr=plotArea, 
        boxPos=boxPos, 
        sepPos=sepPos, 
        labs=names(lrs), 
        cats=nams2labs(software.names),
        mar=mar,
        mgp=mgp
    ))
}

softwareBoxPlotsByIntensityRange = function(species="ECOLI", file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*3, family="Helvetica", pointsize=9)
    }
    par( cfg$par )
    layout( matrix(c(1,2,3,4), ncol=1, byrow = TRUE), heights=c(1,1,1,.4))
    showPartialBoxPlot4Softwares( species=species, fromPart=0, toPart=1/3 )
    showPartialBoxPlot4Softwares( species=species, fromPart=1/3, toPart=2/3 )
    pa = showPartialBoxPlot4Softwares( species=species, fromPart=2/3, toPart=1 )
    par( cfg$par )
    par(  mar = c( 1, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = 0:1, axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0, labels=pa$cats, at=pa$boxPos[(1:5)*3-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsByIntensityRange("ECOLI", file.path(outDir, "software_boxplots_by_range_ecoli.pdf"))
softwareBoxPlotsByIntensityRange("YEAST", file.path(outDir, "software_boxplots_by_range_yeast.pdf"))
softwareBoxPlotsByIntensityRange("HUMAN", file.path(outDir, "software_boxplots_by_range_human.pdf"))


softwareBoxPlotsLowIntensityRange = function(species="ECOLI", file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*3, family="Helvetica", pointsize=9)
    }
    par( cfg$par )
    layout( matrix(c(1,2,3,4), ncol=1, byrow = TRUE), heights=c(1,1,1,.4))
    pa = showPartialBoxPlot4Softwares( species=species, fromPart=0, toPart=1/3 )
    par( cfg$par )
    par(  mar = c( 5.5, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = c(0,1.1), axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0.6, labels=pa$cats, at=pa$boxPos[(1:5)*3-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsLowIntensityRange("ECOLI", file.path(outDir, "software_boxplots_low_intensity_range_ecoli.pdf") )


### Table 1 
library(dplyr)

getIdsMetrics <- function(met){
    #met = metrics[[3]]
    ids = met$identification
    quants = met$quantification[[1]][4] + met$quantification[[1]][10] + met$quantification[[1]][16]
    sepYH = mean(met$separation[[1]][1:6], na.rm=T)
    sepEH = mean(met$separation[[1]][7:12], na.rm=T)
    
    return(c(ids, quants, sepYH, sepEH ))
}

metrics = lapply( ResultSets, getMetrics )
#explainMetrics()
#x = lapply(metrics, showMetrics)
table1 = as.data.frame(lapply(metrics, getIdsMetrics))

table1 = as.data.frame(t(table1))
names(table1) = c("Identifications", "Valid_quant_ratios", "separationYeastHuman", "separationEColiHuman")

write.table(x=table1, file.path(outDir, "table1.tsv"), sep="\t", row.names=T)


## Supp. Table B: shift log2Ratios applied

getAdjustmentValue <- function(softwarename=1,round=1){
    rsIndex = grep(paste0(softwarename, "_peptides_r", round), experimentNames)
    if(length(rsIndex) == 0) return (NA)
    adjustment = ResultSets[[rsIndex]]$result[[1]]$adjustment
    return(adjustment)
}

adjustments = as.data.frame(
    rbind(
        sapply(software.names, getAdjustmentValue, 1),
        sapply(software.names, getAdjustmentValue, 2)
    )
)

row.names(adjustments) = c("iteration.1", "iteration.2")
write.table(x=adjustments, file.path(outDir, "suppTable_Log2adjustments.tsv"), sep="\t", row.names=T)


## Supp. Table C: dynamic ranges of all softwares (all iterations, peptide and protein level)

getDynamicRange <- function(softwarename=1, round=1, species="HUMAN", peptides=T, thresholdMin = 0.01){
    rsIndex = grep(paste0(softwarename, "_r", round), experimentNames)
    if(peptides) rsIndex = grep(paste0(softwarename, "_peptides_r", round), experimentNames)
    if(length(rsIndex) == 0) return (NA)
    
    minInt = quantile(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$x, probs = 0.01, na.rm=T)
    maxInt = max(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$x, na.rm=T)
    dynRange = maxInt - minInt
    return(dynRange)
}

dynamicranges = as.data.frame(
    rbind(
        sapply(software.names, getDynamicRange, 1),
        sapply(software.names, getDynamicRange, 2),
        sapply(software.names, getDynamicRange, 1, peptides=F),
        sapply(software.names, getDynamicRange, 2, peptides=F)
    )
)

row.names(dynamicranges) = c("peptides.iteration.1", 
                             "peptides.iteration.2", 
                             "proteins.iteration.1", 
                             "proteins.iteration.2")

colnames(dynamicranges) = software.names

write.table(x=dynamicranges, file.path(outDir, "Supplementary.Table.C.tsv"), sep="\t", 
            row.names=T, col.names = T)


