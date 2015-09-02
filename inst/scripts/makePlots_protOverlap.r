rm( list=ls() )


DEBUG=T

library(LFQbench)
#source('lfqbench.com.r')
source('lfqbench.plot.r')
source('lfqbench.config.r')
source('lfqbench.defs.r')
#source('lfqbench.access.r')


#if(DEBUG) cat( "R working directory: " + getwd() + "\n")

working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w"
evalCommandLineArguments()
datafile <- file.path(working_dir, "log","ResultSets.rda")
outDir= file.path(working_dir, "output_figures")
load(file = datafile)
mkdir(outDir)


experimentNames = names(ResultSets)

####################################################################################################
plotWithCfg = function(func, ...)
{
    par(cfg$par)
    func(...)
    par(cfg$parBackup)
}
####################################################################################################
getFirstSamplePair = function(rsName, xlims) 
{
    sp = ResultSets[[rsName]]$result[[1]]
    sp$xlim= xlims #c(13, 26)
    return(sp)
}
####################################################################################################
plotRS = function(rsName, scatter=T, box=F, kde=F, showRegLines=F, xlims = c(13,26))
{
    sp = getFirstSamplePair(rsName, xlims)
    if(scatter) showScatterPlot( sp, showRegLines=showRegLines )
    if(box) showScatterAndBoxPlot( sp, showRegLines=showRegLines )
    if(kde) showScatterAndDensityPlot(sp, showRegLines=showRegLines)
}
####################################################################################################
plotRsToFile = function(rsName, ...)
{
    cat(paste("plotting dataset", rsName, " ... "))
    pdf(file.path(outDir, paste0(rsName,".pdf")) , family = "Helvetica", pointsize = 9, width = cfg$PlotWidth, height = cfg$PlotHeight)

    xlims = c(13,26)
    if(grepl("_peptides", rsName)){
        xlims = c(9,26)
    }
    plotRS(rsName, xlims = xlims, ...)
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

protein.experiments.r1.nobuiltin = c(
    "DIAumpire_r1",
    "OpenSWATH_r1", 
    "PeakView_r1",
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

# figure 2: missing values in round 1
mvs = sapply(ResultSets[peptide.experiments.r1], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_peptides_r1.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

mvs = sapply(ResultSets[protein.experiments.r1], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_proteins_r1.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

mvs = sapply(ResultSets[protein.experiments.r1.nobuiltin], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_proteins_r1_nobuiltin.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)


# Supp. to figure 2: missing values in round 2 
mvs = sapply(ResultSets[peptide.experiments.r2], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_peptides_r2.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)

mvs = sapply(ResultSets[protein.experiments.r2], mv4rs)
plotMissingValues( mvs, file = file.path(outDir, "missing_values_proteins_r2.pdf"), labXSize=cfg$AxisAnnotationSize*.70, catNameShift=0)


# make scatterplot files for all datasets
sapply( experimentNames,  plotRsToFile,  scatter=F, box=T, kde=F, showRegLines=T )
speciesLegendH(file.path(outDir, "species_legend_h.pdf"))
speciesLegendV(file.path(outDir, "species_legend_v.pdf"))

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

getPartialProteinLogRatiosNoBuiltIn = function( software="DIAumpire", species="ECOLI", fromPart=0, toPart=1 )
{
    cat(paste0("selecting log ratios for software: ", software, " ...\n"))    
    # r1
    r1LR = getPartialLogRatios( grep(paste0(software, "_r1"), experimentNames), species=species, fromPart, toPart)
    # r2
    r2LR = getPartialLogRatios( grep(paste0(software, "_r2"), experimentNames), species=species, fromPart, toPart)
    return( list("Iteration 1"=r1LR, "Iteration 2"=r2LR) )
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

getPartialPeptideLogRatios = function( software="DIAumpire", species="ECOLI", fromPart=0, toPart=1 )
{
    cat(paste0("selecting log ratios for software: ", software, " ...\n"))    
    # r1
    r1LR = getPartialLogRatios( grep(paste0(software, "_peptides_r1"), experimentNames), species=species, fromPart, toPart)
    # r2
    r2LR = getPartialLogRatios( grep(paste0(software, "_peptides_r2"), experimentNames), species=species, fromPart, toPart)
    return( list("Iteration 1"=r1LR, "Iteration 2"=r2LR) )
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


showPartialBoxPlot4SoftwaresNoBuiltIn = function( species="ECOLI", fromPart=0, toPart=1, file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*2, height=cfg$PlotHeight*1.5, family="Helvetica", pointsize=9)
    }
    mar = c( 0.5, 5.5, .5, 0.5 )
    mgp=c( 3.6, 1.0, 0 )
    par( cfg$par )
    par( mar=mar, mgp=mgp )
    
    lrs = unlist( lapply( software.names, getPartialProteinLogRatiosNoBuiltIn, species, fromPart, toPart ), recursive = F )
    
    
    boxPos = c( 1:2, 4:5, 7:8, 10:11, 13:14 )
    sepPos = c( 3, 6, 9, 12 )
    
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




# boxplots for low/mid/high intensity ranges
showPartialBoxPlot4Softwares = function( species="ECOLI", fromPart=0, toPart=1, file=NULL, usebuiltin = T  )
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


# BoxPlots by intensity range for Peptides

# boxplots for low/mid/high intensity ranges
showPartialBoxPlot4SoftwaresPeptides = function( species="ECOLI", fromPart=0, toPart=1, file=NULL  )
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*2, height=cfg$PlotHeight*1.5, family="Helvetica", pointsize=9)
    }
    mar = c( 0.5, 5.5, .5, 0.5 )
    mgp=c( 3.6, 1.0, 0 )
    par( cfg$par )
    par( mar=mar, mgp=mgp )
    lrs = unlist( lapply( software.names, getPartialPeptideLogRatios, species, fromPart, toPart ), recursive = F )
    
    boxPos = c( 1:2, 4:5, 7:8, 10:11, 13:14 )
    sepPos = c( 3, 6, 9, 12 )
    
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
        axis(1, pos=min(yAxRange)-0.36*yAxSize , labels=nams2labs(software.names), at=boxPos[(1:5)*2-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)    
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

softwareBoxPlotsByIntensityRangePeptides = function(species="ECOLI", file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*3, family="Helvetica", pointsize=9)
    }
    par( cfg$par )
    layout( matrix(c(1,2,3,4), ncol=1, byrow = TRUE), heights=c(1,1,1,.4))
    showPartialBoxPlot4SoftwaresPeptides( species=species, fromPart=0, toPart=1/3 )
    showPartialBoxPlot4SoftwaresPeptides( species=species, fromPart=1/3, toPart=2/3 )
    pa = showPartialBoxPlot4SoftwaresPeptides( species=species, fromPart=2/3, toPart=1 )
    par( cfg$par )
    par(  mar = c( 1, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = 0:1, axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0, labels=pa$cats, at=pa$boxPos[(1:5)*2-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsByIntensityRangePeptides("ECOLI", file.path(outDir, "software_boxplots_by_range_peptides_ecoli.pdf"))
softwareBoxPlotsByIntensityRangePeptides("YEAST", file.path(outDir, "software_boxplots_by_range_peptides_yeast.pdf"))
softwareBoxPlotsByIntensityRangePeptides("HUMAN", file.path(outDir, "software_boxplots_by_range_peptides_human.pdf"))



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
    #par(  mar = c( 5.5, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    par(  mar = c( 10.5, 5.5, .5, 0.5 ), mgp=c( 3.6, 1.0, 0 ) )
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = c(0,1.1), axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0.6, labels=pa$cats, at=pa$boxPos[(1:5)*3-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsLowIntensityRange("ECOLI", file.path(outDir, "software_boxplots_low_intensity_range_ecoli.pdf") )

softwareBoxPlotsLowIntensityRangeNoBuiltIn = function(species="ECOLI", file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*3, family="Helvetica", pointsize=9)
    }
    par( cfg$par )
    layout( matrix(c(1,2,3,4), ncol=1, byrow = TRUE), heights=c(1,1,1,.4))
    pa = showPartialBoxPlot4SoftwaresNoBuiltIn( species=species, fromPart=0, toPart=1/3 )
    par( cfg$par )
    #par(  mar = c( 5.5, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    par(  mar = c( 10.5, 5.5, .5, 0.5 ), mgp=c( 3.6, 1.0, 0 ) )
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = c(0,1.1), axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0.6, labels=pa$cats, at=pa$boxPos[(1:5)*2-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsLowIntensityRangeNoBuiltIn("ECOLI", file.path(outDir, "software_boxplots_low_intensity_range_ecoli_nobuiltin.pdf") )

softwareBoxPlotsLowIntensityRangePeptidesNoBuiltIn = function(species="ECOLI", file=NULL)
{
    if(!is.null(file)) 
    {
        pdf(file=file, onefile=T, width=cfg$PlotWidth*1.5, height=cfg$PlotHeight*3, family="Helvetica", pointsize=9)
    }
    par( cfg$par )
    layout( matrix(c(1,2,3,4), ncol=1, byrow = TRUE), heights=c(1,1,1,.4))
    pa = showPartialBoxPlot4SoftwaresPeptides( species=species, fromPart=0, toPart=1/3 )
    par( cfg$par )
    #par(  mar = c( 5.5, 5.5, .5, 0.5 ), mgp=c( 0, 0, 0 ) )
    par(  mar = c( 10.5, 5.5, .5, 0.5 ), mgp=c( 3.6, 1.0, 0 ) )
    emptyPlot(xRange = c(min(pa$boxPos)-0.5, max(pa$boxPos)+0.5), yRange = c(0,1.1), axes = F, grid=F)
    axis(1, pos=1, labels=pa$labs, at=pa$boxPos, tick=F, cex.axis=cfg$AxisAnnotationSize, las=2)
    axis(1, pos=0.6, labels=pa$cats, at=pa$boxPos[(1:5)*2-1], tick=F, cex.axis=cfg$AxisAnnotationSize, las=1)
    abline(v=pa$sepPos, lty="dotted", lwd=cfg$PlotCurveLineWidth, col="gray" )
    if(!is.null(file)) dev.off()
    par(cfg$parBackup)
}

softwareBoxPlotsLowIntensityRangePeptidesNoBuiltIn("ECOLI", file.path(outDir, "software_boxplots_low_intensity_range_peptides_ecoli.pdf") )


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
    #transform to log10-scale
    dynRange = dynRange / log(10,base = 2)
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


# Table with standard deviations and deviations from expectation of each population (Yeast, Human, E.Coli) and dataset


get_SD <- function(softwarename=1, round=1, species="HUMAN", peptides=T, from=0, to=1){
    rsIndex = grep(paste0(softwarename, "_r", round), experimentNames)
    print(paste("SD for:", rsIndex, "range:", from, "-", to))
    if(peptides) rsIndex = grep(paste0(softwarename, "_peptides_r", round), experimentNames)
    if(length(rsIndex) == 0) return (NA)
    
    lrs = unlist(getPartialLogRatios( rsIndex, species, from, to))
    return(sd(lrs))
    #stdev = sd(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$y)
    #return(stdev)
}

get_Deviation <- function(softwarename=1, round=1, species="HUMAN", peptides=T, from=0, to=1){
    rsIndex = grep(paste0(softwarename, "_r", round), experimentNames)
    if(peptides) rsIndex = grep(paste0(softwarename, "_peptides_r", round), experimentNames)
    if(length(rsIndex) == 0) return (NA)

    lrs = unlist(getPartialLogRatios( rsIndex, species, from, to))
    return(abs(median(lrs) - ResultSets[[rsIndex]]$result[[1]]$data[[species]]$expectation))
    
    #devfromexpectation = abs(ResultSets[[rsIndex]]$result[[1]]$data[[species]]$median - 
    #                             ResultSets[[rsIndex]]$result[[1]]$data[[species]]$expectation)
    
    #return(devfromexpectation)
}

SD = as.data.frame(
    rbind(
        #softwarename=1, round=1, species="HUMAN", peptides=T, from=0, to=1
        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=T, from=2/3, to=1),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=T, from=2/3, to=1),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=T, from=2/3, to=1),

        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="HUMAN", peptides=F, from=2/3, to=1),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="YEAST", peptides=F, from=2/3, to=1),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=1, species="ECOLI", peptides=F, from=2/3, to=1),
        
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=T, from=2/3, to=1),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=T, from=2/3, to=1),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=T, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=T, from=2/3, to=1),
        
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="HUMAN", peptides=F, from=2/3, to=1),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="YEAST", peptides=F, from=2/3, to=1),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=F, from=0, to=1/3),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_SD, round=2, species="ECOLI", peptides=F, from=2/3, to=1)
        
        )
)


deviation = as.data.frame(
    rbind(
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=T, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=T, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=T, from=2/3, to=1),
           
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="HUMAN", peptides=F, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="YEAST", peptides=F, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=1, species="ECOLI", peptides=F, from=2/3, to=1),
       
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=T, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=T, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=T, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=T, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=T, from=2/3, to=1),
        
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="HUMAN", peptides=F, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="YEAST", peptides=F, from=2/3, to=1),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=F, from=0, to=1/3),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=F, from=1/3, to=2/3),
        sapply(software.names, get_Deviation, round=2, species="ECOLI", peptides=F, from=2/3, to=1)        
    )
)


quantile_names = c("1st.quantile", "2nd.quantile", "3rd.quantile")

names.for.SD.and.dev = c(
    paste(rep("HUMAN", 3), rep("it1", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("YEAST", 3), rep("it1", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("ECOLI", 3), rep("it1", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("HUMAN", 3), rep("it1", 3), rep("proteins",3), quantile_names, sep="_"),
    paste(rep("YEAST", 3), rep("it1", 3), rep("proteins",3), quantile_names, sep="_"),
    paste(rep("ECOLI", 3), rep("it1", 3), rep("proteins",3), quantile_names, sep="_"),

    paste(rep("HUMAN", 3), rep("it2", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("YEAST", 3), rep("it2", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("ECOLI", 3), rep("it2", 3), rep("peptides",3), quantile_names, sep="_"),
    paste(rep("HUMAN", 3), rep("it2", 3), rep("proteins",3), quantile_names, sep="_"),
    paste(rep("YEAST", 3), rep("it2", 3), rep("proteins",3), quantile_names, sep="_"),
    paste(rep("ECOLI", 3), rep("it2", 3), rep("proteins",3), quantile_names, sep="_")
    )

row.names(SD) = names.for.SD.and.dev
row.names(deviation) = names.for.SD.and.dev


write.table(x=SD, file.path(outDir, "Standard_deviations.tsv"), sep="\t", 
            row.names=T, col.names = T)
write.table(x=deviation, file.path(outDir, "Deviations_from_expectation.tsv"), sep="\t", 
            row.names=T, col.names = T)

