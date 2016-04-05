#' LFQbench.plotResultSet
#' 
#' Function to generate the scatter plot
#' @param resSet a result set
#' @param showScatterPlot
#' @param showScatterAndDensityPlot
#' @param showScatterAndBoxPlot
#' @param showBoxPlot 
#' @param showDensityPlot
#' @export
LFQbench.plotResultSet = function(resSet,
                         showScatterPlot=T,
                         showScatterAndDensityPlot=T,
                         showScatterAndBoxPlot=T,
                         showBoxPlot=T,
                         showDensityPlot=T){
    
    pdfFile = resSet$docSet$pdfFile
    
    if(!is.null(pdfFile)){
      cat(paste("creating ", pdfFile, " ...", sep = ""))
      pdf(file = pdfFile, onefile = T, width = LFQbench.Config$PlotWidth, height = LFQbench.Config$PlotHeight, family = "Helvetica", pointsize = 9)      
    }
    
    # accuracy and precision
    if(showScatterPlot) sapply( resSet$result, LFQbench.showScatterPlot, showRegLines = T )
    if(showScatterAndDensityPlot) sapply( resSet$result, LFQbench.showScatterAndDensityPlot, showRegLines = T )
    if(showScatterAndBoxPlot) sapply( resSet$result, LFQbench.showScatterAndBoxPlot, showRegLines = T )
    if(showBoxPlot) sapply( resSet$result, LFQbench.showLogRatioBoxPlot )
    if(showDensityPlot) sapply( resSet$result, LFQbench.showDistributionDensityPlot )
    
    if(!is.null(pdfFile)) dev.off()
    cat("[done]\n")
}

#' saveMetrics 
#' 
#' This function save all the metrics related with the experiment 
#' @param resultSets
#' @param File 
saveMetrics = function(resultSets = ResultSets, File = paste(LFQbench.Config$LogFilesLocation, "/metrics.txt", sep = "")){
  metrics = lapply( resultSets, LFQbench.getMetrics )
  sink(file=File, split=T)
  LFQbench.explainMetrics()
  x = lapply(metrics, LFQbench.showMetrics)
  sink()
}

#' LFQbench.explainMetrics 
#' 
#' This function verbose the metrics of the experiment 
#' @export
LFQbench.explainMetrics = function(){
    cat("--------------------------------------------------------------------------------\n")
    cat("# Metrics description: \n\n")
    
    cat("Identification statistics:\n")
    cat("  the number of identifications for benchmark species\n\n")
    
    cat("Quantification statistics:\n")
    cat("   valid ratios - total number of valid log-ratios that meet following conditions:\n")
    cat("    a) peptide/protein could be quantified in both samples.\n")
    cat("    b) calulcated log-ratio falls within the expectation range\n")
    cat("       defined by the sample composition and the spread of spread of the log-ratios observed for each species.\n")
    cat(paste0(
        "    The expectation range is defined as a miximum difference of ",
        LFQbench.Config$LogRatioValidityRangeSDFactor," standard deviations from the average ratio.\n"))
    
    cat("   in plot range - numbers of valid log-ratios inside of the user defined plot range\n")
    cat("   out of plot range - numbers of valid log-ratios outside of the user defined plot range\n")
    cat("   invalid ratios - total numbers of invalid log-ratios\n")
    cat("   invalid ratios, out of validity range - proteins/peptides outside of validity range (see condition b)\n")
    cat("   invalid ratios, missing value - proteins/peptides exclusively quantified in one of the samples (see condition a)\n\n")
    
    cat("Technical variance:\n")
    cat("  the median CV for the background species among replicate runs\n\n")
    
    cat("Global accuracy:\n")
    cat("  accuracy of relative quantification,\n")
    cat("  the median deviation of log-ratios to the expected value\n\n")
    
    cat("Global precision:\n")
    cat("  precision of quantification,\n")
    cat("  the standard deviation of log-ratios\n\n")
    
    cat("Global species overlap:\n")
    cat(" species separation ability,\n")
    cat(" the area under ROC curve between a species pair.\n\n")
    
    cat("Local accuracy:\n")
    cat(" local accuracy of relative quantification,\n")
    cat(" the median deviation of log-ratios to the expected value\n")
    cat(" for partial data split in quantiles by intensity in first sample\n\n")
    
    cat("Local precision:\n")
    cat(" precision of quantification,\n")
    cat(" the standard deviation of log-ratios\n")
    cat(" for partial data split in quantiles by intensity in first sample\n\n")
    
    cat("Local species overlap:\n")
    cat(" species separation ability,\n")
    cat(" the area under ROC curve between a species pair\n")
    cat(" for partial data split in quantiles by intensity in first sample\n")
    cat("--------------------------------------------------------------------------------\n")
}

#' LFQbench.showMetrics 
#' 
#' verbose the metrics of the experiment
#' 
#' @param m the metrics object
#' @export
LFQbench.showMetrics = function(m){
  cat("--------------------------------------------\n")
  cat(paste("\n", m$name, "\n", sep=""))
  sm = function(mn){
      cat( paste(mn,": ", sep="" ))
      v = m[[mn]]
      if( is.atomic(v) ) 
      {
        cat(v)
      }
      else 
      {
        cat("\n")
        show(v)
      }
      cat("\n")
  }
  nix = sapply(names(m), sm)
  cat("\n")
}

#' logIdStatistics
#' 
#' store identification statistics from a list of result sets
#' 
#' @param resultSets the result sets
#' @param csvFile the delimiter separated output file
logIdStatistics = function( resultSets, csvFile=paste(LFQbench.Config$LogFilesLocation,"/identification_statistics.csv",sep="") ){
  cat("storing identification statistics ... ")
  idstats = t( sapply( resultSets, function(rs) rs$idstat ) )
  sums = rowSums(idstats)
  labs = sapply( resultSets, function(rs) rs$docSet$fileBase )
  write.table(
      data.frame( "name"=labs, idstats, "sum"=sums ),
      file=csvFile, 
      sep=LFQbench.Config$CsvColumnSeparator, 
      row.names=F, col.names=T
      )
  cat("[done]\n")
}

#' saveLogRatios
#' 
#' extract and store calculated log-ratios from a list of result sets
#' @param rs the result set
saveLogRatios = function(rs){
    listLR = function(pairName)
    { 
        pairRes = rs$result[[pairName]]
        write.table( 
          pairRes$validLogRatios,
          paste(rs$docSet$logPath,"/",rs$docSet$fileBase," LR (", pairRes$name1,"-" ,pairRes$name2, ").csv",sep = ""),
          sep = ",",dec = ".",row.names = F, col.names = T
        )
        
        write.table( 
          pairRes$validLogRatios[pairRes$validLogRatios$logratio < LFQbench.Config$LogRatioPlotRange[1] | pairRes$validLogRatios$logratio > LFQbench.Config$LogRatioPlotRange[2],],
          paste(rs$docSet$logPath,"/", rs$docSet$fileBase ," LR (", pairRes$name1 , "-" , pairRes$name2, ") out of range.csv", sep=""),
          sep = ",",dec = ".",row.names = F, col.names = T
        )    
        
        write.table( 
          pairRes$invalidLogRatios,
          paste(rs$docSet$logPath ,"/" , rs$docSet$fileBase , " invalid LR (", pairRes$name1 , "-" , pairRes$name2, ").csv",sep=""),
          sep = ",",dec = ".",row.names = F, col.names = T
        )
        
        return( pairRes$validLogRatios )
    }
    logRatios = lapply( names(rs$result), listLR )
    return(logRatios)
}

#' saveSampleMeans
#' 
#' save sample means to file
#' @param resultSet the result set
#' 
saveSampleMeans = function( resultSet )
{
    cat("storing sample means ... ")
    d = resultSet$data
    # id, species, a, b, c
    write.table( 
    data.frame(
      id=d$id,
      species=d$species,
      d$mean
    ), file=resultSet$docSet$avgFile, sep=LFQbench.Config$CsvColumnSeparator, dec=LFQbench.Config$CsvDecimalPointChar, col.names=T, row.names=F
    )
    cat("[done]\n")
}
################################################################################

#' saveSampleCVs
#' 
#' save the in sample CVs to a file
#' @param resultSet the result set
saveSampleCVs = function( resultSet)
{
  cat("storing coeffiecients of in-sample variation ... ")
  d = resultSet$data
  # id, species, a, b, c
  write.table( 
    data.frame(
      id=d$id,
      species=d$species,
      d$cv,
      row.names = NULL
    ), file=resultSet$docSet$rsdFile, sep=LFQbench.Config$CsvColumnSeparator, dec=LFQbench.Config$CsvDecimalPointChar, col.names=T, row.names=F
  )
  cat("[done]\n")
}

#' saveIDs
#' 
#' save identified proteins/peptides to a file
#' @param resultSet the result set
saveIDs = function( resultSet ){
  cat("storing identified protein names ... ")
  write.table( 
    resultSet$data$id
    , file = resultSet$docSet$idsFile, sep=LFQbench.Config$CsvColumnSeparator, dec = LFQbench.Config$CsvDecimalPointChar, col.names = F, row.names = F
  )
  cat("[done]\n")
}

#' saveSpeciesSeparation
#' 
#' This function enable to save the species 
#' @param resultSet
saveSpeciesSeparation = function( resultSet ){
  cat("storing species separation scores ... ")
  # species, a:b, a:c, ...
  d = data.frame(
    species=LFQbench.Config$AllSpeciesPairsLabels,
    sapply( resultSet$result, function(d) d$separation)
  )
  names(d) = c("species", LFQbench.Config$SamplePairsLabels)
  write.table( d, file = resultSet$docSet$rocFile, sep = LFQbench.Config$CsvColumnSeparator, dec = LFQbench.Config$CsvDecimalPointChar, col.names = T, row.names = F  )
  cat("[done]\n")
}

#' LFQbench.getMetrics 
#' 
#' calculate metrics for a result set
#' @param resultSet
#' @export
LFQbench.getMetrics = function(resultSet){
  # identification rate
  # number of identified proteins (only benchmark species)
  # PROBLEM? identification = sum( sapply(LFQbench.Config$AllSpeciesNames, function(s) sum( resultSet$data$species == s, na.rm=T) ))
  identification = sum( sapply(LFQbench.Config$AllSpeciesNames, function(s) length( which(resultSet$data$species == s) ) ) )
  # in case we want all IDs just count all
  # identificationRateMetric = length(resultSet$data$id)
  
  # replication variance
  # median CV for background-species
  replication = median( na.exclude( resultSet$data$cv[ resultSet$data$species == LFQbench.Config$BackgroundSpeciesName ] ) )
    
  # area under ROC curve in predefined quantiles of intensity
  globalSeparation = lapply(resultSet$result, function(d) d$separation )
  localSeparation = lapply(resultSet$result, function(d) d$quantileSeparation )  
  
  # median log-ratio deviation to the expectation in predefined quantiles of intensity
  globalAccuracy = lapply(resultSet$result, function(d) d$accuracy )
  localAccuracy = lapply(resultSet$result, function(d) d$quantileAccuracy )
  
  # standard deviation of log-ratios in predefined quantiles of intensity
  globalPrecision = lapply(resultSet$result, function(d) d$precision )
  localPrecision = lapply(resultSet$result, function(d) d$quantilePrecision )
  
  # count log ratio statistics
  lrs4spc = function(spc, sp)
  {
    iLR = sp$invalidLogRatios$logratio[sp$invalidLogRatios$species == spc]
    vLR = sp$validLogRatios$logratio[sp$validLogRatios$species == spc]
    
    invalid = length(iLR)
    missing = ifelse(invalid == 0, 0, length(which(is.na(iLR))) )
    
    valid = length(vLR)
    plotted = length( which( vLR >= min(LFQbench.Config$LogRatioPlotRange) & vLR <= max(LFQbench.Config$LogRatioPlotRange) ) )
    return( c(
        "invalid ratios" = invalid, 
        "invalid, out of validity range" = (invalid - missing),
        "invalid, missing value" = missing,
        "valid ratios" = valid,
        "in plot range" = plotted, 
        "out of plot range" = (valid - plotted)
    ) )
  }
  
  logRatioStats = lapply(resultSet$result, function(sp) sapply( LFQbench.Config$AllSpeciesNames, lrs4spc, sp  ))
    
  return(
    list(
      "name" = resultSet$docSet$fileBase,
      "Identification statistics" = identification,
      "Quantification statistics" = logRatioStats,
      "Technical variance" = replication,
      "Global accuracy" = globalAccuracy,
      "Global precision" = globalPrecision,
      "Global species overlap" = globalSeparation,
      "Local accuracy" = localAccuracy,
      "Local precision" = localPrecision,
      "Local species overlap" = localSeparation
    )
  )
}

#' getCombinedLogRatios
#' 
#' This function get the Combined Log Ratios 
#' 
#' @param ResultSets 
getCombinedLogRatios = function( ResultSets ){
  listLR = function(pairName, rs)
  { 
    fileName = rs$docSet$fileBase
    spcNames = names( rs$result[[pairName]]$data )
    spcLogRatios = sapply( rs$result[[pairName]]$data, function(d) d$y )
    names(spcLogRatios) = paste(fileName , "  " , pairName , "  " , spcNames, sep = "")
    return(spcLogRatios)
  }
  
  logRatios = unlist( 
    sapply( ResultSets, function(rs) lapply(names(rs$result), listLR, rs) )
    , recursive=F 
  )
  
  return(logRatios)
}

#' getCombinedProteinRSDs 
#' 
#' This fucntion return the combined protein RSDs 
#' @param ResultsSets 
getCombinedProteinRSDs = function(ResultSets){
  CVs = lapply( ResultSets, function(rs) as.vector(rs$data$cv) )
  names(CVs) = as.vector( lapply(ResultSets, function(rs) rs$docSet$fileBase ) )
  return(CVs)
}

#' makeDocSet
#' 
#' construct a document set for a file by paths from user configuration
#' @param inFile
makeDocSet = function( inFile, inputPath=LFQbench.Config$InputFilesLocation, plotPath=LFQbench.Config$PlotFilesLocation, logPath=LFQbench.Config$LogFilesLocation ){
  fileBase = sub(LFQbench.Config$InputExtensionPattern, "", basename(inFile) )
  return( 
    list(
      inputPath = inputPath, 
      plotPath = plotPath, 
      logPath  = logPath,
      fileBase = fileBase, 
      csvFile  = ifelse(file.exists(inFile), inFile, paste(inputPath,"/", inFile, sep = "")),
      pdfFile  = paste(plotPath,"/", fileBase,".pdf", sep = ""),
      logFile  = paste(logPath,"/", fileBase,".log", sep = ""),
      avgFile  = paste(logPath, "/", fileBase, "_sample_means.csv", sep = ""),
      rsdFile  = paste(logPath, "/", fileBase, "_cv.csv", sep = ""),
      l2rFile  = paste(logPath, "/", fileBase, "_log2_ratio.csv", sep = ""),
      rocFile  = paste(logPath, "/", fileBase, "_species_separation.csv", sep = ""),
      idsFile  = paste(logPath, "/", fileBase, "_ids.csv", sep = "")
    )
  )
}

#' makeEverythingInOneFolderDocSet
#' 
#' construct a doc set using the input file's path
#' @param inFile
makeEverythingInOneFolderDocSet = function( inFile ){
    tarDir = dirname(inFile)
    docSet = makeDocSet(inFile = inFile, inputPath=tarDir, plotPath=tarDir, logPath=tarDir)
    return(docSet)
}