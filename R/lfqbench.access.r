#' plotResultSet
#' 
#' Function to generate the scatter plot
#' @param resSet
#' @param showScatterPlot
#' @param showScatterAndDensityPlot
#' @param showScatterAndBoxPlot
#' @param showBoxPlot 
#' @param showDensityPlot 

plotResultSet = function( resSet, 
                          showScatterPlot=T,
                          showScatterAndDensityPlot=T,
                          showScatterAndBoxPlot=T,
                          showBoxPlot=T,
                          showDensityPlot=T){
  cat("creating "+resSet$docSet$pdfFile+" ...")
  pdf(file=resSet$docSet$pdfFile, onefile=T, width=cfg$PlotWidth, height=cfg$PlotHeight, family="Helvetica", pointsize=9)
  
  # accuracy and precision
  if(showScatterPlot) sapply( resSet$result, showScatterPlot, showRegLines=T )
  if(showScatterAndDensityPlot) sapply( resSet$result, showScatterAndDensityPlot, showRegLines=T )
  if(showScatterAndBoxPlot) sapply( resSet$result, showScatterAndBoxPlot, showRegLines=T )
  if(showBoxPlot) sapply( resSet$result, showLogRatioBoxPlot )
  if(showDensityPlot) sapply( resSet$result, showDistributionDensityPlot )
    
  dev.off()
  cat("[done]\n")
}

#' saveMetrics 
#' 
#' This function save all the metrics related with the experiment 
#' 
#' @param resultSets
#' @param File 
#' 
saveMetrics = function(resultSets = ResultSets, File = cfg$LogFilesLocation + "/metrics.txt"){
  metrics = lapply( resultSets, getMetrics )
  sink(file=File, split=T)
  explainMetrics()
  x = lapply(metrics, showMetrics)
  sink()
}

#' explainMetrics 
#' 
#' This function verbose the metrics of the experiment 
#' 
explainMetrics = function(){
  # identification = identification rate, number of identified proteins for benchmark species
  # replication = replication variance, median CV for the background species
  # accuracy = accuracy of relative quantification, AUQC for background species
  # precision = precision of quantification, root mean square of ( medianLR - expectation ) for regulated species
  # separation = species separation ability, area under ROC curve between a species pair
  cat("--------------------------------------------\n")
  cat("\n Metrics description: \n")
  
  cat("\tidentification:\n")
  cat("\t\tidentification rate,\n")
  cat("\t\tthe number of identified proteins for benchmark species\n")
  
  cat("\treplication:\n")
  cat("\t\treplication variance,\n")
  cat("\t\tthe median CV for the background species\n")
  
  cat("\taccuracy:\n")
  cat("\t\taccuracy of relative quantification,\n")
  cat("\t\tthe median deviation of log-ratios to the expected value\n")
  
  cat("\tprecision:\n")
  cat("\t\tprecision of quantification,\n")
  cat("\t\tthe standard deviation of log-ratios\n")

  cat("\tseparation:\n")
  cat("\t\tspecies separation ability,\n")
  cat("\t\tthe area under ROC curve between a species pair\n")
}

#' showMetrics 
#' 
#' This function verbose the metrics of the experiment 
#' 
showMetrics = function(m){
  cat("--------------------------------------------\n")
  cat("\n" + m$name + "\n")
  sm = function(mn) 
  {
      cat( mn + ": " )
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
#' This function trace the Statistics related with each ID
#' 
#' @param resultSets

logIdStatistics = function( resultSets ){
  cat("storing identification statistics ... ")
  idstats = t( sapply( resultSets, function(rs) rs$idstat ) )
  sums = rowSums(idstats)
  labs = sapply( resultSets, function(rs) rs$docSet$fileBase )
  write.table(
    data.frame( "name"=labs, idstats, "sum"=sums ) ,
      file=cfg$LogFilesLocation+"/identification_statistics.csv",
      sep=cfg$CsvColumnSeparator, row.names=F, col.names=T)
  cat("[done]\n")
}

#' saveLogRations
#' 
#' This saveLogRations 
#' @param rs 
saveLogRatios = function( rs ){
  listLR = function(pairName)
  { 
    pairRes = rs$result[[pairName]]
    write.table( 
      pairRes$validLogRatios,
      cfg$LogFilesLocation + "/" + rs$docSet$fileBase + " LR ("+pairRes$name1 + "-" + pairRes$name2+").csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )
    
    write.table( 
      pairRes$validLogRatios[pairRes$validLogRatios$logratio < cfg$LogRatioPlotRange[1] | pairRes$validLogRatios$logratio > cfg$LogRatioPlotRange[2],],
      cfg$LogFilesLocation + "/" + rs$docSet$fileBase + " LR ("+pairRes$name1 + "-" + pairRes$name2+") out of range.csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )    
    
    write.table( 
      pairRes$invalidLogRatios,
      cfg$LogFilesLocation + "/" + rs$docSet$fileBase + " invalid LR ("+pairRes$name1 + "-" + pairRes$name2+").csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )

    return( pairRes$validLogRatios )
  }
  logRatios = lapply( names(rs$result), listLR )
  return(logRatios)
}

#' saveSampleMeans
#' 
#' This function save the Sample Means 
#' @param resultSet
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
    ), file=resultSet$docSet$avgFile, sep=cfg$CsvColumnSeparator, dec=cfg$CsvDecimalPointChar, col.names=T, row.names=F
  )
  cat("[done]\n")
}
################################################################################

#' saveSampleRSD
#' 
#' This function save the Sample RSD
#' @param resultSet
#' 
saveSampleRSD = function( resultSet )
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
    ), file=resultSet$docSet$rsdFile, sep=cfg$CsvColumnSeparator, dec=cfg$CsvDecimalPointChar, col.names=T, row.names=F
  )
  cat("[done]\n")
}

#' saveIDs
#' 
#' This function save the IDs 
#' @param resultSet

saveIDs = function( resultSet ){
  cat("storing identified protein names ... ")
  write.table( 
    resultSet$data$id
    , file=resultSet$docSet$idsFile, sep=cfg$CsvColumnSeparator, dec=cfg$CsvDecimalPointChar, col.names=F, row.names=F
  )
  cat("[done]\n")
}

#' saveSpeciesSeparation
#' 
#' This function enable to save the species 
#' @param resultSet
#' 
saveSpeciesSeparation = function( resultSet ){
  cat("storing species separation scores ... ")
  # species, a:b, a:c, ...
  d = data.frame(
    species=cfg$AllSpeciesPairsLabels,
    sapply( resultSet$result, function(d) d$separation)
  )
  names(d) = c("species", cfg$SamplePairsLabels)
  write.table( d, file=resultSet$docSet$rocFile, sep=cfg$CsvColumnSeparator, dec=cfg$CsvDecimalPointChar, col.names=T, row.names=F  )
  cat("[done]\n")
}

#' getMetrics 
#' 
#' This fucntion get single number metrics
#' @param resultSet
#' 
getMetrics = function(resultSet){
  # identification rate
  # number of identified proteins (only benchmark species)
  # PROBLEM? identification = sum( sapply(cfg$AllSpeciesNames, function(s) sum( resultSet$data$species == s, na.rm=T) ))
  identification = sum( sapply(cfg$AllSpeciesNames, function(s) length( which(resultSet$data$species == s) ) ) )
  # in case we want all IDs just count all
  # identificationRateMetric = length(resultSet$data$id)
  
  # replication variance
  # median CV for background-species
  replication = median( na.exclude( resultSet$data$cv[ resultSet$data$species==cfg$BackgroundSpeciesName ] ) )
    
  # area under ROC curve in predefined ranges of intensity
  separation = lapply(resultSet$result, function(d) d$rangedSeparation )
  
  # median log-ratio deviation to the expectation in predefined ranges of intensity
  accuracy = lapply(resultSet$result, function(d) d$rangedAccuracy )
  
  # standard deviation of log-ratios in predefined ranges of intensity
  precision = lapply(resultSet$result, function(d) d$rangedPrecision )
  
  # count log ratio statistics
  lrs4spc = function(spc, sp)
  {
    iLR = sp$invalidLogRatios$logratio[sp$invalidLogRatios$species==spc]
    vLR = sp$validLogRatios$logratio[sp$validLogRatios$species==spc]
    
    invalid = length(iLR)
    missing = ifelse(invalid==0, 0, length(which(is.na(iLR))) )
    
    valid = length(vLR)
    plotted = length( which( vLR >= min(cfg$LogRatioPlotRange) & vLR <= max(cfg$LogRatioPlotRange) ) )
    return( c(
        "invalid ratios"=invalid, 
        "invalid, out of validity range"=(invalid - missing),
        "invalid, missing value"=missing,
        "valid ratios"=valid,
        "in plot range"=plotted, 
        "out of plot range"=(valid - plotted)
    ) )
  }
  
  logRatioStats = lapply(resultSet$result, function(sp) sapply( cfg$AllSpeciesNames, lrs4spc, sp  ))
    
  return(
    list(
      name = resultSet$docSet$fileBase,
      identification = identification,
      replication = replication,
      accuracy = accuracy,
      precision = precision,
      separation = separation,
      quantification = logRatioStats
    )
  )
}

#' getCombinedLogRatios
#' 
#' This function get the Combined Log Ratios 
#' 
#' @param ResultSets 
#' 
getCombinedLogRatios = function( ResultSets ){
  listLR = function(pairName, rs)
  { 
    fileName = rs$docSet$fileBase
    spcNames = names( rs$result[[pairName]]$data )
    spcLogRatios = sapply( rs$result[[pairName]]$data, function(d) d$y )
    names(spcLogRatios) = fileName + "  " + pairName + "  " + spcNames
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
#' 
getCombinedProteinRSDs = function( ResultSets ){
  CVs = lapply( ResultSets, function(rs) as.vector(rs$data$cv) )
  names(CVs) = as.vector( lapply(ResultSets, function(rs) rs$docSet$fileBase ) )
  return(CVs)
}

#' makeDocSet
#' 
#' This function construct the DocSet 
#' @param inFile 
#' 
makeDocSet = function( inFile ){
  fileBase = sub(cfg$InputExtensionPattern, "", inFile)
  return( 
    list(
      inputPath = cfg$InputFilesLocation,
      plotPath = cfg$PlotFilesLocation,
      logPath = cfg$LogFilesLocation,
      fileBase = fileBase,
      csvFile = cfg$InputFilesLocation + "/" + inFile,
      pdfFile = cfg$PlotFilesLocation + "/" + fileBase + ".pdf",
      logFile = cfg$LogFilesLocation + "/" + fileBase + ".log",
      avgFile = cfg$LogFilesLocation + "/" + fileBase + " sample_means.csv",
      rsdFile = cfg$LogFilesLocation + "/" + fileBase + " cv.csv",
      l2rFile = cfg$LogFilesLocation + "/" + fileBase + " log2_ratio.csv",
      rocFile = cfg$LogFilesLocation + "/" + fileBase + " species_separation.csv",
      idsFile = cfg$LogFilesLocation + "/" + fileBase + " ids.csv"
    )
  )
}