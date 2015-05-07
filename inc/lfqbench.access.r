################################################################################
plotResultSet = function( resSet )
{
  cat("creating "+resSet$docSet$pdfFile+" ...")
  pdf(file=resSet$docSet$pdfFile, onefile=T, width=PlotWidth, height=PlotHeight, family="Helvetica", pointsize=9)
  
  # replication variance
  # plotProteinDispersionBySpecies( resSet$data )
  # plotProteinDispersionBySample( resSet$data )
  
  # accuracy and precision
  sapply( resSet$result, showScatterPlot )
  sapply( resSet$result, showScatterAndDensityPlot, showRegLines=T )
  sapply( resSet$result, showScatterAndBoxPlot, showRegLines=T )
  sapply( resSet$result, showLogRatioBoxPlot )
  sapply( resSet$result, showQuantBarPlot )
  sapply( resSet$result, showDistributionDensityPlot )
  showSingleSpeciesQC( resSet$result, BackgroundSpeciesName )
  
  dev.off()
  cat("[done]\n")
}
################################################################################

################################################################################
saveMetrics = function(resultSets = ResultSets, File = LogFilesLocation + "/metrics.txt")
{
  metrics = lapply( resultSets, getMetrics )
  sink(file=File, split=T)
  explainMetrics()
  x = lapply(metrics, showMetrics)
  sink()
}
################################################################################

################################################################################
explainMetrics = function()
{
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
  cat("\t\tthe area under quantification curve for background species\n")
  
  cat("\tprecision deviation:\n")
  cat("\t\tdeviation of precision of quantification,\n")
  cat("\t\tthe root mean square of ( medianLR - expectation ) for regulated species\n")
  
  cat("\tseparation:\n")
  cat("\t\tspecies separation ability,\n")
  cat("\t\tthe area under ROC curve between a species pair\n")
}
################################################################################

################################################################################
showMetrics = function(m)
{
  cat("--------------------------------------------\n")
  cat("\n" + m$name + "\n")
  cat("identification: "+ m$identification +"\n")
  cat("replication: "+ m$replication +"\n")
  cat("accuracy: \n")
  show(m$accuracy)
  cat("precision deviation: "+ m$precision +"\n")
  cat("separation: \n")
  show(m$separation)
  cat("\n")
}
################################################################################

################################################################################
logIdStatistics = function( resultSets )
{
  cat("storing identification statistics ... ")
  idstats = t( sapply( resultSets, function(rs) rs$idstat ) )
  sums = rowSums(idstats)
  labs = sapply( resultSets, function(rs) rs$docSet$fileBase )
  write.table(
    data.frame( "name"=labs, idstats, "sum"=sums ) ,
      file=LogFilesLocation+"/identification_statistics.csv",
      sep=CsvColumnSeparator, row.names=F, col.names=T)
  cat("[done]\n")
}
################################################################################

################################################################################
saveLogRatios = function( rs )
{
  listLR = function(pairName)
  { 
    pairRes = rs$result[[pairName]]
    write.table( 
      pairRes$validLogRatios,
      LogFilesLocation + "/" + rs$docSet$fileBase + " LR ("+pairRes$name1 + "-" + pairRes$name2+").csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )
    
    write.table( 
      pairRes$validLogRatios[pairRes$validLogRatios$logratio < LogRatioPlotRange[1] | pairRes$validLogRatios$logratio > LogRatioPlotRange[2],],
      LogFilesLocation + "/" + rs$docSet$fileBase + " LR ("+pairRes$name1 + "-" + pairRes$name2+") out of range.csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )    
    
    write.table( 
      pairRes$invalidLogRatios,
      LogFilesLocation + "/" + rs$docSet$fileBase + " invalid LR ("+pairRes$name1 + "-" + pairRes$name2+").csv",
      sep = ",",dec = ".",row.names = F, col.names = T
    )

    return( pairRes$validLogRatios )
  }
  logRatios = lapply( names(rs$result), listLR )
  return(logRatios)
}
################################################################################

################################################################################
# store sample means
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
    ), file=resultSet$docSet$avgFile, sep=CsvColumnSeparator, dec=CsvDecimalPointChar, col.names=T, row.names=F
  )
  cat("[done]\n")
}
################################################################################

################################################################################
# store sample cv
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
    ), file=resultSet$docSet$rsdFile, sep=CsvColumnSeparator, dec=CsvDecimalPointChar, col.names=T, row.names=F
  )
  cat("[done]\n")
}
################################################################################

################################################################################
saveIDs = function( resultSet )
{
  cat("storing identified protein names ... ")
  write.table( 
    resultSet$data$id
    , file=resultSet$docSet$idsFile, sep=CsvColumnSeparator, dec=CsvDecimalPointChar, col.names=F, row.names=F
  )
  cat("[done]\n")
}
################################################################################

################################################################################
# store species separation
saveSpeciesSeparation = function( resultSet )
{
  cat("storing species separation scores ... ")
  # species, a:b, a:c, ...
  d = data.frame(
    species=AllSpeciesPairsLabels,
    sapply( resultSet$result, function(d) d$separation)
  )
  names(d) = c("species", SamplePairsLabels)
  write.table( d, file=resultSet$docSet$rocFile, sep=CsvColumnSeparator, dec=CsvDecimalPointChar, col.names=T, row.names=F  )
  cat("[done]\n")
}
################################################################################

################################################################################
# get single number metrics
getMetrics = function(resultSet)
{
  # identification rate
  # number of identified proteins (only benchmark species)
  # PROBLEM? identification = sum( sapply(AllSpeciesNames, function(s) sum( resultSet$data$species == s, na.rm=T) ))
  identification = sum( sapply(AllSpeciesNames, function(s) length( which(resultSet$data$species == s) ) ) )
  # in case we want all IDs just count all
  # identificationRateMetric = length(resultSet$data$id)
  
  # replication variance
  # median CV for background-species
  replication = median( na.exclude( resultSet$data$cv[ resultSet$data$species==BackgroundSpeciesName ] ) )
  
  # accuracy of relative quantification
  # AUQC of background species
  accuracy = sapply(resultSet$result, function(d) d$data[[BackgroundSpeciesName]]$auqc )
  
  # precision of quantification
  # root mean square ( medianLR - expectation ) for regulated species
  idx = which(AllSpeciesNames != BackgroundSpeciesName)
  dat = sapply( resultSet$result, function( sp ) sapply( idx, function(x) sp$data[[x]]$shift ) )
  precision = rms( dat )
  
  # area under ROC curve
  separation = sapply(resultSet$result, function(d) d$separation)
  
  return(
    list(
      name = resultSet$docSet$fileBase,
      identification = identification,
      replication = replication,
      accuracy = accuracy,
      precision = precision,
      separation = separation
    )  
  )
}
################################################################################

################################################################################
getCombinedLogRatios = function( ResultSets )
{
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
################################################################################

################################################################################
getCombinedProteinRSDs = function( ResultSets )
{
  CVs = lapply( ResultSets, function(rs) as.vector(rs$data$cv) )
  names(CVs) = as.vector( lapply(ResultSets, function(rs) rs$docSet$fileBase ) )
  return(CVs)
}
################################################################################
