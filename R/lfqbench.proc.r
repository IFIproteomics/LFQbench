#' LFQbench.processDocSet
#' 
#' This function process the data before generate the charts
#' @param DocSet is the aata to be processed 
LFQbench.processDocSet = function( DocSet ) {
    rs = LFQbench.processFile( file = DocSet$csvFile )
    rs$docSet = DocSet
    return(rs)
}

#' LFQbench.processFile
#' 
#' This function process the data and generates a result set
#' @param file is the data to be processed
#' @export
LFQbench.processFile = function( file ) {
  if(!exists("DEBUG")) DEBUG <<- F
  cat(paste("processing ", basename(file) ," ... \n",sep = ""))
  
  ################################################################################
  # read run based data
  csvDat = read.table(file, sep=LFQbench.Config$CsvColumnSeparator, dec=LFQbench.Config$CsvDecimalPointChar, header=T)
  
  # csv parsing problems lead to everything being in a single column
  if( ncol(csvDat) < 2 )
    stop("wrong input format! (please check column separator settings)")
  
  # ids are always the first column
  csvIDs = as.vector( csvDat[,1] )
  
  # every numeric column belongs to the amount data
  numericColumnIdx = sapply(csvDat, class) == "numeric" | sapply(csvDat, class) == "integer"
  
  # species column name must start with "speci"
  speciesColumnIdx = grep("speci", names(csvDat), ignore.case = T)
  
  # some checks for data fitness
  if( length( which( numericColumnIdx ) ) < LFQbench.Config$NumberOfSamples )
    stop("number of amount columns is smaller than the number of samples")
  
  if( length( speciesColumnIdx ) < 1 ) 
    stop("no ''species'' column found.")
  
  # 1.0 * is a work-around to force the amounts matrix to numeric class
  csvAmounts = 1.0 * as.matrix( csvDat[, numericColumnIdx] )
  csvSpecies = as.vector( csvDat[, speciesColumnIdx[[1]]] )  
  rm( csvDat )
  ################################################################################
  
  ################################################################################
  # disable low amounts
  csvAmounts[csvAmounts < LFQbench.Config$MinProteinAmount] = NA
  ################################################################################
  
  ################################################################################
  # calculate missing value statistics
  calHistNAsWithSpecies <- function(species){
    sp_hist = hist(csvNumNAs[csvSpecies == species],  breaks=seq(-1, ncol(csvAmounts)), include.lowest=T,  plot=F  )
    return(sp_hist$counts)
  }
  csvNumNAs = apply(csvAmounts,1, function(row) sum(is.na(row)))
  numNAs_histograms = as.data.frame(sapply(LFQbench.Config$AllSpeciesNames, calHistNAsWithSpecies)) 
  numNAs_histograms = cbind( seq(0, ncol(csvAmounts) ), numNAs_histograms)
  names(numNAs_histograms)[1] = c("numNAs")
  ################################################################################
  
  ################################################################################
  # check OTHER species
  benchedSpeciesIndices = unlist(sapply( LFQbench.Config$AllSpeciesNames, function(spc) which(csvSpecies==spc) ), use.names = F)
  otherSpeciesIndices = (1:length(csvIDs))[ -benchedSpeciesIndices ]
  otherSpeciesIDs = csvIDs[ otherSpeciesIndices ]
  NumberOfProteinsForOtherSpecies = length(otherSpeciesIndices)
  # limit data to benched species only
  csvIDs = csvIDs[ benchedSpeciesIndices ]
  csvAmounts = csvAmounts[ benchedSpeciesIndices, ]
  csvSpecies = csvSpecies[ benchedSpeciesIndices ]
  ################################################################################
  # normalize protein amounts to ppm
  if(LFQbench.Config$NormalizeAmountsToPPM) csvAmounts = apply( csvAmounts, 2, function(x) x / sum(x, na.rm=T) * 1000000 )
  ################################################################################
  
  ################################################################################
  # extract data for a single sample
  getSampleAverage = function(sampleIndex, nRuns=floor( ncol(csvAmounts) / LFQbench.Config$NumberOfSamples ))
  {
    amount = as.matrix( csvAmounts[, ((sampleIndex-1)*nRuns + 1):(nRuns * sampleIndex)] )
    avg = rowMeans( amount, na.rm=T )
    std = apply( amount, 1, sd, na.rm=T )
    rsd = std / avg
    return( list( mean=avg, cv=rsd ) )
  }
  ################################################################################
  # generate sample average data
  SampleAverageData = lapply( 1:LFQbench.Config$NumberOfSamples, getSampleAverage )
  SampleAverageAmounts = sapply( 1:LFQbench.Config$NumberOfSamples, function(si) SampleAverageData[[si]]$mean )
  SampleAverageCVs = sapply( 1:LFQbench.Config$NumberOfSamples, function(si) SampleAverageData[[si]]$cv )
  SampleAverageSpecies = csvSpecies
  SampleAverageProteinIDs = csvIDs
  rownames(SampleAverageAmounts) = SampleAverageProteinIDs
  rownames(SampleAverageCVs) = SampleAverageProteinIDs
  colnames(SampleAverageAmounts) = LFQbench.Config$AllSampleNames
  colnames(SampleAverageCVs) = LFQbench.Config$AllSampleNames
  ################################################################################
  
  ################################################################################
  # create protein identification statistics
  NumberOfProteinsBySpecies = sapply( LFQbench.Config$AllSpeciesNames, function(x) length( which(SampleAverageSpecies == x) ) )
  NumberOfProteinsBySpecies["OTHER"] = NumberOfProteinsForOtherSpecies
  ################################################################################
  
  ################################################################################
  # handle missing and low amount values
  # SampleAverageAmounts[ is.na(SampleAverageAmounts) ] = LFQbench.Config$MinProteinAmount
  # SampleAverageAmounts[ SampleAverageAmounts < LFQbench.Config$MinProteinAmount ] = LFQbench.Config$MinProteinAmount
  ################################################################################
  
  ################################################################################
  # process a sample pair
  getSamplePairData = function( SamplePairIndex = 1 )
  {
    ################################################################################
    # extract sample pair data
    SamplePairColor = LFQbench.Config$SamplePairsColors[SamplePairIndex]
    Sample1Index = LFQbench.Config$SamplePairsIndices[SamplePairIndex,1]
    Sample2Index = LFQbench.Config$SamplePairsIndices[SamplePairIndex,2]
    Sample1Name = LFQbench.Config$AllSampleNames[Sample1Index]
    Sample2Name = LFQbench.Config$AllSampleNames[Sample2Index]
    Sample1ProteinAmounts = SampleAverageAmounts[,Sample1Index]
    Sample2ProteinAmounts = SampleAverageAmounts[,Sample2Index]
    SpeciesNames = SampleAverageSpecies
    ProteinIDs = SampleAverageProteinIDs
    LogRatioExpectations = log2( LFQbench.Config$AllExpectedAmounts[,Sample1Index] / LFQbench.Config$AllExpectedAmounts[,Sample2Index] )
    ################################################################################
    
    ################################################################################
    # calculate log ratios
    LogRatio = log2( Sample1ProteinAmounts  / Sample2ProteinAmounts )
    ################################################################################
    
    invalidLogRatios = NULL
    allLogRatios = data.frame(
      entry = ProteinIDs,
      species = SpeciesNames,
      logratio = LogRatio,
      sample1means = Sample1ProteinAmounts,
      sample2means = Sample2ProteinAmounts
    )   
    
    ################################################################################
    # drop invalid log ratios
    if(LFQbench.Config$DropInvalidLogRatio)
    {
      emptyLRs = which(is.na(LogRatio))
      smallLRs = which(LogRatio < min(LFQbench.Config$LogRatioValidityRange))
      bigLRs = which(LogRatio > max(LFQbench.Config$LogRatioValidityRange))      
      badLRs = unique(c(emptyLRs, smallLRs, bigLRs))
      if(length(badLRs)>0)
      {
        invalidLogRatios = data.frame(
  	      	entry = ProteinIDs[badLRs],
  	      	species = SpeciesNames[badLRs],
  	      	logratio = LogRatio[badLRs],
  	      	sample1means = Sample1ProteinAmounts[badLRs],
  	      	sample2means = Sample2ProteinAmounts[badLRs],
  	      	row.names = NULL,
            stringsAsFactors=F
  	    )
        LogRatio = LogRatio[-badLRs]
        SpeciesNames = SpeciesNames[-badLRs]
        ProteinIDs = ProteinIDs[-badLRs]
        Sample1ProteinAmounts = Sample1ProteinAmounts[-badLRs]
        Sample2ProteinAmounts = Sample2ProteinAmounts[-badLRs]
      }
    }
    ################################################################################
    
    ################################################################################
    # calculate log-ratio means and medians
    LogRatioMedians = sapply(LFQbench.Config$AllSpeciesNames, function(x) median( LogRatio[SpeciesNames==x], na.rm=T ))
    LogRatioMeans = sapply(LFQbench.Config$AllSpeciesNames, function(x) mean( LogRatio[SpeciesNames==x], na.rm=T ))
    BackgroundSpeciesMedian = LogRatioMedians[LFQbench.Config$AllSpeciesNames==LFQbench.Config$BackgroundSpeciesName]
    BackgroundSpeciesMean = LogRatioMeans[LFQbench.Config$AllSpeciesNames==LFQbench.Config$BackgroundSpeciesName]
    ################################################################################
    
    ################################################################################
    # center by mean, median, or not at all
    if( LFQbench.Config$CenterLogRatioByBackground || regexpr("median", LFQbench.Config$CenterLogRatioByBackground, ignore.case=T) )
      LogRatioAdjustmentValue = BackgroundSpeciesMedian 
    else if(regexpr("mean", LFQbench.Config$CenterLogRatioByBackground, ignore.case=T))
      LogRatioAdjustmentValue = BackgroundSpeciesMean
    else
      LogRatioAdjustmentValue = 0
    
    ################################################################################
    # center log ratios
    if(LogRatioAdjustmentValue != 0)
    {
      LogRatio = LogRatio - LogRatioAdjustmentValue
      LogRatioMedians = LogRatioMedians - LogRatioAdjustmentValue
      LogRatioMeans = LogRatioMeans - LogRatioAdjustmentValue
    }
    ################################################################################
    
    validLogRatios = data.frame(
      entry = ProteinIDs,
      species = SpeciesNames,
      logratio = LogRatio,
      sample1means = Sample1ProteinAmounts,
      sample2means = Sample2ProteinAmounts,
      row.names = NULL, stringsAsFactors=F
    )
    
    ################################################################################
    # define and clip scatterplot x-axis data
    ScatterPlotXAxisData = log2( Sample2ProteinAmounts )
    # xLim = range(ScatterPlotXAxisData[ScatterPlotXAxisData>0])
    xLim = quantile( ScatterPlotXAxisData[ScatterPlotXAxisData>0], probs=c(0.01,0.99), na.rm = T )
    if(!is.null(LFQbench.Config[["LogIntensityPlotRange"]])) xLim = LFQbench.Config$LogIntensityPlotRange
    yLim = LFQbench.Config$LogRatioPlotRange
    # DISABLED: force log-ratio values into x-axis boundaries
    # ScatterPlotXAxisData[ ScatterPlotXAxisData < xLim[1] ] = xLim[1]
    # ScatterPlotXAxisData[ ScatterPlotXAxisData > xLim[2] ] = xLim[2]
    ################################################################################
    
    ################################################################################
    # package data for a species
    packageSpecies = function( TheSpecies )
    {
      ValueIndices = which(SpeciesNames == TheSpecies)
      SpeciesIndex = LFQbench.Config$AllSpeciesNames == TheSpecies
      
      if(DEBUG) cat( paste(TheSpecies, " has ", length( ValueIndices ), " IDs ... \n", sep = ""))
      
      qcFunc = as.function(
        getQCFunction( LogRatio[ValueIndices] - LogRatioMedians[SpeciesIndex], ensureValueRange=c(0, LFQbench.Config$MaxLogRatioForAUQC) )
      )
      auqc = round( integrate( qcFunc, 0, LFQbench.Config$MaxLogRatioForAUQC, stop.on.error=F )$value / LFQbench.Config$MaxLogRatioForAUQC, 3)  
      return( 
        list(
          species = TheSpecies,
          x = ScatterPlotXAxisData[ ValueIndices ],
          y = LogRatio[ ValueIndices ],
          density = density(LogRatio[ ValueIndices ]),
          col = LFQbench.Config$SpeciesColors[ SpeciesIndex ],
          median = LogRatioMedians[ SpeciesIndex ],
          expectation = LogRatioExpectations[ SpeciesIndex ],
          shift = LogRatioMedians[ SpeciesIndex ] - LogRatioExpectations[ SpeciesIndex ],
          qcfunction = qcFunc,
          auqc = auqc
        )
      )
    }
    
    # package data for all species
    dataBySpecies = lapply( LFQbench.Config$AllSpeciesNames, packageSpecies )
    names( dataBySpecies ) = LFQbench.Config$AllSpeciesNames
    ################################################################################
    
    ################################################################################
    # ROC-AUC for all species pairs
    globalSeparation = apply(LFQbench.Config$AllSpeciesPairs, 1, function(sp) getSepRate( dataBySpecies, LFQbench.Config$AllSpeciesNames[sp] ) )
    names(globalSeparation) = LFQbench.Config$AllSpeciesPairsLabels
    ################################################################################
    
    ################################################################################
    localSeparation = apply(
        LFQbench.Config$AllSpeciesPairs, 1, 
        function(sp) {
            getQuantileSeparation( 
                dataBySpecies, 
                LFQbench.Config$AllSpeciesNames[sp],
                LFQbench.Config$NumberOfIntensityQuantiles )
        } 
        )
    colnames(localSeparation) = LFQbench.Config$AllSpeciesPairsLabels
    ################################################################################
    
    globalAccuracy = sapply( dataBySpecies, function(d) median(d$y, na.rm=T) - d$expectation  )
    names(globalAccuracy) = LFQbench.Config$AllSpeciesNames
    globalPrecision = sapply( dataBySpecies, function(d) sd(d$y, na.rm = T) ) 
    localAccuracy = getQuantileAccuracy( dataBySpecies, LFQbench.Config$AllSpeciesNames, LFQbench.Config$NumberOfIntensityQuantiles )
    localPrecision = getQuantilePrecision( dataBySpecies, LFQbench.Config$AllSpeciesNames, LFQbench.Config$NumberOfIntensityQuantiles )
    
    ################################################################################
    # package sample pair data
    SamplePair = list(
      index = SamplePairIndex,
      index1 = Sample1Index,
      index2 = Sample2Index,
      name1 = Sample1Name,
      name2 = Sample2Name,
      xlim = xLim,
      ylim = yLim,
      col = SamplePairColor,
      data = dataBySpecies,
      separation = globalSeparation,
      accuracy = globalAccuracy,
      precision = globalPrecision,
      quantileSeparation = localSeparation,
      quantileAccuracy = localAccuracy,
      quantilePrecision = localPrecision,
      log2IntensityRanges = LFQbench.Config$Log2IntensityRanges,
      adjustment = LogRatioAdjustmentValue,
      qcrange = LFQbench.Config$AUQCRatioRange,
      allLogRatios = allLogRatios,
      validLogRatios = validLogRatios,
      invalidLogRatios = invalidLogRatios,
      medianLogRatios = LogRatioMedians,
      meanLogRatios = LogRatioMeans,
      expectedLogRatios = LogRatioExpectations
    )
    ################################################################################
    
    return(SamplePair)
  }
  ################################################################################
  
  ################################################################################
  # process all sample pairs
  SamplePairsData = lapply(1:LFQbench.Config$NumberOfSamplePairs, getSamplePairData)
  names(SamplePairsData) = sapply( SamplePairsData, function(d) paste(d$name1, ":", d$name2,sep = "" ))
  ################################################################################
  
  ################################################################################
  # generate a result set
  ResultSet = list(
      docSet = makeEverythingInOneFolderDocSet(file),
      file = file,
      idstat = NumberOfProteinsBySpecies,
      result = SamplePairsData,
      data = list (
        species=SampleAverageSpecies,
        id=SampleAverageProteinIDs,
        mean=SampleAverageAmounts,
        cv=SampleAverageCVs,
        missingvalues = numNAs_histograms
      )
  )
  ################################################################################
  
  cat("[done]\n")
  
  return( ResultSet )
}

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
# RANGE BASED METRICS CALCULATION
################################################################################

################################################################################
# removed from LFQbench.Config
Log2IntensityRanges = rbind(
    "<2"=c(0,2),
    "<4"=c(2,4),
    "<6"=c(4,6),
    "<8"=c(6,8),
    "<10"=c(8,10),
    ">10"=c(10,100)
)

#' getRangeAccuracy
#' 
#' This function compute the median deviation for log-ratios for each species in each range
#' @param dataBySpecies
#' @param ranges 
#' @param spcNames
getRangedAccuracy = function(dataBySpecies, ranges=Log2IntensityRanges, spcNames = LFQbench.Config$AllSpeciesNames){
  el2r = sapply(spcNames, function(sn) dataBySpecies[[sn]]$expectation, USE.NAMES = F )  
  l2rs = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$y ) )
  l2is = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$x ) )
  spcs = unlist( lapply( spcNames, function( sn ) rep( sn, length( dataBySpecies[[sn]]$y ) ) ) )
  acc4rs = function(r, s)
  { 
    idx4r = l2is >= r[1] & l2is < r[2]
    idx4s = spcs %in% s
    vals = l2rs[ idx4r & idx4s ]
    dev = median(vals, na.rm = T) - el2r[s]
    return(dev)
  }
  rangedAcc = sapply( spcNames, function(s) apply( ranges,  1, acc4rs, s) )
  return( rangedAcc )
}

#'getRangePrecision
#'
#'This function generate the standard deviation for log-ratios for each species in each range
#'@param dataBySpecies
#'@param ranges
#'@param spcNames
getRangedPrecision = function(dataBySpecies, ranges=LFQbench.Config$Log2IntensityRanges, spcNames = LFQbench.Config$AllSpeciesNames ){
  l2rs = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$y ) )
  l2is = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$x ) )
  spcs = unlist( lapply( spcNames, function( sn ) rep( sn, length( dataBySpecies[[sn]]$y ) ) ) )
  var4rs = function(r, s)
  { 
    idx4r = l2is >= r[1] & l2is < r[2]
    idx4s = spcs %in% s
    vals = l2rs[ idx4r & idx4s ]
    std = sd(vals, na.rm = T)
    return( std )
  }
  rangedSDs = sapply( spcNames, function(s) apply( ranges,  1, var4rs, s) )
  return( rangedSDs )
}

#' getRangedSepRate
#' 
#'  This function calculate species separation ROC-AUC for a species pair
#' @param dataBySpecies
#' @param spcNames 
#' @param ranges
getRangedSepRate = function(dataBySpecies, spcNames, ranges=LFQbench.Config$Log2IntensityRanges){
  l2rs = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$y ) )
  l2is = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$x ) )
  spcs = unlist( lapply( spcNames, function( sn ) rep( sn, length( dataBySpecies[[sn]]$y ) ) ) )
  spcFlags = as.numeric( factor(spcs) ) - 1
  auc4range = function(r) 
  { 
    idx = l2is >= r[1] & l2is < r[2]
    auc(spcFlags[idx], l2rs[idx])
  }
  rangedAUCs = apply( ranges,  1, auc4range )
  return( rangedAUCs )
}
################################################################################


################################################################################
# QUANTILE BASED METRICS CALCULATION
################################################################################

################################################################################
#' getQuantileSeparation
#' 
#'  This function calculates species separation ROC-AUC for a species pair
#'  in equally sized quantiles of intensity in first sample.
#'  
#' @param dataBySpecies data structure as stored in ResultSet$result[[SamplePairIndex]]$data
#' @param spcNames a vector of two species names
#' @param numberOfQuantiles the number of quantiles
getQuantileSeparation = function( dataBySpecies, spcNames, numberOfQuantiles=LFQbench.Config$NumberOfIntensityQuantiles ){
    # get log-ratios
    l2rs = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$y ) )
    # get intensities of first sample
    l2is = unlist( lapply( spcNames, function( sn ) dataBySpecies[[sn]]$x ) )
    # assign species names to data
    spcs = unlist( lapply( spcNames, function( sn ) rep( sn, length( dataBySpecies[[sn]]$y ) ) ) )
    # flag first species as 0 and the second as 1
    spcFlags = as.numeric( factor(spcs) ) - 1
    
    # flag quantiles of log-ratios
    l2rFlags = as.numeric( cut( rank( l2is ), numberOfQuantiles ) )
    
    auc4q = function( q )
    { 
        idx = l2rFlags == q; 
        auc(spcFlags[idx], l2rs[idx])
    }
    res = sapply( 1:numberOfQuantiles ,  auc4q )
    names(res) = paste0("Q", 1:numberOfQuantiles)
    return( res )
}
################################################################################

################################################################################
#' getQuantileAccuracy
#' 
#' This function compute the median deviation for log-ratios for each species in each range
#' @param dataBySpecies
#' @param spcNames
#' @param numberOfQuantiles
getQuantileAccuracy = function(dataBySpecies, spcNames = LFQbench.Config$AllSpeciesNames, numberOfQuantiles=LFQbench.Config$NumberOfIntensityQuantiles){
    qnts = lapply( spcNames, function(sn) as.numeric( cut( rank( dataBySpecies[[sn]]$x ), numberOfQuantiles ) ) )
    names(qnts) = spcNames
    acc4qs = function(q, s)
    { 
        l2rs = dataBySpecies[[s]]$y
        vals = l2rs[ qnts[[s]] == q  ]
        median(vals, na.rm = T) - dataBySpecies[[s]]$expectation
    }
    res = sapply( spcNames, function(s) sapply( 1:numberOfQuantiles, acc4qs, s) )
    rownames(res) = paste0("Q", 1:numberOfQuantiles)
    return( res )
}
################################################################################

################################################################################
#'getQuantilePrecision
#'
#'This function generate the standard deviation for log-ratios for each species in each range
#'@param dataBySpecies
#'@param spcNames
#' @param numberOfQuantiles
getQuantilePrecision = function(dataBySpecies, spcNames = LFQbench.Config$AllSpeciesNames, numberOfQuantiles=LFQbench.Config$NumberOfIntensityQuantiles )
{
    qnts = lapply( spcNames, function(sn) as.numeric( cut( rank( dataBySpecies[[sn]]$x ), numberOfQuantiles ) ) )
    names(qnts) = spcNames
    var4qs = function(q, s)
    { 
        l2rs = dataBySpecies[[s]]$y
        vals = l2rs[ qnts[[s]] == q  ]
        sd(vals, na.rm = T)
    }
    res = sapply( spcNames, function(s) sapply( 1:numberOfQuantiles,  var4qs, s) )
    rownames(res) = paste0("Q", 1:numberOfQuantiles)
    return( res )
}

################################################################################

################################################################################
# OVERALL METRICS CALCULATION
################################################################################

#' getSepRate
#' 
#' This function calculate species separation ROC-AUC for a species pair
#' @param dataBySpecies
#' @param spcNames
#' @param numberOfQuantiles
getSepRate = function(dataBySpecies, spcNames){
  l2rs = unlist( lapply( spcNames, function( x ) dataBySpecies[[x]]$y ) )
  spcs = unlist( lapply( spcNames, function( x ) rep( x, length( dataBySpecies[[x]]$y ) ) ) )
  spcFlags = as.numeric( factor(spcs) ) - 1
  res = auc( spcFlags, l2rs )
  return( res )
}

################################################################################

