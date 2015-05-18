loadLibrary("SDMTools")

processData = function( DocSet )
{
  #DocSet = DocSets[[1]]    
  cat( "processing " + DocSet$fileBase + " ... \n" )
  
  ################################################################################
  # read run based data
  csvDat = read.table(DocSet$csvFile, sep=cfg$CsvColumnSeparator, dec=cfg$CsvDecimalPointChar, header=T)
  
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
  if( length( which( numericColumnIdx ) ) < cfg$NumberOfSamples )
    stop("number of amount columns is smaller than the number of samples")
  
  if( length( speciesColumnIdx ) < 1 ) 
    stop("no ''species'' column found.")
  
  # 1.0 * is a work-around to force the amounts matrix to numeric class
  csvAmounts = 1.0 * as.matrix( csvDat[, numericColumnIdx] )
  csvSpecies = as.vector( csvDat[, speciesColumnIdx[[1]]] )  
  rm( csvDat )
  
  # disable low amounts
  csvAmounts[csvAmounts < cfg$MinProteinAmount] = NA
  
  ################################################################################
  # check OTHER species
  benchedSpeciesIndices = unlist(sapply( cfg$AllSpeciesNames, function(spc) which(csvSpecies==spc) ), use.names = F)
  otherSpeciesIndices = (1:length(csvIDs))[ -benchedSpeciesIndices ]
  otherSpeciesIDs = csvIDs[ otherSpeciesIndices ]
  NumberOfProteinsForOtherSpecies = length(otherSpeciesIndices)
  # limit data to benched species only
  csvIDs = csvIDs[ benchedSpeciesIndices ]
  csvAmounts = csvAmounts[ benchedSpeciesIndices, ]
  csvSpecies = csvSpecies[ benchedSpeciesIndices ]
  ################################################################################
  # normalize protein amounts to ppm
  if(cfg$NormalizeAmountsToPPM) csvAmounts = apply( csvAmounts, 2, function(x) x / sum(x, na.rm=T) * 1000000 )
  ################################################################################
  
  ################################################################################
  # extract data for a single sample
  getSampleAverage = function(sampleIndex, nRuns=floor( ncol(csvAmounts) / cfg$NumberOfSamples ))
  {
    amount = as.matrix( csvAmounts[, ((sampleIndex-1)*nRuns + 1):(nRuns * sampleIndex)] )
    avg = rowMeans( amount, na.rm=T )
    std = apply( amount, 1, sd, na.rm=T )
    rsd = std / avg
    return( list( mean=avg, cv=rsd ) )
  }
  ################################################################################
  # generate sample average data
  SampleAverageData = lapply( 1:cfg$NumberOfSamples, getSampleAverage )
  SampleAverageAmounts = sapply( 1:cfg$NumberOfSamples, function(si) SampleAverageData[[si]]$mean )
  SampleAverageCVs = sapply( 1:cfg$NumberOfSamples, function(si) SampleAverageData[[si]]$cv )
  SampleAverageSpecies = csvSpecies
  SampleAverageProteinIDs = csvIDs
  rownames(SampleAverageAmounts) = SampleAverageProteinIDs
  rownames(SampleAverageCVs) = SampleAverageProteinIDs
  colnames(SampleAverageAmounts) = cfg$AllSampleNames
  colnames(SampleAverageCVs) = cfg$AllSampleNames
  ################################################################################
  
  ################################################################################
  # create protein identification statistics
  NumberOfProteinsBySpecies = sapply( cfg$AllSpeciesNames, function(x) length( which(SampleAverageSpecies == x) ) )
  NumberOfProteinsBySpecies["OTHER"] = NumberOfProteinsForOtherSpecies
  ################################################################################
  
  ################################################################################
  # handle missing and low amount values
  # SampleAverageAmounts[ is.na(SampleAverageAmounts) ] = cfg$MinProteinAmount
  # SampleAverageAmounts[ SampleAverageAmounts < cfg$MinProteinAmount ] = cfg$MinProteinAmount
  ################################################################################
  
  ################################################################################
  # process a sample pair
  getSamplePairData = function( SamplePairIndex = 1 )
  {
    ################################################################################
    # extract sample pair data
    SamplePairColor = cfg$SamplePairsColors[SamplePairIndex]
    Sample1Index = cfg$SamplePairsIndices[SamplePairIndex,1]
    Sample2Index = cfg$SamplePairsIndices[SamplePairIndex,2]
    Sample1Name = cfg$AllSampleNames[Sample1Index]
    Sample2Name = cfg$AllSampleNames[Sample2Index]
    Sample1ProteinAmounts = SampleAverageAmounts[,Sample1Index]
    Sample2ProteinAmounts = SampleAverageAmounts[,Sample2Index]
    SpeciesNames = SampleAverageSpecies
    ProteinIDs = SampleAverageProteinIDs
    LogRatioExpectations = log2( cfg$AllExpectedAmounts[,Sample1Index] / cfg$AllExpectedAmounts[,Sample2Index] )
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
    if(cfg$DropInvalidLogRatio)
    {
      ValidLogRatioIndices = LogRatio >= cfg$LogRatioValidityRange[1] & LogRatio <= cfg$LogRatioValidityRange[2]
      invalidLogRatios = data.frame(
	      	entry = ProteinIDs[!ValidLogRatioIndices],
	      	species = SpeciesNames[!ValidLogRatioIndices],
	      	logratio = LogRatio[!ValidLogRatioIndices],
	      	sample1means = Sample1ProteinAmounts[!ValidLogRatioIndices],
	      	sample2means = Sample2ProteinAmounts[!ValidLogRatioIndices],
	      	row.names = NULL,
          stringsAsFactors=F
	      )
      LogRatio = LogRatio[ValidLogRatioIndices]
      SpeciesNames = SpeciesNames[ValidLogRatioIndices]
      ProteinIDs = ProteinIDs[ValidLogRatioIndices]
      Sample1ProteinAmounts = Sample1ProteinAmounts[ValidLogRatioIndices]
      Sample2ProteinAmounts = Sample2ProteinAmounts[ValidLogRatioIndices]
    }
    ################################################################################
    
    ################################################################################
    # calculate log-ratio means and medians
    LogRatioMedians = sapply(cfg$AllSpeciesNames, function(x) median( LogRatio[SpeciesNames==x], na.rm=T ))
    LogRatioMeans = sapply(cfg$AllSpeciesNames, function(x) mean( LogRatio[SpeciesNames==x], na.rm=T ))
    BackgroundSpeciesMedian = LogRatioMedians[cfg$AllSpeciesNames==cfg$BackgroundSpeciesName]
    BackgroundSpeciesMean = LogRatioMeans[cfg$AllSpeciesNames==cfg$BackgroundSpeciesName]
    ################################################################################
    
    ################################################################################
    # center by mean, median, or not at all
    if( cfg$CenterLogRatioByBackground || regexpr("median", cfg$CenterLogRatioByBackground, ignore.case=T) )
      LogRatioAdjustmentValue = BackgroundSpeciesMedian 
    else if(regexpr("mean", cfg$CenterLogRatioByBackground, ignore.case=T))
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
    yLim = cfg$LogRatioPlotRange
    # ensure x-axis boundaries
    ScatterPlotXAxisData[ ScatterPlotXAxisData < xLim[1] ] = xLim[1]
    ScatterPlotXAxisData[ ScatterPlotXAxisData > xLim[2] ] = xLim[2]
    ################################################################################
    
    ################################################################################
    # package data for a species
    packageSpecies = function( TheSpecies )
    {
      ValueIndices = which(SpeciesNames == TheSpecies)
      SpeciesIndex = cfg$AllSpeciesNames == TheSpecies
      if(DEBUG) cat( TheSpecies + " has " + length( ValueIndices ) + " IDs ... \n")
      qcFunc = as.function(
        getQCFunction( LogRatio[ValueIndices] - LogRatioMedians[SpeciesIndex], ensureValueRange=c(0, cfg$MaxLogRatioForAUQC) )
      )
      auqc = round( integrate( qcFunc, 0, cfg$MaxLogRatioForAUQC, stop.on.error=F )$value / cfg$MaxLogRatioForAUQC, 3)  
      return( 
        list(
          species = TheSpecies,
          x = ScatterPlotXAxisData[ ValueIndices ],
          y = LogRatio[ ValueIndices ],
          density = density(LogRatio[ ValueIndices ]),
          col = cfg$SpeciesColors[ SpeciesIndex ],
          median = LogRatioMedians[ SpeciesIndex ],
          expectation = LogRatioExpectations[ SpeciesIndex ],
          shift = LogRatioMedians[ SpeciesIndex ] - LogRatioExpectations[ SpeciesIndex ],
          qcfunction = qcFunc,
          auqc = auqc
        )
      )
    }
    
    # package data for all species
    dataBySpecies = lapply( cfg$AllSpeciesNames, packageSpecies )
    names( dataBySpecies ) = cfg$AllSpeciesNames
    ################################################################################
    
    ################################################################################
    # calculate species separation ROC-AUC for a species pair
    getSepRate = function(dataBySpecies, spcNames)
    {
      l2rs = unlist( lapply( spcNames, function( x ) dataBySpecies[[x]]$y ) )
      spcs = unlist( lapply( spcNames, function( x ) rep( x, length( dataBySpecies[[x]]$y ) ) ) )
      spcFlags = as.numeric( factor(spcs) ) - 1
      AreaUnderROCCurveByValue = auc(spcFlags, l2rs)
      return( AreaUnderROCCurveByValue )
    }
    # ROC-AUC for all species pairs
    spcPairsSepRates = apply(cfg$AllSpeciesPairs, 1, function(sp) getSepRate( dataBySpecies, cfg$AllSpeciesNames[sp] ) )
    names(spcPairsSepRates) = cfg$AllSpeciesPairsLabels
    ################################################################################
    
    ################################################################################
    # calculate species separation ROC-AUC for a species pair
    getRangedSepRate = function(dataBySpecies, spcNames, ranges=cfg$Log2IntensityRangesForSpeciesSeparation)
    {
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
    
    spcPairsRangedSeparationRates = apply(cfg$AllSpeciesPairs, 1, function(sp) getRangedSepRate( dataBySpecies, cfg$AllSpeciesNames[sp] ) )
    colnames(spcPairsRangedSeparationRates) = cfg$AllSpeciesPairsLabels
    ################################################################################
    
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
      separation = spcPairsSepRates,
      separationRanged = spcPairsRangedSeparationRates,
      adjustment = LogRatioAdjustmentValue,
      qcrange = cfg$AUQCRatioRange,
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
  SamplePairsData = lapply(1:cfg$NumberOfSamplePairs, getSamplePairData)
  names(SamplePairsData) = sapply( SamplePairsData, function(d) d$name1+":"+d$name2 )
  ################################################################################
  
  ################################################################################
  # generate a result set
  ResultSet = list(
      docSet = DocSet,
      idstat = NumberOfProteinsBySpecies,
      result = SamplePairsData,
      data = list (
        species=SampleAverageSpecies,
        id=SampleAverageProteinIDs,
        mean=SampleAverageAmounts,
        cv=SampleAverageCVs
      )
  )
  ################################################################################
  
  cat("[done]\n")
  
  return( ResultSet )
}

