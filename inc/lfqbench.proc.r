loadLibrary("SDMTools")

processData = function( DocSet )
{
  #DocSet = DocSets[[1]]    
  cat( "processing " + DocSet$fileBase + " ... \n" )
  
  ################################################################################
  # read run based data
  csvDat = read.table(DocSet$csvFile, sep=CsvColumnSeparator, dec=CsvDecimalPointChar, header=T)
  csvIDs = as.vector( csvDat[,1] )
  # csvAmounts = as.matrix( csvDat[,-1] )  ## Old version, it does not support additional factor columns 
  csvAmounts = as.matrix( csvDat[, sapply(csvDat, class) == "numeric" | sapply(csvDat, class) == "integer"] )
  if(any(names(csvDat) == "specie")){
      csvSpecies = as.vector( csvDat[, "specie"] ) 
  }else{
      csvSpecies = as.vector( getSpeciesForIds( csvIDs ) )
  }
  rm( csvDat )
  ################################################################################
  # check OTHER species
  benchedSpeciesIndices = unlist(sapply( AllSpeciesNames, function(spc) which(csvSpecies==spc) ), use.names = F)
  otherSpeciesIndices = (1:length(csvIDs))[ -benchedSpeciesIndices ]
  otherSpeciesIDs = csvIDs[ otherSpeciesIndices ]
  NumberOfProteinsForOtherSpecies = length(otherSpeciesIndices)
  # limit data to benched species only
  csvIDs = csvIDs[ benchedSpeciesIndices ]
  csvAmounts = csvAmounts[ benchedSpeciesIndices, ]
  csvSpecies = csvSpecies[ benchedSpeciesIndices ]
  ################################################################################
  # normalize protein amounts to ppm
  if(NormalizeAmountsToPPM) csvAmounts = apply( csvAmounts, 2, function(x) x / sum(x, na.rm=T) * 1000000 )
  ################################################################################
  
  ################################################################################
  # extract data for a single sample
  getSampleAverage = function(sampleIndex, nRuns=floor( ncol(csvAmounts) / NumberOfSamples ))
  {
    amount = as.matrix( csvAmounts[, ((sampleIndex-1)*nRuns + 1):(nRuns * sampleIndex)] )
    amount[amount < MinProteinAmount] = NA
    avg = rowMeans( amount, na.rm=T )
    std = apply( amount, 1, sd, na.rm=T )
    rsd = std / avg
    return( list( mean=avg, cv=rsd ) )
  }
  ################################################################################
  # generate sample average data
  SampleAverageData = lapply( 1:NumberOfSamples, getSampleAverage )
  SampleAverageAmounts = sapply( 1:NumberOfSamples, function(si) SampleAverageData[[si]]$mean )
  SampleAverageCVs = sapply( 1:NumberOfSamples, function(si) SampleAverageData[[si]]$cv )
  SampleAverageSpecies = csvSpecies
  SampleAverageProteinIDs = csvIDs
  rownames(SampleAverageAmounts) = SampleAverageProteinIDs
  rownames(SampleAverageCVs) = SampleAverageProteinIDs
  colnames(SampleAverageAmounts) = AllSampleNames
  colnames(SampleAverageCVs) = AllSampleNames
  ################################################################################
  
  ################################################################################
  # create protein identification statistics
  NumberOfProteinsBySpecies = sapply( AllSpeciesNames, function(x) length( which(SampleAverageSpecies == x) ) )
  NumberOfProteinsBySpecies["OTHER"] = NumberOfProteinsForOtherSpecies
  ################################################################################
  
  ################################################################################
  # handle missing and low amount values
  # SampleAverageAmounts[ is.na(SampleAverageAmounts) ] = MinProteinAmount
  # SampleAverageAmounts[ SampleAverageAmounts < MinProteinAmount ] = MinProteinAmount
  ################################################################################
  
  ################################################################################
  # process a sample pair
  getSamplePairData = function( SamplePairIndex = 1 )
  {
    ################################################################################
    # extract sample pair data
    SamplePairColor = SamplePairsColors[SamplePairIndex]
    Sample1Index = SamplePairsIndices[SamplePairIndex,1]
    Sample2Index = SamplePairsIndices[SamplePairIndex,2]
    Sample1Name = AllSampleNames[Sample1Index]
    Sample2Name = AllSampleNames[Sample2Index]
    Sample1ProteinAmounts = SampleAverageAmounts[,Sample1Index]
    Sample2ProteinAmounts = SampleAverageAmounts[,Sample2Index]
    SpeciesNames = SampleAverageSpecies
    ProteinIDs = SampleAverageProteinIDs
    LogRatioExpectations = log2( AllExpectedAmounts[,Sample1Index] / AllExpectedAmounts[,Sample2Index] )
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
    if(DropInvalidLogRatio)
    {
      ValidLogRatioIndices = LogRatio >= LogRatioValidityRange[1] & LogRatio <= LogRatioValidityRange[2]
      invalidLogRatios = data.frame(
	      	entry = ProteinIDs[!ValidLogRatioIndices],
	      	species = SpeciesNames[!ValidLogRatioIndices],
	      	logratio = LogRatio[!ValidLogRatioIndices],
	      	sample1means = Sample1ProteinAmounts[!ValidLogRatioIndices],
	      	sample2means = Sample2ProteinAmounts[!ValidLogRatioIndices]
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
    LogRatioMedians = sapply(AllSpeciesNames, function(x) median( LogRatio[SpeciesNames==x], na.rm=T ))
    LogRatioMeans = sapply(AllSpeciesNames, function(x) mean( LogRatio[SpeciesNames==x], na.rm=T ))
    BackgroundSpeciesMedian = LogRatioMedians[AllSpeciesNames==BackgroundSpeciesName]
    BackgroundSpeciesMean = LogRatioMeans[AllSpeciesNames==BackgroundSpeciesName]
    ################################################################################
    
    ################################################################################
    # center by mean, median, or not at all
    if( CenterLogRatioByBackground || regexpr("median", CenterLogRatioByBackground, ignore.case=T) )
      LogRatioAdjustmentValue = BackgroundSpeciesMedian 
    else if(regexpr("mean", CenterLogRatioByBackground, ignore.case=T))
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
      sample2means = Sample2ProteinAmounts
    )
    
    ################################################################################
    # define and clip scatterplot x-axis data
    ScatterPlotXAxisData = log2( Sample2ProteinAmounts )
    # xLim = range(ScatterPlotXAxisData[ScatterPlotXAxisData>0])
    xLim = quantile( ScatterPlotXAxisData[ScatterPlotXAxisData>0], probs=c(0.01,0.99), na.rm = T )
    yLim = LogRatioPlotRange
    # ensure x-axis boundaries
    ScatterPlotXAxisData[ ScatterPlotXAxisData < xLim[1] ] = xLim[1]
    ScatterPlotXAxisData[ ScatterPlotXAxisData > xLim[2] ] = xLim[2]
    ################################################################################
    
    ################################################################################
    # package data for a species
    packageSpecies = function( TheSpecies )
    {
      ValueIndices = which(SpeciesNames == TheSpecies)
      SpeciesIndex = AllSpeciesNames == TheSpecies
      if(DEBUG) cat( TheSpecies + " has " + length( ValueIndices ) + " IDs ... \n")
      qcFunc = as.function(
        getQCFunction( LogRatio[ValueIndices] - LogRatioMedians[SpeciesIndex], ensureValueRange=c(0, MaxLogRatioForAUQC) )
      )
      auqc = round( integrate( qcFunc, 0, MaxLogRatioForAUQC, stop.on.error=F )$value / MaxLogRatioForAUQC, 3)  
      return( 
        list(
          species = TheSpecies,
          x = ScatterPlotXAxisData[ ValueIndices ],
          y = LogRatio[ ValueIndices ],
          density = density(LogRatio[ ValueIndices ]),
          col = SpeciesColors[ SpeciesIndex ],
          median = LogRatioMedians[ SpeciesIndex ],
          expectation = LogRatioExpectations[ SpeciesIndex ],
          shift = LogRatioMedians[ SpeciesIndex ] - LogRatioExpectations[ SpeciesIndex ],
          qcfunction = qcFunc,
          auqc = auqc
        )
      )
    }
    
    # package data for all species
    dataBySpecies = lapply( AllSpeciesNames, packageSpecies )
    names( dataBySpecies ) = AllSpeciesNames
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
    spcPairsSepRates = apply(AllSpeciesPairs, 1, function(sp) getSepRate( dataBySpecies, AllSpeciesNames[sp] ) )
    names(spcPairsSepRates) = AllSpeciesPairsLabels
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
      adjustment = LogRatioAdjustmentValue,
      qcrange = AUQCRatioRange,
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
  SamplePairsData = lapply(1:NumberOfSamplePairs, getSamplePairData)
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

