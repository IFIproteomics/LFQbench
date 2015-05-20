# configuration settings calculated from other parameters

loadLibrary("RColorBrewer")

################################################################################
# create input/output paths
sapply( c(cfg$InputFilesLocation, cfg$PlotFilesLocation, cfg$LogFilesLocation), mkdir )
################################################################################

################################################################################
# process composition of samples
cfg$AllSampleNames = as.vector( colnames(cfg$SampleComposition)[-1] )
cfg$AllSpeciesNames = as.vector( cfg$SampleComposition[,1] )
cfg$AllExpectedAmounts = as.matrix( cfg$SampleComposition[,-1] )
rownames(cfg$AllExpectedAmounts) = cfg$AllSpeciesNames
cfg$NumberOfSamples = length(cfg$AllSampleNames)
cfg$NumberOfSpecies = length(cfg$AllSpeciesNames)
cfg$SampleColors = brewer.pal(max(cfg$NumberOfSamples,3),"Dark2")[1:cfg$NumberOfSamples]
cfg$SpeciesColors = brewer.pal(max(cfg$NumberOfSpecies,3),"Dark2")[1:cfg$NumberOfSpecies]
################################################################################

################################################################################
cfg$AUQCRatioRange = c(0, cfg$MaxLogRatioForAUQC)
################################################################################

################################################################################
# create sample index pairs
cfg$SamplePairsIndices   = createNumericPairs( 1, cfg$NumberOfSamples )
cfg$NumberOfSamplePairs = nrow( cfg$SamplePairsIndices )
cfg$SamplePairsLabels = apply(cfg$SamplePairsIndices, 1, function(sp){nms = cfg$AllSampleNames[sp]; return(nms[1]+":"+nms[2])} )
cfg$SamplePairsColors  	= brewer.pal( max(cfg$NumberOfSamplePairs,3), "Dark2" )[1:cfg$NumberOfSamplePairs]
################################################################################

################################################################################
cfg$AllSpeciesPairs = createNumericPairs(1, cfg$NumberOfSpecies)
cfg$AllSpeciesPairsLabels = apply(cfg$AllSpeciesPairs, 1, function(sp){nms = cfg$AllSpeciesNames[sp]; return(nms[1]+"-"+nms[2])} )
################################################################################
