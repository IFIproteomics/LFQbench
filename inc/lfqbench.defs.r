
loadLibrary("RColorBrewer")

################################################################################
# create input/output paths
sapply( c(InputFilesLocation, PlotFilesLocation, LogFilesLocation), mkdir )
################################################################################

################################################################################
# process composition of samples
AllSampleNames = as.vector( colnames(SampleComposition)[-1] )
AllSpeciesNames = as.vector( SampleComposition[,1] )
AllExpectedAmounts = as.matrix( SampleComposition[,-1] )
rownames(AllExpectedAmounts) = AllSpeciesNames
NumberOfSamples = length(AllSampleNames)
NumberOfSpecies = length(AllSpeciesNames)
SampleColors = brewer.pal(max(NumberOfSamples,3),"Dark2")[1:NumberOfSamples]
SpeciesColors = brewer.pal(max(NumberOfSpecies,3),"Dark2")[1:NumberOfSpecies]
################################################################################

################################################################################
AUQCRatioRange = c(0, MaxLogRatioForAUQC)
################################################################################

################################################################################
# create sample index pairs
SamplePairsIndices   = createNumericPairs( 1, NumberOfSamples )
NumberOfSamplePairs = nrow( SamplePairsIndices )
SamplePairsLabels = apply(SamplePairsIndices, 1, function(sp){nms = AllSampleNames[sp]; return(nms[1]+":"+nms[2])} )
SamplePairsColors  	= brewer.pal( max(NumberOfSamplePairs,3), "Dark2" )[1:NumberOfSamplePairs]
################################################################################

################################################################################
AllSpeciesPairs = createNumericPairs(1, NumberOfSpecies)
AllSpeciesPairsLabels = apply(AllSpeciesPairs, 1, function(sp){nms = AllSpeciesNames[sp]; return(nms[1]+"-"+nms[2])} )
################################################################################

################################################################################
# key-accession-entry-species translation
cat( "reading protein to species database ... " )
IdToSpeciesFile = "inc/hye.id2species.csv"
IdToSpeciesTable = NULL
getSpeciesForIds = function( IDs )
{
  if(is.null(IdToSpeciesTable))
  {
    IdToSpeciesTable <<- read.table(file=IdToSpeciesFile, header=T, sep=";")
  }
  as.vector( IdToSpeciesTable$species[ match(IDs, IdToSpeciesTable$id, nomatch=NA, incomparables=NA) ] )
}
cat( "[done]\n" )
################################################################################
