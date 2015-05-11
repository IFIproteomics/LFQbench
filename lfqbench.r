rm( list=ls() )

DEBUG=T

source('inc/lfqbench.com.r')
source('lfqbench.config.r')
source('inc/lfqbench.defs.r')
source('inc/lfqbench.plot.r')
source('inc/lfqbench.proc.r')
source('inc/lfqbench.access.r')

if(DEBUG) cat( "R working directory: " + getwd() + "\n")

################################################################################
AllInputFiles = list.files( path=InputFilesLocation, pattern=InputExtensionPattern, full.names = FALSE )

if( length(AllInputFiles) < 1 ) 
{
  cat("working directory: " + getwd() + "\n")
  cat("input directory:   " + InputFilesLocation + "\n")
  cat("file extension:    " + InputExtensionPattern + "\n")
  stop( "no input files found!" )
}

DocSets = lapply( AllInputFiles, makeDocSet )
names(DocSets)=sapply(DocSets, function(d)d$fileBase)
ResultSets = lapply( DocSets, processData )

save(DocSets, ResultSets, file = LogFilesLocation + "/ResultSets.rda")
# nix=sapply( ResultSets, saveSampleMeans )
# nix=sapply( ResultSets, saveSpeciesSeparation )
# nix=sapply( ResultSets, saveIDs )
nix=sapply( ResultSets, saveLogRatios )
nix=sapply( ResultSets, plotResultSet )
nix=sapply( ResultSets, saveSampleRSD )
 
plotSpeciesLegends()
plotSampleComposition()
saveMetrics(ResultSets)
logIdStatistics( ResultSets )
################################################################################

# END

################################################################################
# FOR DEBUGGING ONLY
doc = DocSets[[1]]
resultSet = ResultSets[[1]]
################################################################################

# end

