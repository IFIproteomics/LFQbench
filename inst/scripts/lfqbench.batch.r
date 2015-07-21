library(LFQbench)

rm( list=ls() )

DEBUG=T

# load static configuration data
source('lfqbench.config.r')
# override parameters by command line (if any)
evalCommandLineArguments()
# calculate dynamic parameters from config
source('lfqbench.defs.r')
source('lfqbench.plot.r')


if(DEBUG) cat(paste("R working directory: ",getwd(), "\n", sep = ""))

################################################################################
AllInputFiles = list.files( path=cfg$InputFilesLocation, pattern=cfg$InputExtensionPattern, full.names = FALSE )

if( length(AllInputFiles) < 1 ) 
{
  cat(paste("working directory: ", getwd(),"\n", sep=""))
  cat(paste("input directory:   ", cfg$InputFilesLocation, "\n", sep = ""))
  cat(paste("file extension:    ", cfg$InputExtensionPattern, "\n", sep=""))
  stop( "no input files found!" )
}

DocSets = lapply( AllInputFiles, makeDocSet )
names(DocSets)=sapply(DocSets, function(d)d$fileBase)
ResultSets = lapply( DocSets, processData )

save(DocSets, ResultSets, file = paste(cfg$LogFilesLocation, "/ResultSets.rda"))
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

