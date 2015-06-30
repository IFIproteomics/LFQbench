rm( list=ls() )

DEBUG=T

source('lfqbench.com.r')
# load static configuration data
source('lfqbench.config.r')
# override parameters by command line (if any)
evalCommandLineArguments()
# calculate dynamic parameters from config
source('lfqbench.defs.r')

source('lfqbench.plot.r')
source('lfqbench.proc.r')
source('lfqbench.access.r')

if(DEBUG) cat( "R working directory: " + getwd() + "\n")

################################################################################
AllInputFiles = list.files( path=cfg$InputFilesLocation, pattern=cfg$InputExtensionPattern, full.names = FALSE )

if( length(AllInputFiles) < 1 ) 
{
  cat("working directory: " + getwd() + "\n")
  cat("input directory:   " + cfg$InputFilesLocation + "\n")
  cat("file extension:    " + cfg$InputExtensionPattern + "\n")
  stop( "no input files found!" )
}

DocSets = lapply( AllInputFiles, makeDocSet )
names(DocSets)=sapply(DocSets, function(d)d$fileBase)
ResultSets = lapply( DocSets, processData )

save(DocSets, ResultSets, file = cfg$LogFilesLocation + "/ResultSets.rda")
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

