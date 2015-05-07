rm( list=ls() )

source('inc/lfqbench.com.r')
source('lfqbench.config.r')
source('inc/lfqbench.defs.r')
source('inc/lfqbench.plot.r')
source('inc/lfqbench.proc.r')
source('inc/lfqbench.access.r')

DEBUG=T

if(DEBUG) cat( "R working directory: " + getwd() + "\n")

################################################################################
makeDocSet = function( inFile )
{
  fileBase = sub(InputExtensionPattern, "", inFile)
  return( 
    list(
      inputPath = InputFilesLocation,
      plotPath = PlotFilesLocation,
      logPath = LogFilesLocation,
      fileBase = fileBase,
      csvFile = InputFilesLocation + "/" + inFile,
      pdfFile = PlotFilesLocation + "/" + fileBase + ".pdf",
      logFile = LogFilesLocation + "/" + fileBase + ".log",
      avgFile = LogFilesLocation + "/" + fileBase + " sample_means.csv",
      rsdFile = LogFilesLocation + "/" + fileBase + " cv.csv",
      l2rFile = LogFilesLocation + "/" + fileBase + " log2_ratio.csv",
      rocFile = LogFilesLocation + "/" + fileBase + " species_separation.csv",
      idsFile = LogFilesLocation + "/" + fileBase + " ids.csv"
    )
  )
}
################################################################################


################################################################################
AllInputFiles = list.files( path=InputFilesLocation, pattern=InputExtensionPattern, full.names = FALSE )
DocSets = lapply( AllInputFiles, makeDocSet )
names(DocSets)=sapply(DocSets, function(d)d$fileBase)
ResultSets = lapply( DocSets, processData )

save(DocSets, ResultSets, file = LogFilesLocation + "/ResultSets.R")
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

