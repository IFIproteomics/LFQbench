#' LFQbench.batchProcessRootFolder
#' 
#' Wrapper function for processing a batch of input files 
#' automatically discovered from a predefined root folder.
#' Expected structure of root folder:
#' 
#'  - subfolder "input" contains automatically discovered input files
#'  - subfolder "log" contains calculated metrics, stored result sets and other exported information.
#'  - subfolder "plot" contains generated plots.
#' 
#' In case the given root folder is empty, all needed file structures will be created automatically.
#' You can use the freshly created structures to put your input files and rerun LFQbench.
#' 
#' @param rootFolder the root folder for batch processing
#' @export
LFQbench.batchProcessRootFolder = function( rootFolder=LFQbench.Config$DataRootFolder, ... )
{
  if(rootFolder != LFQbench.Config$DataRootFolder) LFQbench.setDataRootFolder( rootFolder )
  inputFiles = list.files( path=LFQbench.Config$InputFilesLocation, pattern=LFQbench.Config$InputExtensionPattern, full.names = T )
  return( LFQbench.batchProcessFiles( inputFiles, ... ) )
}

#' LFQbench.batchProcessFiles
#' 
#' Wrapper function for processing a batch of input files.
#' Plots, metrics and log files will be automatically put into 
#' subfolders of a root folder.
#' PLEASE DEFINE THE ROOT FOLDER BEFORE RUNNING THIS FUNCTION!!!
#' Expected structure of root folder:
#' 
#'  - subfolder "input" contains automatically discovered input files
#'  - subfolder "log" contains calculated metrics, stored result sets and other exported information.
#'  - subfolder "plot" contains generated plots.
#' 
#' In case the given root folder is empty, all needed file structures will be created automatically.
#' You can use the freshly created structures to put your input files and rerun LFQbench.
#' 
#' @param inputFiles a list of files to process
#' @param storeResultSetsImage=T,
#' @param storeSampleMeans=F,
#' @param storeSpeciesSeparation=F,
#' @param storeIDs=F,
#' @param storeIdStats=T,
#' @param storeLogRatios=F,
#' @param storePlots=T,
#' @param storeSampleCVs=F,
#' @param storeMetrics=T
#' @export
LFQbench.batchProcessFiles = function( inputFiles, 
                         storeResultSetsImage=T,
                         storeSampleMeans=F,
                         storeSpeciesSeparation=F,
                         storeIDs=F,
                         storeIdStats=T,
                         storeLogRatios=F,
                         storePlots=T,
                         storeSampleCVs=F,
                         storeMetrics=T
)
{
    ResultSets = NULL
    if( length(inputFiles) < 1 ) 
    {
        cat( paste("working directory: ", getwd(),"\n", sep="") )
        cat( paste("input directory:   ", LFQbench.Config$InputFilesLocation, "\n", sep = "") )
        cat( paste("file extension:    ", LFQbench.Config$InputExtensionPattern, "\n", sep="") )
        stop( "no input files found!" )
    }
    else
    {
        LFQbench.setDataRootFolder( LFQbench.Config$DataRootFolder )
        
        DocSets = lapply( inputFiles, makeDocSet )
        names(DocSets) = sapply( DocSets, function(d) d$fileBase )
        ResultSets = lapply( DocSets, LFQbench.processDocSet )

        if(storeResultSetsImage) save(DocSets, ResultSets, file = paste(LFQbench.Config$LogFilesLocation, "/ResultSets.rda", sep = "") )
        if(storeSampleMeans) nix=sapply( ResultSets, saveSampleMeans )
        if(storeSpeciesSeparation) nix=sapply( ResultSets, saveSpeciesSeparation )
        if(storeIDs) nix=sapply( ResultSets, saveIDs )
        if(storeLogRatios) nix=sapply( ResultSets, saveLogRatios )
        if(storePlots) {
            nix=sapply( ResultSets, LFQbench.plotResultSet )
            LFQbench.plotSpeciesLegends()
            LFQbench.plotSampleComposition()
        }
        if(storeSampleCVs) nix=sapply( ResultSets, saveSampleCVs )
        if(storeMetrics) saveMetrics( ResultSets )
        if(storeIdStats) logIdStatistics( ResultSets )
    }
    
    return( ResultSets )
}
