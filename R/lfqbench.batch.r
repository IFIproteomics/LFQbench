#' processRootFolder
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
processRootFolder = function( rootFolder=cfg$DataRootFolder, ... )
{
  if(rootFolder != cfg$DataRootFolder) setBatchDirectory( rootFolder )
  inputFiles = list.files( path=cfg$InputFilesLocation, pattern=cfg$InputExtensionPattern, full.names = T )
  return( processBatch( inputFiles, ... ) )
}

#' processBatch
#' 
#' Wrapper function for processing a batch of input files.
#' Plots, metrics and log files will be automatically put into 
#' subfolders of a root folder.
#' PLEASE DEFINE THE ROOT FOLDER BEVORE RUNNING THIS FUNCTION!!!
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
#' @param storeSampleMeans=T,
#' @param storeSpeciesSeparation=T,
#' @param storeIDs=T,
#' @param storeIdStats=T,
#' @param storeLogRatios=T,
#' @param storePlots=T,
#' @param storeSampleCVs=T,
#' @param storeMetrics=T,
#' @param returnResultSets=F
#' @export
processBatch = function( inputFiles, 
                         storeResultSetsImage=T,
                         storeSampleMeans=T,
                         storeSpeciesSeparation=T,
                         storeIDs=T,
                         storeIdStats=T,
                         storeLogRatios=T,
                         storePlots=T,
                         storeSampleCVs=T,
                         storeMetrics=T,
                         returnResultSets=F
                         )
{
    if( length(inputFiles) < 1 ) 
    {
        cat( paste("working directory: ", getwd(),"\n", sep="") )
        cat( paste("input directory:   ", cfg$InputFilesLocation, "\n", sep = "") )
        cat( paste("file extension:    ", cfg$InputExtensionPattern, "\n", sep="") )
        stop( "no input files found!" )
    }
    else
    {
        setRootFolder( cfg$DataRootFolder )
        
        DocSets = lapply( inputFiles, makeDocSet )
        names(DocSets) = sapply( DocSets, function(d) d$fileBase )
        ResultSets = lapply( DocSets, processDocSet )
        
        if(storeResultSetsImage) save(DocSets, ResultSets, file = paste(cfg$LogFilesLocation, "/ResultSets.rda", sep = "") )
        if(storeSampleMeans) nix=sapply( ResultSets, saveSampleMeans )
        if(storeSpeciesSeparation) nix=sapply( ResultSets, saveSpeciesSeparation )
        if(storeIDs) nix=sapply( ResultSets, saveIDs )
        if(storeLogRatios) nix=sapply( ResultSets, saveLogRatios )
        if(storePlots) {
            nix=sapply( ResultSets, plotResultSet )
            plotSpeciesLegends()
            plotSampleComposition()
        }
        if(storeSampleCVs) nix=sapply( ResultSets, saveSampleCVs )
        if(storeMetrics) saveMetrics( ResultSets )
        if(storeIdStats) logIdStatistics( ResultSets )
        if(returnResultSets) return(ResultSets)
    }
}
