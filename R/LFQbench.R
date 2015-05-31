#########################################################################/**
#' benchLFQExperiment
#'
#' Analyze complete experiment from original folder and copy the results to the output folder
#
#' @param projectDir A string specifying the input directory
#' @param ouptDir    A string specifying the output directory
#' 
benchLFQExperiment <- function(projectDir, outputDir)
{
    cfg <- loadConfig(projectDir)
    
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
    
    #' @export
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
    
}