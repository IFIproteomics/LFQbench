
#' initConfiguration
#' 
#' Initialization of LFQbench configuration settings.
#' Executing the function without parameters will use default values.
#' Use function changeConfig to change settings after initialization.
#' 
#' @param  DataRootFolder the parent folder for batch processing
#' @param  SampleComposition the quantitative composition of hybrid proteome samples
#' @param  BackgroundSpeciesName the name of non regulated background species
#' @param  InputExtensionPattern the pattern of file extension for batch input files
#' @param  CsvColumnSeparator the input file column delimiter character
#' @param  CsvDecimalPointChar the input file decimal point character
#' @param  BoxPlotWhiskerQuantile the quantile to that whiskers of boxplots will extend, t.m. (1 - quantile*2) portion of data will be inside the whiskers
#' @param  MinProteinAmount the minimum valid protein amount in the input data, protein amounts below this threshold are considered as NA
#' @param  DropInvalidLogRatio if TRUE log-ratios outside validity range will be dropped
#' @param  LogRatioValidityRange the log-ratio validity range 
#' @param  LogRatioPlotRange the log-ratio range for plotting
#' @param  LogIntensityPlotRange the log2-intensity range for plotting
#' @param  MaxLogRatioForAUQC the maximum value for AUQC quantification
#' @param  NumberOfIntensityQuantiles the number of parts for splitting the data for metrics calculation
#' @param  CenterLogRatioByBackground should we center all log-ratios by the median log-ratio of the background species
#' @param  NormalizeAmountsToPPM should we normalize amounts to ppms
#' @param  PlotWidth the plot width in inches used in pdf files
#' @param  PlotHeight the plot height in inches used in pdf files
#' @param  PlotCurveLineWidth the line thickness for curves in plots
#' @param  PlotLegendLineWidth the line thickness in legends
#' @param  PlotPointSize the point size
#' @param  ScatterPlotPointType the point type
#' @param  PlotPointMinAlpha the minimum alpha value used for coloring points in sparse regions
#' @param  AxisLabelSize the relative font size for axis labels
#' @param  AxisAnnotationSize the relative font size for axis labels
#' @param  AxisLineThickness the line thickness of axes
#' @param  par the graphical parameters like mar, mgp, las, ... to set for the plot canvases.
#' @export
initConfiguration = function( 
    # data folder for batch processing
    DataRootFolder=getwd()
    # sample composition
    ,SampleComposition = data.frame(
        species =c("HUMAN", "YEAST", "ECOLI"),
        A       =c(  65,       30,     05   ),
        B       =c(  65,       15,     20   )
    )
  # non regulated background species
  # has equal protein amounts in all samples
  ,BackgroundSpeciesName  = "HUMAN"
  # protein quantification input data CSV format
  ,InputExtensionPattern = "\\..sv$"
  ,CsvColumnSeparator = "\t"
  ,CsvDecimalPointChar = "."
  # whiskers of boxplots will extend to given quantile
  # t.m. (1 - quantile*2) portion of data will be inside the whiskers
  ,BoxPlotWhiskerQuantile = 0.02
  # protein amounts below this threshold are considered as NA
  ,MinProteinAmount    = 0.000001
  # log-ratios outside validity range will be dropped
  ,DropInvalidLogRatio = T
  ,LogRatioValidityRange  = c(-10, 10)
  # log-ratio range for plots
  ,LogRatioPlotRange   = c(-4, 4)
  ,LogIntensityPlotRange = c(9, 26)
  # value used as maximum for AUQC quantification
  ,MaxLogRatioForAUQC  = 2
  # data will be split into this number of quantiles for calculating ranged metrics
  ,NumberOfIntensityQuantiles = 3
  # if TRUE then all log-ratios will be centered
  # by median log-ratio of background species
  ,CenterLogRatioByBackground = T
  # if TRUE then all protein amounts will be translated to ppm values
  ,NormalizeAmountsToPPM = F
  # graphics settings
  # pdf canvas size in inches
  ,PlotWidth = 6
  ,PlotHeight = 5
  # line thickness
  ,PlotCurveLineWidth = 2
  ,PlotLegendLineWidth = 4
  # point size
  ,PlotPointSize = 1.5
  # point type
  ,ScatterPlotPointType = 20
  # minimum alpha threshold for point transparency ramping
  ,PlotPointMinAlpha = 0.1
  # size magnification for axis labels (cex.lab)
  ,AxisLabelSize = 2
  # size magnification for axis annotations (cex.axis)
  ,AxisAnnotationSize = 2
  # line thinkes for plotting axes
  ,AxisLineThickness = 2
  # canvas settings
  ,par = list(
    # plot area margins: c(bottom, left, top, right)
    mar = c( 4.5, 6, 0.5, 0.5 ),
    # axes layout: c(title, label, line)
    mgp = c( 3.7, 1.5, 0 ),
    # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
    las = 1
  )
  )
{
  # list argument names
  argNames = ls()
  # collect arguments and their values in a list
  cfg <<- sapply( argNames, function(n) get(n) )
  # process configuration
  setRootFolder(cfg$DataRootFolder, createSubfolders=F)
  ################################################################################
  processConfig()
}

#' recalculate some parameters from user defined configuration
processConfig = function()
{
  ################################################################################
  # process composition of samples
  cfg$AllSampleNames <<- as.vector( colnames(cfg$SampleComposition)[-1] )
  cfg$AllSpeciesNames <<- as.vector( cfg$SampleComposition[,1] )
  cfg$AllExpectedAmounts <<- as.matrix( cfg$SampleComposition[,-1] )
  rownames(cfg$AllExpectedAmounts) <<- cfg$AllSpeciesNames
  cfg$NumberOfSamples <<- length(cfg$AllSampleNames)
  cfg$NumberOfSpecies <<- length(cfg$AllSpeciesNames)
  cfg$SampleColors <<- brewer.pal(max(cfg$NumberOfSamples,3),"Dark2")[1:cfg$NumberOfSamples]
  cfg$SpeciesColors <<- brewer.pal(max(cfg$NumberOfSpecies,3),"Dark2")[1:cfg$NumberOfSpecies]
  ################################################################################
  
  ################################################################################
  cfg$AUQCRatioRange <<- c(0, cfg$MaxLogRatioForAUQC)
  ################################################################################
  
  ################################################################################
  # create sample index pairs
  cfg$SamplePairsIndices <<- createNumericPairs( 1, cfg$NumberOfSamples )
  cfg$NumberOfSamplePairs <<- nrow( cfg$SamplePairsIndices )
  cfg$SamplePairsLabels <<- apply(cfg$SamplePairsIndices, 1, function(sp){nms = cfg$AllSampleNames[sp]; return(paste(nms[1], ":" ,nms[2], sep = ""))} )
  cfg$SamplePairsColors <<- brewer.pal( max(cfg$NumberOfSamplePairs,3), "Dark2" )[1:cfg$NumberOfSamplePairs]
  ################################################################################
  
  ################################################################################
  cfg$AllSpeciesPairs <<- createNumericPairs(1, cfg$NumberOfSpecies)
  cfg$AllSpeciesPairsLabels <<- apply(cfg$AllSpeciesPairs, 1, function(sp){nms = cfg$AllSpeciesNames[sp]; return(paste(nms[1], "-", nms[2], sep=""))} )
  ################################################################################
  
  # backup graphical parameters
  cfg$parBackup <<- par()[ names(cfg$par) ]
}

#' setRootFolder
#' 
#' change the location of root folder for batch processing
#' @param rootFolder the path to the root folder
#' @export
setRootFolder = function( rootFolder=ifelse(file.exists(cfg$DataRootFolder), cfg$DataRootFolder, getwd()), createSubfolders=T )
{
  cfg$DataRootFolder <<- rootFolder
  ################################################################################
  # path to folder with protein quantification files
  cfg$InputFilesLocation <<- file.path(cfg$DataRootFolder, "input")
  # target location for plot files 
  cfg$PlotFilesLocation <<- file.path(cfg$DataRootFolder, "plot")
  # target location for log files
  cfg$LogFilesLocation  <<- file.path(cfg$DataRootFolder,"log")
  ################################################################################
  
  ################################################################################
  # create input/output paths
  if( file.exists(cfg$DataRootFolder) )
  {
    if(createSubfolders)
    {
        nix = sapply( c(cfg$InputFilesLocation, cfg$PlotFilesLocation, cfg$LogFilesLocation), mkdir )   
    }
  }
  else
  {
    cat("Folder ''"+cfg$DataRootFolder+"'' does not exist.\n")
    cat("No input/output folders were created!\n")
    cat("Please define a valid root folder.")
  }
  ################################################################################
}

#' changeConfig
#' 
#' change some of the initialized LFQbench configuration settings.
#' changeConfig allows to change some configuration settings without loading default values.
#' 
#' @param  DataRootFolder the parent folder for batch processing
#' @param  SampleComposition the quantitative composition of hybrid proteome samples
#' @param  BackgroundSpeciesName the name of non regulated background species
#' @param  InputExtensionPattern the pattern of file extension for batch input files
#' @param  CsvColumnSeparator the input file column delimiter character
#' @param  CsvDecimalPointChar the input file decimal point character
#' @param  BoxPlotWhiskerQuantile the quantile to that whiskers of boxplots will extend, t.m. (1 - quantile*2) portion of data will be inside the whiskers
#' @param  MinProteinAmount the minimum valid protein amount in the input data, protein amounts below this threshold are considered as NA
#' @param  DropInvalidLogRatio if TRUE log-ratios outside validity range will be dropped
#' @param  LogRatioValidityRange the log-ratio validity range 
#' @param  LogRatioPlotRange the log-ratio range for plotting
#' @param  LogIntensityPlotRange the log2-intensity range for plotting
#' @param  MaxLogRatioForAUQC the maximum value for AUQC quantification
#' @param  NumberOfIntensityQuantiles the log2-intensity ranges for range based plots 
#' @param  CenterLogRatioByBackground should we center all log-ratios by the median log-ratio of the background species
#' @param  NormalizeAmountsToPPM should we normalize amounts to ppms
#' @param  PlotWidth the plot width in inches used in pdf files
#' @param  PlotHeight the plot height in inches used in pdf files
#' @param  PlotCurveLineWidth the line thickness for curves in plots
#' @param  PlotLegendLineWidth the line thickness in legends
#' @param  PlotPointSize the point size
#' @param  ScatterPlotPointType the point type
#' @param  PlotPointMinAlpha the minimum alpha value used for coloring points in sparse regions
#' @param  AxisLabelSize the relative font size for axis labels
#' @param  AxisAnnotationSize the relative font size for axis labels
#' @param  AxisLineThickness the line thickness of axes
#' @param  par the graphical parameters like mar, mgp, las, ... to set for the plot canvases.
#' @export
changeConfig = function( 
  DataRootFolder= cfg$DataRootFolder
  ,SampleComposition = cfg$SampleComposition
  ,BackgroundSpeciesName  = cfg$BackgroundSpeciesName
  ,InputExtensionPattern = cfg$InputExtensionPattern
  ,CsvColumnSeparator = cfg$CsvColumnSeparator
  ,CsvDecimalPointChar = cfg$CsvDecimalPointChar
  ,BoxPlotWhiskerQuantile = cfg$BoxPlotWhiskerQuantile
  ,MinProteinAmount    = cfg$MinProteinAmount
  ,DropInvalidLogRatio = cfg$DropInvalidLogRatio
  ,LogRatioValidityRange  = cfg$LogRatioValidityRange
  ,LogRatioPlotRange   = cfg$LogRatioPlotRange
  ,LogIntensityPlotRange = cfg$LogIntensityPlotRange
  ,MaxLogRatioForAUQC  = cfg$MaxLogRatioForAUQC
  ,NumberOfIntensityQuantiles = cfg$NumberOfIntensityQuantiles
  ,CenterLogRatioByBackground = cfg$CenterLogRatioByBackground
  ,NormalizeAmountsToPPM = cfg$NormalizeAmountsToPPM
  ,PlotWidth = cfg$PlotWidth
  ,PlotHeight = cfg$PlotHeight
  ,PlotCurveLineWidth = cfg$PlotCurveLineWidth
  ,PlotLegendLineWidth = cfg$PlotLegendLineWidth
  ,PlotPointSize = cfg$PlotPointSize
  ,ScatterPlotPointType = cfg$ScatterPlotPointType
  ,PlotPointMinAlpha = cfg$PlotPointMinAlpha
  ,AxisLabelSize = cfg$AxisLabelSize
  ,AxisAnnotationSize = cfg$AxisAnnotationSize
  ,AxisLineThickness = cfg$AxisLineThickness
  ,par = cfg$par
)
{
  # list argument names
  argNames = ls()
  # collect arguments and their values in a list
  cfg <<- sapply( argNames, function(n) get(n) )
  # process configuration
  setRootFolder( DataRootFolder )
  processConfig()
}