# data folder
DataRootFolder = "data"

# protein quantification input data CSV format
InputExtensionPattern = "\\.tsv$"
CsvColumnSeparator = "\t"
CsvDecimalPointChar = "."

# quantitative composition of samples
SampleComposition = data.frame(
  species =c("HUMAN", "YEAST", "ECOLI"),
  A       =c(  65,       30,     05   ),
  B       =c(  65,       15,     20   )
)

# non regulated background species
# has equal protein amounts in all samples
BackgroundSpeciesName	= "HUMAN"

# whiskers of boxplots will extend to given quantile
# t.m. (1 - quantile*2) portion of data will be inside the whiskers
BoxPlotWhiskerQuantile = 0.025

# missing protein amounts or amounts below this threshold 
# will be approximated to it
# value in ppm
MinProteinAmount    = 0.0001

# log-ratios outside validity range will be droped
DropInvalidLogRatio = T
LogRatioValidityRange  = c(-10, 10)

# log-ratio range for plots
LogRatioPlotRange   = c(-4, 4)

# value used as maximum for AUQC quantification
MaxLogRatioForAUQC  = 2

# if TRUE then all log-ratios will be centered
# by median log-ratio of background species
CenterLogRatioByBackground = T

# if TRUE then all protein amounts will be translated to ppm values
NormalizeAmountsToPPM = F

# path to folder with protein quantification files
InputFilesLocation = file.path(DataRootFolder, "input")

# target location for plot files 
PlotFilesLocation = file.path(DataRootFolder, "plot")

# target location for log files
LogFilesLocation   = file.path(DataRootFolder,"log")

# graphics settings
# pdf canvas size in inches
PlotWidth	= 6
PlotHeight	= 4
# line thickness
PlotCurveLineWidth = 2
PlotLegendLineWidth = 4
# point size
PlotPointSize = 1
# point type
ScatterPlotPointType = 20
# point transparency
PlotPointAlpha = 0.5
# canvas settings
GraphicsParameters = list(
  # plot area margins: c(bottom, left, top, right)
  mar = c( 3, 3.2, 0.5, 0.5 ),
  # plot axis: c(title, label, line)
  mgp = c( 2, 0.6, 0 ),
  # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
  las = 1
)
