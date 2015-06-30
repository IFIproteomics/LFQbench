cfg = list()

# data folder
cfg$DataRootFolder = "../../../lfqbench.testdata/"

# protein quantification input data CSV format
cfg$InputExtensionPattern = "\\..sv$"
cfg$CsvColumnSeparator = "\t"
cfg$CsvDecimalPointChar = "."

# quantitative composition of samples
cfg$SampleComposition = data.frame(
  species =c("HUMAN", "YEAST", "ECOLI"),
  A       =c(  65,       30,     05   ),
  B       =c(  65,       15,     20   )
)

# non regulated background species
# has equal protein amounts in all samples
cfg$BackgroundSpeciesName	= "HUMAN"

# whiskers of boxplots will extend to given quantile
# t.m. (1 - quantile*2) portion of data will be inside the whiskers
cfg$BoxPlotWhiskerQuantile = 0.01

# missing protein amounts or amounts below this threshold 
# will be approximated to it
# value in ppm
cfg$MinProteinAmount    = 0.000001

# log-ratios outside validity range will be droped
cfg$DropInvalidLogRatio = T
cfg$LogRatioValidityRange  = c(-10, 10)

# log-ratio range for plots
cfg$LogRatioPlotRange   = c(-4, 4)

# value used as maximum for AUQC quantification
cfg$MaxLogRatioForAUQC  = 2

# split positions
cfg$Log2IntensityRanges = rbind(
  "<2"=c(0,2),
  "<4"=c(2,4),
  "<6"=c(4,6),
  "<8"=c(6,8),
  "<10"=c(8,10),
  ">10"=c(10,100)
)

# if TRUE then all log-ratios will be centered
# by median log-ratio of background species
cfg$CenterLogRatioByBackground = T

# if TRUE then all protein amounts will be translated to ppm values
cfg$NormalizeAmountsToPPM = T

# graphics settings
# pdf canvas size in inches
cfg$PlotWidth	= 6
cfg$PlotHeight = 4
# line thickness
cfg$PlotCurveLineWidth = 2
cfg$PlotLegendLineWidth = 4
# point size
cfg$PlotPointSize = 1.5
# point type
cfg$ScatterPlotPointType = 20
# point transparency
cfg$PlotPointMinAlpha = 0.1
# size magnification for axis labels (cex.lab)
cfg$AxisLabelSize = 2
# size magnification for axis annotations (cex.axis)
cfg$AxisAnnotationSize = 2
# line thinkes for plotting axes
cfg$AxisLineThickness = 2
# canvas settings
cfg$par = list(
  # plot area margins: c(bottom, left, top, right)
  mar = c( 4.5, 6, 0.5, 0.5 ),
  # plot axis: c(title, label, line)
  mgp = c( 3.7, 1.5, 0 ),
  # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
  las = 1
)

# set flag for config initialization completeness
cfg$initialized = TRUE
