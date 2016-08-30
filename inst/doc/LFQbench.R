## ----include=FALSE-------------------------------------------------------
# define here global settings
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, results="hide", fig.width = 7, fig.height = 5)

## ----defineSampleSetComosition-------------------------------------------
sampleComposition = data.frame( 
    species = c("HUMAN","YEAST", "ECOLI"), 
    A       = c(  67,     30,       3   ), 
    B       = c(  67,      3,      30   )
)

## ----defineDatasets------------------------------------------------------

dataSets = data.frame(
    "HYE110_SynaptG2S" = c( 
        paste(rep("HYE110_A"), 1:3, sep = "."), 
        paste(rep("HYE110_B"), 1:3, sep = ".")
    ),
    row.names = c( "A1", "A2", "A3", "B1", "B2", "B3" )
)


## ----defineSpeciesTags---------------------------------------------------
speciesTags = list(
    HUMAN = "_HUMAN", 
    YEAST = "_YEAS", 
    ECOLI = "_ECOLI"
)

## ----initLFQbench--------------------------------------------------------
library(LFQbench)

LFQbench.initConfiguration(
    SampleComposition = sampleComposition
)

FSWE.initConfiguration( 
    injectionNames = dataSets,
    speciesTags = speciesTags
)

## ----sourceDir-----------------------------------------------------------
srcDir = "../ext/data/vignette_examples/hye110"

LFQbench.setDataRootFolder( 
    rootFolder = srcDir, 
    createSubfolders = T
)

## ----listSoftwareConfig, results="asis"----------------------------------

print( paste( FSWE.softwareNames, collapse="," ) )


## ----addSoftwareConfig---------------------------------------------------

FSWE.addSoftwareConfiguration(
    # Software configuration name      
    softwareName = "ISOQuant_pep",
    
    # input_format can be wide or long. 
    # Wide contains all quantitative values (all samples and replicates) 
    # for a peptide in a single row, 
    # whereas long contains a single quantitative value (just one replicate) in a row.
    input_format = "wide",
    
    # it is important to know that LFQbench honours the extension: 
    # csv are COMMA separated values, 
    # tsv are TAB separated values
    input.extension = "*.csv$",
    
    # how NA (not available) values are reported
    nastrings = " ",
    
    # in long formats, how the quantitative value column is named
    quantitative.var = make.names("intensity in"),
    
    # in wide formats, how quantitative values are tagged 
    #(they also should include the injection names reported at the datasets object)
    quantitative.var.tag = make.names("intensity in"),
    
    # name of the protein name variable. 
    # Remember: protein names should include species information (speciesTags)
    protein.var = "entry",
    
    # variable name of sequence 
    # (including modifications as defined in FSWE.modificationsToUniMod)
    sequence.mod.var = "sequence",
    
    # variable name of the precursar charge state.
    charge.var = make.names("signal_charge")
)


## ----addModification-----------------------------------------------------

FSWE.addModification( 
    modificationRegExps = "\\[Oxi\\]",
    UniModStrings = "\\(UniMod:35\\)"
)

# list modifications available
print( FSWE.modificationsToUniMod )


## ----FSWE----------------------------------------------------------------

inputFiles = list.files(
    path = LFQbench.Config$DataRootFolder, 
    pattern = "\\..+"
)

nix = sapply(
    inputFiles, 
    FSWE.generateReports,
        softwareSource = "guess",
        keep_original_names = T,
        singleHits = F,
        plotHistogram = T, 
        plotHistNAs = T, 
        reportSequences = F
)


## ----LFQbench_analysis---------------------------------------------------

# some configuration changes for beautifying plots
LFQbench.changeConfiguration(
    LogIntensityPlotRange = c(9,21),
    LogRatioPlotRange = c(-7,7)
)

# run batch analysis and keep result set
res = LFQbench.batchProcessRootFolder()


## ----displayMetrics------------------------------------------------------

# getting the result set of the first benchmarked file
rs = res[[1]]

m = LFQbench.getMetrics( 
    resultSet = rs 
)

# get local accuracy and precision (by intensity tertiles)
acc = m$`Local accuracy`$`A:B`
prec = m$`Local precision`$`A:B`

## ------------------------------------------------------------------------

# get the benchmark result for the first sample pair of the recently used result set
samplePairRes = rs$result[[1]]

# display the scatter plot
LFQbench.showScatterAndBoxPlot( 
    samplePair = samplePairRes, 
    showLegend = T 
)

# display the distributions of log ratios
LFQbench.showDistributionDensityPlot(
    samplePair = samplePairRes, 
    showLegend = F
)


## ----completeExample-----------------------------------------------------

sampleComposition = data.frame(
  species = c("HUMAN", "YEAST", "ECOLI"),
  A       = c(  65,       30,     05   ),
  B       = c(  65,       15,     20   )
)

dataSets = data.frame(
  "HYE124_TTOF6600_64var" = c(
    "lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012", # A
    "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013"  # B
  ),
  row.names = c( "A1", "A2", "A3", "B1", "B2", "B3" )
)

speciesTags = list(
  HUMAN = "_HUMAN", 
  YEAST = "_YEAS", 
  ECOLI = "_ECOLI"
)

LFQbench.initConfiguration(
  SampleComposition = sampleComposition
)

FSWE.initConfiguration( 
  injectionNames = dataSets,
  speciesTags = speciesTags
)

# we don't need to define new software report format in this example
# because Spectronaut and PeakView (SWATH 2.0) report formats are predefined in FSWE

srcDir = "../ext/data/vignette_examples/hye124"

LFQbench.setDataRootFolder( 
  rootFolder = srcDir, 
  createSubfolders = T
)

inputFiles = list.files(
  path = LFQbench.Config$DataRootFolder, 
  pattern = "\\..+"
)

nix = sapply(
  inputFiles, 
  FSWE.generateReports,
  softwareSource = "guess",
  keep_original_names = T,
  singleHits = F,
  plotHistogram = T, 
  plotHistNAs = T, 
  reportSequences = F
)

hye124.res = LFQbench.batchProcessRootFolder()


