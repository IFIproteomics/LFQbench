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


