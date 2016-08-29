## ----loadLibraries, echo=TRUE, warning=FALSE, message=FALSE--------------

library(LFQbench)


## ----defineDatasets------------------------------------------------------

# Define the data sets you want to analyse. 
dataSets = data.frame(
    "HYE124_TTOF6600_64var"=c("lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012", # A
                               "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013") # B
    ,"HYE110_SynaptG2S"=c( paste(rep("HYE110_A"), 1:3, sep = "."), paste(rep("HYE110_B"), 1:3, sep = "."))
    ,row.names=c("A1", "A2", "A3", "B1", "B2", "B3")
)


## ----speciesTags---------------------------------------------------------

speciesTags = list(HUMAN = "_HUMAN", YEAST = "_YEAS", ECOLI = "_ECOLI") 


## ----initLFQbench--------------------------------------------------------

LFQbench.initConfiguration()
FSWE.initConfiguration( dataSets, speciesTags )


## ----addSoftwareConfig---------------------------------------------------

##### Add configuration for ISOQuant (peptides) 
FSWE.addSoftwareConfiguration("ISOQuant_pep",               
                              # Software configuration name      
                              input_format = "wide",        
                              # input_format can be wide or long. Wide contains all quantitative values (all samples and replicates) 
                              # for a peptide in a single row, whereas long contains a single quantitative value (just one replicate)
                              # in a row.
                              input.extension = "*.csv$",
                              # it is important to know that LFQbench honours the extension: csv are COMMA separated values, tsv are TAB separated values
                              nastrings = " ",
                              # how NA (not available) values are reported
                              quantitative.var = make.names("intensity in"),
                              # in long formats, how the quantitative value column is named
                              quantitative.var.tag = make.names("intensity in"),
                              # in wide formats, how quantitative values are tagged 
                              #(they also should include the injection names reported at the datasets object)
                              protein.var = "entry",
                              # name of the protein name variable. Remember: protein names should include species information (speciesTags)
                              sequence.mod.var = "sequence",
                              # variable name of sequence (including modifications as defined in FSWE.modificationsToUniMod)
                              charge.var = make.names("signal_charge")
                              # variable name of the precursar charge state.
                              )
                            


## ----addModification-----------------------------------------------------

# Add a modification
FSWE.addModification( "\\[Oxi\\]", "\\(UniMod:35\\)")

# See modifications available
print(FSWE.modificationsToUniMod)



## ----sourceDir-----------------------------------------------------------

## Source directory and input files
srcDir = "../ext/data/vignette_examples"

##### Set the data root folder and create the subfolder structure
LFQbench.setDataRootFolder( srcDir, createSubfolders = T )


## ----LFQbench_analysis---------------------------------------------------

# LFQbench analysis

#LFQbench.changeConfiguration(SampleComposition = )
LFQbench.Config$LogIntensityPlotRange = c(9,20)

res <- LFQbench.batchProcessRootFolder()
myr <- res[1]

#LFQbench.showAllSpeciesQC(myr)
#LFQbench.plotResultSet(res[1])
#LFQbench.plotResultSet(res[2], showScatterPlot = F, showScatterAndDensityPlot = T, showBoxPlot = F, showDensityPlot = F)

LFQbench.explainMetrics()



