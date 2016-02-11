library(LFQbench)
library(readr)

set.seed(1234)

# For the simulator ===============================================================
numReplicates = 9
numProteinsSpecies = c(2000, 1500, 1000)
peptidesPerProtein = c(3, 3, 3)
species = c("HUMAN", "YEAST", "ECOLI")
speciesRatios = c(0.0, 1.0, -2.0)
stdDevRatios = c(0.001, 0.001, 0.001)
intDistribution = c(13.0, 1.5)
stdDeviationFactorMS = 0.03
BackgroundSignalLevel = 2^14 #2^14
MissingValuesFactor = 0.01
SignalNoiseFactor = 5.0
# General LFQbench settings =======================================================
srcDir = "/Users/napedro/tmp_data/LFQbench/simulator"
expFile = "test_PeptideSummary.tsv"

replicateNames = c(paste0("A", 1:numReplicates), paste0("B", 1:numReplicates))

dataSets = data.frame(
    "Test"= replicateNames 
    ,row.names=replicateNames
)

# Define the species tags of the protein entries in your experiment
speciesTags = list(HUMAN = "_HUMAN", YEAST = "_YEAS", ECOLI = "_ECOLI") 

# =================================================================================

mySimExp <- FSWE.simExperiment(numReplicates = numReplicates,  
                   numProteinsSpecies = numProteinsSpecies, 
                   species = species, 
                   speciesRatios = speciesRatios, 
                   stdDevRatios = stdDevRatios, 
                   peptidesPerProtein = peptidesPerProtein, 
                   intDistribution = intDistribution,
                   stdDeviationFactorMS = stdDeviationFactorMS,
                   BackgroundSignalLevel = BackgroundSignalLevel,
                   MissingValuesFactor = MissingValuesFactor,
                   SignalNoiseFactor = SignalNoiseFactor)

#write.csv2(mySimExp, file = file.path(srcDir, expFile) , sep = "\t", na = "NA", row.names = F, col.names = T)

write_tsv(mySimExp, path = file.path(srcDir, expFile), na = "NA", col_names = T)

LFQbench.initConfiguration()
FSWE.initConfiguration( dataSets, speciesTags )
FSWE.addSoftwareConfiguration("test", 
                              input_format = "wide", 
                              input.extension = "*.tsv$", 
                              nastrings = "NA", 
                              protein_input = F, 
                              quantitative.var = "quant",
                              quantitative.var.tag = "quant", 
                              protein.var = "protein", 
                              sequence.mod.var = "peptide", 
                              charge.var = "charge")
FSWE.switchSoftwareConfiguration("test")

LFQbench.setDataRootFolder( srcDir, createSubfolders = T )

FSWE.generateReports(experimentFile = expFile, 
                     plotHistogram = T, 
                     plotHistNAs = F, 
                     reportSequences = F, 
                     keep_original_names = T)

nix = LFQbench.batchProcessRootFolder()
