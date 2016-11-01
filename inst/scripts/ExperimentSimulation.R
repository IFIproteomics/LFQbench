library(LFQbench)
library(readr)

set.seed(1234)

# For the simulator ===============================================================
numReplicates = 3
numProteinsSpecies = c(2000, 1500, 1000)
peptidesPerProtein = c(3, 3, 3)
species = c("HUMAN", "YEAST", "ECOLI")
speciesRatios = c(0.0, 1.0, -2.0)
stdDevRatios = c(0.001, 0.001, 0.001)
proteinAbundanceDistribution = c(6.0, 1.5)
stdDeviationFactorMS = 0.03
BackgroundSignalLevel = 0.2 #2^14
NMARFactor = 0.05 # Not Missing at Random missing values factor
MARFactor = 0.01 # Missing at Random missing values factor
ProteinAbundanceErrorFactor = 1.5 # 5.0
# General LFQbench settings =======================================================
srcDir = "ext/simulations/simulator"
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
                               species = species, 
                               speciesRatios = speciesRatios, 
                               stdDevRatios = stdDevRatios, 
                               numProteinsSpecies = numProteinsSpecies, 
                               peptidesPerProtein = peptidesPerProtein, 
                               proteinAbundanceDistribution = proteinAbundanceDistribution,
                               stdDeviationFactorMS = stdDeviationFactorMS,
                               BackgroundSignalLevel = BackgroundSignalLevel,
                               NMARFactor = NMARFactor,
                               MARFactor = MARFactor,
                               ProteinAbundanceErrorFactor = ProteinAbundanceErrorFactor)

#write.csv2(mySimExp, file = file.path(srcDir, expFile) , sep = "\t", na = "NA", row.names = F, col.names = T)

mkdir(srcDir)
write_tsv(mySimExp, path = file.path(srcDir, expFile), na = "NA", col_names = T)

LFQbench.initConfiguration()
FSWE.initConfiguration( dataSets, speciesTags )
FSWE.addSoftwareConfiguration("test", 
                              input_format = "wide", 
                              input.extension = "*.tsv$", 
                              nastrings = "NA", 
                              protein_input = F, 
                              quantitative.var = "quant",
                              quantitative.var.tag = "^quant", 
                              protein.var = "protein", 
                              sequence.mod.var = "peptide", 
                              charge.var = "charge")
FSWE.switchSoftwareConfiguration("test")

LFQbench.setDataRootFolder( srcDir, createSubfolders = T )



FSWE.generateReports(experimentFile = expFile, 
                     plotHistogram = T, 
                     plotHistNAs = T, 
                     reportSequences = F, 
                     keep_original_names = T)

nix = LFQbench.batchProcessRootFolder()
