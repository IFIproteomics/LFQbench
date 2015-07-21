## Select Software-depending variable names #########
quantitative.var <- NA
quantitative.var.tag <- NA
protein.var <- NA
filename.var <- NA
sequence.mod.var <- NA
charge.var <- NA
q_filter_threshold <- NA
decoy.var <- NA
sheet.fdr <- NULL
protein_input <- F
input.extension <- "*.csv$"
nastrings = "NA"
decoy.tags <- c("DECOY_","reverse")
intensity.scale <- 1.0

if(software_source == "Spectronaut"){
    qvalue.var <- "EG.Qvalue"
    quantitative.var <- "FG.NormalizedTotalPeakArea"
    protein.var <- "EG.ProteinId"
    filename.var <- "R.FileName"
    sequence.mod.var <- "EG.ModifiedSequence"
    charge.var <- "FG.Charge" 
    nastrings <- "NaN"
    input.extension <- "*.tsv$"
    input_format <- "long"  # Options: "long", "wide"
    intensity.scale <- 393.5628154
}

if(software_source == "DIAumpire"){
    quantitative.var.tag <- "top6"
    quantitative.var <- "top6Area"
    protein.var <- "Protein"
    filename.var <- "ReplicateName"
    sequence.mod.var <- "ModSeq"
    charge.var <- "Charge" 
    nastrings <- "NA"
    input.extension <- "*.tsv$"
    input_format <- "wide"  # Options: "long", "wide"
    intensity.scale <- 65.0507483
}

if(software_source == "Skyline"){
    q_filter_threshold <- 0.01
    qvalue.var <- "annotation_QValue"
    quantitative.var <- "TotalArea"
    protein.var <- "ProteinName"
    filename.var <- "FileName"
    sequence.mod.var <- "ModifiedSequence"
    charge.var <- "PrecursorCharge"  
    decoy.var <- "IsDecoy"
    nastrings <- "#N/A"
    input.extension <- "*.tsv$"
    input_format <- "long"  # Options: "long", "wide"
    intensity.scale <- 0.8087744    
}

if(software_source == "PeakView"){
    q_filter_threshold <- 1.0
    quantitative.var.tag <- "Sample"
    fdr.var.tag <- "FDR"
    quantitative.var <- "TotalAreaFragment"
    protein.var <- "Protein"
    filename.var <- "R.FileName"
    sequence.mod.var <- "Peptide"
    charge.var <- "Precursor Charge"
    #input.extension <- "*.tsv$"
    input.extension <- "*.xls*"
    sheet.data <- "Area - peptides"
    sheet.fdr <- "FDR"
    input_format <- "wide"  # Options: "long", "wide"  
    intensity.scale <- 1.0
}

if(software_source == "openSWATH"){
    q_filter_threshold <- 0.01
    qvalue.var <- "m_score"
    quantitative.var <- "Intensity"
    protein.var <- "ProteinName"
    filename.var <- "filename"
    sequence.mod.var <- "FullPeptideName"
    charge.var <- "Charge"  
    decoy.var <- "decoy"
    input.extension <- "*.tsv.gz$"
    input_format <- "long"  # Options: "long", "wide"
    intensity.scale <- 3.6483363
    ### Wide format
    # q_filter_threshold <- 0.05
    # qvalue.var <- "score"
    # quantitative.var.tag <- "Intensity_"
    # quantitative.var <- "Intensity"
    # protein.var <- "Protein"
    # filename.var <- "FileName"
    # sequence.mod.var <- "Peptide"
    # charge.var <- "Charge"  
    # decoy.var <- "IsDecoy"
    # input.extension <- "*.tsv$"
}

if(software_source == "test"){
    quantitative.var <- "sumArea"
    protein.var <- "ProteinName"
    filename.var <- "ReplicateName"
    sequence.mod.var <- "PeptideSequence"
    charge.var <- "PrecursorCharge" 
    nastrings <- "NA"
    input.extension <- "*.tsv$"
    input_format <- "long"  # Options: "long", "wide"
    
}

if(software_source == "PeakView_builtin_proteins"){
    quantitative.var.tag <- "Sample"
    quantitative.var <- "TotalAreaFragment"
    sheet.data <- "Area - proteins"
    protein.var <- "Protein"
    filename.var <- "R.FileName"
    input.extension <- "*.xls*"
    input_format <- "wide"  # Options: "long", "wide"
    protein_input <- T
}

if(software_source == "DIAumpire_builtin_proteins"){
    quantitative.var.tag <- "Top6Top6Freq"
    quantitative.var <- "top6Area"
    protein.var <- "Protein.Key"
    filename.var <- "ReplicateName"
    sequence.mod.var <- "ModSeq"
    charge.var <- "Charge" 
    nastrings <- "NA"
    input.extension <- "*.tsv$"
    input_format <- "wide"  # Options: "long", "wide"
    protein_input <- T
}

quantitative.var <- gsub(" ", ".", quantitative.var)
protein.var <- gsub(" ", ".", protein.var)
filename.var <- gsub(" ", ".", filename.var)
sequence.mod.var <- gsub(" ", ".", sequence.mod.var)
charge.var <- gsub(" ", ".", charge.var)
if(!is.na(decoy.var)){
    decoy.var <- gsub(" ", ".", decoy.var)
}

#####################################################

