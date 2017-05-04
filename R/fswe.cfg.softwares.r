#' 
#' to be called from init
FSWE.init.defaultSoftwares = function()
{
    FSWE.configFunctionForSoftware <<- list(
        
        "DIAumpire" = function(){
            FSWE.setSoftwareConfiguration(
                quantitative.var.tag = "top6"
                ,quantitative.var = "top6Area"
                ,protein.var = c("Protein", "Proteins")
                ,filename.var = "ReplicateName"
                ,sequence.mod.var = "ModSeq"
                ,charge.var = "Charge" 
                ,nastrings = "NA"
                ,input.extension = "*.tsv$"
                ,input_format = "wide"  # Options: "long", "wide"
            )
        }
        
        ,"DIAumpBuiltinProteins" = function(){
            FSWE.setSoftwareConfiguration(
                quantitative.var.tag = "Top6Top6Freq"
                ,quantitative.var = "top6Area"
                ,protein.var = make.names("Protein Key")
                ,filename.var = "ReplicateName"
                ,sequence.mod.var = "ModSeq"
                ,charge.var = "Charge" 
                ,nastrings = "NA"
                ,input.extension = "*.tsv$"
                ,input_format = "wide"  # Options: "long", "wide"
                ,protein_input = T
            )
        }
        
        ,"OpenSWATH" = function(){
            FSWE.setSoftwareConfiguration(
                q_filter_threshold = 0.01
                ,qvalue.var = "m_score"
                ,quantitative.var = "Intensity"
                ,protein.var = "ProteinName"
                ,filename.var = "filename"
                ,sequence.mod.var = "FullPeptideName"
                ,charge.var = "Charge"  
                ,decoy.var = "decoy"
                ,input.extension = "*.tsv.gz$"
                ,input_format = "long"  # Options: "long", "wide"
            )
        }
        
        ,"PeakView" = function(){
            FSWE.setSoftwareConfiguration(
                q_filter_threshold = 0.01
                ,quantitative.var.tag = "Sample"
                ,fdr.var.tag = "FDR"
                ,quantitative.var = "TotalAreaFragment"
                ,protein.var = "Protein"
                ,filename.var = "R.FileName"
                ,sequence.mod.var = "Peptide"
                ,charge.var = "Precursor Charge"
                ,input.extension = "*.xls*"
                ,sheet.data = "Area - peptides"
                ,sheet.fdr = "FDR"
                ,input_format = "wide"
            )
        }
        
        ,"PViewNoFilter" = function(){
            FSWE.setSoftwareConfiguration(
                q_filter_threshold = 1.0
                ,quantitative.var.tag = "Sample"
                ,fdr.var.tag = "FDR"
                ,quantitative.var = "TotalAreaFragment"
                ,protein.var = "Protein"
                ,filename.var = "R.FileName"
                ,sequence.mod.var = "Peptide"
                ,charge.var = "Precursor Charge"
                ,input.extension = "*.xls*"
                ,sheet.data = "Area - peptides"
                ,sheet.fdr = "FDR"
                ,input_format = "wide"  
            )
        }
        
        ,"PViewBuiltinProteins" = function(){
            FSWE.setSoftwareConfiguration(
                quantitative.var.tag = "Sample"
                ,quantitative.var = "TotalAreaFragment"
                ,sheet.data = "Area - proteins"
                ,protein.var = "Protein"
                ,filename.var = "R.FileName"
                ,input.extension = "*.xls*"
                ,input_format = "wide"  # Options: "long", "wide"
                ,protein_input = T
            )
        }

        ,"Skyline" = function(){
            FSWE.setSoftwareConfiguration(
                q_filter_threshold = 0.01
                ,qvalue.var = "annotation_QValue"
                ,quantitative.var = "TotalArea"
                ,protein.var = "ProteinName"
                ,filename.var = "FileName"
                ,sequence.mod.var = "ModifiedSequence"
                ,charge.var = "PrecursorCharge"  
                ,decoy.var = "IsDecoy"
                ,nastrings = "#N/A"
                ,input.extension = "*.tsv$"
                ,input_format = "long"  # Options: "long", "wide"
            )
        }
        
        ,"Spectronaut" = function(){
            FSWE.setSoftwareConfiguration(
                qvalue.var = "EG.Qvalue"
                ,decoy.var = "EG.IsDecoy"
                ,quantitative.var = "FG.NormalizedTotalPeakArea"
                ,protein.var = "EG.ProteinId"
                ,filename.var = "R.FileName"
                ,sequence.mod.var = "EG.ModifiedSequence"
                ,charge.var = "FG.Charge" 
                ,nastrings = "NaN"
                ,input.extension = "*.tsv$"
                ,input_format = "long"  # Options: "long", "wide"
            )
        }
        
        ,"ISOQuant_peptide" = function(){
            FSWE.setSoftwareConfiguration(
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
        }
    ) 
    
    FSWE.softwareNames <<- names( FSWE.configFunctionForSoftware )
}