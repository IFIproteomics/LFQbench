#' FSWE.init
#' 
#' define software list and the list of initialization function for each software
#' 
#' @export
FSWE.init = function()
{
  FSWE.init.defaultSoftwares()
}


#' FSWE.switchSoftwareConfiguration
#' 
#' change the FSWE configuration to values for a specific software
#' 
#' @param softwareName
#' @export
FSWE.switchSoftwareConfiguration = function(softwareName)
{
  if( is.null(softwareName) | is.na(softwareName) | length(softwareName)==0 ) 
    stop("specifiy a valid software name!")
  if( !any( softwareName == names(FSWE.configFunctionForSoftware) ) )
  {
    stop(paste0( "software \'", softwareName, "\' is undefined.\n" ,
      "please specifiy a valid software name!") )
  }
  FSWE.configFunctionForSoftware[[softwareName]]()
}

#' FSWE.addSoftwareConfiguration
#' 
#' add a new software and its parameters to the list of softwares
#' 
#' @param softwareName
#' @param quantitative.var = NA
#' @param quantitative.var.tag = NA
#' @param protein.var = NA
#' @param filename.var = NA
#' @param sequence.mod.var = NA
#' @param charge.var = NA
#' @param q_filter_threshold = NA
#' @param decoy.var = NA
#' @param sheet.fdr = NULL
#' @param protein_input = F
#' @param input.extension = "*.csv$"
#' @param nastrings = "NA"
#' @param decoy.tags = c("DECOY_""reverse")
#' @param fdr.var.tag = NA
#' @param input_format = "wide"
#' @param qvalue.var = NA
#' @param sheet.data = NA
#' @export
FSWE.addSoftwareConfiguration = function(
  softwareName
  ,quantitative.var = NA
  ,quantitative.var.tag = NA
  ,protein.var = NA
  ,filename.var = NA
  ,sequence.mod.var = NA
  ,charge.var = NA
  ,q_filter_threshold = NA
  ,decoy.var = NA
  ,sheet.fdr = NULL
  ,protein_input = F
  ,input.extension = "*.csv$"
  ,nastrings = "NA"
  ,decoy.tags = c("DECOY_","reverse")
  ,fdr.var.tag = NA
  ,input_format = "wide"
  ,qvalue.var = NA
  ,sheet.data = NA
  )
{
  FSWE.softwareNames <<- c(FSWE.softwareNames, softwareName)
  FSWE.configFunctionForSoftware[[softwareName]] <<- function(){
    FSWE.setSoftwareConfiguration(
      quantitative.var,quantitative.var.tag,protein.var,filename.var 
      ,sequence.mod.var,charge.var,q_filter_threshold,decoy.var
      ,sheet.fdr,protein_input,input.extension,nastrings,decoy.tags,fdr.var.tag,input_format,qvalue.var,sheet.data
    )
  }    
}

#' to be called from init
FSWE.init.defaultSoftwares = function()
{
  FSWE.softwareNames <<- c(
    "DIAumpire",  "DIAumpBuiltinProteins",
    "OpenSWATH",
    "PeakView", "PViewNoFilter", "PViewBuiltinProteins",
    "Skyline", 
    "Spectronaut"
  )
  
  FSWE.configFunctionForSoftware <<- list(
    "Spectronaut" = function(){
      FSWE.setSoftwareConfiguration(
        qvalue.var = "EG.Qvalue"
        ,quantitative.var = "FG.NormalizedTotalPeakArea"
        ,protein.var = "EG.ProteinId"
        ,filename.var = "R.FileName"
        ,sequence.mod.var = "EG.ModifiedSequence"
        ,charge.var = "FG.Charge" 
        ,nastrings = "NaN"
        ,input.extension = "*.tsv$"
        ,input_format = "long"  # Options: "long", "wide"
      )
    },
    "DIAumpire" = function(){
      FSWE.setSoftwareConfiguration(
        quantitative.var.tag = "top6"
        ,quantitative.var = "top6Area"
        ,protein.var = "Protein"
        ,filename.var = "ReplicateName"
        ,sequence.mod.var = "ModSeq"
        ,charge.var = "Charge" 
        ,nastrings = "NA"
        ,input.extension = "*.tsv$"
        ,input_format = "wide"  # Options: "long", "wide"
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
    ,"DIAumpBuiltinProteins" = function(){
      FSWE.setSoftwareConfiguration(
        quantitative.var.tag = "Top6Top6Freq"
        ,quantitative.var = "top6Area"
        ,protein.var = "Protein.Key"
        ,filename.var = "ReplicateName"
        ,sequence.mod.var = "ModSeq"
        ,charge.var = "Charge" 
        ,nastrings = "NA"
        ,input.extension = "*.tsv$"
        ,input_format = "wide"  # Options: "long", "wide"
        ,protein_input = T
      )
    }
  )  
}

FSWE.setSoftwareConfiguration = function(
  quantitative.var = NA
  ,quantitative.var.tag = NA
  ,protein.var = NA
  ,filename.var = NA
  ,sequence.mod.var = NA
  ,charge.var = NA
  ,q_filter_threshold = NA
  ,decoy.var = NA
  ,sheet.fdr = NULL
  ,protein_input = F
  ,input.extension = "*.csv$"
  ,nastrings = "NA"
  ,decoy.tags = c("DECOY_","reverse")
  ,fdr.var.tag = NA
  ,input_format = "wide"
  ,qvalue.var = NA
  ,sheet.data = NA )
{
  # list argument names
  argNames = ls()
  # collect arguments and their values in a list
  FSWE.Config <<- sapply( argNames, function(n) get(n) )
  FSWE.postprocessSoftwareConfiguration()
}

FSWE.postprocessSoftwareConfiguration = function()
{
  FSWE.Config$quantitative.var <<- gsub(" ", ".", FSWE.Config$quantitative.var)
  FSWE.Config$protein.var <<- gsub(" ", ".", FSWE.Config$protein.var)
  FSWE.Config$filename.var <<- gsub(" ", ".", FSWE.Config$filename.var)
  FSWE.Config$sequence.mod.var <<- gsub(" ", ".", FSWE.Config$sequence.mod.var)
  FSWE.Config$charge.var <<- gsub(" ", ".", FSWE.Config$charge.var)
  if(!is.na(FSWE.Config$decoy.var)){
    FSWE.Config$decoy.var <<- gsub(" ", ".", FSWE.Config$decoy.var)
  }
}

# expandList = function(theList)
# {
#   nix = sapply( names(theList), 
#           function(pn) assign( pn,theList[[pn]], envir = parent.env(parent.env(environment()) ) )
#   )
# }
