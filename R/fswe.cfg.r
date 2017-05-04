#' FSWE.initConfiguration
#' 
#' define software list and the list of initialization function for each software
#' 
#' @param injectionNames : a data frame the injection names of each experiment considered. You may process several experiments together (each experiment
#'                         is a column of the data frame. The row names of the data frame are used as column names in the FSWE output.).
#'                         Example: dataSets = data.frame( "TTOF5600_32w"=c("lgillet_L150206_001", "lgillet_L150206_003", "lgillet_L150206_005",
#'                                                                          "lgillet_L150206_002", "lgillet_L150206_013", "lgillet_L150206_014"), 
#'                                                        ,"TTOF5600_64w"=c("lgillet_L150206_007", "lgillet_L150206_009", "lgillet_L150206_011", 
#'                                                                          "lgillet_L150206_008", "lgillet_L150206_010", "lgillet_L150206_012") 
#'                                                        ,"TTOF6600_32w"=c("lgillet_I150211_002", "lgillet_I150211_004", "lgillet_I150211_006",
#'                                                                          "lgillet_I150211_003", "lgillet_I150211_005", "lgillet_I150211_007") 
#'                                                        ,"TTOF6600x64w"=c("lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012",
#'                                                                          "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013")
#'                                                        ,row.names=c("A1", "A2", "A3", "B1", "B2", "B3"))
#'                                                         
#' @param speciesTags : Namesa named vector containing tags to locate the corresponding species at the protein name reported by the software. 
#'                      In SwissProt you may use something like _SPECIES. 
#'                      Names of elements must contain all LFQbench.Config$AllSpecies
#'                      Example: speciesTags = list(HUMAN = "_HUMAN", YEAST = "_YEAS", ECOLI = "_ECOLI") 
#'           
#' @export
FSWE.initConfiguration = function( injectionNames, speciesTags )
{
  FSWE.init.defaultSoftwares()
  FSWE.init.defaultModifications()
  FSWE.setDatasetInjectionNames( injectionNames )
  FSWE.setSpeciesTags( speciesTags )
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
#' @param softwareName          :   a name for this configuration set.
#' @param input_format          :   data format of the file. Options: "wide", "long"
#' @param input.extension       :   file extension. 
#'                                  Warning: FSWE interprets file extensions to read the files. 
#'                                  Ensure that files named i.e. ".tsv" are actually tab separated 
#'                                  values files. Do not use the extension ".xls" in tab separated 
#'                                  values files.
#' @param nastrings             :   string used at the file to label not available / not a number data. 
#'                                  Examples: "NA", "#N/A", "NaN".  
#' @param protein_input         :   Boolean. Set to TRUE if the files contain only protein information 
#'                                  (no sequence/precursors reports).
#' @param filename.var          :   when using "long" report formats, column head name of the injection 
#'                                  file names. 
#' @param quantitative.var      :   column head name of the quantitative values.
#' @param quantitative.var.tag  :   when using "wide" report formats, indicate a common tag for all 
#'                                  quantitative columns. 
#' @param protein.var           :   column head name of the protein names. It must contain a species tag.
#' @param sequence.mod.var      :   column head name of the sequence (including modifications).
#' @param charge.var            :   column head name of the precursor charge state.
#' @param decoy.var             :   (optional) column head name of is_decoy (yes/no) 
#' @param decoy.tags            :   vector of tags for decoy proteins/peptides. 
#'                                  Example: c("DECOY_", "reverse")
#' @param fdr.var.tag           :   when using excel formats (SWATH 2.0), indicate a tag to locate 
#'                                  the FDR columns. 
#' @param qvalue.var            :   (optional) If filtering by Q-value is required, column head name 
#'                                  of the Q-Value.
#' @param q_filter_threshold    :   (optional) Q-value threshold, if you want to filter results by Q-value.
#'                                  Note that you must specify then the qvalue.var (see above).
#' @param sheet.data            :   when using excel formats, name of the sheet containing the 
#'                                  quantitative values.
#' @param sheet.fdr             :   when using excel formats, name of the sheet containing the 
#'                                  False Discovery Rate values.
#'                                  
#' @export
FSWE.addSoftwareConfiguration = function(
  softwareName
  ,input_format = "wide"
  ,input.extension = "*.csv$"
  ,nastrings = "NA"
  ,protein_input = F
  ,filename.var = NA
  ,quantitative.var = NA
  ,quantitative.var.tag = NA
  ,protein.var = NA
  ,sequence.mod.var = NA
  ,charge.var = NA
  ,decoy.var = NA
  ,decoy.tags = c("DECOY_","reverse")
  ,fdr.var.tag = NA
  ,qvalue.var = NA
  ,q_filter_threshold = NA
  ,sheet.data = NA
  ,sheet.fdr = NULL
)
{
  FSWE.softwareNames <<- c(FSWE.softwareNames, softwareName)
  FSWE.configFunctionForSoftware[[softwareName]] <<- function(){
    FSWE.setSoftwareConfiguration(
        input_format, input.extension, nastrings, protein_input, filename.var, 
        quantitative.var, quantitative.var.tag, protein.var, sequence.mod.var, charge.var, 
        decoy.var, decoy.tags, qvalue.var, q_filter_threshold, sheet.data,
        sheet.fdr,fdr.var.tag
    )
  }    
}

FSWE.setSoftwareConfiguration = function(
    input_format = "wide"
    ,input.extension = "*.csv$"
    ,nastrings = "NA"
    ,protein_input = F
    ,filename.var = NA
    ,quantitative.var = NA
    ,quantitative.var.tag = NA
    ,protein.var = NA
    ,sequence.mod.var = NA
    ,charge.var = NA
    ,decoy.var = NA
    ,decoy.tags = c("DECOY_","reverse")
    ,qvalue.var = NA
    ,q_filter_threshold = NA
    ,sheet.data = NA
    ,sheet.fdr = NULL
    ,fdr.var.tag = NA)
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
  if(!is.na(FSWE.Config$qvalue.var)){
      FSWE.Config$qvalue.var <<- gsub(" ", ".", FSWE.Config$qvalue.var)
  }
  
}

#' FSWE.setDatasetInjectionNames
#' 
#' set injections names used in the input data
#'
#' @param theData a data.frame having the injection names of each experiment considered. You may process several experiments together (each experiment
#'                         is a column of the data frame. The row names of the data frame are used as column names in the FSWE output.).
#'                         Example: dataSets = data.frame( "TTOF5600_32w"=c("lgillet_L150206_001", "lgillet_L150206_003", "lgillet_L150206_005",
#'                                                                          "lgillet_L150206_002", "lgillet_L150206_013", "lgillet_L150206_014"), 
#'                                                        ,"TTOF5600_64w"=c("lgillet_L150206_007", "lgillet_L150206_009", "lgillet_L150206_011", 
#'                                                                          "lgillet_L150206_008", "lgillet_L150206_010", "lgillet_L150206_012") 
#'                                                        ,"TTOF6600_32w"=c("lgillet_I150211_002", "lgillet_I150211_004", "lgillet_I150211_006",
#'                                                                          "lgillet_I150211_003", "lgillet_I150211_005", "lgillet_I150211_007") 
#'                                                        ,"TTOF6600x64w"=c("lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012",
#'                                                                          "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013")
#'                                                        ,row.names=c("A1", "A2", "A3", "B1", "B2", "B3")) 
#'                                                        
#' @export
FSWE.setDatasetInjectionNames = function( theData )
{
    # TODO: check the format of theData
    FSWE.dataSets <<- theData
}

#' FSWE.setSpeciesTags
#' 
#' set expressions for identifying species in protein names, e.g. _HUMAN, _YEAST, ...
#'
#' @param speciesTags named list of tags. Names of elements must contain all LFQbench.Config$AllSpeciesNames
#' @export
FSWE.setSpeciesTags = function( speciesTags )
{
    FSWE.speciesTags <<- speciesTags
    if( !all(LFQbench.Config$AllSpeciesNames %in% names(speciesTags)) )
       stop( paste("please define species tags for: ", paste0( LFQbench.Config$AllSpeciesNames, collapse = "," ) , "!") )
}

# expandList = function(theList)
# {
#   nix = sapply( names(theList), 
#           function(pn) assign( pn,theList[[pn]], envir = parent.env(parent.env(environment()) ) )
#   )
# }



FSWE.init.defaultModifications = function()
{
    modsToUniMod <- hash()
    
    # DIA-Umpire
    .set(modsToUniMod, "\\[15\\.995\\(M\\)\\]M" = "M\\(UniMod:35\\)") 
    .set(modsToUniMod, "\\[57\\.021\\(C\\)\\]C" = "C\\(UniMod:4\\)")
    .set(modsToUniMod, "\\[-17\\.027\\(Q\\)\\]Q" = "Q\\(UniMod:28\\)")
    .set(modsToUniMod, "\\[-18\\.011\\(E\\)\\]E" = "E\\(UniMod:27\\)")
    .set(modsToUniMod, "\\[39\\.995\\(C\\)\\]C" = "C\\(UniMod:26\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]A" = "A\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]M" = "M\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]C" = "C\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]D" = "D\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]E" = "E\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]G" = "G\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]P" = "P\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]S" = "S\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]T" = "T\\(UniMod:1\\)")
    .set(modsToUniMod, "\\[42\\.011\\(N-term\\)\\]V" = "V\\(UniMod:1\\)")
    
    # Skyline
    .set(modsToUniMod, "\\[\\+16\\]" = "\\(UniMod:35\\)")
    .set(modsToUniMod, "\\[\\+57\\]" = "\\(UniMod:4\\)")
    
    # SWATH2.0
    .set(modsToUniMod, "\\[Oxi\\]" = "\\(UniMod:35\\)")
    .set(modsToUniMod, "\\[CAM\\]" = "\\(UniMod:4\\)")
    
    # Spectronaut --> Modification format depends on the library used
    # OpenSWATH --> Modifications are already reported in UniMod format
    FSWE.modificationsToUniMod <<- copy(modsToUniMod)
}

#' FSWE.addModification
#'
#' Adds a new modification string to the modifications parser. 
#'
#' @param modificationRegExps    Modification regular expressions as it may be found at the software tool report (input of FSWE). Use here the whole string you would like to substitute by a UniMod string
#' @param UniModStrings      UniMod string that will substitute modificationString at FSWE outputs.  
#'
#' @export
#'
#' @examples FSWE.addModification( "\\[Oxi\\]", "\\(UniMod:35\\)")
FSWE.addModification <- function(modificationRegExps, UniModStrings){
    
    if(!exists("FSWE.modificationsToUniMod")){
        FSWE.modificationsToUniMod <- hash()
    }
    
    .set(FSWE.modificationsToUniMod, c(modificationRegExps), c(UniModStrings))

}

#' FSWE.listModifications
#'
#' @return list of modifications used by FSWE
#' @export
#'
FSWE.listModifications <- function(){
    
    if(!exists("FSWE.modificationsToUniMod")){
        warning( "FSWE modifications parser has not been initialized yet: please use FSWE.initConfiguration first!" )
        return(NULL)
    }
    
    for(mod in ls(FSWE.modificationsToUniMod)){
        print(paste(mod, FSWE.modificationsToUniMod[[mod]] , sep = " : "))
    }
}
