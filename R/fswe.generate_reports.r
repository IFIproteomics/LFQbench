#####################################

#' FSWE.generateReports
#' 
#' It generates peptide and protein reports compatible with LFQbench, by following the data processing
#' suggested in (#Ref_to_article). 
#' 
#' @param experimentFile
#' @export
FSWE.generateReports <- function(
                            experimentFile, 
                            plotHistogram = F, 
                            plotHistNAs = F, 
                            reportSequences = F, 
                            peptide.list = NULL,
                            protein.list = NULL,
                            singleHits = F,
                            softwareSource="guess",
                            working_dir = LFQbench.Config$DataRootFolder,
                            keep_original_names = FALSE,
                            outputFileNameSuffix = NULL,
                            minimumHitsPerSample = 2
                            )
{
    dataSets = as.list(FSWE.dataSets)
    speciesTags = FSWE.speciesTags
    sample.names = names(FSWE.dataSets)
    # TODO: check if datasets are defined
    
    topN.sort_method = "sum"  # "sum", "mean", "idrate_mean"  ## TODO: NOT YET IMPLEMENTED 
    topNindividual = T
    restrictNA = F
    #topN.allowNAs = T  ## Redundant
    top.N = 3 
    top.N.min = 2
    
    results_dir <- file.path( working_dir, "input" )
    supplementary_dir <- file.path( working_dir, "supplementary" )
    LFQbench.setDataRootFolder(working_dir, createSubfolders = F)
    mkdir( results_dir )
    mkdir( supplementary_dir )
    
  if(softwareSource == "guess") {
    softwareSource <- guessSoftwareSource(experimentFile, FSWE.softwareNames)
  }
  
  FSWE.switchSoftwareConfiguration(softwareSource)
  
  # expand config items as variables to the environment of current function
  nix = sapply( names(FSWE.Config), 
        function(pn) assign( pn,FSWE.Config[[pn]], envir = parent.env(environment()) ) 
  )

  ## Ops related to the quantitation variable ##
  sumquant <- paste0("sum(", quantitative.var, ")")
  medianquant <- paste0("median(", quantitative.var,")")
  
  qvalue.filtered = FALSE
  # Read file
  original.experimentFile <- experimentFile
  experimentFile <- file.path(working_dir, experimentFile)
  cat(paste0("Generating peptide report for ", experimentFile, "\n"))
  
  if(grepl(".xls", input.extension)){
    df <- read_excel(experimentFile, sheet=sheet.data, col_names=TRUE)
    names(df) <- gsub(" ", ".", names(df))
    if(!is.na(q_filter_threshold)){
      df$sequence_z <- paste(df[[sequence.mod.var]], df[[charge.var]], sep="_")
      df$sequence_z <- as.character(df$sequence_z)
      df.fdr <- read_excel(experimentFile, sheet=sheet.fdr, col_names=TRUE)
      names(df.fdr) <- gsub(" ", ".", names(df.fdr))
      df.fdr$sequence_z <- paste(df.fdr[[sequence.mod.var]], df.fdr[[charge.var]], sep="_")
      
      #remove possible duplicates
      df <- df %>% filter(!duplicated(sequence_z))
      df.fdr <- df.fdr %>% filter(!duplicated(sequence_z))
      
      tmp1 <- data.frame(sequence_z = as.character(df.fdr$sequence_z))
      tmp1$sequence_z <- as.character(tmp1$sequence_z)
      valid_peptides <- df.fdr[, grepl(fdr.var.tag, colnames(df.fdr), ignore.case=T)] <= q_filter_threshold
      tmp1 <- cbind(tmp1, valid_peptides)
      df.fdr <- tmp1
      rm(tmp1)
      df2 <- left_join(df, df.fdr, by= "sequence_z" )
      df2 <- unique(df2)
      df2.flags <- df2[,grepl(fdr.var.tag, colnames(df2), ignore.case=T)]
      df2.values <-df[,grepl(quantitative.var.tag, colnames(df2), ignore.case=T) & 
                        !grepl(fdr.var.tag, colnames(df2), ignore.case=T)] 
      df2.values[!df2.flags] <- NA
      rm(df2); rm(df2.flags); rm(valid_peptides); rm(df.fdr)
      df <- df[ ,!grepl(quantitative.var.tag, colnames(df), ignore.case=T)]
      df <- cbind(df, df2.values)
      
      qvalue.filtered = TRUE
    }
  }else{
    #df <- read.table(experimentFile, na.strings= nastrings, header=T, 
    #                 sep=guessSep(experimentFile), stringsAsFactors =F, fill = T)

    if(guessSep(experimentFile) == ","){
        df <- read_csv(experimentFile, na = nastrings)
    }
    if(guessSep(experimentFile) == "\t"){
        df <- read_tsv(experimentFile, na = nastrings)
    }
        
  }
  
  # read_csv/read_tsv will produce spaces in column names
  # we need to replace them by . for compatibility with read.table
  names(df) <- make.names(names(df))
  
  # match column names if a column name parameter is a list
  theEnvir = environment()
  matchColumnNames <- function(paramName = "protein.var")
  {
      # find which element matches the real column name 
      # set the var to this value
      paramValues = get0(paramName)
      
      # don't do anything if the value is NULL OR NA
      if( is.null(paramValues) || is.na(paramValues) ) return(F)
      
      # for any kind of vector/array/list/etc.
      if( length(paramValues)>1 ) 
      {
          columnNames = names(df)
          matchingValues = paramValues %in% columnNames
          # the parameter must be there!!!
          if(!any(matchingValues)) 
          {
              stop( paste0("software configuration does not match the input file!\n",
                           "file: ", experimentFile, "\n",
                           paramName, ": ", paste(paramValues, collapse=", "), "\n",
                           "column names: ", paste(columnNames, collapse=", "), "\n")
              )
          }
          assign(paramName, paramValues[matchingValues][1], envir = theEnvir)
      }
  }
  
  # list all variables
  colVars = ls()
  # extract variables' names identifiying data columns
  colVars = colVars[grep("\\.var$", colVars)]
  # check if every parameter is present in the header of input file
  nix = sapply( colVars, matchColumnNames )
  
  # FOR DEBUGGING ONLY
  # sapply( colVars, function(n) cat(paste0("\t", n, " = ", get0(n), "\n")) )
  
  if(protein_input){
    # If the input is already a peptide report, we need to "fake" some columns as a temporary solution 
    # to process the files
    df <- df %>% 
      mutate_( sequence = protein.var, charge = 1)
    sequence.mod.var = "sequence"
    charge.var = "charge"
  }
  
  if(!is.na(q_filter_threshold) & !qvalue.filtered){
    df <- eval( substitute(filter(df, var < q_filter_threshold), list( var = as.name(qvalue.var)) ) )
    qvalue.filtered = TRUE
  }    
  
  # Remove rows containing decoy tags
  for(dectag in decoy.tags){
    df <- df[!grepl(dectag, df[[protein.var]], ignore.case = T),]
  }
  
  # If there is a "decoy" column with boolean values, filter to no decoys
  if(!is.na(decoy.var)){
      #Evaluate first how decoys are reported in the data matrix (numeric (0 or 1) or boolean (false or true))
      is_decoytype_int = typeof(as.matrix(df[, decoy.var])) == "integer"
      if(is_decoytype_int){
          df <- df[tolower(df[[decoy.var]]) == 0,]
      }else{
          df <- df[tolower(df[[decoy.var]]) == "false",]    
      }
  }
  
  # Attach species and remove peptides belonging to multiple species, and not considered species ("NA"s)
  df <- df %>% rowwise()
  df <- eval( substitute( mutate(df, "species" = guessOrganism(var, speciesTags)), list(var = as.name(protein.var) ) ) ) 
  df <- filter(df, species != "NA", species != "multiple")
  
  experiment <- NA
  if(input_format == "wide"){
    # At this point, it is better to transform the data frame into a long-format data frame, and
    # continue the analysis commonly for wide and long inputs.
    
    # find the columns containing the quantitative values (pivoted values). 
    # They will be use to gather the key-value pairs
    experiment <- which( sapply(dataSets, guessExperiment_wide, colnames(df) ) )
    
    # WARNING: the quantitative variables should be the only ones containing the injection names at this point.
    # We filter by the quantitative.var.tag in order to remove any other columns like score, ret time...
    #tmp1 <- df[, grepl(protein.var, colnames(df))]
    #tmp1 <- cbind( tmp1, df[, grepl(sequence.mod.var, colnames(df))] )
    #tmp1 <- cbind( tmp1, df[, grepl(charge.var, colnames(df))] )
    #tmp1 <- cbind( tmp1, df[, grepl("species", colnames(df))])
    tmp1 <- df[, protein.var]
    tmp1 <- cbind( tmp1, df[, sequence.mod.var] )
    tmp1 <- cbind( tmp1, df[, charge.var] )
    tmp1 <- cbind( tmp1, df[, "species"])
    
    # For the quantitative data, we filter by columns that have: the quantitative.var.tag AND any of the injection names
    quant.columns <- grepl(quantitative.var.tag, colnames(df), ignore.case=T) 
    injection.columns <- Reduce(`|`, lapply(dataSets[[experiment]], grepl, colnames(df)))
    tmp1 <- cbind( tmp1, df[, quant.columns & injection.columns])
    
    df <- tmp1
    rm(tmp1)
    
    quantvar.range <- which( grepl( quantitative.var.tag, colnames(df), ignore.case=T ) )
    quantvar.colnames <- names(df)[quantvar.range]
    if(is.na(filename.var))
        filename.var <- "filename.var"
    df <- df %>% gather_(filename.var, quantitative.var, quantvar.colnames) %>%
                 arrange_(protein.var, sequence.mod.var)
    
  }else if(input_format == "long"){
    # If filename.var contains the complete path -> remove it
    df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
    # Do a second check, in case there were extensions of compressed files like .mzxml.gz
    df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
    
    # Guess the experiment by the filename.var column
    injections <- distinct(select_(df, filename.var))  
    experiment <- which(sapply(dataSets, guessExperiment, injections))
  }
  
  # Remove any duplicate: (based on: injectionfile, ModifiedSequence, ChargeState)
  # (These cases are likely to be due to several entries in the library)
  data <- df %>% distinct_(filename.var, sequence.mod.var, charge.var)  # I am not sure I removed all duplicates, or there is still one value per duplicate
  
  # Change zeroes to NA values at the quantitative variable
  data[[quantitative.var]][data[[quantitative.var]] == 0] <- NA
  
  #Remove NAs of the quant variable  
  # TODO: this does not work properly, it may remove everything if there is an empty (NA) variable
  #data <- data %>% na.omit() 
  
  # For each peptide: Sum quantitative values of charge states
  data <- data %>% group_by_(filename.var, sequence.mod.var, protein.var , "species") %>% 
    summarise_( quant_value = sumquant ) # proteinID = protein.var, species = "species" , 
  
  peptides_wide <- spread_(data, filename.var, "quant_value") 
  
  # Change variable names of injections by their right name 
  # (removing additions from software to the name, etc)
  common_names <- names(peptides_wide)[c(1:3)]
  inj_names <- names(peptides_wide)[-c(1:3)]
  inj_names <- as.character( sapply(inj_names, guessInjection, dataSets[[experiment]]) )
  names(peptides_wide) <- c( common_names, inj_names )
  
  experiment.order <- match(dataSets[[experiment]], names(peptides_wide[-c(1:3)])) + 3
  peptides_wide <- peptides_wide[, c(c(1:3), experiment.order)]
  
  #Rename the samples 
  names(peptides_wide) <- c("sequenceID", "proteinID", "species", rownames(FSWE.dataSets) )
  # names(peptides_wide) <- c("sequenceID", "proteinID", "species", inj_names )
  
  # add a sequence column (just to remove it after reporting naked sequences)
  peptides_wide$sequence <- gsub( "*\\[.*?\\]", "", peptides_wide$sequenceID )
  peptides_wide$sequence <- gsub( "*\\(.*?\\)", "", peptides_wide$sequence )
  
  
  # TODO: Move these hard-coded modifications to a config file
  # Convert modifications to UniMod
  peptides_wide$sequenceID <- gsub("\\[15\\.995\\(M\\)\\]M",      "M\\(UniMod:35\\)", peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[57\\.021\\(C\\)\\]C",      "C\\(UniMod:4\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[-17\\.027\\(Q\\)\\]Q",     "Q\\(UniMod:28\\)", peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[-18\\.011\\(E\\)\\]E",     "E\\(UniMod:27\\)", peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[39\\.995\\(C\\)\\]C",      "C\\(UniMod:26\\)", peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]A", "A\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]M", "M\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]C", "C\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]D", "D\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]E", "E\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]G", "G\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]P", "P\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]S", "S\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]T", "T\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[42\\.011\\(N-term\\)\\]V", "V\\(UniMod:1\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[\\+16\\]",                 "\\(UniMod:35\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[\\+57\\]",                 "\\(UniMod:4\\)",   peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[Oxi\\]",                   "\\(UniMod:35\\)",  peptides_wide$sequenceID)
  peptides_wide$sequenceID <- gsub("\\[CAM\\]",                   "\\(UniMod:4\\)",   peptides_wide$sequenceID)
  
  # nums <- sapply(peptides_wide, is.numeric)
    
  #expfile_noext <- file_path_sans_ext(basename(experimentFile))
  expfile_noext <- paste(softwareSource, names(experiment)[1], sep="_")
  if(!is.null(outputFileNameSuffix)) 
  {
    expfile_noext <- paste(expfile_noext, outputFileNameSuffix, sep="_")
  }
  peptidereportname <- paste0(expfile_noext, "_peptides.tsv")
  proteinreportname <- paste0(expfile_noext, "_proteins.tsv")
  sequencereportname <- paste0(expfile_noext, "_sequencelist.csv")
  histPepProteinreportname <- paste0(expfile_noext, "_HistogramPeptidesProtein.pdf")
  histNAsProteinsreportname <- paste0(expfile_noext, "_proteins_HistogramNAs.pdf")
  histNAsPeptidesreportname <- paste0(expfile_noext, "_peptides_HistogramNAs.pdf")
    if(keep_original_names) 
    {
        file_base_name <- original.experimentFile
        if(!is.null(outputFileNameSuffix)) 
        {
          file_base_name <- paste(file_base_name, outputFileNameSuffix, sep="_")
        }
          
        peptidereportname  <- paste(file_base_name, "peptides.tsv", sep="_") 
        proteinreportname  <- paste(file_base_name, "proteins.tsv", sep="_") 
        sequencereportname <- paste(file_base_name, "sequencelist.csv", sep="_")
        histPepProteinreportname <- paste(file_base_name, "HistogramPeptidesProtein.pdf", sep="_")
        histNAsProteinsreportname <- paste(file_base_name, "_proteins_HistogramNAs.pdf", sep="_")
        histNAsPeptidesreportname <- paste(file_base_name, "_peptides_HistogramNAs.pdf", sep="_")
    }
  
  if(reportSequences){
    sequence_list <- as.data.frame(unique(peptides_wide$sequence))
    write.table(sequence_list, 
                file=file.path(supplementary_dir, sequencereportname),
                sep=",", row.names=F, col.names=F)
  }
  
  if(!is.null(peptide.list)){
    peptides_wide <- peptides_wide %>% filter(sequenceID %in% peptide.list)
  }
  
  if(!is.null(protein.list)){
    filterProtein <- function(myprotein, listtosearch){
      which(grepl(myprotein, listtosearch))
    }
    
    matches <- unique(unlist(sapply(protein.list, filterProtein, peptides_wide$proteinID)))
    peptides_wide <- peptides_wide[matches, ] 
  }
  
  attr(peptides_wide, "vars") <- NULL
  peptides_wide <- peptides_wide %>% select(-sequence) 
  # Remove "empty" peptides (all values are NAs).
  common_names <- !sapply(peptides_wide, is.numeric)
  nums <- sapply(peptides_wide, is.numeric)
  peptides_wide <- peptides_wide %>% ungroup() %>%
    mutate(totalIntensity = rowSums(.[nums], na.rm=T)) %>%
    filter(totalIntensity > 0) %>%
    select(-totalIntensity)
  
  
  
  if(protein_input){
    # If the input was already a protein report, we can finish here
    protein_report <- peptides_wide %>% select(-sequenceID)
    write.table(protein_report, file=file.path(results_dir , proteinreportname), 
                sep="\t", row.names=F, col.names=T)
    
    cat("Protein report as input -- No peptide report generated.\n")
    
    return(NA)
    
  }
  
  
  write.table(peptides_wide, file=file.path(results_dir , peptidereportname), 
              sep="\t", row.names=F, col.names=T)
  
  
  ## PROTEIN REPORT
  cat(paste0("Generating protein report for ", experimentFile,"\n"))
  
  if(singleHits){
    print("Summarising protein single hits...")
    proteins_wide <- peptides_wide %>% 
      arrange(proteinID, species) %>%
      group_by(proteinID, species) %>%
      filter(n_distinct(sequenceID) == 1) %>%
      select(-sequenceID) %>% 
      summarise_each(funs(single_hits(.))) 
  }else{
    print("Summarising proteins...")
    if(topNindividual){
      print("using TOP3 individual for each run")
      proteins_wide <- peptides_wide %>% 
        select(-sequenceID) %>% 
        arrange(proteinID, species) %>%
        group_by(proteinID, species) %>%  
        summarise_each(funs(avg_top_n(., top.N, top.N.min))) 
    }else{
      print("using consensus TOP3")
      nums <- sapply(peptides_wide, is.numeric)
      proteins_wide <- peptides_wide %>% ungroup() %>%
        mutate(totalIntensity = rowSums(.[nums], na.rm=T))  %>% 
        group_by(proteinID) %>%
        arrange(desc(totalIntensity)) %>%
        filter(row_number() <= top.N & n() >= top.N.min) %>%
        select(-sequenceID, -totalIntensity) %>%
        group_by(proteinID, species)
      
      if(restrictNA){
        proteins_wide <- proteins_wide %>%
          summarise_each(funs(mean))
      }else{
        proteins_wide <- proteins_wide %>%
          summarise_each(funs(avgNA))
      }
      
    }
  }  
  
  # Remove "empty" proteins (all values are NAs). TODO: I wish I could find a more elegant way to do it. I am tired.
  common_names <- !sapply(proteins_wide, is.numeric)
  nums <- sapply(proteins_wide, is.numeric)
  proteins_wide <- proteins_wide %>% ungroup() %>%
    mutate(totalIntensity = rowSums(.[nums], na.rm=T)) %>%
    filter(totalIntensity > 0) %>%
    select(-totalIntensity)
  
  
  #Remove proteins, which have less than x hits per sample (among the two samples)
  cols_sampleA <- grep("^A", names(proteins_wide))
  cols_sampleB <- grep("^B", names(proteins_wide))
  
  proteins_wide$numHitsA <- apply(proteins_wide, 1, function(elt, cols) length(cols) - sum(is.na(elt[cols])), cols_sampleA)
  proteins_wide$numHitsB <- apply(proteins_wide, 1, function(elt, cols) length(cols) - sum(is.na(elt[cols])), cols_sampleB)
  proteins_wide <- proteins_wide %>% rowwise() %>% 
                                mutate(maxHitsPerSample = max(numHitsA, numHitsB)) %>% 
                                ungroup()
  
  proteins_wide <- proteins_wide %>% 
                        filter(maxHitsPerSample >= minimumHitsPerSample) %>%
                        select(-numHitsA, -numHitsB, -maxHitsPerSample)
  
  write.table(proteins_wide, file=file.path(results_dir , proteinreportname), 
              sep="\t", row.names=F, col.names=T)
  
  ## Histogram: Peptides per protein
  
  peptides_per_protein <- peptides_wide %>%
    group_by(proteinID, species) %>%
    summarise(pep.protein = n_distinct(sequenceID))
  
  if(plotHistogram){
    pdf(file=file.path(supplementary_dir, histPepProteinreportname), width=6, height=4)
    h <- hist(x=peptides_per_protein$pep.protein, breaks=c(0:20,30,40,50,100,1000), main=NULL, xlab="Peptides per protein",
              xlim = c(0,20), plot=T,freq=T)
    dev.off()
  }
  
  ## Histograms: missing values across samples (peptides and proteins)
  peptides_wide$numNAs <- apply(peptides_wide, 1, function(elt) sum(is.na(elt)))
  proteins_wide$numNAs <- apply(proteins_wide, 1, function(elt) sum(is.na(elt)))
  peptides_wide$numNAs <- factor(peptides_wide$numNAs, levels = seq(0, length(which(nums))))
  proteins_wide$numNAs <- factor(proteins_wide$numNAs, levels = seq(0, length(which(nums))))
  
  if(plotHistNAs){
    p <- ggplot(proteins_wide, aes(x = numNAs))
    p <- p + geom_bar() #geom_histogram()
    p <- p + aes() + ylab("# proteins")
    if(exists("histNAs.proteins.scale")){
      p <- p + scale_y_continuous(limits = histNAs.proteins.scale)
    }
    hNAprot <- p + facet_wrap( ~ species, ncol = 3)  + theme_classic() +
      scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
    
    p <- ggplot(peptides_wide, aes(x = numNAs))
    p <- p + geom_bar() # geom_histogram()
    p <- p + aes() + ylab("# peptides")
    if(exists("histNAs.peptides.scale")){
      p <- p + scale_y_continuous(limits = histNAs.peptides.scale)
    }
    hNApep <- p + facet_wrap( ~ species, ncol = 3)  + theme_classic() +
      scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
    
    ggsave(filename = file.path(supplementary_dir, histNAsProteinsreportname), plot = hNAprot , width=6, height=4)
    ggsave(filename = file.path(supplementary_dir, histNAsPeptidesreportname), plot = hNApep,  width=6, height=4)
  }
}