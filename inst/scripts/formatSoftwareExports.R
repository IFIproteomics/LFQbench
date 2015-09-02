rm(list = ls())

library(LFQbench)

DEBUG = T

#working_dir <- "../../ext/data/example_spectronaut"
working_dir <- "../../ext/data/example_peakview"
working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w_commonproteins"

# Options: see fswe.variables.R file
## With the option "guess", input files must start with the software_source name. Then input files from different software sources can be analysed together.
software_source <- "guess"    
keep_original_names <- FALSE
suffix <- "r1"
results_dir <- "input"
supplementary <- "supplementary"

# Use sequencelist when you want to analyse a subset of peptides
#sequencelist <- read.csv("/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/common_peptides/commonPeptides_Spectronaut_DIAumpire.csv", stringsAsFactors=F)$V1
#sequencelist <- read.csv("/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/common_peptides/commonPeptides_allSoftwares_6600_64w.csv", stringsAsFactors=F)$V1
sequencelist <- NULL
proteinlist  <- NULL
#proteinlist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt"
    
plotHistogram = T 
plotHistNAs = T
reportSequences = T
singleHits = F



topN.sort_method = "sum"  # "sum", "mean", "idrate_mean"  ## TODO: NOT YET IMPLEMENTED 
topNindividual = T
restrictNA = F
#topN.allowNAs = T  ## Redundant
top.N = 3 
top.N.min = 2
##

#histNAs.peptides.scale = c(0,2000)
#histNAs.proteins.scale = c(0,700)

#Override config parameters with command line parameters
LFQbench::evalCommandLineArguments()

if(!is.null(sequencelist)){
    sequencelist <- read.table(sequencelist, quote = "", header = F, sep = "\t")$V1
}
if(!is.null(proteinlist)){
    proteinlist <- read.table(proteinlist, quote = "", header = F, sep = "\t")$V1
}

source("fswe.variables.R")
setup_softwaresource_variables(software_source)
source("fswe.functions.R")
source("fswe.datasets.R")

if( !file.exists(file.path(working_dir, results_dir))) { dir.create(file.path(working_dir, results_dir)) }
if( !file.exists(file.path(working_dir, results_dir, supplementary))) { dir.create(file.path(working_dir, results_dir, supplementary)) }
AllInputFiles = list.files( path = working_dir, pattern = input.extension, full.names= FALSE )
if( software_source == "guess") {
    AllInputFiles = list.files( path=working_dir, pattern = ".[\\.].", full.names = FALSE, include.dirs = F, recursive = F, all.files = F)
}

parameter.software_source <- software_source
generateReports <- function(experiment_file, 
                            plotHistogram = F, 
                            plotHistNAs = F, 
                            reportSequences = F, 
                            sequence.list = NULL,
                            protein.list = NULL,
                            singleHits = F){
    
    #  experiment_file <- AllInputFiles[6]
    
    if(parameter.software_source == "guess") {
        software_source <- guessSoftwareSource(experiment_file, software_sources)
    }
    
    setup_softwaresource_variables(software_source)
    
    ## Ops related to the quantitation variable ##
    sumquant <- paste0("sum(", quantitative.var, ")")
    medianquant <- paste0("median(", quantitative.var,")")
    
    qvalue.filtered = FALSE
    # Read file
    original.experiment_file <- experiment_file
    experiment_file <- file.path(working_dir, experiment_file)
    cat(paste0("Generating peptide report for ", experiment_file, "\n"))
    
    if(grepl(".xls", input.extension)){
        df <- read_excel(experiment_file, sheet=sheet.data, col_names=TRUE)
        names(df) <- gsub(" ", ".", names(df))
        if(!is.na(q_filter_threshold)){
            df$sequence_z <- paste(df[[sequence.mod.var]], df[[charge.var]], sep="_")
            df$sequence_z <- as.character(df$sequence_z)
            df.fdr <- read_excel(experiment_file, sheet=sheet.fdr, col_names=TRUE)
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
        df <- read.table(experiment_file, na.strings= nastrings, header=T, 
                         sep=guessSep(experiment_file), stringsAsFactors =F, fill = T)
    }
    
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
        
    
    # Attach specie and remove peptides belonging to multiple species, and not considered species ("NA"s)
    df <- df %>% rowwise()
    df <- eval( substitute(mutate(df, "specie" = guessOrganism(var, species)), list(var = as.name(protein.var)) ) ) 
    df <- filter(df, specie != "NA", specie != "multiple")
    
    experiment <- NA
    if(input_format == "wide"){
        # At this point, it is better to transform the data frame into a long-format data frame, and
        # continue the analysis commonly for wide and long inputs.
                
        # find the columns containing the quantitative values (pivoted values). 
        # They will be use to gather the key-value pairs
        experiment <- which(sapply(experiments, guessExperiment_wide, colnames(df) ))

        # WARNING: the quantitative variables should be the only ones containing the injection names at this point.
        # We filter by the quantitative.var.tag in order to remove any other columns like score, ret time...
        tmp1 <- df[, grepl(protein.var, colnames(df))]
        tmp1 <- cbind( tmp1, df[, grepl(sequence.mod.var, colnames(df))] )
        tmp1 <- cbind( tmp1, df[, grepl(charge.var, colnames(df))] )
        tmp1 <- cbind( tmp1, df[, grepl("specie", colnames(df))])
        tmp1 <- cbind( tmp1, df[, grepl(quantitative.var.tag, colnames(df), ignore.case=T)])
        df <- tmp1
        rm(tmp1)
        
        quantvar.range <- which(grepl(quantitative.var.tag, colnames(df), ignore.case=T))

        df <- df %>% gather_(filename.var, quantitative.var, quantvar.range) %>%  
                    arrange_(protein.var, sequence.mod.var)
        
    }else if(input_format == "long"){
        # If filename.var contains the complete path -> remove it
        df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
        # Do a second check, in case there were extensions of compressed files like .mzxml.gz
        df[[filename.var]] <- basename(file_path_sans_ext(df[[filename.var]]))
        
        # Guess the experiment by the filename.var column
        injections <- distinct(select_(df, filename.var))  
        experiment <- which(sapply(experiments, guessExperiment, injections))
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
    data <- data %>% group_by_(filename.var, sequence.mod.var, protein.var , "specie") %>% 
        summarise_( quant_value = sumquant ) # proteinID = protein.var, specie = "specie" , 
    
    peptides_wide <- spread_(data, filename.var, "quant_value") 
    
    # Change variable names of injections by their right name 
    # (removing additions from software to the name, etc)
    common_names <- names(peptides_wide)[c(1:3)]
    inj_names <- names(peptides_wide)[-c(1:3)]
    inj_names <- as.character(sapply(inj_names, guessInjection, experiment))
    names(peptides_wide) <- c( common_names, inj_names )
    
    experiment.order <- match(experiments[[experiment]], names(peptides_wide[-c(1:3)])) + 3
    peptides_wide <- peptides_wide[, c(c(1:3), experiment.order)]

    #Rename the samples 
    names(peptides_wide) <- c("sequenceID", "proteinID", "specie", sample.names) 
    
    # add a sequence column (just to remove it afterwards)
    peptides_wide$sequence <- gsub( "*\\[.*?\\]", "", peptides_wide$sequenceID )
    peptides_wide$sequence <- gsub( "*\\(.*?\\)", "", peptides_wide$sequence )
    
    # Scale intensities to a common intensity scale (CIS)
    nums <- sapply(peptides_wide, is.numeric)
    peptides_wide[, nums] <- peptides_wide[, nums] * intensity.scale 
    
    #expfile_noext <- file_path_sans_ext(basename(experiment_file))
    expfile_noext <- paste(software_source, names(experiment)[1], suffix, sep="_")
    peptidereportname <- paste0(expfile_noext, "_peptides.tsv")
    proteinreportname <- paste0(expfile_noext, "_proteins.tsv")
    sequencereportname <- paste0(expfile_noext, "_sequencelist.csv")
    histPepProteinreportname <- paste0(expfile_noext, "_HistogramPeptidesProtein.pdf")
    histNAsProteinsreportname <- paste0(expfile_noext, "_proteins_HistogramNAs.pdf")
    histNAsPeptidesreportname <- paste0(expfile_noext, "_peptides_HistogramNAs.pdf")
    if(keep_original_names) {
        peptidereportname  <- paste(original.experiment_file, suffix, "peptides.tsv", sep="_") 
        proteinreportname  <- paste(original.experiment_file, suffix, "proteins.tsv", sep="_") 
        sequencereportname <- paste(original.experiment_file, suffix, "sequencelist.csv", sep="_")
        histPepProteinreportname <- paste(original.experiment_file, suffix, "HistogramPeptidesProtein.pdf", sep="_")
        histNAsProteinsreportname <- paste(original.experiment_file, suffix, "_proteins_HistogramNAs.pdf", sep="_")
        histNAsPeptidesreportname <- paste(original.experiment_file, suffix, "_peptides_HistogramNAs.pdf", sep="_")
    }
    
    if(reportSequences){
        sequence_list <- as.data.frame(unique(peptides_wide$sequence))
        write.table(sequence_list, 
                    file=file.path(working_dir, results_dir, supplementary, 
                                   sequencereportname),
                    sep=",", row.names=F, col.names=F)
    }
    
    if(!is.null(sequence.list)){
        peptides_wide <- peptides_wide %>% filter(sequenceID %in% sequence.list)
    }

    if(!is.null(protein.list)){
        filterProtein <- function(myprotein, listtosearch){
            which(grepl(myprotein, listtosearch))
        }
        
        matches <- unique(unlist(sapply(protein.list, filterProtein, peptides_wide$proteinID)))
        peptides_wide <- peptides_wide[matches, ] 
    }
    
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
        write.table(protein_report, file=file.path(working_dir, results_dir , proteinreportname), 
                    sep="\t", row.names=F, col.names=T)
        
        cat("Protein report as input -- No peptide report generated.\n")
        
        return(NA)
        
    }
    

    write.table(peptides_wide, file=file.path(working_dir, results_dir , peptidereportname), 
                sep="\t", row.names=F, col.names=T)
    
    
    ## PROTEIN REPORT
    cat(paste0("Generating protein report for ", experiment_file,"\n"))

    if(singleHits){
        print("Summarising protein single hits...")
        proteins_wide <- peptides_wide %>% 
            arrange(proteinID, specie) %>%
            group_by(proteinID, specie) %>%
            filter(n_distinct(sequenceID) == 1) %>%
            select(-sequenceID) %>% 
            summarise_each(funs(single_hits(.))) 
    }else{
        print("Summarising proteins...")
        if(topNindividual){
            print("using TOP3 individual for each run")
            proteins_wide <- peptides_wide %>% 
                select(-sequenceID) %>% 
                arrange(proteinID, specie) %>%
                group_by(proteinID, specie) %>%  
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
                group_by(proteinID, specie)
            
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
    
    write.table(proteins_wide, file=file.path(working_dir, results_dir , proteinreportname), 
                sep="\t", row.names=F, col.names=T)
    
    ## Histogram: Peptides per protein
    
    peptides_per_protein <- peptides_wide %>%
        group_by(proteinID, specie) %>%
        summarise(pep.protein = n_distinct(sequenceID))
    
    if(plotHistogram){
        pdf(file=file.path(working_dir, results_dir, supplementary, histPepProteinreportname), width=6, height=4)
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
        p <- p + geom_histogram()
        p <- p + aes() + ylab("# proteins")
        if(exists("histNAs.proteins.scale")){
            p <- p + scale_y_continuous(limits = histNAs.proteins.scale)
        }
        hNAprot <- p + facet_wrap( ~ specie, ncol = 3)  + theme_classic() +
                        scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
 
        p <- ggplot(peptides_wide, aes(x = numNAs))
        p <- p + geom_histogram()
        p <- p + aes() + ylab("# peptides")
        if(exists("histNAs.peptides.scale")){
            p <- p + scale_y_continuous(limits = histNAs.peptides.scale)
        }
        hNApep <- p + facet_wrap( ~ specie, ncol = 3)  + theme_classic() +
                       scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
        
        ggsave(filename = file.path(working_dir, results_dir, supplementary, histNAsProteinsreportname), plot = hNAprot , width=6, height=4)
        ggsave(filename = file.path(working_dir, results_dir, supplementary, histNAsPeptidesreportname), plot = hNApep,  width=6, height=4)
        
    }
}


nix <- sapply(AllInputFiles, generateReports, 
              plotHistogram = plotHistogram, 
              plotHistNAs = plotHistNAs, 
              reportSequences = reportSequences, 
              sequence.list = sequencelist,
              protein.list = proteinlist,
              singleHits = singleHits)

