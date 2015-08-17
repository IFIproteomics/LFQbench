# This script re-scales all peptides/proteins output files from formatSoftwareExports in order to have a common scale.

# Assumption: since we lose  at this point the peak retention times, we will assume in this script that most of the peaks detected by different software
#             tools match. 

library(LFQbench)

# 0. Parameters
working_dir = "./"
software_scale_base = ""
software_sources = c("PViewNoFilter")


# software_scale_base = "PViewNoFilter"
# software_sources = c("PViewNoFilter", "Spectronaut", "OpenSWATH", "DIAumpire", "Skyline")
# working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/ttof5600_64w/input"

evalCommandLineArguments()

working_dir = file.path(working_dir, "input")

# 1. Read all peptide files in folder (scaling will be based on these files).
AllInputFiles = list.files( path=working_dir, pattern=".[\\.].", full.names= FALSE, include.dirs = F, recursive = F, all.files = F)

PeptidesFiles = list.files( path=working_dir, pattern="._peptides.tsv$", 
                            full.names= FALSE, include.dirs = F, recursive = F, all.files = F)


validFiles <- sapply(PeptidesFiles, guessSoftwareSource, software_sources, TRUE)
validFiles2 <- validFiles[nchar(validFiles) > 0]

# 2. Read also all protein files in folder (in order to be scaled too).
ProteinsFiles <- list.files( path=working_dir, pattern="._proteins.tsv$", 
                             full.names= FALSE, include.dirs = F, recursive = F, all.files = F)

# 3. Create a data frame with common peptides (no modifications) for each different peptide file (each software tool). 
#    Identify each column with its corresponding software tool.

## sum all quantitative values
readAndsumQuantFile <- function(qfile, working_dir, software_sources){
    df <- read.table(file =file.path(working_dir, qfile), header = T, sep = "\t")
    quant.cols <- sapply(df, is.numeric)
    quant.val <- rowSums(df[,quant.cols], na.rm = T)
    ids <- as.data.frame(df$sequenceID)
    names(ids) <- c("sequenceID")
    ids$softsource <- guessSoftwareSource(qfile, software_sources, F)
    ids$quantitation <- quant.val
    ids <- ids %>% filter(!duplicated(sequenceID))
    return(ids)
}

peptides.all <- do.call("rbind", lapply(PeptidesFiles, readAndsumQuantFile, working_dir, software_sources))
peptides.all$softsource <- factor(peptides.all$softsource)
peptides.all <- peptides.all %>% spread(key = "softsource", value = "quantitation")

## We use only common peptides to all software sources
peptides.all <- peptides.all[complete.cases(peptides.all),]

# 4. Model each tool against the software_base.
callm <- function(y, x, thedata){
    if(!(x == y)){
        slope = lm(formula =  get(x) ~ get(y) - 1,data = thedata)    
        return(as.vector(slope$coefficients))
    }else{
        return(1.0)
    }
}

software_sources_nobase <- software_sources[software_sources != software_scale_base]
models <- sapply(software_sources, callm, software_scale_base, peptides.all)

# 5. Re-scale all files (peptides and proteins).
rescalefile <- function(filename, working_dir, scaling){
    df <- read.table(file.path(working_dir, filename), header=T, sep="\t")
    soft <- guessSoftwareSource(filename, software_sources, T)
    if(length(soft)==0) {
        return(NA)
    }
    quant.cols <- sapply(df, is.numeric)
    df[, quant.cols] <- df[, quant.cols] * scaling[soft]
    write.table(df, file.path(working_dir, filename), row.names=F, col.names = T, sep="\t")
}

nix <- sapply(PeptidesFiles, rescalefile, working_dir, models)
nix <- sapply(ProteinsFiles, rescalefile, working_dir, models)

write(rbind(names(models),models), file.path(working_dir, "CIS.models.txt"), sep = "\t", ncolumns = 2)
