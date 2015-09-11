# This script re-scales all peptides/proteins output files from formatSoftwareExports in order to have a common scale.

# Assumption: since we lose  at this point the peak retention times, we will assume in this script that most of the peaks detected by different software
#             tools match. 

# library(LFQbench)

# 0. Parameters
working_dir = "./"
subfolder = "input"
software_scale_base = ""
software_sources = c("PViewNoFilter")


# software_scale_base = "PViewNoFilter"
# software_sources = c("PViewNoFilter", "Spectronaut", "OpenSWATH", "DIAumpire", "Skyline")
# working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/ttof5600_64w/input"

evalCommandLineArguments()

# working_dir = file.path(working_dir, "input")
working_dir = file.path(working_dir, subfolder)


#' FSWE.scaleIntensities
#' 
#' It scales FSWE-processed files to a common intensity scale defined by one software reference.
#' 
#' @param working_dir directory containing files processed by FSWE that you want to scale
#' @param fileModelingPattern pattern for files you want to use for modeling
#' @param filesToModelPattern pattern for files to apply the scaling models 
#' @param softwareReference software name for scaling reference
#' 
#' @export
FSWE.scaleIntensities <- function(
                                    working_dir = file.path(LFQbench.Config$DataRootFolder, "input"),
                                    fileModelingPattern = "_peptides",
                                    filesToModelPattern = ".[\\.].",
                                    softwareReference   = "PeakView"
                                    
){
    # 1. Read all peptide files in folder (scaling will be based on these files).
    filesToModel = list.files( path=working_dir, pattern=filesToModelPattern, full.names= FALSE, include.dirs = F, recursive = F, all.files = F)
    
    modelingFiles = list.files( path=working_dir, pattern=fileModelingPattern, 
                                full.names= FALSE, include.dirs = F, recursive = F, all.files = F)
    
    
    validFiles <- sapply(modelingFiles, guessSoftwareSource, software_sources, TRUE)
    validFiles2 <- validFiles[nchar(validFiles) > 0]
    
    # 2. Create a data frame with common peptides (no modifications) for each different peptide file (each software tool). 
    #    Identify each column with its corresponding software tool.
    readAndsumQuantFile <- function(qfile){
        df <- read.table(file =file.path(working_dir, qfile), header = T, sep = "\t")
        quant.cols <- sapply(df, is.numeric)
        quant.val <- rowSums(df[,quant.cols], na.rm = T)
        ids <- as.data.frame(df$sequenceID)
        names(ids) <- c("sequenceID")
        ids$softsource <- guessSoftwareSource(qfile, FSWE.softwareNames, F)
        ids$quantitation <- quant.val
        ids <- ids %>% filter(!duplicated(sequenceID))
        return(ids)
    }
    
    peptides.all <- do.call("rbind", lapply(modelingFiles, readAndsumQuantFile))
    peptides.all$softsource <- factor(peptides.all$softsource)
    peptides.all <- peptides.all %>% spread(key = "softsource", value = "quantitation")
    
    ## We use only common peptides to all software sources
    peptides.all <- peptides.all[complete.cases(peptides.all),]
    
    # 4. Model each tool against the software_base.
    callm <- function(y, x, thedata, percentile = 0.98){
        if(!(x == y)){
            data_percentile = thedata[thedata[, y] < quantile(thedata[, y], percentile), ]
            slope = lm(formula =  get(x) ~ get(y) - 1,data = data_percentile) 
            pdf(file = file.path(working_dir, paste("CIS", x, y, ".pdf", sep="_")), width = 6, height = 6)
            plot(thedata[,y], thedata[,x], xlab = y, ylab = x, pch=19)
            abline(slope, col="red")
            midx = (max(thedata[, y]) - min(thedata[, x])) / 2
            lowy = (max(thedata[, y]) - min(thedata[, x])) / 5
            sumlm = summary(slope)
            rsquared = sumlm$adj.r.squared 
            text(midx, lowy, adj=c(0,0), labels = paste("slope =" , round(slope$coefficients, 2), 
                                                        "adj. R^2 =", round(rsquared,2),  sep=" ") )
            dev.off()
            return(as.vector(slope$coefficients))
        }else{
            return(1.0)
        }
    }
    
    software_sources_nobase <- software_sources[software_sources != software_scale_base]
    models <- sapply(software_sources, callm, software_scale_base, peptides.all)
    
    # 5. Re-scale all files (peptides and proteins).
    rescalefile <- function(filename, scaling){
        df <- read.table(file.path(working_dir, filename), header=T, sep="\t")
        soft <- guessSoftwareSource(filename, software_sources, T)
        if(length(soft)==0) {
            return(NA)
        }
        quant.cols <- sapply(df, is.numeric)
        df[, quant.cols] <- df[, quant.cols] * scaling[soft]
        write.table(df, file.path(working_dir, filename), row.names=F, col.names = T, sep="\t")
    }
    
    nix <- sapply(filesToModel, rescalefile, models)
    
    write(rbind(names(models),models), file.path(working_dir, "CIS.models.txt"), sep = "\t", ncolumns = 2)
    
}


