library(LFQbench)

# 0. Parameters
working_dir = "./"
subfolder = "input"

#software_sources <- c("PViewNoFilter", "PeakView", "Skyline", "Spectronaut", "OpenSWATH", "DIAumpire", "Skyline")
#slopes <- c(1, 1, 1, 1, 1, 1, 1)

evalCommandLineArgumentsSafe = function(){
    args=commandArgs(trailingOnly = T)
    sapply(args[grep("<-", args)], 
           function(arg) eval(expr = parse(text=arg), envir = globalenv())
    ) 
}

evalCommandLineArguments()

working_dir = file.path(working_dir, subfolder)

names(slopes) <- software_sources
print(slopes)

# 1. Read all peptide files in folder (scaling will be based on these files).
PeptidesFiles = list.files( path=working_dir, pattern="._peptides.tsv$", 
                            full.names= FALSE, include.dirs = F, recursive = F, all.files = F)


# 2. Read also all protein files in folder (in order to be scaled too).
ProteinsFiles <- list.files( path=working_dir, pattern="._proteins.tsv$", 
                             full.names= FALSE, include.dirs = F, recursive = F, all.files = F)

# 3. Re-scale all files (peptides and proteins).
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

nix <- sapply(PeptidesFiles, rescalefile, working_dir, slopes)
nix <- sapply(ProteinsFiles, rescalefile, working_dir, slopes)
