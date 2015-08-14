# This script re-scales all peptides/proteins output files from formatSoftwareExports in order to have a common scale.

# Assumption: since we lose  at this point the peak retention times, we will assume in this script that most of the peaks detected by different software
#             tools match. 

library(LFQbench)

loadLibrary("VennDiagram")
loadLibrary("tidyr")
loadLibrary("dplyr")
loadLibrary("ggplot2")
loadLibrary("readxl")
loadLibrary("classInt")
loadLibrary("scales")
loadLibrary("pander")
loadLibrary("GGally")

panel.txt <- function(x, y, labels, cex, font, ...)
{
    lab <- labels
    lab <- gsub(pattern="_int", "", lab)
    text(0.5, 0.5, lab, cex=cex, font=font)
}

panel.plot <- function(x, y) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    ct <- cor.test(x,y)
    sig <- symnum(ct$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
    r <- ct$estimate
    rt <- format(r, digits=2)[1]
    cex <- 0.3/strwidth(rt)
    
    text(.5, .5, rt, cex=cex * abs(r))
    text(.8, .8, sig, cex=cex*0.8, col='blue')
}

panel.smooth <- function (x, y) {
    points(x, y, pch=20)
    abline(lm(y~x -1), col="red")
    lines(stats::lowess(y~x), col="blue")
}
# 0. Parameters
working_dir = "./"
software_scale_base = "PViewNoFilter"
software_sources = c("PViewNoFilter", "Spectronaut", "OpenSWATH", "DIAumpire", "Skyline")

working_dir = "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/ttof5600_32w/input"

# 1. Read all peptide files in folder (scaling will be based on these files).
AllInputFiles = list.files( path=working_dir, pattern=".[\\.].", full.names= FALSE, include.dirs = F, recursive = F, all.files = F)

PeptidesFiles = list.files( path=working_dir, pattern="._peptides.tsv$", full.names= FALSE, include.dirs = F, recursive = F, all.files = F)


validFiles <- sapply(PeptidesFiles, guessSoftwareSource, software_sources, TRUE)
validFiles2 <- validFiles[nchar(validFiles) > 0]

# 2. Read also all protein files in folder (in order to be scaled too).

# 3. Create a data frame with common peptides (no modifications) for each different peptide file (each software tool). 
#    Identify each column with its corresponding software tool.

# 4. Model each tool against the software_base.

# 5. Re-scale all files (peptides and proteins).
