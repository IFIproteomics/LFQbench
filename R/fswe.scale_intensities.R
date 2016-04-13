#' FSWE.scaleIntensities
#' 
#' It scales FSWE-processed files to a common intensity scale defined by one software reference.
#' 
#' Assumption: since we lose  at this point the peak retention times, 
#' we will assume in this script that most of the peaks detected by different software tools match. 
#' 
#' @param working_dir LFQbench-root directory with the subfolder input containing files processed by FSWE that you want to scale
#' @param fileModelingPattern pattern for files you want to use for modeling
#' @param filesToModelPattern pattern for files to apply the scaling models 
#' @param softwareReference software name for scaling reference
#' 
#' @export
FSWE.scaleIntensities <- function(
                                    working_dir = LFQbench.Config$DataRootFolder,
                                    fileModelingPattern = "_peptides",
                                    filesToModelPattern = ".[\\.].",
                                    filesNotToModelPattern = "FILE_I_DONT_WANT_TO_HAVE",
                                    softwareReference = "PeakView",
                                    softToModelMap = list(),
                                    softwareLabels = NULL,
                                    outputFileNameSuffix = NULL
)
{
    results_dir = file.path(working_dir, "input")
    supplementary_dir <- file.path( working_dir, "supplementary" )
    mkdir( results_dir )
    mkdir( supplementary_dir )
    
    # 1. Read all peptide files in folder (scaling will be based on these files).
    filesToModel = list.files( path=results_dir, pattern=filesToModelPattern, full.names= FALSE, include.dirs = F, recursive = F, all.files = F)
    
    # exclude files which should not be modelled
    filesToModel = filesToModel[!grepl(filesNotToModelPattern, filesToModel)]
    
    modelingFiles = list.files( path=results_dir, pattern=fileModelingPattern, 
                                full.names= FALSE, include.dirs = F, recursive = F, all.files = F)
    
    # collect present software names
    softwareNamesPresent <- sapply(filesToModel, guessModelSoftware, softToModelMap)
    softwareNamesPresent <- unique( softwareNamesPresent[!is.na(softwareNamesPresent)] )
    
    if(!any(softwareNamesPresent==softwareReference))
    {
        stop(
            paste0(
                "no files for reference software found!\n",
                "software names present: ", paste(softwareNamesPresent), "\n",
                "reference software: ", softwareReference, "\n"
            )
        )
    }
        
    
    # 2. Create a data frame with common peptides (no modifications) for each different peptide file (each software tool). 
    #    Identify each column with its corresponding software tool.
    readAndsumQuantFile <- function(qfile){
        df <- read.table(file =file.path(results_dir, qfile), header = T, sep = "\t")
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
    callm <- function(y, x, thedata, percentile = 0.98)
    {
        if(!(x == y))
        {
            outputFileBase = paste("CIS", x, y, sep="_")
            if(!is.null(outputFileNameSuffix))
                outputFileBase = paste("CIS", x, y, outputFileNameSuffix, sep="_")    
                
            data_percentile = thedata[thedata[, y] < quantile(thedata[, y], percentile), ]
            slope = lm(formula =  get(x) ~ get(y) - 1, data = data_percentile)
            
            save( thedata, 
                  LFQbench.Config, 
                  softwareLabels, 
                  y, x, slope, 
                  file=file.path(supplementary_dir, paste0(outputFileBase, ".rda") ) 
            )
            
        pdf(
            file=file.path(supplementary_dir, paste0( outputFileBase, ".pdf") ),
            width=LFQbench.Config$PlotWidth, height=LFQbench.Config$PlotHeight, family="Helvetica", pointsize=9)
        par(
            # plot area margins: c(bottom, left, top, right)
            mar = c( 5, 9, 0.5, 3 ),
            # axes layout: c(title, label, line)
            mgp = c( 3.7, 1.5, 0 ),
            # axis labels orientation: 0: parallel, 1: horizontal, 2: perpendicular, 3: vertical
            las = 1
        )
        
        #### plot model
        yLab = x
        xLab = y
        if(!is.null(softwareLabels))
        {
        if(any(names(softwareLabels)==yLab)) yLab = softwareLabels[[yLab]]
        if(any(names(softwareLabels)==xLab)) xLab = softwareLabels[[xLab]]
        }
        plot(thedata[,y], thedata[,x], xlab = "", ylab = "", pch=19,
             lwd = LFQbench.Config$AxisLineThickness, 
             cex.axis = LFQbench.Config$AxisAnnotationSize
        )
        title(xlab=xLab, mgp = c( 3.5, 1.5, 0 ), cex.lab=LFQbench.Config$AxisLabelSize)
        title(ylab=yLab, mgp = c( 7, 1.5, 0 ), cex.lab=LFQbench.Config$AxisLabelSize)
        abline(slope, col="red")
        midx = (max(thedata[, y]) - min(thedata[, x])) / 2
        lowy = (max(thedata[, y]) - min(thedata[, x])) / 5
        sumlm = summary(slope)
        rsquared = sumlm$adj.r.squared 
        text(midx, lowy, adj=c(0,0), labels = paste("slope =" , round(slope$coefficients, 2), 
                                                "adj. R^2 =", round(rsquared,2),  sep=" ") )
        dev.off()
        par(LFQbench.Config$parBackup)
            return(as.vector(slope$coefficients))
        }else{
            return(1.0)
        }
    }
    
    models <- sapply(softwareNamesPresent, callm, softwareReference, peptides.all)
    
    # 5. Re-scale all files (peptides and proteins).
    rescalefile <- function(filename, scaling){
        df <- read.table(file.path(results_dir, filename), header=T, sep="\t")
        soft <- guessModelSoftware(filename, softToModelMap)
        if(length(soft)==0) {
            return(NA)
        }
        quant.cols <- sapply(df, is.numeric)
        # skip files of the reference software
        if(scaling[soft]!=1.0) 
        {
            df[, quant.cols] <- df[, quant.cols] * scaling[soft]
            write.table(df, file.path(results_dir, filename), row.names=F, col.names = T, sep="\t")   
        }
    }
    
    nix <- sapply(filesToModel, rescalefile, models)
    
    modelsFile = "CIS.models.txt"
    if(!is.null(outputFileNameSuffix)) modelsFile = paste0("CIS.models.",outputFileNameSuffix,".txt")
    write(rbind(names(models),models), file.path(supplementary_dir, modelsFile), sep = "\t", ncolumns = 2)
    
}

# find model software for a file and replace softwares without a model by the right ones
guessModelSoftware = function(filename, softToModelMap)
{
    src = guessSoftwareSource(filename, FSWE.softwareNames, T)
    if(any(names(softToModelMap)==src))
        src = softToModelMap[[src]]
    return(src)
}   

