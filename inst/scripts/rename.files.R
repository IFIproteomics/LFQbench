library(LFQbench)

working_dir <- "./"


#working_dir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w"
evalCommandLineArguments()

################################################

working_dir = file.path(working_dir, "input")

AllInputFiles = list.files( path=working_dir, pattern=".[\\.].", full.names= FALSE, include.dirs = F, recursive = F, all.files = F)

renamefile <- function(myfile){
    #DIA-Umpire
    check <- grepl("(?=^DIAumpire)", myfile, ignore.case = T, perl = TRUE) * 
             grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
             grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_r1.tsv"))
    check <- grepl("(?=^DIAumpire)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_r2.tsv"))
    
    check <- grepl("(?=^DIAumpire)", myfile, ignore.case = T, perl = TRUE) * 
             grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
             grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_peptides_r1.tsv"))
    check <- grepl("(?=^DIAumpire)", myfile, ignore.case = T, perl = TRUE) * 
             grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
             grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_peptides_r2.tsv"))

    check <- grepl("(?=^DIAumpBuiltinProteins)", myfile, ignore.case = T, perl = TRUE) * 
             grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_builtin_r1.tsv"))
    check <- grepl("(?=^DIAumpBuiltinProteins)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"DIAumpire_builtin_r2.tsv"))
    
    #OpenSWATH
    check <- grepl("(?=^OpenSWATH)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"OpenSWATH_r1.tsv"))
    check <- grepl("(?=^OpenSWATH)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"OpenSWATH_r2.tsv"))
    
    check <- grepl("(?=^OpenSWATH)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"OpenSWATH_peptides_r1.tsv"))
    check <- grepl("(?=^OpenSWATH)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"OpenSWATH_peptides_r2.tsv"))
    
    #PeakView
    check <- grepl("(?=PViewNoFilter)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_r1.tsv"))
    check <- grepl("(?=PViewNoFilter)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_r2.tsv"))
    
    check <- grepl("(?=PViewNoFilter)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_peptides_r1.tsv"))
    check <- grepl("(?=PViewNoFilter)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_peptides_r2.tsv"))
    
    check <- grepl("(?=^PeakView)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_r1.tsv"))
    check <- grepl("(?=^PeakView)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_r2.tsv"))
    
    check <- grepl("(?=^PeakView)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_peptides_r1.tsv"))
    check <- grepl("(?=^PeakView)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_peptides_r2.tsv"))
    
    check <- grepl("(?=^PViewBuiltinProteins)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_builtin_r1.tsv"))
    check <- grepl("(?=^PViewBuiltinProteins)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"PeakView_builtin_r2.tsv"))
    
    #Skyline
    check <- grepl("(?=^Skyline)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Skyline_r1.tsv"))
    check <- grepl("(?=^Skyline)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Skyline_r2.tsv"))
    
    check <- grepl("(?=^Skyline)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Skyline_peptides_r1.tsv"))
    check <- grepl("(?=^Skyline)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Skyline_peptides_r2.tsv"))   
    
    #Spectronaut
    check <- grepl("(?=^Spectronaut)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Spectronaut_r1.tsv"))
    check <- grepl("(?=^Spectronaut)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=proteins\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Spectronaut_r2.tsv"))
    
    check <- grepl("(?=^Spectronaut)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it1_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Spectronaut_peptides_r1.tsv"))
    check <- grepl("(?=^Spectronaut)", myfile, ignore.case = T, perl = TRUE) * 
        grepl("(?=peptides\\.tsv$)", myfile, ignore.case = T, perl = TRUE) *
        grepl("(?=_it2_)", myfile, ignore.case = T, perl = TRUE)
    if(check>0) file.rename(file.path(working_dir, myfile), file.path(working_dir,"Spectronaut_peptides_r2.tsv"))      
}

nix = sapply(AllInputFiles, renamefile)
