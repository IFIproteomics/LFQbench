#!/bin/bash

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_source="guess" suffix="it1"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_source="guess" suffix="it1"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_source="guess" suffix="it1"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_source="guess" suffix="it1"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_source="guess" suffix="it2" q_filter_threshold <- 0.01' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_source="guess" suffix="it2" q_filter_threshold <- 0.01' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_source="guess" suffix="it2" q_filter_threshold <- 0.01' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_source="guess" suffix="it2" q_filter_threshold <- 0.01' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w"' rename.files.R


# Copy files of iteration 1 and 2 to a common folder
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w/input

cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w/input/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w/input/

# Run the LFQbench batch
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w"' lfqbench.batch.r

# Generate the plots
#R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w"' makePlots_protOverlap_nobuiltin.r
#R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w"' makePlots_protOverlap_nobuiltin.r
#R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w"' makePlots_protOverlap_nobuiltin.r
#R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w"' makePlots_protOverlap_nobuiltin.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w"' makePlots_protOverlap.r

# Build main figures
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w/output_figures/

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w/output_figures/
pdflatex figure2_new
pdflatex figure2_new_peptides
pdflatex Supp.Figure.H
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w/output_figures/
pdflatex figure2_new
pdflatex figure2_new_peptides
pdflatex Supp.Figure.H
rm *.log *aux  *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w/output_figures/
pdflatex figure2_new
pdflatex figure2_new_peptides
pdflatex Supp.Figure.H
rm *.log *aux  *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w/output_figures/
pdflatex figure2_new
pdflatex figure2_new_peptides
pdflatex Supp.Figure.H
rm *.log *aux *.tex
