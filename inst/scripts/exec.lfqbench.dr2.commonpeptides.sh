#!/bin/bash

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_source="guess" suffix="it1" peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_source="guess" suffix="it1" peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_source="guess" suffix="it1" peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_source="guess" suffix="it1" peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_source="guess" suffix="it2"  peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" q_filter_threshold=0.01 results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_source="guess" suffix="it2"  peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" q_filter_threshold=0.01 results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_source="guess" suffix="it2"  peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" q_filter_threshold=0.01 results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"  subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" subfolder="input_commonpeptides"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_source="guess" suffix="it2"  peptidelist="/Users/napedro/github_prj/LFQbench_package/LFQbench/inst/scripts/peptide_overlap/peptides_in_all_tools.tsv" q_filter_threshold=0.01 results_dir="input_commonpeptides"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonpeptides"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" subfolder="input_commonpeptides"' rename.files.R


# Copy files of iteration 1 and 2 to a common folder
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/input

cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w/input_commonpeptides/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/input/

# Remove the empty built-in files created
rm /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/input/*_builtin*
rm /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/input/*_builtin*
rm /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/input/*_builtin*
rm /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/input/*_builtin*


# Run the LFQbench batch
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides"' lfqbench.batch.r

# Generate the plots
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides"' makePlots_protOverlap.r

# Build main figures
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/output_figures/

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonpeptides/output_figures/
pdflatex figure2_new 
pdflatex figure2_new_peptides 
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonpeptides/output_figures/
pdflatex figure2_new 
pdflatex figure2_new_peptides 
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonpeptides/output_figures/
pdflatex figure2_new 
pdflatex figure2_new_peptides 
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonpeptides/output_figures/
pdflatex figure2_new 
pdflatex figure2_new_peptides 
rm *.log *aux *.tex
