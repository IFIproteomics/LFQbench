#!/bin/bash

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_source="guess" suffix="it1" proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_source="guess" suffix="it1" proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_source="guess" suffix="it1" proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_source="guess" suffix="it1" proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") software_scale_base="PViewNoFilter" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_source="guess" suffix="it2"  proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" q_filter_threshold=0.01 results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_source="guess" suffix="it2"  proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" q_filter_threshold=0.01 results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_source="guess" suffix="it2"  proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" q_filter_threshold=0.01 results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView"  subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" subfolder="input_commonproteins"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_source="guess" suffix="it2"  proteinlist="/Users/napedro/git_repos/LFQbench/inst/scripts/protein_overlap/list_common_proteins.txt" q_filter_threshold=0.01 results_dir="input_commonproteins"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") software_scale_base="PeakView" subfolder="input_commonproteins"' fswe.CIS.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" subfolder="input_commonproteins"' rename.files.R


# Copy files of iteration 1 and 2 to a common folder
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins/input

cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w_commonproteins/input_commonproteins/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins/input/

# Run the LFQbench batch
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins"' lfqbench.batch.r

# Generate the plots
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins"' makePlots_protOverlap.r

# Build main figures
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins/output_figures/

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_commonproteins/output_figures/
pdflatex figure2_new 
pdflatex figure4_commonProteins
pdflatex figure4_new_commonProteins
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_commonproteins/output_figures/
pdflatex figure2_new 
pdflatex figure4_commonProteins
pdflatex figure4_new_commonProteins
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_commonproteins/output_figures/
pdflatex figure2_new 
pdflatex figure4_commonProteins
pdflatex figure4_new_commonProteins
rm *.log *aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_commonproteins/output_figures/
pdflatex figure2_new 
pdflatex figure4_commonProteins
pdflatex figure4_new_commonProteins
rm *.log *aux *.tex
