#!/bin/bash

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_source="guess" suffix="it1" singleHits=T results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") slopes=c(1,512.1616,3.7966,94.2101,1.0956) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_source="guess" suffix="it1" singleHits=T results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") slopes=c(1,478.3153,3.8189,89.3611,1.0873) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_source="guess" suffix="it1" singleHits=T results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") slopes=c(1,296.371,3.8115,40.0234,1.0262) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_source="guess" suffix="it1" singleHits=T results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","DIAumpire","Skyline") slopes=c(1,315.0351,4.1612,44.3178,1.1174) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_source="guess" suffix="it2" singleHits=T q_filter_threshold=0.01 results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" software_sources=c("PViewNoFilter","Spectronaut","OpenSWATH","Skyline") slopes=c(1,498.2715,4.1569,161.348) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_source="guess" suffix="it2" singleHits=T q_filter_threshold=0.01 results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") slopes=c(1,495.4002,4.1311,182.1853) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_source="guess" suffix="it2" singleHits=T q_filter_threshold=0.01 results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") slopes=c(1,277.4423,3.837,150.1638) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w" subfolder="input_singleHits"' rename.files.R

R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_source="guess" suffix="it2" singleHits=T q_filter_threshold=0.01 results_dir="input_singleHits"' formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" software_sources=c("PeakView","Spectronaut","OpenSWATH","Skyline") slopes=c(1,311.1526,4.3588,174.278) subfolder="input_singleHits"' fswe.CIS.justscaling.R
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w" subfolder="input_singleHits"' rename.files.R


# Copy files of iteration 1 and 2 to a common folder
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits/input
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits
mkdir /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits/input

cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_32w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF5600_64w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_32w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration1/TTOF6600_64w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_32w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF5600_64w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_32w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits/input/
cp /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/iteration2/TTOF6600_64w/input_singleHits/*.tsv /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits/input/

# Run the LFQbench batch
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits"' lfqbench.batch.r
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits"' lfqbench.batch.r

# Generate the plots
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits"' makePlots_protOverlap.r
R CMD BATCH --no-save --no-restore '--args working_dir="/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits"' makePlots_protOverlap.r

# Build main figures
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits/output_figures/
cp ./layout/* /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits/output_figures/

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_32w_singleHits/output_figures/
pdflatex figure2_new
#rm *.log *.aux *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF5600_64w_singleHits/output_figures/
pdflatex figure2_new
#rm *.log *.aux  *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_32w_singleHits/output_figures/
pdflatex figure2_new
#rm *.log *.aux  *.tex

cd /Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/output.from.softwares/draft.v2/dataanalysis_TOF6600_64w_singleHits/output_figures/
pdflatex figure2_new
#rm *.log *.aux *.tex
