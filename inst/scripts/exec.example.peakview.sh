R CMD BATCH --no-save --no-restore formatSoftwareExports.R
R CMD BATCH --no-save --no-restore '--args cfg$DataRootFolder="../../ext/data/example_peakview"' lfqbench.batch.r
