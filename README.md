# epiGBS-benchmarking Github

This github contains all the scripts used for the benchmarking of the epiGBS data as described in (Paper).

It requires the epiGBS output using the Arabidopsis data in the paper.
All packages and programs can be install using conda:
`conda env create -f src/env/epiGBS-benchmarking.yaml`

To generate the figures in the manuscript the scripts `src/methylation_calling_benchmarking.sh` and `src/snp_calling_benchmarking.sh` should be run.
These bash scripts contain all the commands and call all Rscripts / Snakemake scripts used in the analysis.
