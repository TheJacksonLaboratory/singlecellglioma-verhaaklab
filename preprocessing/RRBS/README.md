Preprocessing steps:

1) Run bismark_meth_extract.sh for each BAM file, generating a bismark coverage file
2) Run bulk-calcEpiallele.sh for each BAM file, creating a directory of epiallele output files
3) Run bulk-calcDNAme_disorder.sh for each subdirectory of epiallele output files, creating per-CpG PDR
4) Run bulk-merge_filter_CpG_PDR.R for per-chromosome per-CpG PDR files for one sample, merging and filtering per-CpG PDR

