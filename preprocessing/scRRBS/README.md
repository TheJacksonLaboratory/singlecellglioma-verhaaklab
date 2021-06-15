Preprocessing steps:

1) Run bismark_meth_extract.sh for each BAM file, generating a bismark coverage file
2) Run bulk-calcEpiallele.sh for a set of BAM files corresponding to a given tumor, creating a directory of epiallele output files
3) Run calc_DNAme_disorder.sh for each directory containing cell output files for one tumor.
