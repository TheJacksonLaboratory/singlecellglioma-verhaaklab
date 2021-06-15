#!/bin/bash


# This script runs bismark methylation extractor for an input BAM file 
# Format for calling script:
#	bash ./bismark_meth_extract.sh ${SLURM_ARRAY_TASK_ID}

# Define the fastq BAM file listed in the .txt file as well as the sample-specific barcode.
BAM_FILE="$1"

# Set up output filenames and locations.
OUT=${OUTDIR}

# Location to genome reference FASTA file. GRCh37/hg19. Prepared with Bowtie2.
GENOME=/projects/verhaak-lab/scgp/reference/genomes/hsapiens/ENSEMBL/grch37/

# Run the Bismark Methylation Extractor for bam file aligned to reference genome.
bismark_methylation_extractor -p --multicore 2 --merge_non_CpG --no_header --gzip --bedGraph --comprehensive --buffer_size 10G --cytosine_report --genome_folder ${GENOME} ${OUT}/aligned_bam/${BAM_FILE}.bam -o ${OUT}/bedGraph
