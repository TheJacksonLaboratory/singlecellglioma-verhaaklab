#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=bulk-calcDNAme_disorder-%j
#SBATCH --chdir=${WORKDIR}
#SBATCH --error=${LOG_DIR}/calcDNAme_disorder-%j.err
#SBATCH --output=${LOG_DIR}/calcDNAme_disorder-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --qos=batch


# This script generates chromosome level iput files for DNAme disorder calculation and calculates per-CpG DNAme disorder
# Format for calling script: 
#	 sbatch bulk-calcDNAme_disorder.sh ${linker_row_num}

sample_dir=$(awk -v line=${ind} 'NR==line' ${RRBS_DNAme_disorder_output_dirs.txt})

# Output the sample name and its corresponding JobID to output file
echo "${sample_dir}	${SLURM_JOB_ID}" >> ${LOG_DIR}/bulk-calcDNAme_disorder_jobIDs.txt


# Generate DNAme disorder data files for each chromosome
Rscript --verbose bulk-prepare_calc_CpG_PDR.R \
	-i ${sample_dir} \
	-o ${sample_dir}

# Run PDR calculations on each set of input data files
sbatch --job-name=bulk-calcDNAme_disorder-${sample_dir}-chr_${SLURM_ARRAY_TASK_ID} --export=ALL,sample_dir=${sample_dir} bulk-calc_CpG_PDR.sh 


