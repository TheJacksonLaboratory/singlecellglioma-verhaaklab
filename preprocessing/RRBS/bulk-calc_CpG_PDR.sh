#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=calc_CpG_PDR-%j-chr_%a
#SBATCH --chdir=${WORKDIR}
#SBATCH --error=${LOG_DIR}/calc_CpG_PDR-%j-chr_%a.err
#SBATCH --output=${LOG_DIR}/calc_CpG_PDR-%j-chr_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=72:00:00
#SBATCH --array=1-22
#SBATCH --partition=compute
#SBATCH --qos=batch


# This script submits a job array for running bulk-calc_chr_CpG_PDR.R for each autosomal chromosome
# Format for calling script: 
#	sbatch --export=ALL,sample_dir=${sample_dir} bulk-calc_CpG_PDR.sh 


# Output the sample name and its corresponding JobID to output file
echo "${sample_dir}	${SLURM_ARRAY_JOB_ID}" >> ${LOG_DIR}/bulk-calc_CpG_PDR_jobIDs.txt

# Run calc_chr_CpG.PDR.R for each chromosome
Rscript --verbose bulk-calc_chr_CpG_PDR.R \
	-i ${sample_dir}/CpG_PDR \
	-n "${SLURM_ARRAY_TASK_ID}" \
	-o ${sample_dir}



