#! /bin/sh

#SBATCH --get-user-env
#SBATCH --job-name=calc_context_PDR-%j
#SBATCH --chdir=${WORKDIR}
#SBATCH --error=${LOG_DIR}/calc_context_PDR-%j.err
#SBATCH --output=${LOG_DIR}/calc_context_PDR-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --qos=batch


# Format for calling script:
#	sbatch --export=ALL,INDIR=<filename> calcDNAme_disorder.sh


# Output the input directory name and its corresponding PBS JobID to output file
echo "${INDIR}	${SLURM_JOBID}" >> ${LOG_DIR}/calc_context_PDR_jobIDs.txt


Rscript --verbose calc_global_and_context_PDR.R \
-i ${INDIR} \
-a /projects/verhaak-lab/scgp/reference/genome_annotations/hg19 \
-o ${OUTDIR}
