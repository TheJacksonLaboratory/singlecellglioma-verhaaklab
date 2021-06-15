#!/bin/sh


# This script determines if the input methylation bam files are single or paired-end
# and then runs bam2epiallele.R to identify concordant and discordant reads

# The input is a text file containing the full file pathname of desired bam file(s)

# Format for calling script:
#	qsub -v BAMS=<filename> calcEpimutation_v3.sh

# Format for calling script:
#	sbatch calcEpiallele.sh <BAMS>


BAMS="$1"

# Output the BAM file list name and its corresponding Job ID to log file
echo "${BAMS}	${SLURM_JOBID}" >> ${LOG_DIR}/calcEpiallele_jobIDs.txt

# Minimum number of CpGs for a read to be considered for analysis
count=4

# Output file subdirectory. Will be created if it doesn't exist
outdir=${PWD}/$(basename $BAMS .txt)
mkdir -p $outdir

### Loop through input file and process bams in parallel
# The function bam_processing uses bam2epiallele.R to calculate the proportion of reads with more than "count"
# CpGs that contain all methylated, all unmethylated, or heterogeneously methylated CpGs
bam_processing() {
	echo "Analyzing $(basename ${1} .bam) at: " `date`
	# Determine if bam file is single or paired end
	if [ "$(samtools flagstat ${1} | awk -F'[ ]' 'NR==6' -)" == "0 + 0 paired in sequencing" ]
	then 
		# For single-end samples; convert from BAM to text format if not already done
		if [ ! -e "${2}/$(basename ${1} .bam).txt" ]
		then
			samtools view -F 512 -q 20 ${1} > ${2}/$(basename ${1} .bam).txt
		fi
	else 
		# For paired-end samples; convert from BAM to text format if not already done
		if [ ! -e "${2}/$(basename ${1} .bam).txt" ]
		then
			samtools view -f 66 -F 512 -q 20 ${1} > ${2}/$(basename ${1} .bam).txt
		fi
	fi

	# Run bam2PDR.R for text-format input bam
	Rscript bam2epiallele.R -b ${2}/$(basename ${1} .bam).txt -c ${3} -o ${2}

	# Coordinate sort output files and remove unsorted versions
	sort -k1,1V -k2,2n ${2}/$(basename ${1} .bam)_epiallele_CpG_loci-unsorted.bed > ${2}/$(basename ${1} .bam)_epiallele_CpG_loci.bed
	sort -k1,1V ${2}/$(basename ${1} .bam)_epiallele_CpG_site_depth-unsorted.txt > ${2}/$(basename ${1} .bam)_epiallele_CpG_site_depth.txt
	sort -k1,1V ${2}/$(basename ${1} .bam)_epiallele_CpG_siteID-unsorted.txt > ${2}/$(basename ${1} .bam)_epiallele_CpG_siteID.txt
	sort -k1,1n -k2,2n ${2}/$(basename ${1} .bam)_epiallele_CpG_status-unsorted.txt > ${2}/$(basename ${1} .bam)_epiallele_CpG_status.txt

	rm ${2}/$(basename ${1} .bam)_epiallele_CpG_loci-unsorted.bed
	rm ${2}/$(basename ${1} .bam)_epiallele_CpG_site_depth-unsorted.txt
	rm ${2}/$(basename ${1} .bam)_epiallele_CpG_siteID-unsorted.txt
	rm ${2}/$(basename ${1} .bam)_epiallele_CpG_status-unsorted.txt

	echo "Done analyzing $(basename ${1} .bam) at: " `date`
}

export -f bam_processing

parallel 'bam_processing {1} {2} {3}' ::: "$(cat ${BAMS})" ::: "$outdir" ::: $count


# Combine the cell epiallele metrics tables into one sample table
echo "Combining ${BAMS} PDR tables at: " `date`
echo -e "Sample\tTotal_Reads\tAll_Meth_Reads\tAll_Meth_Pct\tMixed_Meth_Reads\tMixed_Meth_Pct\tAll_Unmeth_Reads\tAll_Unmeth_Pct\tDiscarded_Reads\tDiscarded_Pct\tUnique_CpGs" \
| cat - ${outdir}/*_read_distribution.txt > $(basename $BAMS .txt)_read_distribution.txt
echo "Done combining ${BAMS} PDR tables at: " `date`

# Combine the cell epiallele CpG IDs into one sample table
echo "Combining ${BAMS} epiallele CpG ID tables at: " `date`
cat ${outdir}/*_epiallele_CpG_siteID.txt > $(basename $BAMS .txt)_epiallele_CpG_siteID.txt
echo "Done combining ${BAMS} epiallele CpG ID tables at: " `date`
