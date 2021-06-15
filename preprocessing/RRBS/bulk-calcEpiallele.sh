#!/bin/sh


# This script determines if the input methylation bam file is single or paired-end
# and then runs bulk-bam2epiallele.R to identify concordant and discordant reads

# The input is a text file containing the full file pathname of desired bam file(s)

# Format for calling script:
#	sbatch bulk-calcEpiallele.sh <BAM>


BAM="$1"

# Output the BAM file name and its corresponding JobID to output file
echo "${BAM}	${SLURM_JOBID}" >> ${LOG_DIR}/bulk-calcEpiallele_jobIDs.txt

# Minimum number of CpGs for a read to be considered for analysis
count=4

# Directory containing bed files of sample high coverage CpGs; files are assumed to have same prefix as bam file
coverage_list_dir=${COVERAGE_BED_FILE_DIR}

# Output file subdirectory. Will be created if it doesn't exist
outdir=${PWD}/$(basename $BAM .bam)
mkdir -p $outdir


### Use bulk-bam2epiallele.R to calculate the proportion of reads with more than "count"
# CpGs that contain all methylated, all unmethylated, or heterogeneously methylated CpGs

echo "Analyzing $(basename ${BAM} .bam) at: " `date`
# Determine if bam file is single or paired end
if [ "$(samtools flagstat ${BAM} | awk -F'[ ]' 'NR==6' -)" == "0 + 0 paired in sequencing" ]
then 
	# For single-end samples; convert from BAM to text format if not already done
	if [ ! -e "${outdir}/$(basename ${BAM} .bam).txt" ]
	then
		samtools view -F 512 -q 20 ${BAM} > ${outdir}/$(basename ${BAM} .bam).txt
	fi
else 
	# For paired-end samples; convert from BAM to text format if not already done
	if [ ! -e "${outdir}/$(basename ${BAM} .bam).txt" ]
	then
		samtools view -f 66 -F 512 -q 20 ${BAM} > ${outdir}/$(basename ${BAM} .bam).txt
	fi
fi

# Run bulk-bam2epiallele.R for text-format input bam
Rscript bulk-bam2epiallele.R \
-b ${outdir}/$(basename ${BAM} .bam).txt \
-l ${coverage_list_dir}/$(basename ${BAM} .bam)_10x_coverage_CpGs.bed \
-c ${count} \
-o ${outdir}

# Coordinate sort output files and remove unsorted versions
sort -k1,1V -k2,2n ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_loci-unsorted.bed > ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_loci.bed
sort -k1,1V ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_site_depth-unsorted.txt > ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_site_depth.txt
sort -k1,1V ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_siteID-unsorted.txt > ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_siteID.txt
sort -k1,1n -k2,2n ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_status-unsorted.txt > ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_status.txt

rm ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_loci-unsorted.bed
rm ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_site_depth-unsorted.txt
rm ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_siteID-unsorted.txt
rm ${outdir}/$(basename ${BAM} .bam)_epiallele_CpG_status-unsorted.txt

echo "Done analyzing $(basename ${BAM} .bam) at: " `date`
###

