## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Revert GDC-aligned legacy BAM to unaligned SAM file
## Clear BAM attributes, such as re-calibrated base quality scores, alignment information
## Moreover, BAM file is split into per-readgroup BAM files
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6484
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## NOTE THAT
## --MAX_DISCARD_FRACTION=0.05
## SHOULD BE USED FOR TESTING PURPOSES ONLY!!
## This has now been added to "config[revertsam_extra_args]"
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## 05/29/18
## Because of a snakemake limitation with dynamic output (as this rule outputs an variable)
## number of files depending on the number of readgroups, we have added python code to
## create empty files ("touch") for those readgroups not relevant for this analysis
## Issue: https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
## Unlikely this will get fixed any time soon
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## 08/28/18
## TCGA BAMs often have uninformative readgroup tags. Added an additional step
## "bam2ubam" to rename these tags. This includes the RGID tag. Added a variable
## ALIQUOT_TO_LEGACY_RGID that maps to legacy RGID tags. The variable ALIQUOT_TO_RGID
## maps to the new (sometiomes renamed) RGIDs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

from python.glassfunc import touch_file

rule revertsam:
    input:
        lambda wildcards: ancient(manifest.getSourceBAM(wildcards.aliquot_barcode))
    output:
        map = "results/align/revertsam/{aliquot_barcode}/{aliquot_barcode}.output_map.txt",
        bams = temp(expand("results/align/revertsam/{{aliquot_barcode}}/{{aliquot_barcode}}.{rg}.revertsam.bam", rg = manifest.getAllReadgroups(limit_bam = True)))
    params:
        dir = "results/align/revertsam/{aliquot_barcode}",
        mem = CLUSTER_META["revertsam"]["mem"],
        walltime = CLUSTER_META["revertsam"]["walltime"]
    threads:
        CLUSTER_META["revertsam"]["ppn"]
    log: 
        "logs/align/revertsam/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/revertsam/{aliquot_barcode}.txt"
    message:
        "Reverting sample back to unaligned BAM file, stripping any previous "
        "pre-processing and restoring original base quality scores. Output files are split "
        "by readgroup.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Input BAM: {input}"
    run:
        old_rg = manifest.getLegacyRGIDs(wildcards["aliquot_barcode"]) # if wildcards["aliquot_barcode"] in ALIQUOT_TO_LEGACY_RGID else ""
        new_rg = manifest.getRGIDs(wildcards["aliquot_barcode"])
        ## Create a readgroup name / filename mapping file
        rgmap = pd.DataFrame(
            {
                "READ_GROUP_ID": old_rg,
                "OUTPUT": ["results/align/revertsam/{aliquot_barcode}/{aliquot_barcode}.{rg}.revertsam.bam".format(aliquot_barcode = wildcards["aliquot_barcode"], rg = rg) for rg in new_rg]
            },
            columns = ["READ_GROUP_ID", "OUTPUT"]
        )
        rgmap.to_csv(output["map"], sep="\t", index=False)

        ## Create empty files ("touch") for readgroups not in this BAM file
        ## Workaround for issue documented here: https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
        other_rg_f = ["results/align/revertsam/{aliquot_barcode}/{aliquot_barcode}.{rg}.revertsam.bam".format(aliquot_barcode = wildcards["aliquot_barcode"], rg = rg) for rg in manifest.getRGIDsNotInAliquot(wildcards["aliquot_barcode"])] #manifest.getAllReadgroups(limit_bam = True)] #manifest.getRGIDsNotInAliquot(wildcards["aliquot_barcode"])]
        for f in other_rg_f:
            touch_file(f)

        shell("gatk --java-options -Xmx{params.mem}g RevertSam \
            --INPUT={input} \
            --OUTPUT_BY_READGROUP=true \
            --OUTPUT_BY_READGROUP_FILE_FORMAT=bam \
            --OUTPUT_MAP={output.map} \
            --RESTORE_ORIGINAL_QUALITIES=true \
            --VALIDATION_STRINGENCY=SILENT \
            --ATTRIBUTE_TO_CLEAR=AS \
            --ATTRIBUTE_TO_CLEAR=FT \
            --ATTRIBUTE_TO_CLEAR=CO \
            --ATTRIBUTE_TO_CLEAR=XT \
            --ATTRIBUTE_TO_CLEAR=XN \
            --ATTRIBUTE_TO_CLEAR=OC \
            --ATTRIBUTE_TO_CLEAR=OP \
            --SANITIZE=true \
            --SORT_ORDER=queryname {config[revertsam_extra_args]}\
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Replaces existing readgroup tag with a more informative one
# Removing RGDT because of bug
# See https://gatkforums.broadinstitute.org/gatk/discussion/12085/fastqtosam-no-value-found-for-tagged-argument-for-iso-8601-run-date-parameter
#             --RGDT=\"{params.RGDT}\" \
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule bam2ubam:
    input:
        bam = lambda wildcards: ancient(manifest.getSourceBAM(wildcards.aliquot_barcode)), ### Input BAM is not actually used, but needs to check if present to prioritize this rule over others
        rgbam = "results/align/revertsam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.revertsam.bam"
    output:
        temp("results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam")
    params:
        RGID = lambda wildcards: wildcards.readgroup,                                                     ## ID
        RGPL = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform"), ## Platform
        RGPU = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform_unit"), ## Platform unit
        RGLB = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_library"), ## Library
        RGDT = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_timestamp"), ## Date
        RGSM = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_sample_id"), ## Sample name
        RGCN = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_center"), ## Center
        mem = CLUSTER_META["bam2ubam"]["mem"],
        walltime = CLUSTER_META["bam2ubam"]["walltime"]
    threads:
        CLUSTER_META["bam2ubam"]["ppn"]
    log:
        "logs/align/bam2ubam/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/bam2ubam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Updating RG tags\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g AddOrReplaceReadGroups \
            --INPUT={input.rgbam} \
            --OUTPUT={output} \
            --RGID=\"{params.RGID}\" \
            --RGPU=\"{params.RGPU}\" \
            --RGSM=\"{params.RGSM}\" \
            --RGPL=\"{params.RGPL}\" \
            --RGLB=\"{params.RGLB}\" \
            --RGCN=\"{params.RGCN}\" \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1" 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert from FASTQ pair to uBAM
## This step eases follow up steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fq2ubam:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        temp("results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam")
    params:
        RGID = lambda wildcards: wildcards.readgroup,                                                     ## ID
        RGPL = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform"), ## Platform
        RGPU = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform_unit"), ## Platform unit
        RGLB = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_library"), ## Library
        RGDT = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_timestamp"), ## Date
        RGSM = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_sample_id"), ## Sample ID
        RGCN = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_center"), ## Center
        mem = CLUSTER_META["fq2ubam"]["mem"]
    threads:
        CLUSTER_META["fq2ubam"]["ppn"]
    log:
        "logs/align/fq2ubam/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/fq2ubam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Converting FASTQ file to uBAM format\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g FastqToSam \
            --FASTQ={input.R1} \
            --FASTQ2={input.R2} \
            --OUTPUT={output} \
            --READ_GROUP_NAME=\"{params.RGID}\" \
            --PLATFORM_UNIT=\"{params.RGPU}\" \
            --SAMPLE_NAME=\"{params.RGSM}\" \
            --PLATFORM=\"{params.RGPL}\" \
            --LIBRARY_NAME=\"{params.RGLB}\" \
            --SEQUENCING_CENTER=\"{params.RGCN}\" \
            --SORT_ORDER=queryname \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

ruleorder: fq2ubam > bam2ubam

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on uBAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## 6/11/18
## Added --extract parameter, which includes a "summary.txt" file which gives WARN, FAIL, PASS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fastqc:
    input:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam"
    output:
        "results/align/fastqc/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned_fastqc.html"
    params:
        dir = "results/align/fastqc/{aliquot_barcode}",
        mem = CLUSTER_META["fastqc"]["mem"],
        walltime = CLUSTER_META["fastqc"]["walltime"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["fastqc"]["ppn"]
    log:
        "logs/align/fastqc/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/fastqc/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Running FASTQC\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Illumina Adapters
## Add XT tag to read records to mark the 5' start position of adapter sequences
## Adapter sequences are then removed by subsequent steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markadapters:
    input:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam"
    output:
        bam = temp("results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.bam"),
        metric = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.metrics.txt"
    params:
        mem = CLUSTER_META["markadapters"]["mem"],
        walltime = CLUSTER_META["revertsam"]["walltime"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["markadapters"]["ppn"]
    params:
        mem = CLUSTER_META["markadapters"]["mem"]
    log: 
        dynamic("logs/align/markadapters/{aliquot_barcode}.{readgroup}.log")
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MarkIlluminaAdapters \
            --INPUT={input} \
            --OUTPUT={output.bam} \
            --METRICS={output.metric} \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## (1) BAM to FASTQ
## Converts cleaned-up uBAM to FASTQ format, one FASTQ pair per readgroup
##
## (2) Align reads using BWA-MEM
## This is optimzed for >70 bp PE reads and this pipeline will need update if we intend
## to use it with shorter reads
##
## (3) Merge BAM Alignment
## Restore altered data and apply and adjust meta information lost during alignment
## ie. restore all the original readgroup tags
##
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
##
## Update 06/01: Added markadapter metrics as input even though not required. Metrics file
## is saved all the way at the end of markadapters step, and adding it makes sure that the
## input data is good
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin I am new to this step. Are read groups, including RGID per read incorporated in final merge bam
# I am more familiar with bwa mem -R '@RG\tID:foo\tSM:bar' format. 
# https://github.com/TheJacksonLaboratory/glass_wgs_alignment/blob/d72fb20659bd20fddf952d331533b9ffd88d446e/runner/preprocess_fqs.R#L103
# @sbamin If original bam file (and so coverted fastq) are using ref genome with chr prefix for chromosomes, 
# does workflow take care of removing chr prefix (or vice versa) while aligning to ref genome from GATK b37 legacy bundle?

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.bam",
        metric = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.metrics.txt"
    output:
        bam = temp("results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.realn.bam"),
        bai = temp("results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.realn.bai")
    threads:
        CLUSTER_META["samtofastq_bwa_mergebamalignment"]["ppn"]
    conda:
        "../envs/align.yaml"
    params:
        mem = CLUSTER_META["samtofastq_bwa_mergebamalignment"]["mem"],
        walltime = CLUSTER_META["samtofastq_bwa_mergebamalignment"]["walltime"]
    log: 
        "logs/align/samtofastq_bwa_mergebamalignment/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/revertsam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
            --INPUT={input.bam} \
            --FASTQ=/dev/stdout \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=true \
            --NON_PF=true \
            --TMP_DIR={config[tempdir]} | \
         bwa mem -M -t {threads} -p {config[reference_fasta]} /dev/stdin | \
         gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
            --ALIGNED_BAM=/dev/stdin \
            --UNMAPPED_BAM={input.bam} \
            --OUTPUT={output.bam} \
            --REFERENCE_SEQUENCE={config[reference_fasta]} \
            --CREATE_INDEX=true \
            --ADD_MATE_CIGAR=true \
            --CLIP_ADAPTERS=false \
            --CLIP_OVERLAPPING_READS=true \
            --INCLUDE_SECONDARY_ALIGNMENTS=true \
            --MAX_INSERTIONS_OR_DELETIONS=-1 \
            --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            --ATTRIBUTES_TO_RETAIN=XS \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Duplicates & merge readgroups
## This step marks duplicate reads
## See: https://gatkforums.broadinstitute.org/gatk/discussion/2799
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markduplicates:
    input:
        lambda wildcards: expand("results/align/bwa/{sample}/{sample}.{rg}.realn.bam", sample=wildcards.aliquot_barcode, rg=manifest.getRGIDs(wildcards.aliquot_barcode))
    output:
        bam = temp("results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam"),
        bai = temp("results/align/markduplicates/{aliquot_barcode}.realn.mdup.bai"),
        metrics = "results/align/markduplicates/{aliquot_barcode}.metrics.txt"
    params:
        max_records = 6000000,
        walltime = lambda wildcards: CLUSTER_META["markduplicates"]["walltime"],
        mem = lambda wildcards: CLUSTER_META["markduplicates"]["mem"]
    #resources:
    # 	mem = lambda wildcards, attempt: CLUSTER_META["markduplicates"]["mem"] if attempt == 1 else CLUSTER_META["markduplicates"]["mem_if_fail"]
    threads:
        CLUSTER_META["markduplicates"]["ppn"]
    log:
        "logs/align/markduplicates/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/markduplicates/{aliquot_barcode}.txt"
    message:
        "Readgroup-specific BAM files are combined into a single BAM. "
        "Potential PCR duplicates are marked.\n"
        "Sample: {wildcards.aliquot_barcode}"
    run:
        multi_input = " ".join(["--INPUT=" + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            --TMP_DIR={config[tempdir]} \
            --MAX_RECORDS_IN_RAM={params.max_records} \
            > {log} 2>&1")

# @sbamin A few notes on IndelRealignment step at annotated link: https://hyp.is/8_20bK-aEeerk1MduBFv6w/gatkforums.broadinstitute.org/gatk/discussion/7847 
# IndelRealinger adds OC:Z tag to realigned reads (GATK Doc # 7156), shifts MAPQ by +10. 
# For high-quality data (>75bp read length, 30x plus, PCR-free lib), and MAPQ > 20, there seems to have minimal impact of +MAPQ10. 
# Perhaps, we should get clarification from GATK team based on GATK Doc # 7847 on pros and cons of keeping IndelRealigner step in the workflow. 
# If there is minimal impact, we should add IndelRealigner. I have not yet seen but OC:Z tag may be required by other indel callers.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Recalibrate base quality scores
## This steps computes a bare recalibration table
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule baserecalibrator:
    input:
        "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam"
    output:
        "results/align/bqsr/{aliquot_barcode}.bqsr.txt"
    params:
        mem = CLUSTER_META["baserecalibrator"]["mem"],
        walltime = CLUSTER_META["baserecalibrator"]["walltime"]
    threads:
        CLUSTER_META["baserecalibrator"]["ppn"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/align/bqsr/{aliquot_barcode}.recal.log"
    benchmark:
        "benchmarks/align/bqsr/{aliquot_barcode}.recal.txt"
    message:
        "Calculating base recalibration scores.\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g BaseRecalibrator \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --known-sites {config[gnomad_vcf]} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Apply BQSR
## This step applies a base recalibration table to an input BAM file
## Formerly "PrintReads"
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule applybqsr:
    input:
        bam = "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam",
        bqsr = "results/align/bqsr/{aliquot_barcode}.bqsr.txt"
    output:
        protected("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    params:
        mem = CLUSTER_META["applybqsr"]["mem"],
        walltime = CLUSTER_META["applybqsr"]["walltime"]
    threads:
        CLUSTER_META["applybqsr"]["ppn"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/align/bqsr/{aliquot_barcode}.apply.log"
    benchmark:
        "benchmarks/align/bqsr/{aliquot_barcode}.apply.txt"
    message:
        "Applying base recalibration scores and generating final BAM file\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ApplyBQSR \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -OQ true \
            -O {output} \
            -bqsr {input.bqsr} \
            --create-output-bam-md5 true \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## WGS Metrics
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calculate coverage metrics for the WGS BAM file. These metrics are shown by MultiQC
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule wgsmetrics:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        "results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt"
    params:
        mem = CLUSTER_META["wgsmetrics"]["mem"],
        walltime = CLUSTER_META["wgsmetrics"]["walltime"]
    threads:
        CLUSTER_META["wgsmetrics"]["ppn"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.log"
    benchmark:
        "benchmarks/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt"
    message:
        "Computing WGS Metrics\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectWgsMetrics \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --USE_FAST_ALGORITHM true \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Validate BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Final check to ensure no errors in final analysis-ready BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Jun 28: Added "|| true" to for exit code zero. Snakemake deletes output if exit code
## != zero, and ValidateSamFile returns exit code 2 (errors) or 3 (warnings) if notable
## events are found
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule validatebam:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        "results/align/validatebam/{aliquot_barcode}.ValidateSamFile.txt"
    params:
        mem = CLUSTER_META["validatebam"]["mem"],
        walltime = CLUSTER_META["validatebam"]["walltime"]
    threads:
        CLUSTER_META["validatebam"]["ppn"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/align/validatebam/{aliquot_barcode}.ValidateSamFile.log"
    benchmark:
        "benchmarks/align/validatebam/{aliquot_barcode}.ValidateSamFile.txt"
    message:
        "Validating BAM file\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ValidateSamFile \
            -I {input} \
            -O {output} \
            -M SUMMARY \
            > {log} 2>&1 \
            || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MultiQC
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MultiQC (and its Click dependency) can have problwms with Python locale settings
## Fix these by specifying the following options in .bash_profile
## export LC_ALL=en_US.UTF-8
## export LANG=en_US.UTF-8
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule multiqc:
    input:
        expand("results/align/validatebam/{sample}.ValidateSamFile.txt", sample = manifest.getSelectedAliquots()),
        expand("results/align/wgsmetrics/{sample}.WgsMetrics.txt", sample = manifest.getSelectedAliquots()),
        lambda wildcards: ["results/align/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample=sample, rg=readgroup)
          for readgroup in manifest.getAllReadgroups()]
    output:
        "results/align/multiqc/multiqc_report.html"
    params:
        dir = "results/align/multiqc",
        mem = CLUSTER_META["multiqc"]["mem"],
        walltime = CLUSTER_META["multiqc"]["walltime"]
    threads:
        CLUSTER_META["multiqc"]["ppn"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/align/multiqc/multiqc.log"
    benchmark:
        "benchmarks/align/multiqc/multiqc.txt"
    message:
        "Running MultiQC"
    shell:
        "multiqc \
            --interactive \
            -o {params.dir} {config[workdir]}/results/align \
            > {log} 2>&1; \
            cp -R {params.dir}/* {config[html_dir]}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on aligned BAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fastqc_bam:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/align/fastqc/{aliquot_barcode}/{aliquot_barcode}.aligned_fastqc.html"
    params:
        dir = "results/align/fastqc/{aliquot_barcode}",
        mem = CLUSTER_META["fastqc_bam"]["mem"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["fastqc_bam"]["ppn"]
    log:
        "logs/align/fastqc-bam/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/fastqc-bam/{aliquot_barcode}.txt"
    message:
        "Running FASTQC\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Get genic coverage given a BAM file
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule gencode_coverage:
#     input:
#         "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
#     output:
#         "results/align/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv"
#     params:
#         mem = CLUSTER_META["gencode_coverage"]["mem"]
#     conda:
#         "../envs/align.yaml"
#     threads:
#         CLUSTER_META["gencode_coverage"]["ppn"]
#     log:
#         "logs/align/gencode-coverage/{aliquot_barcode}.log"
#     benchmark:
#         "benchmarks/align/gencode-coverage/{aliquot_barcode}.txt"
#     message:
#         "Computing coverage using flattened gencode GTF\n"
#         "Sample: {wildcards.aliquot_barcode}"
#     shell:"""
#         module load samtools;
#         module load bedtools;
#         samtools view -q 10 -b {input} | 
#             bedtools coverage -a {config[gencode_gtf_flat]} -b stdin -d | 
#             bedtools groupby -i stdin -g 1,2,3,4,5 -c 7 -o sum | 
#             bedtools groupby -i stdin -g 5 -c 4,6 -o sum,sum
#             1> {output}
#             2> {log}
#         """

## END ##