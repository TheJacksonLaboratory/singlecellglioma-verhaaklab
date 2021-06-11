## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Kallisto pipeline for quantifying transcript expression in TPM
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run fastp to clean up RNAseq fastq:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule fastp:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        R1 = temp("results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq"),
        R2 = temp("results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq"),
        json = "results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.json",
        html = "results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.html"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/RNAseq/fastp/{aliquot_barcode}.{readgroup}.log"
    message:
        "Cleaning up fastq with fastp \n"
        "Sample: {wildcards.aliquot_barcode}.{wildcards.readgroup}"
    shell:
    	"(fastp \
		-i {input.R1} \
		-I {input.R2}  \
		-o {output.R1} \
		-O {output.R2}  \
		-j {output.json} \
		-h {output.html}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run kallisto to quantify transcripts:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule kallisto:
    input:
        lambda wildcards: expand("results/kallisto/fastp/{{aliquot_barcode}}/{{aliquot_barcode}}.{readgroup}.{read}.fastp.fastq", readgroup=manifest.getRGIDs(wildcards.aliquot_barcode),read=['R1','R2'])
    output:
        tsv = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.tsv",
        h5 = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.h5",
        json = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/run_info.json",
        fusion = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/fusion.txt"
    params:
    	output_dir = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/"
    conda:
        "../envs/kallisto.yaml"
    log:
        "logs/RNAseq/kallisto/aliquot/{aliquot_barcode}.log"
    message:
        "Quantifying transcript with kallisto \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"(kallisto quant \
		-i {config[kallisto_idx]} \
		--fusion \
		-o {params.output_dir} \
		{input}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge tpms together:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule mergetpm:
    input:
        expand("results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        protected("results/kallisto/kallisto/final/transcript_tpms_all_samples.tsv")
    log:
        "logs/RNAseq/kallisto/mergetpm.log"
    message:
        "Merging aliquot TPMs into one file and uploading to database"
    shell:
    	"""
    	set +o pipefail; 
    	cat {input} | head -1 | sed 's/^/aliquot_barcode\t/' > {output}    	
    	for f in {input}
    	do
    		al=$(echo $f | cut -c 35-63)				#Get aliquot barcode from file name
    		sed "s/^/$al\t/g" $f | tail -n+2 -q >> {output}
    	done 
    	
    	Rscript R/snakemake/tpm2db.R {output}
    	2>{log}
    	"""

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create a tpm matrix to store locally:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule tpm_matrix:
    input:
        expand("results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        transcript = protected("results/kallisto/kallisto/final/transcript_tpm_matrix_all_samples.tsv"),
        gene = protected("results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv")
    params:
    	dir_str = "results\/kallisto\/kallisto\/aliquot\/",
    	head_tmp = "results/kallisto/kallisto/final/header.tsv",
    	mat_tmp = "results/kallisto/kallisto/final/transcript_tpm_matrix_all_samples.tsv2"
    log:
        "logs/RNAseq/post/tpm_matrix.log"
    message:
        "Merging aliquot TPMs into a transcript and gene-level matrix"
    shell:
    	"""
		num=$(ls -f1 {input} | wc -l)
		upper=$(echo "$((5 * $num))")
		myseq=$(seq 5 5 $upper | sed 's/^\|$//g' | paste -sd,)
		myseq=$(echo "1,2,"$myseq)
		paste {input} | cut -f $myseq > {output.transcript}

		ls -f1 {input} | sed 's/{params.dir_str}//g' | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){{print "\t$1"}}' | perl -ne 'print "target_id\tlength$_\n"' > {params.head_tmp}
		cat {params.head_tmp} {output.transcript} | grep -v "tpm" > {params.mat_tmp}
		mv {params.mat_tmp} {output.transcript}
		rm -f {params.head_tmp} 
		
		Rscript R/snakemake/transcript2gene.R {output.transcript} {config[ensembl_transcript_mapping]} {output.gene}
		2>{log}
    	"""	
			
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create a tpm gct file for certain applications:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule tpm_gct:
    input:
        "results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
    output:
        "results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.gct"
    log:
        "logs/RNAseq/post/tpm_gct.log"
    message:
        "Creating gct file"
    shell:
    	"""
		Rscript R/snakemake/gct_create.R {input} {output}
		2>{log}
    	"""			

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run transcriptional subtype classifier and upload to db:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule transcript_class:
    input:
        "results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.gct"
    output:
        protected("results/kallisto/kallisto/final/p_result_gene_tpm_matrix_all_samples.gct.txt")
    conda:
        "../envs/transcript_class.yaml"
    log:
        "logs/RNAseq/post/transcript_class.log"
    message:
        "Running transcriptional subtype classifier and uploading to db"
    shell:
    	"""
		Rscript R/snakemake/transcriptclass2db.R {input} {output}
		2>{log}
    	"""				

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run pizzly to call transcript fusions:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule pizzly:
    input:
        kfus = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/fusion.txt",
        h5 = "results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.h5"
    output:
        cache = "results/kallisto/pizzly/{aliquot_barcode}/index.cache.txt",
        fusions = "results/kallisto/pizzly/{aliquot_barcode}/{aliquot_barcode}.fusions.fasta",
        json = "results/kallisto/pizzly/{aliquot_barcode}/{aliquot_barcode}.json",
        tsv = "results/kallisto/pizzly/{aliquot_barcode}/{aliquot_barcode}.fusions.tsv"
    params:
    	output = "results/kallisto/pizzly/{aliquot_barcode}/{aliquot_barcode}",
    	align_score = 2
    conda:
        "../envs/pizzly.yaml"
    log:
        "logs/RNAseq/pizzly/{aliquot_barcode}.log"
    message:
        "Calling fusions with pizzly \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
    	fr=$(pizzly_get_fragment_length.py {input.h5})
    	echo $fr
    	
		pizzly -k 31 \
		--gtf {config[kallisto_gtf]} \
		--cache {output.cache} \
		--align-score {params.align_score} \
        --insert-size $fr \
        --fasta {config[kallisto_ref]} \
        --output {params.output} \
        {input.kfus}
 		
 		pizzly_flatten_json.py {output.json} > {output.tsv}
 		2>{log}
	   	"""					

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge fusions together into a single table
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule mergefusions:
    input:
        expand("results/kallisto/pizzly/{aliquot_barcode}/{aliquot_barcode}.fusions.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        protected("results/kallisto/pizzly/final/fusions_all_samples.tsv")
    log:
        "logs/RNAseq/pizzly/mergefusions.log"
    message:
        "Merging all fusions into one file and uploading to database"
    shell:
    	"""
    	set +o pipefail; 
    	cat {input} | head -1 | sed "s/^/aliquot_barcode\t/" > {output}    	
    	for f in {input}
    	do
    		al=$(echo $f | cut -c 25-53)				#Get aliquot barcode from file name
    		sed "s/^/$al\t/g" $f | tail -n+2 -q >> {output}
    	done 
    	
    	Rscript R/snakemake/pizzly2db.R {output}
    	2>{log}
    	"""
