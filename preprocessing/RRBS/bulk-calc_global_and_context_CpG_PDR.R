# bulk-calc_global_and_context_CpG_PDR.R - This script takes the proportion of discordant reads 
# (PDR) per CpG for one sample and calculates global and context-specific DNAme disorder

require(optparse)
require(dplyr)
require(genomation)
require(GenomicRanges)


# Input arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing DNAme disorder output files for one sample."),
  make_option(c("-a","--annotationDir"), action="store", default="/projects/verhaak-lab/scgp/reference/genome_annotations/hg19", type="character",
              help="Directory containing genomic context annotations."),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory.")
  )

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
annotationDir <- args$annotationDir
outdir <- args$outdir


##### Input arguments #####
# Create output directory if not existing
if(!dir.exists(outdir)) {dir.create(outdir)}

# Change directory to working directory
setwd(outdir)

# Load annotations
exon_anno <- readBed(paste0(annotationDir,"/","Ensembl_hg19_exon-sorted.bed"))
gene_body_anno <- readBed(paste0(annotationDir,"/","Ensembl_hg19_gene_body-sorted.bed"))
intergenic_anno <- readBed(paste0(annotationDir,"/","Ensembl_hg19_intergenic-sorted.bed"))
intron_anno <- readBed(paste0(annotationDir,"/","Ensembl_hg19_intron-sorted.bed"))
enhancer_anno <- readBed(paste0(annotationDir,"/","FANTOM5_hg19_human_permissive_enhancers_phase_1_and_2-sorted.bed"))
promoter_anno <- readBed(paste0(annotationDir,"/","FANTOM5_hg19_gene_matched_TSS.padded_1500u_500d.bed"))
cgi_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_CGI-sorted.bed"))
cgi_shore_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_CGI_shore-sorted.bed"))
dnaseI_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_DNaseI_hypersensitive_sites-sorted.bed"))
alu_repeat_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_RepeatMasker-Alu_family-sorted.bed"))
l1_repeat_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_RepeatMasker-L1_family-sorted.bed"))
l2_repeat_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_RepeatMasker-L2_family-sorted.bed"))
mir_repeat_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_RepeatMasker-MIR_family-sorted.bed"))
gliobla_ctcf_anno <- readBed(paste0(annotationDir,"/","ENCFF001USC_Gliobla_CTCF-sorted.bed"))
h1hesc_ctcf_1_anno <- readBed(paste0(annotationDir,"/","ENCFF001UBA_H1HESC_CTCF_1-sorted.bed"))
h1hesc_ctcf_2_anno <- readBed(paste0(annotationDir,"/","ENCFF001UBB_H1HESC_CTCF_2-sorted.bed"))
h1hesc_ezh2_anno <- readBed(paste0(annotationDir,"/","ENCFF001SUU_H1HESC_EZH2-sorted.bed"))
nha_ctcf_anno <- readBed(paste0(annotationDir,"/","ENCFF001TAT_NHA_CTCF-sorted.bed"))
nha_ezh2_anno <- readBed(paste0(annotationDir,"/","ENCFF001TAU_NHA_EZH2-sorted.bed"))

annos <- list(
  exon=exon_anno,
  gene_body=gene_body_anno,
  intergenic=intergenic_anno,
  intron=intron_anno,
  enhancer=enhancer_anno,
  promoter=promoter_anno,
  cgi=cgi_anno,
  cgi_shore=cgi_shore_anno,
  dnaseI=dnaseI_anno,
  alu_repeat=alu_repeat_anno,
  l1_repeat=l1_repeat_anno,
  l2_repeat=l2_repeat_anno,
  mir_repeat=mir_repeat_anno,
  gliobla_ctcf=gliobla_ctcf_anno,
  h1hesc_ctcf_1=h1hesc_ctcf_1_anno,
  h1hesc_ctcf_2=h1hesc_ctcf_2_anno,
  h1hesc_ezh2=h1hesc_ezh2_anno,
  nha_ctcf=nha_ctcf_anno,
  nha_ezh2=nha_ezh2_anno
)

# Load merged CpG PDR file
merged_file <- paste0(indir,"/",basename(indir),"_individual_CpG_PDR-merged.txt")

sample_CpG_PDR <- read.delim(merged_file, header = FALSE)
names(sample_CpG_PDR) <- c("cpg_id","site_depth","PDR")
###########################


##### Calculate global and context-specific PDR for all input CpGs #####
# Convert CpG site IDs to a GRanges object
sample_CpG_ranges <- GRanges(seqnames=paste0("chr",unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',1))),
                             IRanges(start = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2))),
                                     end = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2)))))

# Initialize genomic context output table (assuming 19 genomic annotations)
# for calculating PDR by taking weighted average of per-CpG PDRs
sample_context_PDR <- data.frame(matrix(nrow = 1,ncol = 40))
sample_context_PDR[,1] <- unlist(strsplit(basename(merged_file),"_individual_CpG_PDR-merged.txt"))
names(sample_context_PDR) <- c("sample_barcode","global_PDR","exon_epiallele_count","exon_PDR","gene_body_epiallele_count",
                               "gene_body_PDR","intergenic_epiallele_count","intergenic_PDR","intron_epiallele_count",
                               "intron_PDR","enhancer_epiallele_count","enhancer_PDR","promoter_epiallele_count",
                               "promoter_PDR","cgi_epiallele_count","cgi_PDR","cgi_shore_epiallele_count","cgi_shore_PDR",
                               "dnaseI_epiallele_count","dnaseI_PDR","alu_repeat_epiallele_count","alu_repeat_PDR",
                               "l1_repeat_epiallele_count","l1_repeat_PDR","l2_repeat_epiallele_count","l2_repeat_PDR","mir_repeat_epiallele_count",
                               "mir_repeat_PDR","gliobla_ctcf_epiallele_count","gliobla_ctcf_PDR","h1hesc_ctcf_1_epiallele_count",
                               "h1hesc_ctcf_1_PDR","h1hesc_ctcf_2_epiallele_count","h1hesc_ctcf_2_PDR","h1hesc_ezh2_epiallele_count",
                               "h1hesc_ezh2_PDR","nha_ctcf_epiallele_count","nha_ctcf_PDR","nha_ezh2_epiallele_count","nha_ezh2_PDR")

# Initialize genomic context output table for calculating PDR by taking weighted average
sample_context_PDR.wm <- sample_context_PDR

# Calculate global PDR by taking weighted average of all CpG PDRs (using site depth for weights)
sample_context_PDR.wm$global_PDR <- weighted.mean(sample_CpG_PDR$PDR,sample_CpG_PDR$site_depth/sum(sample_CpG_PDR$site_depth))


### Loop through target genomic contexts and compute context-specific PDR
### as the unweighted or weighted average of CpGs overlapping genomic context
# Partition CpGs by genomic context
eCpG_loci.context <- vector(mode = "list", length = length(annos))
names(eCpG_loci.context) <- names(annos)
for (i in 1:length(annos)) {
  eCpG_loci.context[[i]] <- unique(findOverlaps(sample_CpG_ranges,annos[[i]])@from)
}

# For each annotation, record the number of overlapping epialleles
# and the weighted average of overlapping CpG PDRs
eCpG_PDR.wm <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
context_fraction_ind <- seq(1,length(annos)*2,2)
context_PDR_ind <- seq(2,length(annos)*2,2)

for (i in 1:length(annos)) {
  # Calculate number of reads overlapping genomic contexts for weighted mean PDR table
  eCpG_PDR.wm[context_fraction_ind[i]] <- length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges[eCpG_loci.context[[i]]])@from))
  names(eCpG_PDR.wm)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_epiallele_count")

  # Calculate context-specific PDR as weighted mean of PDR values for eCpGs overlapping target genomic context
  eCpG_PDR.wm[context_PDR_ind[i]] <- weighted.mean(sample_CpG_PDR$PDR[eCpG_loci.context[[i]]],
                                                sample_CpG_PDR$site_depth[eCpG_loci.context[[i]]]/sum(sample_CpG_PDR$site_depth[eCpG_loci.context[[i]]]))
  names(eCpG_PDR.wm)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")
}

# Fill sample output table with cell data
sample_context_PDR.wm[-c(1:2)] <- eCpG_PDR.wm

# Save output table to file
write.table(sample_context_PDR.wm, file = paste0(outdir,"/",basename(indir),"_global_and_context-specific_PDR.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
#######################################################################

