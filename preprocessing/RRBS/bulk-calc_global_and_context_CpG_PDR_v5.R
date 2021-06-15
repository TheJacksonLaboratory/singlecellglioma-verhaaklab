# bulk-calc_global_and_context_CpG_PDR-v5.R - This script takes the proportion of discordant reads 
# (PDR) per CpG for one sample and calculates global and context-specific epimutation burden 

require(optparse)
require(dplyr)
require(genomation)
require(GenomicRanges)


# Input arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing epimutation output files for one sample."),
  make_option(c("-a","--annotationDir"), action="store", default="/projects/verhaak-lab/scgp/reference/genome_annotations/hg19", type="character",
              help="Directory containing genomic context annotations."),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory.")
  )

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
annotationDir <- args$annotationDir
outdir <- args$outdir


# ##### For troubleshooting #####
# indir <- "~/Documents/scgp/epimutation/results/hypoxia_experiment_RRBS-ge20x_coverage/SCGP-HF-2354-01D-5MC-7201_pe_sorted"
# indir <- "/projects/verhaak-lab/scgp/results/epimutation/hypoxia_experiment_RRBS-ge10x_coverage/SCGP-HF-2354-05D-5MC-0921_pe_sorted"
# annotationDir <- "/projects/verhaak-lab/scgp/reference/genome_annotations/hg19"
# outdir <- "/projects/verhaak-lab/scgp/results/epimutation/hypoxia_experiment_RRBS-ge20x_coverage/SCGP-HF-2354-05D-5MC-0921_pe_sorted"
# sample_cov <- read.delim("~/sumner/verhaak-lab/scgp/results/RRBS/human/bed_graph/SCGP-HF-2354-05D-5MC-0921_pe_sorted.bismark.cov.gz", header = FALSE)
# consensus_CpGs <- readBed("~/Documents/scgp/reference/hypoxia_experiment_RRBS_annotations/HF2354_consensus_10x_coverage_CpGs.bed")
# ###############################


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


# # Load UCSC bed file of DNA sequences with interspersed repeats and low complexity; 
# # Methylation observations overlapping these regions may be erroneous
# rmsk_anno <- readBed(paste0(annotationDir,"/","UCSC_hg19_RepeatMasker-sorted.bed"))

# Load sample bismark coverage file
sample_cov <- read.delim(paste0("/projects/verhaak-lab/scgp/results/RRBS/human/bed_graph","/",
                                gsub("_sorted","",basename(indir)),".bismark.cov.gz"), header = FALSE)
names(sample_cov) <- c("chr","start","end","meth_percentage","count_methylated","count_unmethylated")

# Sort coverage file and remove CpGs aligned to sex chromosomes or scaffolds
sample_cov <- sample_cov[which(sample_cov$chr %in% c(1:22)),]
sample_cov <- sample_cov[order(as.integer(sample_cov$chr),as.integer(sample_cov$start)),]

# Calculate site depth in coverage file as the total number of CpGs at each position
sample_cov$site_depth <- sample_cov$count_methylated + sample_cov$count_unmethylated

# Convert coverage data frame to a GRanges object
sample_cov.ranges <- GRanges(seqnames=paste0("chr",sample_cov$chr),
                             IRanges(start = as.integer(sample_cov$start),
                                     end = as.integer(sample_cov$end)))

# Load sample coordinates for epialleles
epiallele_loci <- readBed(paste0(indir,"/",basename(indir),"_epiallele_CpG_loci.bed"))

# Check if a merged PDR table exists; if not, merge individual chromosome PDR tables
# NOTE: this assumes individual files have the suffix "_individual_CpG_PDR-chrX.txt" 
# and merged file has suffix "_individual_CpG_PDR-merged.txt"
merged_file <- paste0(indir,"/",basename(indir),"_individual_CpG_PDR-merged.txt")

if (!checkmate::testFileExists(merged_file)) {
  # List individual chromosome files
  chr_files <- list.files(path = indir, pattern = "*_individual_CpG_PDR-chr*", full.names = TRUE)
  
  # Stop program execution if any chromosome files are missing
  try(if(length(chr_files) < 22) stop("Error: chromosome file(s) missing"))
  
  ### Order and concatenate chromosome files
  # Extract chromosome number from file name
  chr_num <- as.integer(gsub("chr","",gsub(".txt","",unlist(lapply(strsplit(basename(chr_files),"-"),'[',2)))))
  
  # Order chromosome files by number
  chr_files.ordered <- chr_files[order(chr_num)]
  
  # NOTE: Currently chromosome files are being printed with column names (NEED TO FIX); remove header from files before concatenating them
  for (i in 1:length(chr_files.ordered)) {
    # Check if header has been previously removed; if not remove first line of file
    check <- system(paste0("head -n 1 ",chr_files.ordered[i]), intern = TRUE)
    if (pracma::strcmp(check, "cpg_id\tcpg_site_depth\tPDR")) {
      system(paste0("sed '1d' ",chr_files.ordered[i]," > ",indir,"/tmpfile; mv ",indir,"/tmpfile ",chr_files.ordered[i]))
    }
  }
  
  # Concatenate chromosome files
  system(paste0("cat ",paste(chr_files.ordered, collapse = " "), "> ", indir,"/",basename(indir),"_individual_CpG_PDR-merged.txt"))
  
  # Remove tmpfile
  system("rm tmpfile")
  ###
}

# Load merged CpG PDR file
sample_CpG_PDR <- read.delim(merged_file, header = FALSE)
names(sample_CpG_PDR) <- c("cpg_id","site_depth","PDR")

# Remove any rows with missing data (NOTE: should be no missing data with updated calculation format)
if (any(is.na(sample_CpG_PDR$PDR))) {sample_CpG_PDR <- sample_CpG_PDR[-which(is.na(sample_CpG_PDR$PDR)),]}

# # Load bed file of consensus CpGs (with methylation observations in all samples in cell line)
# if (strsplit(basename(indir),"-")[[1]][3] == "2354") { # HF-2354 cell line
#   consensus_CpGs <- readBed("/projects/verhaak-lab/scgp/reference/hypoxia_experiment_RRBS_annotations/HF2354_consensus_10x_coverage_CpGs.bed")
# }
# if (strsplit(basename(indir),"-")[[1]][3] == "3016") { # HF-3016 cell line
#   consensus_CpGs <- readBed("/projects/verhaak-lab/scgp/reference/hypoxia_experiment_RRBS_annotations/HF3016_consensus_10x_coverage_CpGs.bed")
# }
# 
# # Filter consensus CpGs to autosomes
# seqlevels(consensus_CpGs, pruning.mode="coarse") <- as.character(c(1:22))
# 
# # Format consensus CpGs for consistency with input data
# seqlevels(consensus_CpGs) <- paste0("chr",seqlevels(consensus_CpGs))
###########################


##### Calculate global and context-specific PDR for all input CpGs #####
# Convert CpG site IDs to a GRanges object
sample_CpG_ranges <- GRanges(seqnames=paste0("chr",unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',1))),
                             IRanges(start = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2))),
                                     end = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2)))))


### Additional CpG filtering
# For reference, record the total number of CpGs with methylation observations
CpG_count.total <- nrow(sample_cov)

# # For reference, record the number of CpGs not overlapping 
# # regions with interspersed repeats and low complexity DNA sequences
# CpG_count.nonoverlapping_rmsk <- length(which(!seq(1,length(sample_cov.ranges)) %in% 
#                                                 unique(findOverlaps(sample_cov.ranges,rmsk_anno)@from)))

# For reference, record the number of CpGs with >= 10x coverage
CpG_count.ge10x <- length(which(sample_cov$site_depth >= 10))

# For reference, record the number of CpGs with >= 20x coverage
CpG_count.ge20x <- length(which(sample_cov$site_depth >= 20))

# For reference, record the number of CpGs overlapping epialleles
sample_cov.epiallele_ind <- unique(findOverlaps(sample_cov.ranges,unique(epiallele_loci))@from)
CpG_count.epiallele <- length(sample_cov.epiallele_ind)

# For reference, record the number of epiallele CpGs with >= 10x coverage
CpG_count.epiallele.ge10x <- length(which(sample_cov$site_depth[sample_cov.epiallele_ind] >= 10))

# For reference, record the number of epiallele CpGs with >= 20x coverage
CpG_count.epiallele.ge20x <- length(which(sample_cov$site_depth[sample_cov.epiallele_ind] >= 20))

# For reference, record number of CpGs considered for epimutation analysis
# (epiallele CpGs with >= 10x coverage in cell line)
CpG_count.epimutation <- nrow(sample_CpG_PDR)

# Determine CpGs considered for epimutation analysis with coverage >= 20x
sample_CpG_PDR.ge10x_ind <- which(sample_CpG_PDR$site_depth >= 10)

# For reference, record number of CpGs considered for epimutation analysis with coverage >= 10x 
CpG_count.epimutation.ge10x <- length(which(sample_CpG_PDR$site_depth >= 10))

# Calculate percentiles for site depth
sample_CpG_PDR <- sample_CpG_PDR %>% mutate(depth_percentile = ntile(site_depth,100))

# Determine threshold for site depth (= 100 * 95th percentile of coverage)
max_site_depth <- median(sample_CpG_PDR$site_depth[which(sample_CpG_PDR$depth_percentile == 95)]) * 100

# Determine CpGs considered for epimutation analysis with coverage <= max_site_depth
sample_CpG_PDR.coverage_filtering_ind <- which(sample_CpG_PDR$site_depth <= max_site_depth)

# Record number of CpGs considered for epimutation analysis with coverage <= max_site_depth
CpG_count.epimutation.DP_filtered <- length(sample_CpG_PDR.coverage_filtering_ind)

# Determine CpGs considered for epimutation analysis with coverage >= 20x
sample_CpG_PDR.ge20x_ind <- which(sample_CpG_PDR$site_depth >= 20)

# Record number of CpGs considered for epimutation analysis with coverage >= 20x
CpG_count.epimutation.ge20x <- length(sample_CpG_PDR.ge20x_ind)

# # Determine CpGs considered for epimutation analysis not overlapping 
# # regions with interspersed repeats and low complexity DNA sequences
# sample_CpG_PDR.nonoverlapping_rmsk_ind <- which(!seq(1,length(sample_CpG_ranges)) %in% 
#                                                   unique(findOverlaps(sample_CpG_ranges,rmsk_anno)@from))
# 
# # Record number of CpGs considered for epimutation analysis not overlapping 
# # regions with interspersed repeats and low complexity DNA sequences
# CpG_count.epimutation.nonoverlapping_rmsk <- length(sample_CpG_PDR.nonoverlapping_rmsk_ind)

# Filter CpGs considered for epimutation analysis:
#   coverage >= 20x and <= 100 times the 95th percentile of coverage
#   not overlapping regions with interspersed repeats and low complexity DNA sequences
sample_CpG_PDR <- sample_CpG_PDR[Reduce(intersect, list(sample_CpG_PDR.coverage_filtering_ind,
                                                        sample_CpG_PDR.ge20x_ind)),]

# Write table to file
write.table(sample_CpG_PDR[,c(1:3)], file = paste0(outdir,"/",basename(indir),"_individual_CpG_PDR-merged.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

sample_CpG_ranges <- sample_CpG_ranges[Reduce(intersect, list(sample_CpG_PDR.coverage_filtering_ind,
                                                              sample_CpG_PDR.ge20x_ind))]

# Aggregate metrics for various CpG counts
CpG_count <- data.frame(sample=basename(indir),
                        totalCpGs = CpG_count.total,
                        totalCpGs.ge10x = CpG_count.ge10x,
                        totalCpGs.ge20x = CpG_count.ge20x,
                        epialleleCpGs = CpG_count.epiallele,
                        epialleleCpGs.ge10x = CpG_count.epiallele.ge10x,
                        epialleleCpGs.ge20x = CpG_count.epiallele.ge20x,
                        epimutationCpGs = CpG_count.epimutation,
                        epimutationCpGs.le_high_coverage = CpG_count.epimutation.DP_filtered,
                        epimutationCpGs.ge10x = CpG_count.epimutation.ge10x,
                        epimutationCpGs.ge20x = CpG_count.epimutation.ge20x,
                        epimutationCpGs.filtered = nrow(sample_CpG_PDR))

# Write table to file
write.table(CpG_count, file = paste0(outdir,"/",basename(indir),"_CpG_filtering_counts.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
###


# Initialize genomic context output table (assuming 19 genomic annotations)
# for calculating PDR by taking unweighted or weighted average of CpG PDRs
sample_context_PDR <- data.frame(matrix(nrow = 1,ncol = 40))
sample_context_PDR[,1] <- unlist(strsplit(basename(merged_file),"_individual_CpG_PDR-merged.txt"))
names(sample_context_PDR) <- c("sample","global","exon_epiallele_count","exon_PDR","gene_body_epiallele_count",
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

# Calculate global PDR by taking average of all CpG PDRs
sample_context_PDR$global <- weighted.mean(sample_CpG_PDR$PDR)

# Calculate global PDR by taking weighted average of all CpG PDRs (using site depth for weights)
sample_context_PDR.wm$global <- weighted.mean(sample_CpG_PDR$PDR,sample_CpG_PDR$site_depth/sum(sample_CpG_PDR$site_depth))


### Loop through target genomic contexts and compute context-specific PDR
### as the unweighted or weighted average of CpGs overlapping genomic context
# Partition CpGs by genomic context
eCpG_loci.context <- vector(mode = "list", length = length(annos))
names(eCpG_loci.context) <- names(annos)
for (i in 1:length(annos)) {
  eCpG_loci.context[[i]] <- unique(findOverlaps(sample_CpG_ranges,annos[[i]])@from)
}

# For each annotation, record the number of overlapping epialleles
# and the unweighted and weighted average of overlapping CpG PDRs
eCpG_PDR <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
eCpG_PDR.wm <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
context_fraction_ind <- seq(1,length(annos)*2,2)
context_PDR_ind <- seq(2,length(annos)*2,2)

for (i in 1:length(annos)) {
  # Calculate number of reads overlapping genomic contexts for mean PDR table
  eCpG_PDR[context_fraction_ind[i]] <- length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges[eCpG_loci.context[[i]]])@from))
  names(eCpG_PDR)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_epiallele_count")

  # Calculate number of reads overlapping genomic contexts for weighted mean PDR table
  eCpG_PDR.wm[context_fraction_ind[i]] <- length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges[eCpG_loci.context[[i]]])@from))
  names(eCpG_PDR.wm)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_epiallele_count")

  # Calculate context-specific PDR as mean of PDR values for eCpGs overlapping target genomic context
  eCpG_PDR[context_PDR_ind[i]] <- mean(sample_CpG_PDR$PDR[eCpG_loci.context[[i]]])
  names(eCpG_PDR)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")

  # Calculate context-specific PDR as weighted mean of PDR values for eCpGs overlapping target genomic context
  eCpG_PDR.wm[context_PDR_ind[i]] <- weighted.mean(sample_CpG_PDR$PDR[eCpG_loci.context[[i]]],
                                                sample_CpG_PDR$site_depth[eCpG_loci.context[[i]]]/sum(sample_CpG_PDR$site_depth[eCpG_loci.context[[i]]]))
  names(eCpG_PDR.wm)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")
}

# Fill sample output table with cell data
sample_context_PDR[-c(1:2)] <- eCpG_PDR
sample_context_PDR.wm[-c(1:2)] <- eCpG_PDR.wm

# Save output tables to file
write.table(sample_context_PDR, file = paste0(outdir,"/",basename(indir),"_ge20x_global_and_context_set1and2-specific_unweighted_mean_PDR.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

write.table(sample_context_PDR.wm, file = paste0(outdir,"/",basename(indir),"_ge20x_global_and_context_set1and2-specific_weighted_mean_PDR.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
#######################################################################


# ##### Calculate global and context-specific PDR including only CpGs with observations in all cell-line samples #####
# # Filter input data to consensus CpGs
# consensus_ind <- findOverlaps(sample_CpG_ranges,consensus_CpGs)@from
# sample_CpG_ranges.consensus <- sample_CpG_ranges[consensus_ind]
# sample_CpG_PDR.consensus <- sample_CpG_PDR[consensus_ind,]
# 
# # Initialize genomic context output table (assuming 9 genomic annotations)
# # for calculating PDR by taking average of consensus CpG PDRs
# sample_context_PDR.consensus <- data.frame(matrix(nrow = 1,ncol = 18))
# sample_context_PDR.consensus[,1] <- unlist(strsplit(basename(merged_file),"_individual_CpG_PDR-merged.txt"))
# names(sample_context_PDR.consensus) <- c("sample","global","exon_fraction","exon_PDR","gene_body_fraction",
#                                          "gene_body_PDR","intergenic_fraction","intergenic_PDR","intron_fraction",
#                                          "intron_PDR","enhancer_fraction","enhancer_PDR","promoter_fraction",
#                                          "promoter_PDR","cgi_fraction","cgi_PDR","dnaseI_fraction","dnaseI_PDR")
# 
# # Initialize genomic context output table for calculating PDR by taking weighted average
# sample_context_PDR.wm.consensus <- sample_context_PDR.consensus
# 
# # Calculate global PDR by taking average of consensus CpG PDRs (using site depth for weights)
# sample_context_PDR.consensus$global <- mean(sample_CpG_PDR.consensus$PDR)
# 
# # Calculate global PDR by taking weighted average of consensus CpG PDRs (using site depth for weights)
# sample_context_PDR.wm.consensus$global <- weighted.mean(sample_CpG_PDR.consensus$PDR,sample_CpG_PDR.consensus$site_depth/sum(sample_CpG_PDR.consensus$site_depth))
# 
# 
# ### Loop through target genomic contexts and compute context-specific PDR
# ### as the unweighted or weighted average of CpGs overlapping genomic context
# # Partition CpGs by genomic context
# eCpG.consensus_loci.context <- vector(mode = "list", length = length(annos))
# names(eCpG.consensus_loci.context) <- names(annos)
# for (i in 1:length(annos)) {
#   eCpG.consensus_loci.context[[i]] <- unique(findOverlaps(sample_CpG_ranges.consensus,annos[[i]])@from)
# }
# 
# # For each annotation, record the fraction of overlapping epialleles
# # and the unweighted and weighted average of overlapping CpG PDRs
# eCpG.consensus_PDR <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
# eCpG.consensus_PDR.wm <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
# context_fraction_ind <- seq(1,length(annos)*2,2)
# context_PDR_ind <- seq(2,length(annos)*2,2)
# 
# for (i in 1:length(annos)) {
#   # Calculate fraction of reads with consensus CpGs overlapping genomic contexts for mean PDR table
#   eCpG.consensus_PDR[context_fraction_ind[i]] <- length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges.consensus[eCpG.consensus_loci.context[[i]]])@from))/
#     length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges.consensus)@from))
#   names(eCpG.consensus_PDR)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_fraction")
#   
#   # Calculate fraction of reads with consensus CpGs overlapping genomic contexts for weighted mean PDR table
#   eCpG.consensus_PDR.wm[context_fraction_ind[i]] <- length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges.consensus[eCpG.consensus_loci.context[[i]]])@from))/
#     length(unique(findOverlaps(epiallele_loci,sample_CpG_ranges.consensus)@from))
#   names(eCpG.consensus_PDR.wm)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_fraction")
#   
#   # Calculate context-specific PDR as mean of PDR values for eCpGs overlapping target genomic context
#   eCpG.consensus_PDR[context_PDR_ind[i]] <- mean(sample_CpG_PDR.consensus$PDR[eCpG.consensus_loci.context[[i]]])
#   names(eCpG.consensus_PDR)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")
#   
#   # Calculate context-specific PDR as weighted mean of PDR values for eCpGs overlapping target genomic context
#   eCpG.consensus_PDR.wm[context_PDR_ind[i]] <- weighted.mean(sample_CpG_PDR.consensus$PDR[eCpG.consensus_loci.context[[i]]],
#                                                 sample_CpG_PDR.consensus$site_depth[eCpG.consensus_loci.context[[i]]]/sum(sample_CpG_PDR.consensus$site_depth[eCpG.consensus_loci.context[[i]]]))
#   names(eCpG.consensus_PDR.wm)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")
# }
# 
# # Fill sample output table with cell data
# sample_context_PDR.consensus[-c(1:2)] <- eCpG.consensus_PDR
# sample_context_PDR.wm.consensus[-c(1:2)] <- eCpG.consensus_PDR.wm
# 
# # Save output tables to file
# write.table(sample_context_PDR.consensus, file = paste0(outdir,"/",basename(indir),"_ge20x_global_and_context-specific_unweighted_mean_PDR-consensus_CpGs.txt"), 
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# 
# write.table(sample_context_PDR.wm.consensus, file = paste0(outdir,"/",basename(indir),"_ge20x_global_and_context-specific_weighted_mean_PDR-consensus_CpGs.txt"), 
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# ####################################################################################################################
