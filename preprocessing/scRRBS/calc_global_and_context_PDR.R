require(optparse)
require(genomation)
require(GenomicRanges)

# Input arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing cell output files for one sample."),
  make_option(c("-a","--annotationDir"), action="store", default="/Users/anderk/Documents/SCGP/reference/genome_annotations", type="character",
              help="Directory containing genomic context annotations."),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory.")
  )

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
annotationDir <- args$annotationDir
outdir <- args$outdir


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



# Extract the methylation status and loci files
status_files <- list.files(path = indir, pattern = "_epiallele_CpG_status.txt")
loci_files <- list.files(path = indir, pattern = "_epiallele_CpG_loci.bed")

# Initialize output table (assuming 9 genomic annotations)
# sample_context_PDR <- data.frame(matrix(nrow = length(status_files),ncol = 23))
sample_context_PDR <- data.frame(matrix(nrow = length(status_files),ncol = 40))
sample_context_PDR[,1] <- unlist(strsplit(status_files,"_epiallele_CpG_status.txt"))
names(sample_context_PDR) <- c("cell_barcode","PDR","exon_epiallele_count","exon_PDR","gene_body_epiallele_count","gene_body_PDR",
                               "intergenic_epiallele_count","intergenic_PDR","intron_epiallele_count","intron_PDR",
                               "enhancer_epiallele_count","enhancer_PDR","promoter_epiallele_count","promoter_PDR",
                               "cgi_epiallele_count","cgi_PDR","dnaseI_epiallele_count","dnaseI_PDR","cgi_shore_epiallele_count",
                               "cgi_shore_PDR","alu_repeat_epiallele_count","alu_repeat_PDR","l1_repeat_epiallele_count",
                               "l1_repeat_PDR","l2_repeat_epiallele_count","l2_repeat_PDR","mir_repeat_epiallele_count",
                               "mir_repeat_PDR","gliobla_ctcf_epiallele_count","gliobla_ctcf_PDR","h1hesc_ctcf_1_epiallele_count",
                               "h1hesc_ctcf_1_PDR","h1hesc_ctcf_2_epiallele_count","h1hesc_ctcf_2_PDR","h1hesc_ezh2_epiallele_count",
                               "h1hesc_ezh2_PDR","nha_ctcf_epiallele_count","nha_ctcf_PDR","nha_ezh2_epiallele_count","nha_ezh2_PDR")


# Loop through methylation status and loci files
# and compute context-specific PDR for each cell
for (file in 1:length(status_files)) {
  # Load sample data
  eCpG_status <- read.delim(paste0(indir,"/",status_files[file]), header = FALSE)
  eCpG_loci <- readBed(paste0(indir,"/",loci_files[file]))
  
  # Partition reads by genomic context
  eCpG_loci.context <- vector(mode = "list", length = length(annos))
  for (i in 1:length(annos)) {
    eCpG_loci.context[[i]] <- unique(findOverlaps(eCpG_loci,annos[[i]])@from)
  }
  
  # Extract methylation calls from eCpG_status, and 
  # identify concordant and discordant epialleles
  eCpG_status.calls <- eCpG_status[, seq(3, dim(eCpG_status)[2]-1, 2)]
  eCpG_status.eSum <- rowSums(eCpG_status.calls, na.rm = TRUE)
  
  # Calculate the proportion of discordant epialleles for entire dataset
  sample_context_PDR[file,2] <- length(which(eCpG_status.eSum != 0 & 
                                               eCpG_status.eSum != 
                                               apply(eCpG_status.calls, 1, 
                                                     function(x) length(which(!is.na(x))))))/length(eCpG_status.eSum)
  
  # Count the proportion of discordant epialleles per genomic context
  eCpG_PDR <- data.frame(matrix(nrow = 1,ncol = length(annos)*2))
  context_fraction_ind <- seq(1,length(annos)*2,2)
  context_PDR_ind <- seq(2,length(annos)*2,2)
  for (i in 1:length(annos)) {
    # For each annotation, record the total number of epialleles
    # and the percentage of epialleles with mixed methylation status
    eCpG_PDR[context_fraction_ind[i]] <- length(eCpG_status.eSum[eCpG_loci.context[[i]]])
    names(eCpG_PDR)[context_fraction_ind[i]] <- paste0(names(annos[i]),"_epiallele_count")
    eCpG_PDR[context_PDR_ind[i]] <- length(which(eCpG_status.eSum[eCpG_loci.context[[i]]] != 0 & 
                                                 eCpG_status.eSum[eCpG_loci.context[[i]]] != 
                                                   apply(eCpG_status.calls[eCpG_loci.context[[i]],], 1, 
                                                         function(x) length(which(!is.na(x))))))/length(eCpG_status.eSum[eCpG_loci.context[[i]]])
    names(eCpG_PDR)[context_PDR_ind[i]] <- paste0(names(annos[i]),"_PDR")
  }

  # Fill sample output table with cell data
  sample_context_PDR[file,-c(1:2)] <- eCpG_PDR
}

# Extract case barcode from cell barcode and add as separate column
sample_context_PDR <- data.frame(cell_barcode=sample_context_PDR[,1],
                                 case_barcode=gsub("-","",substr(sample_context_PDR[,1],6,11)),
                                 sample_context_PDR[,-1])

# Save output table to file
write.table(sample_context_PDR, file = paste0(outdir,"/",basename(indir),"_global_and_context-specific_PDR.txt"), 
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

