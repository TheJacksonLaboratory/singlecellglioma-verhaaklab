# bulk-filter_CpG_PDR.R - This script takes the proportion of discordant reads 
# (PDR) per CpG for one sample and performs data filtering

require(optparse)
require(dplyr)
require(genomation)
require(GenomicRanges)


# Input arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing epiallele output files for one sample."),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory.")
  )

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
outdir <- args$outdir


##### Input arguments #####
# Create output directory if not existing
if(!dir.exists(outdir)) {dir.create(outdir)}

# Change directory to working directory
setwd(outdir)


# Load sample bismark coverage file
sample_cov <- read.delim(paste0(BISMARK_COVERAGE_DIR,"/",
                                basename(indir),".bismark.cov.gz"), header = FALSE)
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
  
  # Concatenate chromosome files
  system(paste0("cat ",paste(chr_files.ordered, collapse = " "), "> ", indir,"/",basename(indir),"_individual_CpG_PDR-merged.txt"))
  
  # Remove tmpfile
  system("rm tmpfile")
  ###
}

# Load merged CpG PDR file
sample_CpG_PDR <- read.delim(merged_file, header = FALSE)
names(sample_CpG_PDR) <- c("cpg_id","site_depth","PDR")
###########################


##### Calculate global and context-specific PDR for all input CpGs #####
# Convert CpG site IDs to a GRanges object
sample_CpG_ranges <- GRanges(seqnames=paste0("chr",unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',1))),
                             IRanges(start = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2))),
                                     end = as.integer(unlist(lapply(strsplit(as.character(sample_CpG_PDR$cpg_id),"_"),'[',2)))))


### Additional CpG filtering
# For reference, record the total number of CpGs with methylation observations
CpG_count.total <- nrow(sample_cov)

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