# bulk-calc_chr_CpG_PDR.R - This script takes input epiallele data for one chromosome of one sample 
# and calculates the proportion of discordant reads (PDR) per CpG

require(optparse)
require(genomation)
require(GenomicRanges)
require(tictoc)


# Input arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing chromosome-specific input DNAme disorder calculation files."),
  make_option(c("-n","--chr"), action="store", default=NULL, type="character",
              help="Chromosome for which data will be processed."),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory.")
  )

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
chr <- args$chr
outdir <- args$outdir


print(paste0("Analyzing ",strsplit(indir,"/")[[1]][grep("SCGP-HF",strsplit(indir,"/")[[1]])]," chr",chr," at ",date()))


##### Input arguments #####
# Identify the methylation and CpG site depth files for given chromosme
eCpG_depth_file <- list.files(path = indir, pattern = paste0("_eCpG_depth-chr",chr,".txt"), full.names = TRUE)
eCpG_status.eSum_file <- list.files(path = indir, pattern = paste0("_eCpG_status.eSum-chr",chr,".txt"), full.names = TRUE)
eCpG_status.calls_file <- list.files(path = indir, pattern = paste0("_eCpG_status.calls-chr",chr,".txt"), full.names = TRUE)
eCpG_status.siteID_file <- list.files(path = indir, pattern = paste0("_eCpG_status.siteID-chr",chr,".txt"), full.names = TRUE)

# Load CpG site depth table
eCpG_depth <- read.delim(eCpG_depth_file, header = FALSE)

# Load sum of methylation calls per read
eCpG_status.eSum <- as.numeric(read.delim(eCpG_status.eSum_file, header = FALSE)$V1)

# Load methylation call matrix
eCpG_status.calls <- read.delim(eCpG_status.calls_file, header = FALSE)

# Load CpG ID matrix
eCpG_status.siteID <- as.matrix(read.delim(eCpG_status.siteID_file, header = FALSE))
###########################


##### Determine the number of reads overlapping each CpG and calculate PDR #####
# Calculate the sum of methylation calls per read assuming the read is concordantly methylated
eCpG_status.eNum <- apply(eCpG_status.calls, 1, function(x) length(which(!is.na(x))))

# Flag read as discordant if eCpG_status.eSum != 0 or eCpG_status.eNum
eCpG_status.discordance <- ifelse(eCpG_status.eSum != 0 & eCpG_status.eSum != eCpG_status.eNum,TRUE,FALSE)

# Initialize output table of CpG-level PDR by copying eCpG_depth
sample_CpG_PDR <- eCpG_depth
names(sample_CpG_PDR) <- c("cpg_id", "cpg_site_depth")
sample_CpG_PDR$PDR <- as.numeric(NA)


print(paste0("Calculating per-CpG PDR for chr ",chr," at ",date()))

### Calculate PDR for each CpG
# For each CpG, determine the indices for all overlapping reads
ind <- lapply(eCpG_depth$V1, function(x) which(eCpG_status.siteID == as.character(x), arr.ind = TRUE))

# For each CpG, calculate PDR as the fraction of overlapping reads with discordant flag
for (i in 1:length(ind)) {
  sample_CpG_PDR$PDR[i] <- length(which(eCpG_status.discordance[ind[[i]][,1]]))/length(ind[[i]][,1])
}
###


# Save CpG PDR table to file
write.table(sample_CpG_PDR, file = paste0(outdir,"/",strsplit(indir,"/")[[1]][grep("SCGP-HF",strsplit(indir,"/")[[1]])],"_individual_CpG_PDR-chr",chr,".txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

print(paste0("Finished calculating per-CpG PDR for chr ",chr," at ",date()))
################################################################################

