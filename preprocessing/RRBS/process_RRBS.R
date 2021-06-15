# process_RRBS.R - this script assesses quality of RRBS data generated using bismark methylation extractor

require(methylKit)
require(genomation)
require(GenomicRanges)
require(GenomicFeatures)


##### Input arguments #####
# Sample linker file
#linker <- read.delim("~/sumner/verhaak-lab/scgp/data/RRBS/20-verhaak-006/20-verhaak-006-linker.txt", header = FALSE)
linker <- read.delim("/projects/verhaak-lab/scgp/data/RRBS/20-verhaak-006/20-verhaak-006-linker.txt", header = FALSE)

# Directory containing sample bismark methylation coverage data
#meth_dir <- "~/sumner/verhaak-lab/scgp/results/RRBS/human/bed_graph"
meth_dir <- "/projects/verhaak-lab/scgp/results/RRBS/human/bed_graph"

# Working directory
#workdir <- "~/sumner/verhaak-lab/scgp/results/RRBS/human/methylKit"
workdir <- "/projects/verhaak-lab/scgp/results/RRBS/human/methylKit"
###########################


# Load bismark methylation coverage files
meth_files <- list.files(path = meth_dir, pattern = "*.cov.gz", full.names = TRUE)

# Extract sample barcode from coverage files
sample_id <- substr(basename(meth_files),1,25)

# Separate coverage files and associated barcode into one group per cell line
meth_files.HF2354 <- meth_files[grep("HF-2354",meth_files)]
meth_files.HF3016 <- meth_files[grep("HF-3016",meth_files)]

sample_id.HF2354 <- sample_id[grep("HF-2354", sample_id)]
sample_id.HF3016 <- sample_id[grep("HF-3016", sample_id)]

# Define treatment groups for each cell line, with 21% oxygen as control (treatment = 0)
oxygen_concentration.HF2354 <- substr(sample_id.HF2354,24,25)
treatment.HF2354 <- ifelse(oxygen_concentration.HF2354 == "21",0,1)
oxygen_concentration.HF3016 <- substr(sample_id.HF3016,24,25)
treatment.HF3016 <- ifelse(oxygen_concentration.HF3016 == "21",0,1)


### Processing samples: minimum 10x coverage
### NOTE: manipulating methylkit object columns converts the object from methylkit to regular list, data.frame, etc.
## HF2354
# Process the HF2354 samples as a methylRawList object, retaining CpGs with 10x coverage
meth_dat.HF2354.10x <- methRead(as.list(meth_files.HF2354),
                                sample.id=as.list(sample_id.HF2354),
                                assembly="hg19",
                                pipeline = "bismarkCoverage",
                                treatment=treatment.HF2354,
                                context="CpG",
                                mincov = 10
)

# Save methylRawList object
saveRDS(meth_dat.HF2354.10x, file = paste0(workdir,"/","meth_dat.HF2354.10x.Rds"))

# Unite sample objects to one table containing CpGs present in all samples
meth_dat.merged.HF2354.10x <- unite(meth_dat.HF2354.10x)

# Save merged methylRaw object
saveRDS(meth_dat.merged.HF2354.10x, file = paste0(workdir,"/","meth_dat.merged.HF2354.10x.Rds"))


# Filter CpGs in individual tables to those overlapping autosomal chromosomes, chr X
meth_dat.HF2354.10x <- lapply(meth_dat.HF2354.10x, function(x) x[which(x$chr %in% c(1:22,"X")),])

# Sort ranges in individual tables
meth_dat.HF2354.10x <- lapply(meth_dat.HF2354.10x, function(x) x[order(x$chr,x$start),])

# Generate BED files of 10x coverage CpGs for each sample
for (i in 1:length(meth_dat.HF2354.10x)) {
  # Create bed-format data frame of sample meth, with score corresponding to base coverage
  temp <- data.frame(chr=meth_dat.HF2354.10x[[i]]$chr,
                     start=meth_dat.HF2354.10x[[i]]$start - 1,
                     end=meth_dat.HF2354.10x[[i]]$end,
                     name=0,
                     score=meth_dat.HF2354.10x[[i]]$coverage,
                     strand=".")
  
  # Save data frame as bed file
  write.table(temp, file = paste0(workdir,"/",meth_dat.HF2354.10x[[i]]@sample.id,"_10x_coverage_CpGs.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


# Sort ranges in merged tables and convert from methylRaw to data frame
meth_dat.merged.HF2354.10x <- meth_dat.merged.HF2354.10x[order(meth_dat.merged.HF2354.10x$chr,meth_dat.merged.HF2354.10x$start),]
meth_dat.merged.HF2354.10x <- getData(meth_dat.merged.HF2354.10x)

# Create bed-format data frame of merged meth, with score corresponding to median sample base coverage
median_meth_coverage <- apply(meth_dat.merged.HF2354.10x[,grep("coverage",names(meth_dat.merged.HF2354.10x))],1,median)
temp.merged <- data.frame(chr=meth_dat.merged.HF2354.10x$chr,
                          start=meth_dat.merged.HF2354.10x$start - 1,
                          end=meth_dat.merged.HF2354.10x$end,
                          name=0,
                          score=median_meth_coverage,
                          strand=".")

# Save data frame as bed file
write.table(temp.merged, file = paste0(workdir,"/","HF2354_consensus_10x_coverage_CpGs.bed"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


## HF3016
# Process the HF3016 samples as a methylRawList object, retaining CpGs with 10x coverage
meth_dat.HF3016.10x <- methRead(as.list(meth_files.HF3016),
                                sample.id=as.list(sample_id.HF3016),
                                assembly="hg19",
                                pipeline = "bismarkCoverage",
                                treatment=treatment.HF3016,
                                context="CpG",
                                mincov = 10
)

# Save methylRawList object
saveRDS(meth_dat.HF3016.10x, file = paste0(workdir,"/","meth_dat.HF3016.10x.Rds"))

# Unite sample objects to one table containing CpGs present in all samples
meth_dat.merged.HF3016.10x <- unite(meth_dat.HF3016.10x)

# Save merged methylRaw object
saveRDS(meth_dat.merged.HF2354.10x, file = paste0(workdir,"/","meth_dat.merged.HF2354.10x.Rds"))


# Filter CpGs in individual tables to those overlapping autosomal chromosomes, chr X
meth_dat.HF3016.10x <- lapply(meth_dat.HF3016.10x, function(x) x[which(x$chr %in% c(1:22,"X")),])

# Sort ranges in individual tables
meth_dat.HF3016.10x <- lapply(meth_dat.HF3016.10x, function(x) x[order(x$chr,x$start),])

# Generate BED files of 10x coverage CpGs for each sample
for (i in 1:length(meth_dat.HF3016.10x)) {
  # Create bed-format data frame of sample meth, with score corresponding to base coverage
  temp <- data.frame(chr=meth_dat.HF3016.10x[[i]]$chr,
                     start=meth_dat.HF3016.10x[[i]]$start - 1,
                     end=meth_dat.HF3016.10x[[i]]$end,
                     name=0,
                     score=meth_dat.HF3016.10x[[i]]$coverage,
                     strand=".")
  
  # Save data frame as bed file
  write.table(temp, file = paste0(workdir,"/",meth_dat.HF3016.10x[[i]]@sample.id,"_10x_coverage_CpGs.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


# Sort ranges in merged tables and convert from methylRaw to data frame
meth_dat.merged.HF3016.10x <- meth_dat.merged.HF3016.10x[order(meth_dat.merged.HF3016.10x$chr,meth_dat.merged.HF3016.10x$start),]
meth_dat.merged.HF3016.10x <- getData(meth_dat.merged.HF3016.10x)

# Create bed-format data frame of merged meth, with score corresponding to median sample base coverage
median_meth_coverage <- apply(meth_dat.merged.HF3016.10x[,grep("coverage",names(meth_dat.merged.HF3016.10x))],1,median)
temp.merged <- data.frame(chr=meth_dat.merged.HF3016.10x$chr,
                          start=meth_dat.merged.HF3016.10x$start - 1,
                          end=meth_dat.merged.HF3016.10x$end,
                          name=0,
                          score=median_meth_coverage,
                          strand=".")

# Save data frame as bed file
write.table(temp.merged, file = paste0(workdir,"/","HF3016_consensus_10x_coverage_CpGs.bed"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
###


### Processing samples: minimum 5x coverage
### NOTE: manipulating methylkit object columns converts the object from methylkit to regular list, data.frame, etc.
## HF2354
# Process the HF2354 samples as a methylRawList object, retaining CpGs with 5x coverage
meth_dat.HF2354.5x <- methRead(as.list(meth_files.HF2354),
                                sample.id=as.list(sample_id.HF2354),
                                assembly="hg19",
                                pipeline = "bismarkCoverage",
                                treatment=treatment.HF2354,
                                context="CpG",
                                mincov = 5
)

# Save methylRawList object
saveRDS(meth_dat.HF2354.5x, file = paste0(workdir,"/","meth_dat.HF2354.5x.Rds"))

# Unite sample objects to one table containing CpGs present in all samples
meth_dat.merged.HF2354.5x <- unite(meth_dat.HF2354.5x)

# Save merged methylRaw object
saveRDS(meth_dat.merged.HF2354.5x, file = paste0(workdir,"/","meth_dat.merged.HF2354.5x.Rds"))


# Filter CpGs in individual tables to those overlapping autosomal chromosomes, chr X
meth_dat.HF2354.5x <- lapply(meth_dat.HF2354.5x, function(x) x[which(x$chr %in% c(1:22,"X")),])

# Sort ranges in individual tables
meth_dat.HF2354.5x <- lapply(meth_dat.HF2354.5x, function(x) x[order(x$chr,x$start),])

# Generate BED files of 5x coverage CpGs for each sample
for (i in 1:length(meth_dat.HF2354.5x)) {
  # Create bed-format data frame of sample meth, with score corresponding to base coverage
  temp <- data.frame(chr=meth_dat.HF2354.5x[[i]]$chr,
                     start=meth_dat.HF2354.5x[[i]]$start - 1,
                     end=meth_dat.HF2354.5x[[i]]$end,
                     name=0,
                     score=meth_dat.HF2354.5x[[i]]$coverage,
                     strand=".")
  
  # Save data frame as bed file
  write.table(temp, file = paste0(workdir,"/",meth_dat.HF2354.5x[[i]]@sample.id,"_5x_coverage_CpGs.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


# Sort ranges in merged tables and convert from methylRaw to data frame
meth_dat.merged.HF2354.5x <- meth_dat.merged.HF2354.5x[order(meth_dat.merged.HF2354.5x$chr,meth_dat.merged.HF2354.5x$start),]
meth_dat.merged.HF2354.5x <- getData(meth_dat.merged.HF2354.5x)

# Create bed-format data frame of merged meth, with score corresponding to median sample base coverage
median_meth_coverage <- apply(meth_dat.merged.HF2354.5x[,grep("coverage",names(meth_dat.merged.HF2354.5x))],1,median)
temp.merged <- data.frame(chr=meth_dat.merged.HF2354.5x$chr,
                          start=meth_dat.merged.HF2354.5x$start - 1,
                          end=meth_dat.merged.HF2354.5x$end,
                          name=0,
                          score=median_meth_coverage,
                          strand=".")

# Save data frame as bed file
write.table(temp.merged, file = paste0(workdir,"/","HF2354_consensus_5x_coverage_CpGs.bed"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


## HF3016
# Process the HF3016 samples as a methylRawList object, retaining CpGs with 10x coverage
meth_dat.HF3016.5x <- methRead(as.list(meth_files.HF3016),
                                sample.id=as.list(sample_id.HF3016),
                                assembly="hg19",
                                pipeline = "bismarkCoverage",
                                treatment=treatment.HF3016,
                                context="CpG",
                                mincov = 5
)

# Save methylRawList object
saveRDS(meth_dat.HF3016.5x, file = paste0(workdir,"/","meth_dat.HF3016.5x.Rds"))

# Unite sample objects to one table containing CpGs present in all samples
meth_dat.merged.HF3016.5x <- unite(meth_dat.HF3016.5x)

# Save merged methylRaw object
saveRDS(meth_dat.merged.HF3016.5x, file = paste0(workdir,"/","meth_dat.merged.HF3016.5x.Rds"))


# Filter CpGs in individual tables to those overlapping autosomal chromosomes, chr X
meth_dat.HF3016.5x <- lapply(meth_dat.HF3016.5x, function(x) x[which(x$chr %in% c(1:22,"X")),])

# Sort ranges in individual tables
meth_dat.HF3016.5x <- lapply(meth_dat.HF3016.5x, function(x) x[order(x$chr,x$start),])

# Generate BED files of 5x coverage CpGs for each sample
for (i in 1:length(meth_dat.HF3016.5x)) {
  # Create bed-format data frame of sample meth, with score corresponding to base coverage
  temp <- data.frame(chr=meth_dat.HF3016.5x[[i]]$chr,
                     start=meth_dat.HF3016.5x[[i]]$start - 1,
                     end=meth_dat.HF3016.5x[[i]]$end,
                     name=0,
                     score=meth_dat.HF3016.5x[[i]]$coverage,
                     strand=".")
  
  # Save data frame as bed file
  write.table(temp, file = paste0(workdir,"/",meth_dat.HF3016.5x[[i]]@sample.id,"_5x_coverage_CpGs.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


# Sort ranges in merged tables and convert from methylRaw to data frame
meth_dat.merged.HF3016.5x <- meth_dat.merged.HF3016.5x[order(meth_dat.merged.HF3016.5x$chr,meth_dat.merged.HF3016.5x$start),]
meth_dat.merged.HF3016.5x <- getData(meth_dat.merged.HF3016.5x)

# Create bed-format data frame of merged meth, with score corresponding to median sample base coverage
median_meth_coverage <- apply(meth_dat.merged.HF3016.5x[,grep("coverage",names(meth_dat.merged.HF3016.5x))],1,median)
temp.merged <- data.frame(chr=meth_dat.merged.HF3016.5x$chr,
                          start=meth_dat.merged.HF3016.5x$start - 1,
                          end=meth_dat.merged.HF3016.5x$end,
                          name=0,
                          score=median_meth_coverage,
                          strand=".")

# Save data frame as bed file
write.table(temp.merged, file = paste0(workdir,"/","HF3016_consensus_5x_coverage_CpGs.bed"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
###


# ##### Visualize QC metrics #####
# # Calculate
# ### Plot number unique CpGs for each cohort
# pdf(file = paste0(workdir,"/","hypoxia_experiment_RRBS-%_methylation_per_base.pdf"), width = 10, height = 10)
# par(mfrow = c(2, 2))
# # HF-2354 5x
# 
# # HF-2354 10x
# # HF-3016 5x
# # HF-3016 10x
# def.off()
# ###
# ################################
