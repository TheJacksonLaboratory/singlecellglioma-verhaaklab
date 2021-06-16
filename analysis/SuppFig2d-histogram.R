##############################################
# Generate a karyogram for a single-cell example
# Updated: 2020.05.29
# Author: Kevin J.
##################################################

# Working directory for this analysis on Helix.
mybasedir = "/Users/johnsk/mnt/verhaak-lab/scgp/"
datadir  = "/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/rerun-dedup_bed_graph"
pattern   = '_pe.deduplicated.bismark.cov.gz$'
setwd(mybasedir)

###################################
library(tidyverse)
library(data.table)
library(openxlsx)
library(parallel)
library(GenomicRanges)
library(ggbio)
library(RColorBrewer)
###################################
## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())


## Data are available on Synapse for single cells and 50-cell populations under the bismark coverage files.
## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="/Users/johnsk/github/data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("/Users/johnsk/github/data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(meta_data, by="case_barcode")  %>% 
  mutate(file_path = paste("/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/rerun-dedup_bed_graph/", cell_barcode, "_pe.deduplicated.bismark.cov.gz", sep=""))


## Read in an example coverage files.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)
files <- gsub("\\//", "\\/", files)

# Restrict to the *single* cells that passed our QC metrics.
files = files[files%in%rrbs_qc_pass$file_path]

## Sample SM019 with a cell that has around the average CpG count: 143,839
## SCGP-SM-019-01D-S5M-78S5
files[827]

dat = data.table::fread(cmd = sprintf("zcat < %s", file_to_read), verbose = FALSE, showProgress = FALSE)
colnames(dat) <- c("chr", "start", "end", "meth", "meth_read", "unmeth_read")

pdf(file = "github/results/Fig1/SFigd-example-histogram.pdf", width = 5, height = 3.5)
ggplot(dat, aes(x=meth)) + geom_histogram( fill="#377eb8") +
  labs(y="CpG frequency", x="DNA methylation percentage (%)") +
  plot_theme
dev.off()  

dat_filt <- dat %>%
  filter(!grepl("GL|X|Y|MT", chr)) %>% 
  mutate(chr = paste("chr", chr, sep="")) %>% 
  mutate(meth = ifelse(meth > 1 & meth <49, 0, meth)) %>% 
  mutate(meth = ifelse(meth > 51 & meth < 99, 100, meth)) %>% 
  mutate(meth = as.factor(meth))
colnames(dat_filt)[4] <- "DNAme"

# Create GRanges object.
cpg_distr <- makeGRangesFromDataFrame(dat_filt, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="start", end.field="end")
seqlengths(cpg_distr)

pdf(file = "github/results/Fig1/SuppFig2e-cpg-distribution-karyogram.pdf", width=5, height=9, useDingbats = FALSE)
autoplot(cpg_distr, layout="karyogram", aes(color = DNAme)) + 
  scale_color_manual(values = c("0" = "#fc8d59", "50" = "#ffffbf", "100" = "#91bfdb")) + 
  plot_theme +
  null_y
dev.off()


### END ###