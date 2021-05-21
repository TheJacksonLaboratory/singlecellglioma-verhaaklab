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

## Load in clinical data (subtype, age, treatment, hypermutation_status).
# Supply metadata so that 10X filenames and samples can be linked together.
metadata = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
metadata$`10X_id_short` <- as.character(metadata$`10X_id_short`)
clin_data = metadata %>% 
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11))) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  select(case_barcode,idh_status, idh_codel_subtype = subtype, grade = who_grade, timepoint = initial_recurrence,  is_hypermutator = hypermutation)


## Load in the epimutation data:
epimut_cpg <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples_pdr_score_table-all_single_cells.txt", sep="\t", header=T, stringsAsFactors = F)

## Load in the quality control data for these samples.
rrbs_qc <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/scgp_cnv_status.txt", sep="\t", header=T, stringsAsFactors = F)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number.
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cell_num == 1, cpg_unique > 40000, conversion_rate > 95, tumor_cnv == 1) %>% 
  # remove the single-end libraries as these may have founding effects.
  filter(!library_id %in% c("SCGP-UC-917-1-1", "SCGP-UC-917-1-2", "SCGP-UC-917-1-3", "SCGP-UC-917-2-2")) %>% 
  filter(!grepl("SCGP-HF-", sample_barcode)) %>%
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11))) %>% 
  left_join(clin_data, by="case_barcode") %>%  
  mutate(file_path = paste("/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/rerun-dedup_bed_graph/", sample_barcode, "_pe.deduplicated.bismark.cov.gz", sep=""))

## Calculate the mean epimutation for each tumor.
epimut_pass_qc <- epimut_cpg %>% 
  filter(Sample%in%rrbs_qc_pass$sample_barcode)

## Read in an example coverage files.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)
files <- gsub("\\//", "\\/", files)
# Restrict to the *single* cells that passed our QC metrics.
files = files[files%in%rrbs_qc_pass$file_path]

## Sample UC917 with a cell that has around the average CpG count: 143,839
## SCGP-UC-917-01D-S5M-49S5
files[827]

dat = data.table::fread(cmd = sprintf("zcat < %s", files[827]), verbose = FALSE, showProgress = FALSE)
colnames(dat) <- c("chr", "start", "end", "meth", "meth_read", "unmeth_read")

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SFigd-example-histogram.pdf", width = 5, height = 3.5)
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

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig2e-cpg-distribution-karyogram.pdf", width=5, height=9, useDingbats = FALSE)
autoplot(cpg_distr, layout="karyogram", aes(color = DNAme)) + 
  scale_color_manual(values = c("0" = "#fc8d59", "50" = "#ffffbf", "100" = "#91bfdb")) + 
  plot_theme +
  null_y
dev.off()

png(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig2e-cpg-distribution-karyogramv2.png", width=5, height=9, res = 300, units ="in")
autoplot(cpg_distr, layout="karyogram", aes(color = DNAme)) + 
  scale_color_manual(values = c("0" = "#fc8d59", "50" = "#ffffbf", "100" = "#91bfdb")) + 
  plot_theme +
  null_y
dev.off()

### END ###