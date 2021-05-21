##################################
# Determine the DNAme differences between cell states (separately for IDHmut & IDHwt)
# Updated: 2021.05.16
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
# Load the essential R-packages.
library(tidyverse)
library(nlme)
library(ggpubr)
library(openxlsx)
library(RColorBrewer)
library(LOLA)
library(GenomicRanges)
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

## Read in the liger data for each patient.
sm001_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM001/SM001-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm002_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM002/SM002-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm004_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM004/SM004-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm006_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM006/SM006-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm008_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM008/SM008-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm011_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM011/SM011-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm012_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM012/SM012-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm015_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM015/SM015-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm017_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM017/SM017-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm018_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/SM018/SM018-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
uc917_classes = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/UC917/UC917-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
## Restrict to the dna dataset.
scgp_classes <- bind_rows(sm001_classes, sm002_classes, sm004_classes, sm006_classes, sm008_classes,
                          sm011_classes, sm012_classes, sm015_classes, sm017_classes, sm018_classes,
                          uc917_classes)

## Summarized DNA methylation data across the 914 cells passing QC.
tiles_10kb = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(tiles_10kb) <- gsub("\\.", "-",  colnames(tiles_10kb))

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/scgp_cnv_status.txt", sep="\t", header=T, stringsAsFactors = F)

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number.
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cell_num == 1, cpg_unique > 40000, conversion_rate > 95, tumor_cnv == 1) %>% 
  left_join(meta_data, by=c("sample"="subject_id")) 

## The list of IDHmutated samples.
idh_mut <- c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917")
scgp_classes_meth <- scgp_classes %>% 
  filter(dataset=="dna") %>% 
  mutate(binary_liger_state = ifelse(liger_state%in%c("stemcell", "prolif_stemcell"), "stemcell", "differentiated"),
         subtype = ifelse(sample_names%in%idh_mut, "IDHmut", "IDHwt")) %>% 
  select(sample_barcode = cell_name, binary_liger_state, subtype)

## Properly arrange QC data.
rrbs_qc_pass_ordered <- rrbs_qc_pass %>% 
  inner_join(scgp_classes_meth, by="sample_barcode") %>% 
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11))) %>% 
  arrange(sample_barcode) %>% 
  select(sample_barcode, case_barcode, binary_liger_state, subtype = subtype.y)

## Separate into the two major glioma groups.
sc_idh_mut = rrbs_qc_pass_ordered %>%
  filter(subtype == "IDHmut")
sc_idh_wt = rrbs_qc_pass_ordered %>%
  filter(subtype == "IDHwt")

##########################
### IDHmut cells
##########################
## Limit to only those **IDHmut** samples.
tiles_10kb_beta_mut <- as.matrix(tiles_10kb[ ,colnames(tiles_10kb)%in%sc_idh_mut$sample_barcode])
all(colnames(tiles_10kb_beta_mut)==sc_idh_mut$sample_barcode)

## This function will normalize the beta-distribution of the tiled-regions.
logit2 = function(x) {
  log2(x/(1-x))
}
## Add small error value so that "-Inf" is not generated.
tiles_10kb_mut_mval = logit2(tiles_10kb_beta_mut+0.001)

## Prepare output table for linear mixed effects model.
mut_10kb_results = matrix(NA, nrow = dim(tiles_10kb_mut_mval)[1], ncol = 5)
colnames(mut_10kb_results) =  c("region", "celltype_coef", "celltype_pval", "differentiated_mean", "stemcell_mean")
mut_10kb_results = as.data.frame(mut_10kb_results)
mut_10kb_results$region = paste(tiles_10kb$Chromosome, tiles_10kb$Start, tiles_10kb$End, sep="-")

## Run the lme on 10kb tiles.
system.time(
  for ( i in 1:dim(tiles_10kb_mut_mval)[1]) { 
    Yi = tiles_10kb_mut_mval[i, ]
    fit = try(lme(Yi ~ binary_liger_state, random = ~1|case_barcode, data = sc_idh_mut, na.action = na.omit),  silent = T)
    if(!inherits(fit, "try-error")) {
      mut_10kb_results[i, "celltype_coef"] = summary(fit)$tTable[2,1]
      mut_10kb_results[i, "celltype_pval"] = summary(fit)$tTable[2,5]
      mut_10kb_results[i, "differentiated_mean"] = mean(tiles_10kb_beta_mut[i, which(sc_idh_mut$binary_liger_state=="differentiated")], na.rm = T)
      mut_10kb_results[i, "stemcell_mean"] = mean(tiles_10kb_beta_mut[i, which(sc_idh_mut$binary_liger_state=="stemcell")], na.rm = T)
    
    }
  }
)

## Process the results for IDHmut (10-kb tiles).
mut_10kb_results_set <- mut_10kb_results %>% 
  separate(region, c("chr", "start", "end"), sep="-")
## Perform multiple-hypotheses testing correction.
mut_10kb_results_set$celltype_pval_adj = p.adjust(mut_10kb_results_set$celltype_pval, method = "fdr", n = length(mut_10kb_results_set$celltype_pval))
## Calculate the mean difference between the two cell populations.
mut_10kb_results_set$mean_diff <- mut_10kb_results_set$stemcell_mean-mut_10kb_results_set$differentiated_mean
## Does it appear correct?
ggplot(mut_10kb_results_set, aes(mean_diff)) + geom_histogram() + geom_vline(xintercept = c(-0.25, 0.25), col="red")

## Define sets of results based on HYPERmethylation in stem cells and HYPOmethylation in stem cells.
stem_cell_sig_hyper_mut = mut_10kb_results_set %>% 
  filter(celltype_pval_adj < 0.05, mean_diff > 0.25)
stem_cell_sig_hypo_mut = mut_10kb_results_set %>% 
  filter(celltype_pval_adj < 0.05, mean_diff < -0.25)

# Create a genomic ranges object for HYPERmethylation.
stem_cell_sig_hyper_mut <- stem_cell_sig_hyper_mut[, c('chr', 'start', 'end')]
stem_cell_sig_hyper_mut_gr <- makeGRangesFromDataFrame(stem_cell_sig_hyper_mut, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                   seqnames.field = "chr", start.field="start", end.field="end")

# Create a genomic ranges object for HYPOmethylation.
stem_cell_sig_hypo_mut <- stem_cell_sig_hypo_mut[, c('chr', 'start', 'end')]
stem_cell_sig_hypo_mut_gr <- makeGRangesFromDataFrame(stem_cell_sig_hypo_mut, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                  seqnames.field = "chr", start.field="start", end.field="end")

# Create a genomic ranges object for the background set.
mut_10kb_results_set <- mut_10kb_results_set[, c('chr', 'start', 'end')]
mut_10kb_results_set_gr <- makeGRangesFromDataFrame(mut_10kb_results_set, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                    seqnames.field = "chr", start.field="start", end.field="end")

# Set Universe the total tiles that were covered in 25% of samples.
universe <- GRanges(mut_10kb_results_set_gr)
geneset_stem_hyper_mut <- GRanges(stem_cell_sig_hyper_mut_gr)
geneset_stem_hypo_mut <- GRanges(stem_cell_sig_hypo_mut_gr)

## Load the LOLA set that is specific to the JASPAR database.
lolaDB <- loadRegionDB("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/LOLAJaspar/hg19/")
## Perform association with 10kb tiles and HYPERmethylation patterns in stem cells.
locResults_stem_hyper_mut_10kb <- runLOLA(geneset_stem_hyper_mut, universe, lolaDB)
## convert -log10(p-value) to p-value and correct for FDR. 
locResults_stem_hyper_mut_10kb$pvalue <- 10^-locResults_stem_hyper_mut_10kb$pValueLog
## Perform multiple-hypotheses testing correction.
locResults_stem_hyper_mut_10kb$qValue = p.adjust(locResults_stem_hyper_mut_10kb$pvalue , method = "fdr", n = length(locResults_stem_hyper_mut_10kb$pvalue))
locResults_stem_hyper_mut_10kb_filt <- locResults_stem_hyper_mut_10kb[which(locResults_stem_hyper_mut_10kb$qValue<0.05) ,]
# Found: SP1, ZEB1, ELF5, FOXC1

## Perform association with 10kb tiles and HYPOmethylation patterns in stem cells.
locResults_stem_hypo_mut_10kb <- runLOLA(geneset_stem_hypo_mut, universe, lolaDB)
## convert -log10(p-value) to p-value and correct for FDR. 
locResults_stem_hypo_mut_10kb$pvalue <- 10^-locResults_stem_hypo_mut_10kb$pValueLog
## Perform multiple-hypotheses testing correction.
locResults_stem_hypo_mut_10kb$qValue = p.adjust(locResults_stem_hypo_mut_10kb$pvalue , method = "fdr", n = length(locResults_stem_hypo_mut_10kb$pvalue))
locResults_stem_hypo_mut_10kb_filt <- locResults_stem_hypo_mut_10kb[which(locResults_stem_hypo_mut_10kb$qValue<0.05) ,]
# Found: FOXD1


##########################
### IDHwt cells
##########################
## Limit to only those **IDHwt** samples.
tiles_10kb_beta_wt <- as.matrix(tiles_10kb[ ,colnames(tiles_10kb)%in%sc_idh_wt$sample_barcode])
all(colnames(tiles_10kb_beta_wt)==sc_idh_wt$sample_barcode)

## This function will normalize the beta-distribution of the tiled-regions.
logit2 = function(x) {
  log2(x/(1-x))
}
## Add small error value so that "-Inf" is not generated.
tiles_10kb_wt_mval = logit2(tiles_10kb_beta_wt+0.001)

## Prepare output table for linear mixed effects model.
wt_10kb_results = matrix(NA, nrow = dim(tiles_10kb_wt_mval)[1], ncol = 5)
colnames(wt_10kb_results) =  c("region", "celltype_coef", "celltype_pval", "differentiated_mean", "stemcell_mean")
wt_10kb_results = as.data.frame(wt_10kb_results)
wt_10kb_results$region = paste(tiles_10kb$Chromosome, tiles_10kb$Start, tiles_10kb$End, sep="-")


## Run the lme on 10kb tiles.
system.time(
  for ( i in 1:dim(tiles_10kb_wt_mval)[1]) { 
    Yi = tiles_10kb_wt_mval[i, ]
    fit = try(lme(Yi ~ binary_liger_state, random = ~1|case_barcode, data = sc_idh_wt, na.action = na.omit),  silent = T)
    if(!inherits(fit, "try-error")) {
      wt_10kb_results[i, "celltype_coef"] = summary(fit)$tTable[2,1]
      wt_10kb_results[i, "celltype_pval"] = summary(fit)$tTable[2,5]
      wt_10kb_results[i, "differentiated_mean"] = mean(tiles_10kb_beta_wt[i, which(sc_idh_wt$binary_liger_state=="differentiated")], na.rm = T)
      wt_10kb_results[i, "stemcell_mean"] = mean(tiles_10kb_beta_wt[i, which(sc_idh_wt$binary_liger_state=="stemcell")], na.rm = T)
      
    }
  }
)


## Process the results for IDHwt (10-kb tiles).
wt_10kb_results_set <- wt_10kb_results %>% 
  separate(region, c("chr", "start", "end"), sep="-")

## Perform multiple-hypotheses testing correction.
wt_10kb_results_set$celltype_pval_adj = p.adjust(wt_10kb_results_set$celltype_pval, method = "fdr", n = length(wt_10kb_results_set$celltype_pval))
## Calculate the mean difference between the two cell populations.
wt_10kb_results_set$mean_diff <- wt_10kb_results_set$stemcell_mean-wt_10kb_results_set$differentiated_mean
## Does it appear correct?
ggplot(wt_10kb_results_set, aes(mean_diff)) + geom_histogram() + geom_vline(xintercept = c(-0.25, 0.25), col="red")

## Define sets of results based on HYPERmethylation in stem cells and HYPOmethylation in stem cells.
stem_cell_sig_hyper_wt = wt_10kb_results_set %>% 
  filter(celltype_pval_adj < 0.05, mean_diff > 0.25)
stem_cell_sig_hypo_wt = wt_10kb_results_set %>% 
  filter(celltype_pval_adj < 0.05, mean_diff < -0.25)

# Create a genomic ranges object for HYPERmethylation.
stem_cell_sig_hyper_wt <- stem_cell_sig_hyper_wt[, c('chr', 'start', 'end')]
stem_cell_sig_hyper_wt_gr <- makeGRangesFromDataFrame(stem_cell_sig_hyper_wt, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                       seqnames.field = "chr", start.field="start", end.field="end")

# Create a genomic ranges object for HYPOmethylation.
stem_cell_sig_hypo_wt <- stem_cell_sig_hypo_wt[, c('chr', 'start', 'end')]
stem_cell_sig_hypo_wt_gr <- makeGRangesFromDataFrame(stem_cell_sig_hypo_wt, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                      seqnames.field = "chr", start.field="start", end.field="end")

# Create a genomic ranges object for the background set.
wt_10kb_results_set <- wt_10kb_results_set[, c('chr', 'start', 'end')]
wt_10kb_results_set_gr <- makeGRangesFromDataFrame(wt_10kb_results_set, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                                   seqnames.field = "chr", start.field="start", end.field="end")

# Set Universe the total tiles that were covered in 25% of samples.
universe <- GRanges(wt_10kb_results_set_gr)
geneset_stem_hyper_wt <- GRanges(stem_cell_sig_hyper_wt_gr)
geneset_stem_hypo_wt <- GRanges(stem_cell_sig_hypo_wt_gr)

## Load the LOLA set that is specific to the JASPAR database.
#lolaDB <- loadRegionDB("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/LOLAJaspar/hg19/")
## Perform association with 10kb tiles and HYPERmethylation patterns in stem cells.
locResults_stem_hyper_wt_10kb <- runLOLA(geneset_stem_hyper_wt, universe, lolaDB)
## convert -log10(p-value) to p-value and correct for FDR. 
locResults_stem_hyper_wt_10kb$pvalue <- 10^-locResults_stem_hyper_wt_10kb$pValueLog
## Perform multiple-hypotheses testing correction.
locResults_stem_hyper_wt_10kb$qValue = p.adjust(locResults_stem_hyper_wt_10kb$pvalue , method = "fdr", n = length(locResults_stem_hyper_wt_10kb$pvalue))
locResults_stem_hyper_wt_10kb_filt <- locResults_stem_hyper_wt_10kb[which(locResults_stem_hyper_wt_10kb$qValue<0.05) ,]
# Found: SP1, ZEB1, HOX5A, GATA3

## Perform association with 10kb tiles and HYPOmethylation patterns in stem cells.
locResults_stem_hypo_10kb_wt <- runLOLA(geneset_stem_hypo_wt, universe, lolaDB)
## convert -log10(p-value) to p-value and correct for FDR. 
locResults_stem_hypo_10kb_wt$pvalue <- 10^-locResults_stem_hypo_10kb_wt$pValueLog
## Perform multiple-hypotheses testing correction.
locResults_stem_hypo_10kb_wt$qValue = p.adjust(locResults_stem_hypo_10kb_wt$pvalue , method = "fdr", n = length(locResults_stem_hypo_10kb_wt$pvalue))
locResults_stem_hypo_10kb_wt_filt <- locResults_stem_hypo_10kb_wt[which(locResults_stem_hypo_10kb_wt$qValue<0.05) ,]
# Found: ELK4


#####################################
### Produce figure detailing results
#####################################
locResults_stem_hyper_mut_10kb_filt$direction <- "hyper"
locResults_stem_hyper_mut_10kb_filt$subtype <- "IDHmut"
locResults_stem_hypo_mut_10kb_filt$direction <- "hypo"
locResults_stem_hypo_mut_10kb_filt$subtype <- "IDHmut"

locResults_stem_hyper_wt_10kb_filt$direction <- "hyper"
locResults_stem_hyper_wt_10kb_filt$subtype <- "IDHwt"
locResults_stem_hypo_10kb_wt_filt$direction <- "hypo"
locResults_stem_hypo_10kb_wt_filt$subtype <- "IDHwt"
lola_results_10kb_mut = bind_rows(locResults_stem_hyper_mut_10kb_filt, locResults_stem_hypo_mut_10kb_filt)
lola_results_10kb_wt = bind_rows(locResults_stem_hyper_wt_10kb_filt, locResults_stem_hypo_10kb_wt_filt)
  
lola_results_10kb_mut = lola_results_10kb_mut %>% 
  mutate(filename = gsub("\\.bed", "", filename)) %>% 
  arrange(filename) %>% 
  ## Remove non-human TFs.
  filter(!grepl("[a-z]", filename)) %>% 
  ## Remove "MZF1_1-4" because there appears to be redundancy with a later version.
  filter(!filename%in%c("MZF1_1-4"))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/IDHmut-methylation-diff-cell-states.pdf", height = 5, width = 7)
ggplot(lola_results_10kb_mut, aes(x=filename, y=-log10(qValue), color=direction, size=oddsRatio)) + 
  geom_point() +
  ylim(0, 4.5) +
  labs(x="", y="LOLA enrichment\n-log10(adjusted p-value)") +
  scale_color_manual(values=c("hyper" = "#fb6a4a", 
                              "hypo" = "#fcbba1")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  plot_theme +
  facet_grid(~direction, scales = "free", space="free_x")
dev.off()

lola_results_10kb_wt = lola_results_10kb_wt %>% 
  mutate(filename = gsub("\\.bed", "", filename)) %>% 
  arrange(filename) %>% 
  ## Remove non-human TFs.
  filter(!grepl("[a-z]", filename)) 
  
## All together.
lola_results_10kb = bind_rows(locResults_stem_hyper_mut_10kb_filt, locResults_stem_hypo_mut_10kb_filt,
                                  locResults_stem_hyper_wt_10kb_filt, locResults_stem_hypo_10kb_wt_filt)
lola_results_10kb_final = lola_results_10kb %>% 
  mutate(filename = gsub("\\.bed", "", filename)) %>% 
  arrange(filename) %>% 
  ## Remove non-human TFs.
  filter(!grepl("[a-z]", filename)) %>% 
  filter(!filename%in%c("MZF1_1-4")) %>% 
  distinct() %>% 
  mutate(tf = sapply(strsplit(filename, "_|__"), "[[", 1))

tf_order <- lola_results_10kb_final %>% 
  arrange(desc(qValue)) %>% 
  select(tf) %>% 
  distinct()
lola_results_10kb_final$tf <-  factor(lola_results_10kb_final$tf, levels = tf_order$tf )

pdf(file = "github/results/Fig3/Fig3f-TFBS-enrich-DNAMe-diff-sorted.pdf", height = 4, width = 8, useDingbats = FALSE)
ggplot(lola_results_10kb_final, aes(x=tf, y=-log10(qValue), color=subtype, size=2)) + 
  ylim(0, 7) +
  geom_point(alpha=0.8) +
  labs(x="", y="LOLA enrichment\n-log10(adjusted p-value)", color="Glioma\nsubtype") +
  scale_color_manual(values=c("IDHwt" = "#7fbf7b", 
                              "IDHmut" = "#af8dc3")) +
  geom_hline(yintercept = -log10(0.05), color="red", linetype = "dashed") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(~direction, scales = "free", space="free_x") +
  guides(size=FALSE)
dev.off()

### END ####

