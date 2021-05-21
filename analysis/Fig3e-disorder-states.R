##################################
# Determine the DNAme and DNAme disorder between cell states.
# Updated: 2020.05.16
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(RColorBrewer)
# Function for plotting densities by group.
source("/Users/johnsk/Documents/Single-Cell-DNAmethylation/scripts/quality-metrics/densityRRBS.R")
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

## Summarized DNA methylation data across the 914 cells passing QC.
tiles_10kb = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(tiles_10kb) <- gsub("\\.", "-",  colnames(tiles_10kb))

## Create a mean methylation value across the 10-kb tiled region.
bin_counts_10kb <- colSums(!is.na(tiles_10kb[, 5:918]))
bin_meth_10kb <- colMeans(tiles_10kb[, 5:918], na.rm = TRUE)
broad_meth <- as.data.frame(t(bind_rows(bin_counts_10kb, bin_meth_10kb)))
colnames(broad_meth) <- c("bins_covered_10kb", "mean_methylation_10kb")
broad_meth$sample_barcode <- rownames(broad_meth)

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/scgp_cnv_status.txt", sep="\t", header=T, stringsAsFactors = F)

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number.
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cell_num == 1, cpg_unique > 40000, conversion_rate > 95, tumor_cnv == 1) %>% 
  left_join(meta_data, by=c("sample"="subject_id")) 

## Order the cells by their sample barcode.
rrbs_qc_pass_ordered <- rrbs_qc_pass %>% 
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11))) %>% 
  arrange(sample_barcode)
tiles_10kb_tumor <- as.matrix(tiles_10kb[ ,colnames(tiles_10kb)%in%rrbs_qc_pass_ordered$sample_barcode])

## Ensure that the cells are positioned in the same order.
all(colnames(tiles_10kb_tumor)==rrbs_qc_pass_ordered$sample_barcode)

## Create a density plot for these variable regions. Seems to follow a beta-distribution. 
densityRRBS(tiles_10kb_tumor, sampGroups = rrbs_qc_pass_ordered$case_barcode, main="10-kb tiled regions", legend = F)

## Load in the DNAme disroder table. 
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)

## Create a new variable to indicate "non_tumor" or shortened case barcode.
gg_epiallele = epiallele_info %>% 
  mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) 

## Read in the liger cell state data for each patient.
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

## Restrict to only the DNA methylation data.
idh_mut <- c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917")

## Promoter averaged DNAme values.
prom_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_FANTOM5_gene_matched_promoter_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(prom_cpgs) <- c("sample_barcode", "promoter_mean_meth", "promoter_cpgs")
prom_cpgs$sample_barcode <- gsub("_pe.deduplicated", "", prom_cpgs$sample_barcode)

scgp_classes_meth <- scgp_classes %>% 
  filter(dataset=="dna") %>% 
  inner_join(gg_epiallele, by=c("cell_name"="sample_barcode")) %>% 
  inner_join(broad_meth, by=c("cell_name"="sample_barcode")) %>% 
  # Collapsing the stem and proliferating stem cells into a single group.
  mutate(binary_liger_state = ifelse(liger_state%in%c("stemcell", "prolif_stemcell"), "Stem-like", "Diff.-like"),
         subtype = ifelse(case_normal_barcode%in%idh_mut, "IDHmut", "IDHwt")) %>% 
  inner_join(prom_cpgs, by=c("cell_name"="sample_barcode"))
#%>% 
  ## These samples do not have multiple populations.
  #filter(!case_normal_barcode%in%c("SM004", "SM018", "SM011"))

## Providing the appropriate case_order.
case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
scgp_classes_meth <- scgp_classes_meth %>% mutate(case_normal_barcode = factor(case_normal_barcode, levels = case_order))

## Define groups based on IDH mutation status.
scgp_classes_meth_idh <- scgp_classes_meth %>% filter(subtype=="IDHmut")
scgp_classes_meth_wt <- scgp_classes_meth %>% filter(subtype=="IDHwt")

## Examine the differences in PDR by cell state.
# Across individual subjects. PDR appears elevated in most samples.
## IDH-mut samples with at least scDNAm two populations.
pdf(file = "github/results/Fig3/IDHmut-promoter-DNAme-disorder.pdf", height = 4, width = 8, useDingbats = FALSE)
ggplot(scgp_classes_meth_idh, aes(x=binary_liger_state,  y=promoter_PDR, fill=binary_liger_state)) + 
  geom_boxplot(outlier.shape = NA)  +
  labs(y="Promoter DNAme disorder", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~case_normal_barcode, scales = "free", space="free_x") + 
  guides(alpha=FALSE) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  stat_n_text()
dev.off()



## IDHwt samples with at least scDNAm two populations.
pdf(file = "github/results/Fig3/IDHwt-promoter-DNAme-disorder.pdf", height = 4, width = 8, useDingbats = FALSE)
ggplot(scgp_classes_meth_wt, aes(x=binary_liger_state,  y=promoter_PDR, fill=binary_liger_state)) + 
  geom_boxplot(outlier.shape = NA)  +
  labs(y="Promoter DNAme disorder", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~case_normal_barcode, scales = "free", space="free_x") + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  stat_n_text()
dev.off()


## What is the mean methylation across these 10kb tiled and promoter regions?
## IDHmut 10kb.
ggplot(scgp_classes_meth_idh, aes(x=binary_liger_state,  y= promoter_mean_meth, fill=binary_liger_state)) + 
  geom_boxplot(outlier.shape = NA)  +
  labs(y="Mean promoter DNAme", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~case_normal_barcode, scales = "free", space="free_x") + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  stat_n_text()



## IDHwt 10-kb tiles.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/liger/tiles-10kb-IDHwt.pdf", height = 5, width = 7)
ggplot(scgp_classes_meth_wt, aes(x=binary_liger_state,  y= promoter_mean_meth, fill=binary_liger_state)) + 
  geom_boxplot(outlier.shape = NA)  +
  labs(y="Mean promoter DNAme", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~case_normal_barcode, scales = "free", space="free_x") + 
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  stat_n_text()
dev.off()


### Averaged across all subjects (epimutation and 10-kb tiles):
## DNAme disorder
pdf(file = "github/results/Fig3/Fig3e-disorder-all.pdf", height = 4, width = 6, useDingbats = FALSE)
ggplot(scgp_classes_meth, aes(x=binary_liger_state,  y=promoter_PDR, fill=binary_liger_state)) + 
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  #geom_boxplot(outlier.shape = NA)  +
  labs(y="Promoter DNAme disorder", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme  +
  facet_grid(~subtype, scales = "free", space="free_x") + 
  guides(alpha=FALSE) + stat_compare_means(method = "wilcox.test") +
  stat_n_text()
dev.off()

## Mean methylation across 10-kb tiles.
pdf(file = "github/results/Fig3/Fig3f-tiles-10kb-all.pdf", height = 4, width = 6, useDingbats = FALSE)
ggplot(scgp_classes_meth, aes(x=binary_liger_state,  y=mean_methylation_10kb, fill=binary_liger_state)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(y="10kb tiled DNAme", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme +
  facet_grid(~subtype, scales = "free", space="free_x") + 
  guides(alpha=FALSE) + stat_compare_means(method = "wilcox.test") +
  stat_n_text()
dev.off()

pdf(file = "github/results/Fig3/Fig3f-promoter-DNAme-all.pdf", height = 4, width = 6, useDingbats = FALSE)
ggplot(scgp_classes_meth, aes(x=binary_liger_state,  y=promoter_mean_meth, fill=binary_liger_state)) + 
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  #geom_boxplot(outlier.shape = NA) +
  labs(y="Mean promoter DNAme", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme +
  facet_grid(~subtype, scales = "free", space="free_x") + 
  guides(alpha=FALSE) + 
  stat_compare_means(method = "wilcox.test") +
  stat_n_text()
dev.off()


## Restrict to samples with at least a 20:80 cell state distribution.
scgp_classes_meth_filt <- scgp_classes_meth %>% 
  filter(case_normal_barcode%in%c("SM001", "SM002", "SM006", "SM008", "SM012", "SM017"))

pdf(file = "github/results/Fig3/Fig3f-promoter-disorder-filtered.pdf", height = 4, width = 8, useDingbats = FALSE)
ggplot(scgp_classes_meth_filt, aes(x=case_normal_barcode,  y=promoter_PDR, fill=binary_liger_state)) + 
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y="Mean promoter DNAme disorder", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme +
  facet_grid(~subtype, scales = "free") + 
  guides(alpha=FALSE) + 
  stat_compare_means(method = "wilcox.test", label = "p.format") 
dev.off()

pdf(file = "github/results/Fig3/Fig3f-promoter-DNAme-filtered.pdf", height = 4, width = 8, useDingbats = FALSE)
ggplot(scgp_classes_meth_filt, aes(x=case_normal_barcode,  y=promoter_mean_meth, fill=binary_liger_state)) + 
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y="Mean promoter DNAme", x="LIGER-defined cell state", fill = "Cell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1")) +
  plot_theme +
  facet_grid(~subtype, scales = "free") + 
  guides(alpha=FALSE) + 
  stat_compare_means(method = "wilcox.test", label = "p.format")
dev.off()

#### END ####
