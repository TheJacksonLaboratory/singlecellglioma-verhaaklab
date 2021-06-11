##################################
# Determine the DNAme and DNAme disorder between cell states.
# Updated: 2020.05.16
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(RColorBrewer)
library(EnvStats)
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

## Summarized DNA methylation data across the 914 cells passing QC (wide version of `analysis_scRRBS_10kb_tiled_methylation.tsv.gz`).
tiles_10kb = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(tiles_10kb) <- gsub("\\.", "-",  colnames(tiles_10kb))
## colnames(tiles_10kb) <- gsub("UC-917", "SM-019", colnames(tiles_10kb))

## Create a mean methylation value across the 10-kb tiled region.
bin_counts_10kb <- colSums(!is.na(tiles_10kb[, 5:918]))
bin_meth_10kb <- colMeans(tiles_10kb[, 5:918], na.rm = TRUE)
broad_meth <- as.data.frame(t(bind_rows(bin_counts_10kb, bin_meth_10kb)))
colnames(broad_meth) <- c("bins_covered_10kb", "mean_methylation_10kb")
broad_meth$cell_barcode <- rownames(broad_meth)

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(meta_data, by="case_barcode") 

## Order the cells by their sample barcode.
rrbs_qc_pass_ordered <- rrbs_qc_pass %>% 
  arrange(cell_barcode)
tiles_10kb_tumor <- as.matrix(tiles_10kb[ ,colnames(tiles_10kb)%in%rrbs_qc_pass_ordered$cell_barcode])

## Ensure that the cells are positioned in the same order.
all(colnames(tiles_10kb_tumor)==rrbs_qc_pass_ordered$cell_barcode)

## Load in the DNAme disorder table. 
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epiallele_info$sample_barcode <- gsub("UC-917", "SM-019", epiallele_info$sample_barcode)
disorder_dat <- read.table(file="/Users/johnsk/github/data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)
all(epiallele_info$sample_barcode==disorder_dat$cell_barcode)

epiallele_info$promoter_PDR
disorder_dat$promoter_PDR

## Create a new variable to indicate "non_tumor" or shortened case barcode.
gg_epiallele = disorder_dat %>% 
  inner_join(rrbs_qc, by=c("cell_barcode", "case_barcode")) %>% 
  mutate(case_normal_barcode = ifelse(tumor_status==0, "Non-tumor", case_barcode)) 

## Read in the liger cell state data for each patient.
scgp_classes <- read.table(file="/Users/johnsk/github/data/analysis_scRRBS_liger_classifications.csv", sep = ",", header = TRUE)

## Restrict to only the DNA methylation data.
idh_mut <- c("SM001", "SM002", "SM004", "SM008", "SM015", "SM)19")

## Promoter averaged DNAme values for measured promoters (not available on Synapse, but could be generated with analysis_scRRBS_individual_promoter_methylation.csv).
prom_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_FANTOM5_gene_matched_promoter_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(prom_cpgs) <- c("sample_barcode", "promoter_mean_meth", "promoter_cpgs")
prom_cpgs$cell_barcode <- gsub("_pe.deduplicated", "", prom_cpgs$sample_barcode)

scgp_classes_meth <- scgp_classes %>% 
  filter(dataset=="dna") %>% 
  inner_join(gg_epiallele, by="cell_barcode") %>% 
  inner_join(broad_meth, by="cell_barcode") %>% 
  # Collapsing the stem and proliferating stem cells into a single group.
  mutate(binary_liger_state = ifelse(liger_state%in%c("stemcell", "prolif_stemcell"), "Stem-like", "Diff.-like"),
         subtype = ifelse(case_normal_barcode%in%idh_mut, "IDHmut", "IDHwt")) %>% 
  inner_join(prom_cpgs, by="cell_barcode")

## Providing the appropriate case_order.
case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
scgp_classes_meth <- scgp_classes_meth %>% mutate(case_normal_barcode = factor(case_normal_barcode, levels = case_order))

## Define groups based on IDH mutation status.
scgp_classes_meth_idh <- scgp_classes_meth %>% filter(subtype=="IDHmut")
scgp_classes_meth_wt <- scgp_classes_meth %>% filter(subtype=="IDHwt")

## Examine the differences in PDR by cell state.
# Across individual subjects. PDR appears elevated in most samples.
## IDH-mut samples with at least scDNAm two populations.
pdf(file = "results/Fig3/IDHmut-promoter-DNAme-disorder.pdf", height = 4, width = 8, useDingbats = FALSE)
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
