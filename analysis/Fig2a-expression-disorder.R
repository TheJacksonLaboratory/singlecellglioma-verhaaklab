##############################################
# Investigate DNAme disorder and gene expression levels.
# Updated: 2021.05.13
# Author: Kevin J.
##################################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
library(tidyverse)
library(Matrix)
library(openxlsx)
library(Seurat)
library(EnvStats)
library(ggpubr)
# Method to determine "high" levels based on ranked data.
source("singlecellglioma-verhaaklab/analysis//ROSE-high-delineation.R")
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

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


## Load the SCGP raw counts matrix that has already been filtered for quality metrics.
scgp_obj <- ReadH5AD("data/analysis_scRNAseq_tumor_counts.h5ad")
scgp_obj@meta.data$cell_barcode <- rownames(scgp_obj@meta.data)
avail_genes = rownames(scgp_obj@assays$RNA@data)

## Cell state information.
meta_10x <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

# Limit to IDHmut samples.
meta_10x_mut <- meta_10x %>% 
  filter(case_barcode%in%c("SM019", "SM001", "SM002", "SM004", "SM008", "SM015"))

## Keep only tumor cells.
cells_to_keep = which(meta_10x_mut$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
mut_clust_annot <- meta_10x_mut[cells_to_keep, ]

## Subset to IDHmut tumor cells.
scgp_obj_idh_mut_tumor <- subset(x = scgp_obj, subset = cell_barcode %in% mut_clust_annot$cell_barcode)

## Make a Seurat object with the standard pipeline through PCA.
scgp_obj_idh_mut_tumor <- scgp_obj_idh_mut_tumor %>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

## Create a data.frame with the expression variability metrics.
rna_var_mut = scgp_obj_idh_mut_tumor@assays$RNA@meta.features
rna_var_mut$gene_id <- rownames(rna_var_mut)


######################################
### DNAme disorder / PDR estimates ###
######################################
## Determine promoter-level PDR levels based on categories (i.e., high vs. low):
prom_epimut <- read.table("data/Samples-passQC_single_cells_individual_promoter-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)
## Revise the sample names to match between the promoter methylation data.
colnames(prom_epimut) <- gsub("\\.", "-",  colnames(prom_epimut))

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(meta_data, by="case_barcode") %>% 
  mutate(IDH_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))

## Define subgroups based on IDHmut status.
rrbs_qc_pass_wt = rrbs_qc_pass %>% filter(IDH_status=="IDHwt")
rrbs_qc_pass_mut = rrbs_qc_pass %>% filter(IDH_status=="IDHmut")

## Restrict to only tumor cells.
prom_epimut_tumor = prom_epimut[, colnames(prom_epimut)%in%rrbs_qc_pass$cell_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(prom_epimut_tumor)/dim(prom_epimut_tumor)[2])
## Identify the number of gene promoters measured in at least 100 cells.
sum(missingness<(1-100/844))

## Restrict to genes covered in at least 100 tumor cells samples:
epimut_prom_cov <- as.data.frame(rowSums(!is.na(prom_epimut_tumor)))
epimut_prom_avg <- as.data.frame(rowMeans(prom_epimut_tumor, na.rm = TRUE))
epimut_prom_total = cbind(epimut_prom_avg, epimut_prom_cov)
colnames(epimut_prom_total) <- c("tumor_pdr", "tumor_cov")
epimut_prom_total$gene_id <- rownames(epimut_prom_total)
epimut_prom_total_filt = epimut_prom_total %>% 
  filter(tumor_cov>100)

## The distribution is heavily tilted toward low epimutation rate.
hist(epimut_prom_total_filt$tumor_pdr)
ggplot(epimut_prom_total_filt, aes(x=tumor_cov, y=tumor_pdr)) + geom_point(alpha=0.2) +
  xlim(0, 844) + stat_cor(method="spearman") +
  plot_theme

## Create sort-able object for plotting by increasing gene promoter epimutation.
sort_df <- epimut_prom_total_filt %>%
  arrange(desc(tumor_pdr)) 
gene_order <- rev(unique(sort_df$gene_id))
epimut_prom_total_filt <- epimut_prom_total_filt %>% mutate(gene_id = factor(gene_id, levels = gene_order))

## Quick visualization.
ggplot(epimut_prom_total_filt, aes(x=gene_id)) +
  geom_point(aes(y=tumor_pdr)) +
  plot_theme +
  null_x

## Determine the CUTOFF for high gene promoter epimutation.
epimut_vector = epimut_prom_total_filt$tumor_pdr
calculate_cutoff(epimut_vector) # 0.399 is the cut-off for "high" DNAme disorder (aka epimutation).

## Add the new variable for classification:
epimut_prom_total_filt = epimut_prom_total_filt %>% 
  mutate(tumor_class = ifelse(tumor_pdr > 0.399, "high", "low"))

## Combine RNA data with promoter epimutation data.
rna_var_mut_epi = rna_var_mut %>% 
  inner_join(epimut_prom_total_filt, by="gene_id")

## Break into groups based on levels of DNAme disorder (epimutation; low, intermediate, and high).
## Breaks were chosen to represent easy to follow groups of 0.1 increments with 0.4 and greater being high epimutation.
rna_var_mut_epi$epimut_groups <-cut(rna_var_mut_epi$tumor_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))

## How do these groups look in terms of total number of genes.
table(rna_var_mut_epi$epimut_groups) # Fewer genes with high DNAme disorder.

## Plots for both gene expression and variability that is controlled for expression levels.
## The dispersion method helps control for the relationship between variability and average expression.
my_comparisons <- list( c("(-0.01,0.1]", "(0.1,0.2]"), c("(-0.01,0.1]", "(0.3,0.4]"), c("(0.1,0.2]", "(0.4,1.1]") )

## Calculate the summary values across groups:
rna_var_mut_epi %>% 
  group_by(epimut_groups) %>% 
  summarise(avg_expression = median(mvp.mean, na.rm = T), 
            avg_dispersion = median(mvp.dispersion.scaled, na.rm = T))

ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.dispersion.scaled)) + 
  geom_boxplot(outlier.shape = NA) +
  #ylim(-2.25, 2.25) +  
  labs(x="Promoter epimutation groups", y="Dispersion (mean-scaled)") +
  plot_theme + 
  stat_compare_means(comparisons = my_comparisons)  + 
  stat_n_text(y.pos = -3.25) 

## There are some outliers that squish the graph making it difficult to visually observe the shift.
ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.mean)) + 
  geom_boxplot() +
  stat_compare_means(method="kruskal.test") + 
  labs(x="Promoter epimutation groups", y="Mean expression") +
  plot_theme + 
  stat_n_text(y.pos = 1)


######################
### IDHwt         ####
######################
# Limit to IDHwt samples.
meta_10x_wt <- meta_10x %>% 
  filter(case_barcode%in%c("SM006", "SM012", "SM017", "SM018", "SM011"))

## Keep only tumor cells.
cells_to_keep = which(meta_10x_wt$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
wt_clust_annot <- meta_10x_wt[cells_to_keep, ]

## Subset to IDHwt tumor cells.
scgp_obj_idh_wt_tumor <- subset(x = scgp_obj, subset = cell_barcode %in% wt_clust_annot$cell_barcode)

## Make a Seurat object with the standard pipeline through PCA.
scgp_obj_idh_wt_tumor <- scgp_obj_idh_wt_tumor %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

## Create a data.frame with the expression variability metrics.
rna_var_wt = scgp_obj_idh_wt_tumor@assays$RNA@meta.features
rna_var_wt$gene_id <- rownames(rna_var_wt)

## Combine RNA data with promoter epimutation data.
rna_var_wt_epi = rna_var_wt %>% 
  inner_join(epimut_prom_total_filt, by="gene_id")

## Break into groups based on levels of epimutation (low, intermediate, and high). 
rna_var_wt_epi$epimut_groups <-cut(rna_var_wt_epi$tumor_pdr,c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))

## How do these groups look in terms of total number of genes.
table(rna_var_wt_epi$epimut_groups) # Fewer genes with high epimutation.

## Plots for both gene expression and variability that is controlled for expression levels.
## The dispersion method helps control for the relationship between variability and average expression.
my_comparisons <- list( c("(-0.01,0.1]", "(0.1,0.2]"), c("(-0.01,0.1]", "(0.3,0.4]"), c("(0.1,0.2]", "(0.4,1.1]") )

## Calculate the summary values across groups:
rna_var_wt_epi %>% 
  group_by(epimut_groups) %>% 
  summarise(avg_expression = median(mvp.mean, na.rm = T), 
            avg_dispersion = median(mvp.dispersion.scaled, na.rm = T))

ggplot(rna_var_wt_epi, aes(x=epimut_groups, y=mvp.dispersion.scaled)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x="Promoter epimutation groups", y="Dispersion (mean-scaled)") +
  plot_theme + 
  stat_compare_means(comparisons = my_comparisons)  + 
  stat_n_text(y.pos = -3.25) 

ggplot(rna_var_wt_epi, aes(x=epimut_groups, y=mvp.mean)) + 
  geom_boxplot() +
  stat_compare_means(method="kruskal.test") + 
  labs(x="Promoter epimutation groups", y="Mean expression") +
  plot_theme + 
  stat_n_text(y.pos = 1)


###########################
### Combine plots (color w/subtype designation)
###########################
rna_var_all_epi = rna_var_wt_epi %>%
  inner_join(rna_var_mut_epi, by=c("gene_id", "tumor_pdr", "epimut_groups"))
colnames(rna_var_all_epi) <- gsub("\\.x", ".IDHwt", colnames(rna_var_all_epi))
colnames(rna_var_all_epi) <- gsub("\\.y", ".IDHmut", colnames(rna_var_all_epi))

rna_var_all_epi_express = rna_var_all_epi %>% 
  gather(IDHstatus, mean_expression, c("mvp.mean.IDHwt", "mvp.mean.IDHmut")) %>% 
  mutate(IDHstatus = recode(IDHstatus, "mvp.mean.IDHwt" = "IDHwt", "mvp.mean.IDHmut" = "IDHmut"))


pdf(file = "results/Fig2/Fig2a-violin-expression-disorder-restricted-violin.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
ggplot(rna_var_all_epi_express, aes(x=epimut_groups, y=mean_expression, fill=IDHstatus)) + 
  geom_violin(width=1.1) +
  geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
  ylim(0, 0.9) + 
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Promoter DNAme disorder groups", y="Mean expression log(cpm)", fill="Glioma subtype") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free")
dev.off()

pdf(file = "github/results/Fig2/Fig2a-boxplot-expression-disorder.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(rna_var_all_epi_express, aes(x=epimut_groups, y=mean_expression, fill=IDHstatus)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Promoter DNAme disorder groups", y="Mean expression log(cpm)", fill="Glioma subtype") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free")
dev.off()


### END ####