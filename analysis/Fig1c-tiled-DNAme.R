##################################
# Generate boxplots for tiled 10kb DNA methylation.
# Updated: 2021.05.12
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
tiles_10kb <- read.table("Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
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

## Restrict the average methylation to only those cells passing qc.
idh_mut <- c("SM001", "SM002", "SM004", "SM008", "SM015", "SM019")
codel <- c("SM001", "SM004")
broad_meth_pass = broad_meth %>% 
  filter(cell_barcode%in%rrbs_qc_pass_ordered$cell_barcode) %>% 
  mutate(idh_codel_subtype = ifelse(!case_barcode%in%idh_mut, "IDHwt", ifelse(case_barcode%in%codel, "IDHmut-codel", "IDHmut-noncodel")),
         idh_status = ifelse(!case_barcode%in%idh_mut, "IDHwt", "IDHmut"))
case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
broad_meth_pass <- broad_meth_pass %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

pdf(file = "results/Fig1/Fig1c-tiled-meth.pdf", height = 4, width = 6, useDingbats = FALSE, bg="transparent")
ggplot(broad_meth_pass, aes(x=case_barcode,  y = mean_methylation_10kb, fill= idh_status)) + 
  geom_boxplot(outlier.shape = NA)  +
  labs(y="10-kb tiled DNA methylation", x="", fill="Glioma subtype") +
  scale_fill_manual(values=c("IDHwt" = "#7FBF7B", 
                             "IDHmut" = "#AF8DC3")) +
  ylim(0.40, 0.65) +
  plot_theme  +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(~idh_status, scales = "free", space="free_x")
dev.off()


## 2.398488e-115 between IDHmut and IDHwt tumors.
wilcox.test(broad_meth_pass$mean_methylation_10kb~broad_meth_pass$idh_status, alternative = "two.sided")$p.value

#### END ####
