##################################
# Determine whether there is a difference in DNAme + DNAme disorder based on ecDNA.
# Updated: 2020.05.06
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)


###################################
# Necessary packages:
library(tidyverse)
library(RColorBrewer)
library(openxlsx)
library(gplots)
library(grDevices)
library(GenomicRanges)
library(ggpubr)
library(EnvStats)
###################################
## Plotting theme:
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




# Load in meta data.
meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

### Epimutation
epimut_cpg <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epimut_cpg = epimut_cpg %>% 
  select(sample_barcode, pdr = PDR)

## Additional information about single-cells passing QC.
rrbs_qc_all <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)
rrbs_qc = rrbs_qc_all %>% 
  filter(tumor_status==1) %>% 
inner_join(epimut_cpg, by=c("cell_barcode" ="sample_barcode"))

## Determine which cells possess evidence for ecDNA.
## Samples with ecDNA: SM006, SM012, SM017, and SM018.


#############
## SM006  ##
#############
### 1Mb variable bins as processed by Ginkgo.
SM006 <- read.table("data/SM006-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
SM006 <- SM006[-which(SM006$CHR%in%c("chrX", "chrY")), ]
colnames(SM006) <- gsub("\\.", "-", colnames(SM006))

## Load in the QC data for scRRBS:
rrbs_SM006 = rrbs_qc %>% 
  filter(case_barcode == "SM006", cell_barcode%in%colnames(SM006), bisulfite_conversion_rate > 95, cpg_unique > 40000) %>% 
  arrange(cell_barcode)

# Subset the SM006 CNV data to the cells with passing QC.
index <- colnames(SM006)%in%rrbs_SM006$cell_barcode
index[1:3] <- c(TRUE, TRUE, TRUE)
SM006_filtered <- SM006[ ,index]

# Check to make sure the variables are in the same order.
all(rrbs_SM006$cell_barcode == colnames(SM006_filtered)[4:85])

# What about EGFR copy number??? ecDNA???
egfr_SM006 <- as.numeric(SM006_filtered[1153, 4:85])
hist(egfr_SM006)

# Add EGFR copy number to metadata.
rrbs_SM006$egfr_cn <- egfr_SM006
cor.test(rrbs_SM006$egfr_cn, rrbs_SM006$pdr)
ggplot(rrbs_SM006, aes(x=egfr_cn, y=pdr)) + geom_point()

rrbs_SM006$ecDNA <- ifelse(rrbs_SM006$egfr_cn > 6, "ecDNA+", "ecDNA-")
ggplot(rrbs_SM006, aes(x=ecDNA, y=pdr)) + geom_boxplot()
# There is a significant difference in epimutation rate.
wilcox.test(rrbs_SM006$pdr~rrbs_SM006$ecDNA)

#############
## SM012  ##
#############
### 1Mb variable bins as processed by Ginkgo.
SM012 <- read.table("data/SM012-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
SM012 <- SM012[-which(SM012$CHR%in%c("chrX", "chrY")), ]
colnames(SM012) <- gsub("\\.", "-", colnames(SM012))

rrbs_SM012 = rrbs_qc %>% 
  filter(case_barcode == "SM012", cell_barcode%in%colnames(SM012), bisulfite_conversion_rate > 95, cpg_unique > 40000) %>% 
  arrange(cell_barcode)

# Subset the SM012 CNV data to the cells with passing QC.
index <- colnames(SM012)%in%rrbs_SM012$cell_barcode
index[1:3] <- c(TRUE, TRUE, TRUE)
SM012_filtered <- SM012[ ,index]

# Check to make sure the variables are in the same order.
all(rrbs_SM012$cell_barcode == colnames(SM012_filtered)[4:72])

# What about EGFR copy number??? ecDNA???
egfr_SM012 <- as.numeric(SM012_filtered[1153, 4:72])
hist(egfr_SM012)
myc_SM012 <- as.numeric(SM012_filtered[1358, 4:72])
hist(myc_SM012)
pdgfra_SM012 <- as.numeric(SM012_filtered[655, 4:72])
hist(pdgfra_SM012)

# Add EGFR copy number to metadata.
rrbs_SM012$egfr_cn <- egfr_SM012
cor.test(rrbs_SM012$egfr_cn, rrbs_SM012$pdr)
ggplot(rrbs_SM012, aes(x=egfr_cn, y=pdr)) + geom_point()

rrbs_SM012$ecDNA <- ifelse(rrbs_SM012$egfr_cn > 6, "ecDNA+", "ecDNA-")
ggplot(rrbs_SM012, aes(x=ecDNA, y=pdr)) + geom_boxplot()
# There is a significant difference in epimutation rate.
wilcox.test(rrbs_SM012$pdr~rrbs_SM012$ecDNA)


#############
## SM017  ##
#############
### 1Mb variable bins as processed by Ginkgo.
SM017 <- read.table("data/SM017-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
SM017 <- SM017[-which(SM017$CHR%in%c("chrX", "chrY")), ]
colnames(SM017) <- gsub("\\.", "-", colnames(SM017))

rrbs_SM017 = rrbs_qc %>% 
  filter(case_barcode == "SM017", cell_barcode%in%colnames(SM017), bisulfite_conversion_rate > 95, cpg_unique > 40000) %>% 
  arrange(cell_barcode)

# Subset the SM017 CNV data to the cells with passing QC.
index <- colnames(SM017)%in%rrbs_SM017$cell_barcode
index[1:3] <- c(TRUE, TRUE, TRUE)
SM017_filtered <- SM017[ ,index]

# Check to make sure the variables are in the same order.
all(rrbs_SM017$cell_barcode == colnames(SM017_filtered)[4:94])

# What about EGFR copy number??? ecDNA???
egfr_SM017 <- as.numeric(SM017_filtered[1153, 4:94])
hist(egfr_SM017)

# Add EGFR copy number to metadata.
rrbs_SM017$egfr_cn <- egfr_SM017
cor.test(rrbs_SM017$egfr_cn, rrbs_SM017$pdr)
ggplot(rrbs_SM017, aes(x=egfr_cn, y=pdr)) + geom_point()

rrbs_SM017$ecDNA <- ifelse(rrbs_SM017$egfr_cn > 6, "ecDNA+", "ecDNA-")
ggplot(rrbs_SM017, aes(x=ecDNA, y=pdr)) + geom_boxplot()
# There is a significant difference in epimutation rate.
wilcox.test(rrbs_SM017$pdr~rrbs_SM017$ecDNA)


#############
## SM018  ##
#############
### 1Mb variable bins as processed by Ginkgo.
SM018 <- read.table("data/SM018-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
SM018 <- SM018[-which(SM018$CHR%in%c("chrX", "chrY")), ]
colnames(SM018) <- gsub("\\.", "-", colnames(SM018))

rrbs_SM018 = rrbs_qc %>% 
  filter(case_barcode == "SM018", cell_barcode%in%colnames(SM018), bisulfite_conversion_rate > 95, cpg_unique > 40000) %>% 
  arrange(cell_barcode)

# Subset the SM018 CNV data to the cells with passing QC.
index <- colnames(SM018)%in%rrbs_SM018$cell_barcode
index[1:3] <- c(TRUE, TRUE, TRUE)
SM018_filtered <- SM018[ ,index]

# Check to make sure the variables are in the same order.
all(rrbs_SM018$cell_barcode == colnames(SM018_filtered)[4:59])

# What about EGFR copy number??? ecDNA???
egfr_SM018 <- as.numeric(SM018_filtered[1153, 4:59])
hist(egfr_SM018)
mdm4_SM018 <- as.numeric(SM018_filtered[164, 4:59])
hist(mdm4_SM018)

# Add EGFR copy number to metadata.
rrbs_SM018$egfr_cn <- egfr_SM018
cor.test(rrbs_SM018$egfr_cn, rrbs_SM018$pdr)
ggplot(rrbs_SM018, aes(x=egfr_cn, y=pdr)) + geom_point()
rrbs_SM018$ecDNA <- ifelse(rrbs_SM018$egfr_cn > 6, "ecDNA+", "ecDNA-")
ggplot(rrbs_SM018, aes(x=ecDNA, y=pdr)) + geom_boxplot()
# There is a significant difference in epimutation rate.
wilcox.test(rrbs_SM018$pdr~rrbs_SM018$ecDNA)

#######################################
## Final image of ecDNA count per cell
#######################################
## Histograms and stacked barplots with greater than *6* copies of oncogene (assumption of tetraploidy).
SM006_sub <- rrbs_SM006 %>% 
  mutate(case_barcode = gsub("-","", substr(cell_barcode, 6, 11))) %>% 
  select(cell_barcode, case_barcode, pdr, egfr_cn, ecDNA)
SM012_sub <- rrbs_SM012 %>% 
  mutate(case_barcode = gsub("-","", substr(cell_barcode, 6, 11))) %>% 
  select(cell_barcode, case_barcode, pdr, egfr_cn, ecDNA)
SM017_sub <- rrbs_SM017 %>% 
  mutate(case_barcode = gsub("-","", substr(cell_barcode, 6, 11))) %>% 
  select(cell_barcode, case_barcode, pdr, egfr_cn, ecDNA)
SM018_sub <- rrbs_SM018 %>% 
  mutate(case_barcode = gsub("-","", substr(cell_barcode, 6, 11))) %>% 
  select(cell_barcode, case_barcode, pdr, egfr_cn, ecDNA)


## Summarized DNA methylation data across the 914 cells passing QC.
tiles_10kb = read.table("data/Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(tiles_10kb) <- gsub("\\.", "-",  colnames(tiles_10kb))

## TILED DNA methylation.
## Create a mean methylation value across the 10-kb tiled region.
bin_counts_10kb <- colSums(!is.na(tiles_10kb[, 5:918]))
bin_meth_10kb <- colMeans(tiles_10kb[, 5:918], na.rm = TRUE)
broad_meth <- as.data.frame(t(bind_rows(bin_counts_10kb, bin_meth_10kb)))
colnames(broad_meth) <- c("bins_covered_10kb", "mean_methylation_10kb")
broad_meth$cell_barcode <- rownames(broad_meth)
broad_meth$case_barcode <- gsub("-", "", substr(broad_meth$cell_barcode, 6, 11))

## Combine all data into one plot.
plot_ecDNA <- bind_rows(SM006_sub, SM012_sub, SM017_sub, SM018_sub)
plot_ecDNA = plot_ecDNA %>% 
  inner_join(broad_meth, by=c("cell_barcode", "case_barcode"))

pdf(file = "results/Fig6/SuppFig10f-ecDNA-distribution.pdf", width = 7, height = 4)
ggplot(plot_ecDNA, aes(x = egfr_cn, fill=ecDNA)) + 
  geom_histogram(binwidth = 1) + 
  facet_grid(~case_barcode, scales = "fixed", space="free") + 
  geom_vline(xintercept = 6, alpha=0.8, linetype=2, col="red") +
  labs(y="Cell count", x="EGFR copy number state") +
  scale_fill_manual(name='ecDNA status', values=c('ecDNA+'='#377eb8', "ecDNA-"= "gray")) +
  plot_theme
dev.off()
  
## Create a plot of general DNAme levels across ecDNA+/- cells.
gg_ecDNA_meth = ggplot(plot_ecDNA, aes(x=ecDNA,  y=mean_methylation_10kb, fill=ecDNA)) + 
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  ylab("DNAme (10-kb tiles)") + xlab("cell-specific ecDNA status") +
  scale_fill_manual(name='ecDNA status', values=c('ecDNA+'='#377eb8', "ecDNA-"= "gray")) +
  plot_theme +
  facet_grid(~case_barcode, scales = "free", space="free_x") + 
  guides(alpha=FALSE)

pdf(file = "results/Fig6/SuppFig10g-ecDNA-DNAme.pdf", width = 7, height = 4, useDingbats = FALSE)
gg_ecDNA_meth + stat_compare_means(method = "wilcox.test") +
  stat_n_text()
dev.off()

### END ####