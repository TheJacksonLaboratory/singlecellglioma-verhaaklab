##################################
# Visualize the DNAme disorder and DNAme values across promoters
# Updated: 2020.04.30
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
## Load the essential packages.
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ggpubr)
library(viridis)
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
                       axis.text.y=element_blank())

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank())

## Load the SCGP subject-level metadata.
full_meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Need to extract subtype and IDHmut status.
meta_data = full_meta_data %>% 
  select(case_barcode = subject_id, subtype) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(case_barcode_short = gsub("-", "", substr(case_barcode, 6, 11)))


## Load in the processed promoter methylation data for each sample.
prom2kb <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-final_FANTOM5_gene_matched_promoter_methylation.txt", sep="\t", header=T, stringsAsFactors = F)

## Load in the data with the total number of epialleles as well as context specific PDR.
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epiallele_tumor = epiallele_info %>% 
  filter(tumor_cnv == 1) %>% 
  left_join(meta_data, by=c("sample"="case_barcode"))


## Remove the annotated regions (chromosome, start, end, strand, gene).
prom2kb_mat = as.matrix(prom2kb[ , 6:1127])
rownames(prom2kb_mat) = prom2kb[ ,5]
## The output for the averaged promoter methylation needs to be adjusted.
colnames(prom2kb_mat) <- gsub(".deduplicated", "",  colnames(prom2kb_mat))
colnames(prom2kb_mat) <- gsub("_pe", "",  colnames(prom2kb_mat))
colnames(prom2kb_mat) <- gsub("\\.", "-",  colnames(prom2kb_mat))

## Restrict only to those samples that passed qc.
prom2kb_mat_sc <- prom2kb_mat[ , colnames(prom2kb_mat) %in%epiallele_tumor$sample_barcode]

## Restrict to promoters that are measured in at least 10% of samples.
missingness <- rowSums(is.na(prom2kb_mat_sc)/dim(prom2kb_mat_sc)[2])
sum(missingness < 0.9)

## Which promoters to keep to ensure comparability?
prom_keep = which(missingness < 0.9) # regions with at least 10% of the samples have a measurement.
prom2kb_clust = prom2kb_mat_sc[prom_keep, ]

## Do the rows match?
all(colnames(prom2kb_clust)==epiallele_tumor$sample_barcode)

# Function for plotting densities by group.
source("/Users/johnsk/Documents/Single-Cell-DNAmethylation/scripts/quality-metrics/densityRRBS.R")
densityRRBS(prom2kb_clust, main="promoter methylation")

## Calculate the median, average promoter methylation for each cell.
epiallele_tumor$avg_prom_meth <- colMeans(prom2kb_clust, na.rm=TRUE)

####### Promoter epimutation ##########
## Promoter epimutation rates.
prom_epimut <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-passQC_single_cells_individual_promoter-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)

## Revise the sample names to match between the promoter methylation data.
colnames(prom_epimut) <- gsub("\\.", "-",  colnames(prom_epimut))

## Restrict to only tumor cells.
prom_epimut_tumor = prom_epimut[, colnames(prom_epimut)%in%epiallele_tumor$sample_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(prom_epimut_tumor)/dim(prom_epimut_tumor)[2])
## Identify the number of gene promoters measured in at least 10% of cells.
sum(missingness < 0.9)
## Which promoters to keep to ensure comparability?
epimut_prom_keep = which(missingness < 0.9) # regions with at least 10% of the samples have a measurement.
prom_epimut_tumor_filt = data.matrix(prom_epimut_tumor[epimut_prom_keep, ])

#########################
## Combine two datasets
#########################
## Sanity check:
hist(as.numeric(prom2kb_clust["PANK4",]))
hist(as.numeric(prom_epimut_tumor_filt["PANK4",]))

## Build a data.frame with promoter-level DNA methylation:
prom2kb_clust = as.data.frame(prom2kb_clust)
prom2kb_clust$gene = rownames(prom2kb_clust)
gg_prom_meth = prom2kb_clust %>% 
  gather(sample, methylation, `SCGP-SM-001-01D-S5M-01S5`:`SCGP-UC-917-01D-S5M-95S5`) %>% 
  select(sample, gene, methylation)

## Build a data.frame with promoter-level epimutation:
prom_epimut_tumor_filt = as.data.frame(prom_epimut_tumor_filt)
prom_epimut_tumor_filt$gene = rownames(prom_epimut_tumor_filt)
gg_prom_epimut = prom_epimut_tumor_filt %>% 
  gather(sample, epimutation, `SCGP-SM-001-01D-S5M-01S5`:`SCGP-UC-917-01D-S5M-95S5`) %>% 
  select(sample, gene, epimutation)

## Sanity check number 2:
gg_prom_meth %>% filter(gene=="PANK4")
gg_prom_epimut  %>% filter(gene=="PANK4")

## Combine these two regulatory matrices.
IDHsamples <- c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917")
gg_prom_reg = gg_prom_meth %>% 
  inner_join(gg_prom_epimut, by=c("sample", "gene")) %>% 
  filter(!is.na(epimutation)) %>% 
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11)),
         subtype = ifelse(case_barcode%in%IDHsamples, "IDHmut", "IDHwt"))

gg_prom_reg_sum <- gg_prom_reg %>% 
  group_by(gene, subtype) %>% 
  summarise(epim_avg = mean(epimutation, na.rm=TRUE),
            meth_avg = mean(methylation, na.rm=TRUE))

## Generate a promoter plot for epimutation vs. DNA methylation by subtype.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5b-promoter-disorder-DNAme.pdf", height = 4, width = 6, useDingbats=FALSE)
ggplot(gg_prom_reg_sum, aes(x=epim_avg, y=meth_avg) ) + 
  geom_point(alpha=0.5, size=0.75, aes(colour = meth_avg)) +
  scale_colour_viridis() +
  geom_smooth(method="lm", color="black") +
  plot_theme +
  facet_grid(~subtype) +
  stat_cor(method = "spearman") +
  labs(x="Mean promoter-level DNAme disorder", y="Mean promoter-level DNAme", color = "DNAme\nlevel")
dev.off()

png(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5b-promoter-disorder-DNAme.png", height = 4, width = 6)
p1 = ggplot(gg_prom_reg_sum, aes(x=epim_avg, y=meth_avg) ) + 
  geom_point(alpha=0.5, size=0.75, aes(colour = meth_avg)) +
  scale_colour_viridis() +
  geom_smooth(method="lm", color="black") +
  plot_theme +
  null_x +
  null_y +
  facet_grid(~subtype) +
  labs(x="", y="", color = "DNAme\nlevel") 
dev.off()

ggsave("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5b-promoter-disorder-DNAme.png",
       p1,
       width = 6,
       height = 4,
       dpi=300)

## How many distinct promoter regions?
n_distinct(gg_prom_reg_sum$gene)

## What's the exact p-value?
gg_prom_reg_sub_mut = gg_prom_reg_sum %>% 
  filter(subtype=="IDHmut")
cor.test(gg_prom_reg_sub_mut$meth_avg, gg_prom_reg_sub_mut$epim_avg, method="s")$p.value
gg_prom_reg_sub_wt = gg_prom_reg_sum %>% 
  filter(subtype=="IDHwt")
cor.test(gg_prom_reg_sub_wt$meth_avg, gg_prom_reg_sub_wt$epim_avg, method="s")$p.value

#####################################################
### Gene-body DNA methylation - epimutation analyses
#####################################################
## Load in the processed gene body methylation data for each sample.
genebody_meth <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-final_Ensembl_gene_body_methylation.txt", sep="\t", header=T, stringsAsFactors = F)

# Remove the annotated regions (chromosome, start, end, strand, gene) for clustering.
gene_body_mat = genebody_meth[ , 6:1127]
rownames(gene_body_mat) <- genebody_meth$Name

# Revise the sample names to match between the promoter methylation data.
colnames(gene_body_mat) <- gsub(".deduplicated", "",  colnames(gene_body_mat))
colnames(gene_body_mat) <- gsub("_pe", "",  colnames(gene_body_mat))
colnames(gene_body_mat) <- gsub("\\.", "-",  colnames(gene_body_mat))
gene_body_mat_tumor = gene_body_mat[ ,colnames(gene_body_mat)%in%epiallele_tumor$sample_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(gene_body_mat_tumor)/dim(gene_body_mat_tumor)[2])
## Identify the number of gene promoters measured in at least 10% of cells.
sum(missingness < 0.9)
## Which promoters to keep to ensure comparability?
meth_genes_keep = which(missingness < 0.9) # regions with at least 10% of the samples have a measurement.
genebody_tumor_filt = data.matrix(gene_body_mat_tumor[meth_genes_keep, ])

## Gene body epimutation rates.
gb_epimut <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-passQC_single_cells_individual_gene_body-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)

## Revise the sample names to match between the promoter methylation data.
colnames(gb_epimut) <- gsub("\\.", "-",  colnames(gb_epimut))

## Restrict to only tumor cells.
gb_epimut_tumor = gb_epimut[ ,colnames(gb_epimut)%in%epiallele_tumor$sample_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(gb_epimut_tumor)/dim(gb_epimut_tumor)[2])
## Identify the number of gene promoters measured in at least 10% of cells.
sum(missingness < 0.9)
## Which promoters to keep to ensure comparability?
epimut_genes_keep = which(missingness < 0.9) # regions with at least 10% of the samples have a measurement.
gb_epimut_tumor_filt = data.matrix(gb_epimut_tumor[epimut_genes_keep, ])


##############################
### Combine datasets
##############################
#########################
## Combine two datasets
#########################
## Sanity check:
hist(as.numeric(genebody_tumor_filt["ENSG00000146648",]))
hist(as.numeric(gb_epimut_tumor_filt["ENSG00000146648",]))

## Build a data.frame with promoter-level DNA methylation:
genebody_tumor_filt = as.data.frame(genebody_tumor_filt)
genebody_tumor_filt$gene = rownames(genebody_tumor_filt)
gg_genebody_meth = genebody_tumor_filt %>% 
  gather(sample, methylation, `SCGP-SM-001-01D-S5M-01S5`:`SCGP-UC-917-01D-S5M-95S5`) %>% 
  select(sample, gene, methylation)

## Build a data.frame with promoter-level epimutation:
gb_epimut_tumor_filt = as.data.frame(gb_epimut_tumor_filt)
gb_epimut_tumor_filt$gene = rownames(gb_epimut_tumor_filt)
gg_genebody_epimut = gb_epimut_tumor_filt %>% 
  gather(sample, epimutation, `SCGP-SM-001-01D-S5M-01S5`:`SCGP-UC-917-01D-S5M-95S5`) %>% 
  select(sample, gene, epimutation)

## Sanity check number 2:
gg_genebody_meth %>% filter(gene=="ENSG00000146648")
gg_genebody_epimut  %>% filter(gene=="ENSG00000146648")

## Combine these two regulatory matrices.
IDHsamples <- c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917")
gg_genebody_reg = gg_genebody_meth %>% 
  inner_join(gg_genebody_epimut, by=c("sample", "gene")) %>% 
  filter(!is.na(epimutation)) %>% 
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11)),
         subtype = ifelse(case_barcode%in%IDHsamples, "IDHmut", "IDHwt"))

## Distinct genes: 17132
gg_genebody_reg_sub = gg_genebody_reg %>% 
  group_by(gene, subtype) %>% 
  summarise(meth_avg = mean(methylation),
            epim_avg = mean(epimutation))

## Plot gene-body DNA methylation.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5b-genebody-disorder-DNAme.pdf", height = 4, width = 6, useDingbats=FALSE)
ggplot(gg_genebody_reg_sub, aes(x=epim_avg, y=meth_avg) ) + 
  geom_point(alpha=0.5, size=0.75, aes(colour = meth_avg)) +
  scale_colour_viridis() +
  geom_smooth(method="lm", color="black") +
  plot_theme +
  facet_grid(~subtype) +
  stat_cor(method = "spearman") +
  labs(x="Mean gene body-level DNAme disorder", y="Mean genebody-level DNAme", color = "DNAme\nlevel")
dev.off()


p2 = ggplot(gg_genebody_reg_sub, aes(x=epim_avg, y=meth_avg) ) + 
  geom_point(alpha=0.5, size=0.75, aes(colour = meth_avg)) +
  scale_colour_viridis() +
  geom_smooth(method="lm", color="black") +
  plot_theme +
  null_x +
  null_y +
  facet_grid(~subtype) +
  labs(x="", y="", color = "DNAme\nlevel") 

ggsave("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5b-genebody-disorder-DNAme.png",
       p2,
       width = 6,
       height = 4,
       dpi=300)


## How many distinct gene body regions?
n_distinct(gg_genebody_reg_sub$gene)


### END ###

