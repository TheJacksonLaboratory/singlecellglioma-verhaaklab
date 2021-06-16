##############################################
# Analyze bulk SCGP samples for ssGSEA enrichment scores for stress signalling
# Updated: 2021.05.18
# Author: Kevin J.
###############################################

## Working directory for this analysis in the SCGP-analysis project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###############################################
## Load the necessary packages.
library(tidyverse)
library(openxlsx)
library(DBI)
library(GSVA)
library(ggpubr)
###############################################
## Theme #####
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

## Establish connection with the SCGP database.
con <- DBI::dbConnect(odbc::odbc(), "scgpDB")

## Load in some critical tables from the bulk SCGP RNA-sequencing analysis.
cases <- dbReadTable(con,  Id(schema="clinical", table="cases"))
transcriptional_subtype <- dbReadTable(con,  Id(schema="analysis", table="transcriptional_subtype"))
transcript_tpm <- dbReadTable(con,  Id(schema="analysis", table="transcript_tpm"))
gene_tpm <- dbReadTable(con,  Id(schema="analysis", table="gene_tpm"))

## Create a gene x sample matrix.
gene_df = gene_tpm %>% 
  dplyr::select(-est_counts) %>% 
  pivot_wider(names_from = aliquot_barcode, values_from = tpm)
gene_df <- as.data.frame(gene_df)


rownames(gene_df) <- gene_df$gene_symbol
gene_df$gene_symbol <- NULL
gene_mat <- as.matrix(gene_df)

## Load in the referent gene sets for stress, oxidative stress, hypoxia..
stress_genes <- read.table(file="data/GO_REGULATION_OF_RESPONSE_TO_STRESS_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
stress_genes_list = list(stress_genes$V1)
ox_stress_genes <- read.table(file="data/GO_RESPONSE_TO_OXIDATIVE_STRESS_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
ox_stress_genes_list = list(ox_stress_genes$V1)
harris_hypoxia_genes <- read.table(file="data/HARRIS_HYPOXIA_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
harris_hypoxia_list = list(harris_hypoxia_genes$V1)
ELVIDGE_hypoxia_genes <- read.table(file="data/ELVIDGE_HYPOA_UP-geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
ELVIDGE_hypoxia_list = list(ELVIDGE_hypoxia_genes$V1)

## Try the enrichment with a set of randomly selected genes.
set.seed(43)
random_genes = sample(rownames(gene_mat), 1500)
random_gene_list =  list(random_genes)

## Determine the ssGSEA enrichment scores for each gene_set.
stress_gsea = gsva(gene_mat, stress_genes_list, method="ssgsea")
ox_stress_gsea = gsva(gene_mat, ox_stress_genes_list, method="ssgsea")
harris_gsea = gsva(gene_mat, harris_hypoxia_list, method="ssgsea")
elvidge_gsea = gsva(gene_mat, ELVIDGE_hypoxia_list, method="ssgsea")
random_gsea = gsva(gene_mat, random_gene_list, method="ssgsea")

cor.test(stress_gsea, ox_stress_gsea, method="spearman")
cor.test(stress_gsea, harris_gsea, method="spearman")
cor.test(elvidge_gsea, harris_gsea, method="spearman")
cor.test(random_gsea, harris_gsea, method="spearman")
cor.test(random_gsea, elvidge_gsea, method="spearman")
cor.test(random_gsea, stress_gsea, method="spearman")

#############################
### Global DNAme disorder
#############################
## Load in the DNAme disorder table. 
epiallele_info <- read.table(file="data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) 

## Calculate the mean disorder for each tumor.
disorder_mean <- epiallele_info %>% 
  filter(cell_barcode%in%rrbs_qc_pass$cell_barcode) %>% 
  group_by(case_barcode) %>% 
  summarise(disorder_mean = mean(PDR)) %>% 
  ungroup() %>% 
  select(case_barcode, disorder = disorder_mean) %>% 
  filter(!case_barcode%in%c("SM015", "SM008", "SM019"))

## Sanity check that data are aligned.
gsub("-", "", substr(colnames(stress_gsea), 6, 11))==disorder_mean$case_barcode

## They appear to be tightly correlated. 
cor.test(disorder_mean$disorder, as.numeric(stress_gsea), method="s")
cor.test(disorder_mean$disorder, as.numeric(harris_gsea), method="s")

## Combine into a single data.frame for interpretation of results.
epimut_stress = cbind(disorder_mean, as.numeric(harris_gsea), as.numeric(stress_gsea), as.numeric(ox_stress_gsea), as.numeric(random_gsea))
colnames(epimut_stress)[3:6] = c("HARRIS_HYPOXIA_ssGSEA", "GO_RESPONSE_TO_STRESS_ssGSEA", "GO_RESPONSE_TO_OX_STRESS_ssGSEA", "RANDOM_ssGSEA")
epimut_stress$idh_status <- ifelse(epimut_stress$case_barcode%in%c("SM001", "SM002", "SM004"), "IDHmut", "IDHwt")


epimut_stress_long = epimut_stress %>% 
  gather(pathway, ssGSEA, c(HARRIS_HYPOXIA_ssGSEA, GO_RESPONSE_TO_OX_STRESS_ssGSEA, RANDOM_ssGSEA)) %>% 
  mutate(pathway = recode(pathway, "HARRIS_HYPOXIA_ssGSEA" = "Hypoxia (Harris) genes",
                          "GO_RESPONSE_TO_OX_STRESS_ssGSEA" = "GO ox. stress response genes",
                          "RANDOM_ssGSEA" = "Random genes"))

pdf("results/Fig2/SuppFig5disorder-RNA-ssGSEA.pdf", height=4, width=7.5, useDingbats = FALSE, bg="transparent")
ggplot(epimut_stress_long, aes(x=ssGSEA, y=disorder)) + 
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(method="s") + 
  labs(x="ssGSEA enrichment score", y="Mean scRRBS DNAme disorder") +
  plot_theme +
  theme(panel.spacing = unit(1.5, "lines")) +
  facet_grid(. ~ pathway, scales = "free_x", space = "free")
dev.off()



### END ###