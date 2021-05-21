##################################
# Plot the DNAme and methylation across genomic elements
# Updated: 2021.05.14
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/"
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

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())

## Load the SCGP subject-level metadata.
full_meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Need to extract subtype and IDHmut status.
meta_data = full_meta_data %>% 
  select(case_barcode = subject_id, subtype) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(case_barcode_short = gsub("-", "", substr(case_barcode, 6, 11)))

## Load in the data with the total number of epialleles as well as context specific PDR.
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epiallele_info_set2 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context_set2-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epiallele_info_set2_clean <- epiallele_info_set2 %>% 
  select(sample_barcode, cgi_shore_epiallele_count:nha_ezh2_PDR)

epiallele_tumor = epiallele_info %>% 
  inner_join(epiallele_info_set2_clean, by="sample_barcode") %>% 
  filter(tumor_cnv == 1) %>% 
  left_join(meta_data, by=c("sample"="case_barcode")) %>% 
  select(sample_barcode, subtype:case_barcode_short, adapter:nha_ezh2_PDR)

epiallele_tumor_pdr = epiallele_tumor %>% 
  select(sample_barcode, case_barcode = case_barcode_short, idh_status, ends_with('PDR')) %>% 
  gather(feature_name, epimutation_burden, ends_with('PDR')) 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="PDR"] <- "global" 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="h1hesc_ctcf_2_PDR"] <- "ctcf_esc" 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="nha_ctcf_PDR"] <- "ctcf_nha" 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="h1hesc_ezh2_PDR"] <- "ezh2_esc" 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="nha_ezh2_PDR"] <- "ezh2_nha" 
epiallele_tumor_pdr$feature_name[epiallele_tumor_pdr$feature_name=="alu_repeat_PDR"] <- "alu" 


epiallele_tumor_pdr$feature_name <- gsub("_PDR", "", epiallele_tumor_pdr$feature_name)
epiallele_tumor_pdr_filt = epiallele_tumor_pdr %>% 
  filter(feature_name%in%c("global", "alu", "intergenic", "intron", "cgi_shore", "gene_body",
                           "ezh2_nha", "ezh2_esc", "ctcf_nha", "ctcf_esc", "dnaseI", "exon", "cgi", "promoter", "tss"))

case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
feature_order <- c("global", "alu", "intergenic", "intron", "cgi_shore", "gene_body", 
                   "ezh2_nha", "ezh2_esc", "dnaseI", "ctcf_nha", "ctcf_esc", "exon", "cgi", "promoter", "tss")
epiallele_tumor_pdr_filt$case_barcode <- factor(epiallele_tumor_pdr_filt$case_barcode, levels = case_order)
epiallele_tumor_pdr_filt$feature_name <- factor(epiallele_tumor_pdr_filt$feature_name, levels = feature_order)

epiallele_tumor_pdr_filt %>% 
  group_by(feature_name) %>%
  summarise(avg_epimutation = mean(epimutation_burden)) %>% 
  arrange(desc(avg_epimutation))


##############################
## DNA methylation
#############################

## Load in the processed methylation data for these samples and create aggregation (promoters, gene body, DNase)
## Aggregated measures:
shore_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_CGI_shore_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
shore_cpgs$Sample <- gsub("_pe.deduplicated", "", shore_cpgs$Sample)
shore_cpgs$case_barcode <- gsub("-", "", substr(shore_cpgs$Sample, 6, 11))
shore_cpgs$case_barcode <- factor(shore_cpgs$case_barcode, levels = case_order)


## CGI
cgi_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_CGI_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(cgi_cpgs) <- c("sample_barcode", "cgi_mean_meth", "cgi_cpgs")

## CGI shores
shore_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_CGI_shore_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(shore_cpgs) <- c("sample_barcode", "cgi_shore_mean_meth", "cgi_shore_cpgs")

## Promoters
prom_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_FANTOM5_gene_matched_promoter_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(prom_cpgs) <- c("sample_barcode", "promoter_mean_meth", "promoter_cpgs")

## Gene body
genebody_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_Ensembl_gene_body_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(genebody_cpgs) <- c("sample_barcode", "gene_body_mean_meth", "gene_body_cpgs")

## DNaseI hypersensitivity
dnase_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_DNaseI_hypersensitive_site_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(dnase_cpgs) <- c("sample_barcode", "dnaseI_mean_meth", "dnaseI_cpgs")

## Alu:
alu_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_Repeat-Alu_family_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(alu_cpgs) <- c("sample_barcode", "alu_mean_meth", "alu_cpgs")

## Other repeat elements:
l1_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_Repeat-L1_family_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(l1_cpgs) <- c("sample_barcode", "l1_mean_meth", "l1_cpgs")
l2_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_Repeat-L2_family_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(l2_cpgs) <- c("sample_barcode", "l2_mean_meth", "l2_cpgs")
mir_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_UCSC_Repeat-MIR_family_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(mir_cpgs) <- c("sample_barcode", "mir_mean_meth", "mir_cpgs")

## EZH2:
ezh2_nha_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_ENCODE_NHA_EZH2_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(ezh2_nha_cpgs) <- c("sample_barcode", "ezh2_nha_mean_meth", "ezh2_nha_cpgs")
ezh2_esc_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_ENCODE_H1HESC_EZH2_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(ezh2_esc_cpgs) <- c("sample_barcode", "ezh2_esc_mean_meth", "ezh2_esc_cpgs")

## CTCF:
ctcf_nha_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_ENCODE_NHA_CTCF_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(ctcf_nha_cpgs) <- c("sample_barcode", "ctcf_nha_mean_meth", "ctcf_nha_cpgs")
ctcf_esc_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_ENCODE_H1HESC_CTCF_1_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(ctcf_esc_cpgs) <- c("sample_barcode", "ctcf_esc_mean_meth", "ctcf_esc_cpgs")

## Intergenic:
intergenic_cpgs <- read.table("/Users/johnsk/mnt/verhaak-lab/scgp/results/meth_clustering/rerun-reformatted_deduplicated-final/Samples-final_Ensembl_intergenic_aggregated_methylation.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(intergenic_cpgs) <- c("sample_barcode", "intergenic_mean_meth", "intergenic_cpgs")


## Combine the data into a single data.frame:
feature_cpgs = shore_cpgs %>% 
  inner_join(cgi_cpgs, by="sample_barcode") %>%
  inner_join(prom_cpgs, by="sample_barcode") %>%
  inner_join(genebody_cpgs, by="sample_barcode") %>%
  inner_join(dnase_cpgs, by="sample_barcode") %>% 
  inner_join(intergenic_cpgs, by="sample_barcode") %>% 
  inner_join(alu_cpgs, by="sample_barcode") %>%
  inner_join(l1_cpgs, by="sample_barcode") %>%
  inner_join(l2_cpgs, by="sample_barcode") %>% 
  inner_join(mir_cpgs, by="sample_barcode") %>% 
  inner_join(ezh2_nha_cpgs, by="sample_barcode") %>% 
  inner_join(ezh2_esc_cpgs, by="sample_barcode") %>%
  inner_join(ctcf_nha_cpgs, by="sample_barcode") %>% 
  inner_join(ctcf_esc_cpgs, by="sample_barcode")
  

feature_cpgs$sample_barcode <- gsub("_pe.deduplicated", "", feature_cpgs$sample_barcode)
feature_cpgs$case_barcode <- gsub("-", "", substr(feature_cpgs$sample_barcode, 6, 11))
feature_cpgs$case_barcode <- factor(feature_cpgs$case_barcode, levels = case_order)

feature_cpgs_long <- feature_cpgs %>% 
  select(-ends_with("_mean_meth")) %>% 
  gather(feature_name, cpg_count, ends_with('_cpgs')) %>% 
  filter(sample_barcode%in%epiallele_tumor$sample_barcode)
feature_cpgs_long$feature_name <- gsub("_cpgs", "", feature_cpgs_long$feature_name)

feature_meth_long <- feature_cpgs %>% 
  select(-ends_with("_cpgs")) %>% 
  gather(feature_name, mean_meth, ends_with('_mean_meth'))  %>% 
  filter(sample_barcode%in%epiallele_tumor$sample_barcode)
feature_meth_long$feature_name <- gsub("_mean_meth", "", feature_meth_long$feature_name)

feature_long <- feature_cpgs_long %>% 
  inner_join(feature_meth_long, by=c("sample_barcode", "case_barcode", "feature_name")) %>% 
  mutate(idh_status = ifelse(case_barcode%in%c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917"), "IDHmut", "IDHwt"))
feature_long$case_barcode <- factor(feature_long$case_barcode, levels = case_order)
feature_order <- c("alu", "l1", "l2", "mir","cgi_shore", "gene_body", "intergenic", "dnaseI", "ezh2_nha", "ezh2_esc",
                   "ctcf_nha", "ctcf_esc", "cgi",  "promoter")
feature_long$feature_name <- factor(feature_long$feature_name, levels = feature_order)

###############
## Integrate the DNA methylation and epimutation data at the single cell level.
###############
## Combine the epimutation and DNA methylation values based on sample_barcode + feature_name.
table(epiallele_tumor_pdr_filt$feature_name)
table(feature_long$feature_name)

## Combine DNA methylation and epimutation for common features.
comb_meth_epimut <- feature_long %>% 
  inner_join(epiallele_tumor_pdr_filt, by=c("sample_barcode", "feature_name", "case_barcode", "idh_status")) 

case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
feature_order <- c("alu", "intergenic", "cgi_shore", "gene_body", "ezh2_nha", "ezh2_esc", "dnaseI", 
                   "ctcf_nha", "ctcf_esc", "cgi", "promoter")
comb_meth_epimut$case_barcode <- factor(comb_meth_epimut$case_barcode, levels = case_order)
comb_meth_epimut$feature_name <- factor(comb_meth_epimut$feature_name, levels = feature_order)

comb_meth_epimut_filt <- comb_meth_epimut %>% 
  filter(feature_name%in% c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi", "promoter"))

####
## Is there a region that has a higher MAD methylation per subject?
####
summarized_mad <- comb_meth_epimut %>%
  group_by(case_barcode, feature_name, idh_status) %>% 
  summarise(mad_epimut = mad(epimutation_burden),
            mad_meth = mad(mean_meth),
            mad_cpg = mad(cpg_count)) %>% 
  filter(feature_name%in% c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi"))

feature_order <- c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi")
summarized_mad$feature_name <- factor(summarized_mad$feature_name, levels = feature_order)

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig4-intratumoral-variation.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(summarized_mad, aes(x=feature_name, y=mad_meth, fill=idh_status)) +
  geom_boxplot(outlier.shape=NA) +
  labs(x="Genomic element", y="Median absolute deviation\nsingle-cell DNAme", fill="IDH status") +
  plot_theme +
  ylim(0, 0.08) +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B")) +
  ggtitle("Intra-patient DNAme variation") +
  stat_compare_means(method="wilcox", label = "p.signif") +
  theme(panel.spacing.x = unit(1.5, "lines"), 
        axis.text.x = element_text(angle=45, hjust=1))
dev.off()


###############
### INTER-tumor variation
###############
summarized_mad_inter <- comb_meth_epimut %>%
  group_by(feature_name, idh_status) %>% 
  summarise(mad_epimut = mad(epimutation_burden),
            mad_meth = mad(mean_meth),
            mad_cpg = mad(cpg_count)) %>% 
  filter(feature_name%in% c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi"))

feature_order <- c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi")
summarized_mad_inter$feature_name <- factor(summarized_mad_inter$feature_name, levels = feature_order)

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig4-intertumoral-variation.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(summarized_mad_inter, aes(x=feature_name, y=mad_meth, fill=idh_status)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x="Genomic element", y="Median absolute deviation\nsingle-cell DNAme", fill="IDH status") +
  plot_theme +
  ylim(0, 0.08) +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B")) +
  ggtitle("Inter-patient DNAme variation") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()


## What's the difference between INTER and INTRA:
summarized_mad_idhwt_intra <- summarized_mad %>% 
  filter(idh_status=="IDHwt") %>% 
  group_by(feature_name) %>% 
  summarise(feature_med = median(mad_meth))
summarized_mad_alu_idhmut_inter <- summarized_mad_inter %>% 
  filter(idh_status=="IDHwt") %>% 
  group_by(feature_name)
## Calculate for IDHwt tumors.
summary(summarized_mad_alu_idhmut_inter$mad_meth/median(summarized_mad_idhwt_intra$feature_med))

## Across all different features:
pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/methylation-epimutation-all-contexts.pdf", width = 10, height = 4)
ggplot(comb_meth_epimut_filt, aes(x=epimutation_burden, y=mean_meth)) + 
  geom_point() +
  labs(x="", y="Mean DNA methylation") +
  plot_theme +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B")) +
  facet_grid( ~ feature_name, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  stat_cor(method = "spearman")
dev.off()

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig4-DNAme-disorder-contexts.pdf", width = 10, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(comb_meth_epimut_filt, aes(x=epimutation_burden, y=mean_meth, color=idh_status)) + 
  geom_point(alpha = 0.8) +
  labs(x="DNAme disorder (PDR)", y="Mean DNAme", color="IDHmut status") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position="bottom") +
  ylim(0.1, 0.9) +
  scale_color_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B")) +
  facet_grid(idh_status ~ feature_name, scale = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE, color="gray25") + 
  stat_cor(method = "spearman", size=2)
dev.off()

### Feature count - SuppFig2f ####
pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig2f-feature-cpg-count.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(comb_meth_epimut, aes(x=feature_name, y=cpg_count)) + 
  geom_boxplot(outlier.shape=NA) +
  labs(x="Genomic element", y="CpG count per feature") +
  plot_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/feature-epimutation-value.pdf", width = 9, height = 5)
ggplot(comb_meth_epimut, aes(x=feature_name, y=epimutation_burden, fill=idh_status)) + 
  geom_boxplot(outlier.shape=NA) +
  labs(x="Genomic element", y="CpG count per feature", fill = "IDHmut status") +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B")) +
  plot_theme
dev.off()

comb_meth_epimut_cgi_chi_shores <- comb_meth_epimut %>% 
  filter(feature_name%in%c("cgi", "cgi_shore"))
ggplot(comb_meth_epimut_cgi_chi_shores, aes(x=feature_name, y=epimutation_burden)) + 
  geom_boxplot() +
  labs(x="Genomic element", y="Genomic context-specific\nepimutation burden") +
  plot_theme +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  stat_compare_means(method="wilcox")


#### Are there any differences in cell states by genomic context? ###
comb_meth_epimut_filt_liger <- comb_meth_epimut %>% 
  filter(feature_name%in% c("cgi_shore", "cgi")) %>% 
  filter(!case_barcode%in%c("SM004", "SM018", "SM011"))

ggplot(comb_meth_epimut_filt_liger, aes(x=case_barcode, y=epimutation_burden, fill=liger_cell_state)) +
  geom_boxplot(outlier.shape=NA) +
  labs(x="case", y="Epimutation burden", fill="LIGER cell state") +
  scale_fill_manual(values=c("stemcell" = "#fb6a4a", 
                             "differentiated" = "#fcbba1")) +
  plot_theme +
  guides(fill=FALSE) +
  facet_grid(~ feature_name, scales = "free_y", space = "free") +
  stat_compare_means(method="wilcox", label = "p.signif")

ggplot(comb_meth_epimut_filt_liger, aes(x=liger_cell_state, y=epimutation_burden, fill=liger_cell_state)) +
  geom_boxplot(outlier.shape=NA) +
  labs(x="case", y="Epimutation burden", fill="LIGER cell state") +
  scale_fill_manual(values=c("stemcell" = "#fb6a4a", 
                             "differentiated" = "#fcbba1")) +
  plot_theme +
  guides(fill=FALSE) +
  facet_grid(idh_status~ feature_name, scales = "free_y", space = "free") +
  stat_compare_means(method="wilcox", label = "p.signif")

ggplot(comb_meth_epimut_filt_liger, aes(x=case_barcode, y=mean_meth, fill=liger_cell_state)) +
  geom_boxplot(outlier.shape=NA) +
  labs(x="case", y="DNA methylation", fill="LIGER cell state") +
  scale_fill_manual(values=c("stemcell" = "#fb6a4a", 
                             "differentiated" = "#fcbba1")) +
  plot_theme +
  guides(fill=FALSE) +
  facet_grid( ~ feature_name, scales = "free_y", space = "free") +
  stat_compare_means(method="wilcox", label = "p.signif")



## Restrict to specific features and plot by sample:
## Promoter
comb_meth_epimut_promoter <- comb_meth_epimut %>% 
  filter(feature_name=="promoter")

ggplot(comb_meth_epimut_promoter, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point() +
  labs(x="Epimutation Burden Promoter", y="Mean DNA methylation\nPromoter") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")

## CGI
comb_meth_epimut_cgi <- comb_meth_epimut %>% 
  filter(feature_name=="cgi")

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/methylation-epimutation-cgi-paper.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_cgi, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE, color=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  #stat_cor(method = "spearman") +
  plot_theme 
dev.off()

### CpG island shore
comb_meth_epimut_cgi_shore <- comb_meth_epimut %>% 
  filter(feature_name=="cgi_shore")

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/methylation-epimutation-cgi-shore.pdf", width = 8, height = 6)
ggplot(comb_meth_epimut_cgi_shore, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Epimutation Burden CGI Shores", y="Mean DNA methylation\nCGI shore", color="Subject") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")
dev.off()

### intergenic
comb_meth_epimut_intergenic <- comb_meth_epimut %>% 
  filter(feature_name=="intergenic")

ggplot(comb_meth_epimut_intergenic, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point() +
  labs(x="Epimutation Burden Intergenic", y="Mean DNA methylation\nIntergenic") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")

### alu *****
comb_meth_epimut_alu <- comb_meth_epimut %>% 
  filter(feature_name=="alu")

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/methylation-epimutation-alu-paper.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_alu, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE, color=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  #stat_cor(method = "spearman") +
  plot_theme 
dev.off()

### DNaseI
comb_meth_epimut_dnaseI <- comb_meth_epimut %>% 
  filter(feature_name=="dnaseI")

ggplot(comb_meth_epimut_dnaseI, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point() +
  labs(x="Epimutation Burden DNaseI", y="Mean DNA methylation\nDNaseI") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")

### EZH2 ESC
comb_meth_epimut_ezh2_esc <- comb_meth_epimut %>% 
  filter(feature_name=="ezh2_esc")

ggplot(comb_meth_epimut_ezh2_esc, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point() +
  labs(x="Epimutation Burden EZH2 ESC", y="Mean DNA methylation\nEZH2 ESC") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")

### CTCF ESC
comb_meth_epimut_ctcf_esc <- comb_meth_epimut %>% 
  filter(feature_name=="ctcf_esc")

ggplot(comb_meth_epimut_ctcf_esc, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point() +
  labs(x="Epimutation Burden CTCF ESC", y="Mean DNA methylation\nCTCF ESC") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor(method = "spearman")

#################################################
### Epimutation and intratumoral DNAm variation
#################################################

## Is there a relationship between the level of epimutation burden for a region and DNA methylation variation.
epimut_meth_mad = comb_meth_epimut_filt %>% 
  group_by(case_barcode, feature_name) %>% 
  summarise(epimut_median = median(epimutation_burden),
            meth_mad = mad(mean_meth)) %>% 
  ungroup()

## Not very helpful
ggplot(epimut_meth_mad, aes(x=epimut_median, y=meth_mad, color=case_barcode)) + 
  geom_point() +
  labs(x="Median epimutation burden", y="Median Aboslute Deviation DNA methylation") +
  plot_theme +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid( ~ feature_name, scales = "fixed") 

ggplot(epimut_meth_mad, aes(x=epimut_median, y=meth_mad, color=feature_name)) + 
  geom_point() +
  labs(x="Median epimutation burden", y="Median Aboslute Deviation DNA methylation") +
  plot_theme 

ggplot(epimut_meth_mad, aes(x=epimut_median, y=meth_mad, color=feature_name)) + 
  geom_point() +
  labs(x="Median epimutation burden", y="Median Aboslute Deviation DNA methylation") +
  plot_theme +
  facet_grid( ~ case_barcode, scales = "fixed") +
  stat_cor(method = "spearman")

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/methylation-mad-epimutation.pdf", width = 9, height = 5)
ggplot(epimut_meth_mad, aes(x=epimut_median, y=meth_mad)) + 
  geom_point() +
  labs(x="Median epimutation burden", y="Median Aboslute Deviation DNA methylation") +
  plot_theme +
  facet_grid( ~ feature_name, scales = "fixed") +
  stat_cor(method = "spearman", size=2)
dev.off()

ggplot(epimut_meth_mad, aes(x=case_barcode, y=meth_mad)) + 
  geom_boxplot() +
  labs(x="Subject", y="Median Aboslute Deviation DNA methylation\n(different contexts)") +
  plot_theme  

ggplot(epimut_meth_mad, aes(x=feature_name, y=meth_mad)) + 
  geom_boxplot() +
  labs(x="Subject", y="Median Aboslute Deviation DNA methylation\n(different contexts)") +
  plot_theme + 
  facet_grid( ~ case_barcode, scales = "fixed") 


#####################################
##### Tumor vs. Normal Comparisons
#####################################
epiallele_all = epiallele_info %>% 
  inner_join(epiallele_info_set2_clean, by="sample_barcode") %>% 
  left_join(meta_data, by=c("sample"="case_barcode")) %>% 
  select(sample_barcode, subtype:case_barcode_short, adapter:nha_ezh2_PDR) %>% 
  mutate(normal_tumor = ifelse(tumor_cnv==0, "normal", "tumor"))
epiallele_all$normal_tumor[epiallele_all$normal_tumor=="tumor"&epiallele_all$idh_status=="IDHmut"] <- "IDHmut"
epiallele_all$normal_tumor[epiallele_all$normal_tumor=="tumor"&epiallele_all$idh_status=="IDHwt"] <- "IDHwt"


epiallele_all_pdr = epiallele_all %>% 
  select(sample_barcode, case_barcode = case_barcode_short, idh_status, normal_tumor, ends_with('PDR')) %>% 
  gather(feature_name, epimutation_burden, ends_with('PDR')) 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="PDR"] <- "global" 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="h1hesc_ctcf_2_PDR"] <- "ctcf_esc" 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="nha_ctcf_PDR"] <- "ctcf_nha" 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="h1hesc_ezh2_PDR"] <- "ezh2_esc" 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="nha_ezh2_PDR"] <- "ezh2_nha" 
epiallele_all_pdr$feature_name[epiallele_all_pdr$feature_name=="alu_repeat_PDR"] <- "alu" 


epiallele_all_pdr$feature_name <- gsub("_PDR", "", epiallele_all_pdr$feature_name)
epiallele_all_pdr_filt = epiallele_all_pdr %>% 
  filter(feature_name%in%c("global", "alu", "intergenic", "intron", "cgi_shore", "gene_body",
                           "ezh2_nha", "ezh2_esc", "ctcf_nha", "ctcf_esc", "dnaseI", "exon", "cgi", "promoter", "tss"))

case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
feature_order <- c("global", "alu", "intergenic", "intron", "cgi_shore", "gene_body", 
                   "ezh2_nha", "ezh2_esc", "dnaseI", "ctcf_nha", "ctcf_esc", "exon", "cgi", "promoter", "tss")
epiallele_all_pdr_filt$case_barcode <- factor(epiallele_all_pdr_filt$case_barcode, levels = case_order)
epiallele_all_pdr_filt$feature_name <- factor(epiallele_all_pdr_filt$feature_name, levels = feature_order)

epiallele_all_pdr_filt %>% 
  group_by(feature_name) %>%
  summarise(avg_epimutation = mean(epimutation_burden)) %>% 
  arrange(desc(avg_epimutation))

feature_cpgs_long <- feature_cpgs %>% 
  select(-ends_with("_mean_meth")) %>% 
  gather(feature_name, cpg_count, ends_with('_cpgs')) %>% 
  filter(sample_barcode%in%epiallele_all$sample_barcode)
feature_cpgs_long$feature_name <- gsub("_cpgs", "", feature_cpgs_long$feature_name)

feature_meth_long <- feature_cpgs %>% 
  select(-ends_with("_cpgs")) %>% 
  gather(feature_name, mean_meth, ends_with('_mean_meth'))  %>% 
  filter(sample_barcode%in%epiallele_all$sample_barcode)
feature_meth_long$feature_name <- gsub("_mean_meth", "", feature_meth_long$feature_name)

feature_long <- feature_cpgs_long %>% 
  inner_join(feature_meth_long, by=c("sample_barcode", "case_barcode", "feature_name")) %>% 
  mutate(idh_status = ifelse(case_barcode%in%c("SM001", "SM002", "SM004", "SM008", "SM015", "UC917"), "IDHmut", "IDHwt"))
feature_long$case_barcode <- factor(feature_long$case_barcode, levels = case_order)
feature_order <- c("alu", "l1", "l2", "mir","cgi_shore", "gene_body", "intergenic", "dnaseI", "ezh2_nha", "ezh2_esc",
                   "ctcf_nha", "ctcf_esc", "cgi",  "promoter")
feature_long$feature_name <- factor(feature_long$feature_name, levels = feature_order)


## Combine DNA methylation and epimutation for common features.
comb_meth_epimut_all <- feature_long %>% 
  inner_join(epiallele_all_pdr_filt, by=c("sample_barcode", "feature_name", "case_barcode", "idh_status")) 

feature_order <- c("alu", "intergenic", "cgi_shore", "gene_body", "dnaseI", "ezh2_esc", "ctcf_esc", "ezh2_nha", "ctcf_nha", "cgi",  "promoter")
comb_meth_epimut_all$feature_name <- factor(comb_meth_epimut_all$feature_name, levels = feature_order)
tissue_order <- c("normal", "IDHmut", "IDHwt")
comb_meth_epimut_all$normal_tumor <- factor(comb_meth_epimut_all$normal_tumor, levels = tissue_order)


pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig4b-disorder-normal-vs-tumor.pdf", width = 7, height = 5, useDingbats = FALSE)
ggplot(comb_meth_epimut_all, aes(x=feature_name, y=epimutation_burden, fill=normal_tumor)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x="", y="DNAme disorder (PDR)", fill="cell type") +
  plot_theme +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B", 
                             "normal" = "#808080"))
dev.off()

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig1/SuppFig4b-disorder-normal-vs-tumor-wilcox.pdf", width = 9, height = 5)
my_comparisons <- list( c("normal", "IDHmut"), c("normal", "IDHwt")) 
ggplot(comb_meth_epimut_all, aes(x=normal_tumor, y=epimutation_burden)) + 
  geom_boxplot() +
  labs(x="", y="Epimutation burden", fill="cell class") +
  plot_theme +
  facet_grid( ~ feature_name, scales = "free_y", space = "free") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox") 
dev.off()
#####################################
##### Boxplots of CpGs per feature
#####################################
feature_long_filt_plot = feature_long %>% 
  filter(feature_name%in%c("alu","intergenic", "cgi_shore", "gene_body", "dnaseI", "ezh2_esc", "ctcf_esc", "cgi",  "promoter"))

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/qc/average-cpg-per-feature.pdf", width = 7, height = 5)
ggplot(feature_long_filt_plot, aes(x=feature_name, y=cpg_count)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x="Feature", y="CpG count per feature") +
  plot_theme 
dev.off()


### END ###
