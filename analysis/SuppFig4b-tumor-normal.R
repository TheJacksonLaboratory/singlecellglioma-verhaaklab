##################################
# Plot the DNAme and methylation across genomic elements
# Updated: 2021.05.14
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

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95) %>% 
  left_join(meta_data, by="case_barcode") 


## Load in the DNAme disorder table. 
epiallele_info <- read.table(file="data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)

## Create a new variable to indicate "non_tumor" or shortened case barcode.
epiallele_all = epiallele_info %>% 
  inner_join(rrbs_qc_pass, by=c("cell_barcode", "case_barcode")) %>% 
  mutate(case_normal_barcode = ifelse(tumor_status==0, "Non-tumor", case_barcode),
         idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"),
         normal_tumor = ifelse(tumor_status==0, "non-tumor", "tumor")) 
epiallele_all$normal_tumor[epiallele_all$normal_tumor=="tumor"&epiallele_all$idh_status=="IDHmut"] <- "IDHmut"
epiallele_all$normal_tumor[epiallele_all$normal_tumor=="tumor"&epiallele_all$idh_status=="IDHwt"] <- "IDHwt"

epiallele_all_pdr = epiallele_all %>% 
  select(cell_barcode, case_barcode, idh_status, normal_tumor, ends_with('PDR')) %>% 
  gather(genomic_context, disorder, ends_with('PDR')) 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="PDR"] <- "global" 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="h1hesc_ctcf_2_PDR"] <- "ctcf_esc" 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="nha_ctcf_PDR"] <- "ctcf_nha" 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="h1hesc_ezh2_PDR"] <- "ezh2_esc" 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="nha_ezh2_PDR"] <- "ezh2_nha" 
epiallele_all_pdr$genomic_context[epiallele_all_pdr$genomic_context=="alu_repeat_PDR"] <- "alu" 

epiallele_all_pdr$genomic_context <- gsub("_PDR", "", epiallele_all_pdr$genomic_context)
epiallele_all_pdr_filt = epiallele_all_pdr %>% 
  filter(genomic_context%in%c("alu", "intergenic", "cgi_shore", "gene_body", "dnaseI", "ezh2_esc", "ctcf_esc", "ezh2_nha", "ctcf_nha", "cgi",  "promoter"))


feature_order <- c("alu", "intergenic", "cgi_shore", "gene_body", "dnaseI", "ezh2_esc", "ctcf_esc", "ezh2_nha", "ctcf_nha", "cgi",  "promoter")
epiallele_all_pdr_filt$genomic_context <- factor(epiallele_all_pdr_filt$genomic_context, levels = feature_order)
tissue_order <- c("non-tumor", "IDHmut", "IDHwt")
epiallele_all_pdr_filt$normal_tumor <- factor(epiallele_all_pdr_filt$normal_tumor, levels = tissue_order)

pdf("results/Fig1/SuppFig4b-disorder-normal-vs-tumor.pdf", width = 7, height = 5, useDingbats = FALSE)
ggplot(epiallele_all_pdr_filt, aes(x=genomic_context, y=disorder, fill=normal_tumor)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x="", y="DNAme disorder (PDR)", fill="cell type") +
  plot_theme +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("IDHmut" = "#AF8DC3", 
                             "IDHwt" = "#7FBF7B", 
                             "non-tumor" = "#808080"))
dev.off()


## Load the context-specific epiallele DNA methylation data.
meth <- read.csv("data/analysis_scRRBS_context_specific_methylation.csv", sep = ",", header = TRUE)
meth$genomic_context[meth$genomic_context=="h1hesc_ctcf_2"] <- "ctcf_esc" 
meth$genomic_context[meth$genomic_context=="nha_ctcf"] <- "ctcf_nha" 
meth$genomic_context[meth$genomic_context=="h1hesc_ezh2"] <- "ezh2_esc" 
meth$genomic_context[meth$genomic_context=="nha_ezh2"] <- "ezh2_nha" 
meth$genomic_context[meth$genomic_context=="alu_repeat"] <- "alu" 

meth_filt <- meth %>% 
  filter(genomic_context%in%c("alu","intergenic", "cgi_shore", "ezh2_nha", "ezh2_esc", "dnaseI", "ctcf_nha", "ctcf_esc", "cgi", "promoter")) %>% 
  dplyr::select(-c(case_barcode, num_cpgs, case_barcode)) 

## A data.frame containing both DNAme disorder and DNA methylation quantified across epialleles.
epiallele_meth <- epiallele_all_pdr_filt %>% 
  inner_join(meth_filt, by=c("cell_barcode", "genomic_context")) %>% 
  filter(normal_tumor!="non-tumor")


####
## Is there a region that has a higher MAD methylation per subject?
####
summarized_mad <- epiallele_meth %>%
  group_by(case_barcode, genomic_context, idh_status) %>% 
  summarise(mad_epimut = mad(disorder),
            mad_meth = mad(beta_value)) %>% 
  filter(genomic_context%in% c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi"))

feature_order <- c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi")
summarized_mad$genomic_context <- factor(summarized_mad$genomic_context, levels = feature_order)

pdf("results/Fig1/SuppFig4-intratumoral-variation.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(summarized_mad, aes(x=genomic_context, y=mad_meth, fill=idh_status)) +
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
summarized_mad_inter <- epiallele_meth %>%
  group_by(genomic_context, idh_status) %>% 
  summarise(mad_epimut = mad(disorder),
            mad_meth = mad(beta_value)) %>% 
  filter(genomic_context%in% c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi"))

feature_order <- c("alu","intergenic", "cgi_shore", "ezh2_esc", "ezh2_nha", "dnaseI", "ctcf_esc", "ctcf_nha", "cgi")
summarized_mad_inter$genomic_context <- factor(summarized_mad_inter$genomic_context, levels = feature_order)

pdf("results/Fig1/SuppFig4-intertumoral-variation.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(summarized_mad_inter, aes(x=genomic_context, y=mad_meth, fill=idh_status)) +
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
  group_by(genomic_context) %>% 
  summarise(feature_med = median(mad_meth))
summarized_mad_alu_idhmut_inter <- summarized_mad_inter %>% 
  filter(idh_status=="IDHwt") %>% 
  group_by(genomic_context)

## Calculate for IDHwt tumors.
summary(summarized_mad_alu_idhmut_inter$mad_meth/median(summarized_mad_idhwt_intra$feature_med))

## Across all different features:
feature_meth_order <- c("alu","intergenic", "cgi_shore", "ezh2_nha", "ezh2_esc", "dnaseI", "ctcf_nha", "ctcf_esc", "cgi", "promoter")
epiallele_meth$genomic_context <- factor(epiallele_meth$genomic_context, levels = feature_meth_order)

pdf("results/Fig1/SuppFig4-DNAme-disorder-contexts.pdf", width = 10, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(epiallele_meth, aes(x=disorder, y=beta_value, color=idh_status)) + 
  geom_point(alpha = 0.8) +
  labs(x="DNAme disorder (PDR)", y="Mean DNAme", color="IDHmut status") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position="bottom") +
  ylim(0.1, 0.9) +
  scale_color_manual(values=c("IDHmut" = "#AF8DC3", 
                              "IDHwt" = "#7FBF7B")) +
  facet_grid(idh_status ~ genomic_context, scale = "free") +
  guides(fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE, color="gray25") + 
  stat_cor(method = "spearman", size=2)
dev.off()

### END ###