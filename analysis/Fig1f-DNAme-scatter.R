##################################
# Plot the CpG coverage and methylation across Alu and CpG island elements.
# Updated: 2021.05.12
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


## Restrict to specific features and plot by sample:
## CGI
comb_meth_epimut_cgi <- comb_meth_epimut %>% 
  filter(feature_name=="cgi")

pdf("github/results/Fig1//Fig1f-g-scatter.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_cgi, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free")
dev.off()

### alu *****
comb_meth_epimut_alu <- comb_meth_epimut %>% 
  filter(feature_name=="alu")

pdf("github/results/Fig1/Fig1f-alu-scatter.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_alu, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  guides(color=FALSE, fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free")
dev.off()


### CpG islands
comb_meth_epimut_cgi <- comb_meth_epimut %>% 
  filter(feature_name=="cgi")

pdf("github/results/Fig1/Fig1f-g-legend.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_cgi, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free")
dev.off()

pdf("github/results/Fig1/Fig1g-alu-legend.pdf", width = 5, height = 3.5)
ggplot(comb_meth_epimut_cgi, aes(x=epimutation_burden, y=mean_meth, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  guides(color=FALSE, fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free")
dev.off()



### END ###