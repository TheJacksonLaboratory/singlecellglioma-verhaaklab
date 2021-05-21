##############################################
# Analysis of epimutation at different replication times
# Updated: 2020.06.01
# Author: Kevin J.
###############################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/mnt/verhaak-lab/scgp"
setwd(mybasedir)

###############################################
## Load the necessary packages.
library(tidyverse)
library(ggpubr)
library(openxlsx)
###############################################
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

## Analysis of epimutation at different replication times
## First, generate the context-specific PDR plot for reference.
# Load the SCGP subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
meta_data = meta_data %>%
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11)))

## Final epiallele information.
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)

## Generate a basis for context-specific epimutation. 
epiallele_context_plot = epiallele_info %>% 
  mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) %>% 
  left_join(meta_data, by=c("case_normal_barcode"="case_barcode")) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(idh_codel_subtype = ifelse(case_normal_barcode=="Non-tumor", "Non-tumor", idh_status)) %>% 
  gather(genomic_context, pdr, c(promoter_PDR, gene_body_PDR, intergenic_PDR, dnaseI_PDR)) 
## Formatting for plot.
epiallele_context_plot$genomic_context = gsub("dnaseI_PDR", "DNaseI", epiallele_context_plot$genomic_context)
epiallele_context_plot$genomic_context = gsub("promoter_PDR", "Promoter", epiallele_context_plot$genomic_context)
epiallele_context_plot$genomic_context = gsub("gene_body_PDR", "Gene body", epiallele_context_plot$genomic_context)
epiallele_context_plot$genomic_context = gsub("intergenic_PDR", "Intergenic", epiallele_context_plot$genomic_context)
epiallele_context_plot$genomic_context <- factor(epiallele_context_plot$genomic_context, levels = c("Promoter", "Gene body", "Intergenic", "DNaseI"))
subtype_order <- c("Non-tumor", "IDHmut", "IDHwt")
case_order <- c("Non-tumor","SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
epiallele_context_plot <- epiallele_context_plot %>% mutate(case_normal_barcode = factor(case_normal_barcode, levels = case_order))
epiallele_context_plot <- epiallele_context_plot %>% mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = subtype_order))

## Recreate a supplemental figure of epimutation by context.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/manuscript/figure-drafts/Fig1/SFig-epimutation-context.pdf", width = 9, height = 5)
ggplot(epiallele_context_plot, aes(x=case_normal_barcode, y = pdr, fill=idh_codel_subtype)) +
  geom_boxplot(outlier.shape = NA) +
  plot_theme +
  scale_fill_manual(name='Glioma subtype', values=c('Non-tumor'='gray', 'IDHmut'='#AF8DC3', "IDHwt"="#7FBF7B")) +
  labs(y = "Epimutation burden", x="") + 
  guides(alpha = FALSE, color = FALSE) +
  facet_grid(~genomic_context, scales = "free", space="free_x") 
dev.off()

## Load in the calculated promoter epimutation at different replication times.
reptime_SM001 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-001-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM002 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-002-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM004 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-004-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM006 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-006-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM008 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-008-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM011 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-011-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM012 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-012-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM015 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-015-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM017 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-017-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM018 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-018-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_UC917 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-UC-917-01D-S5M_cells.deduplicated-full_list-final_RT-specific_PDR.txt", sep="\t", header=T, stringsAsFactors = F)

## Create a single data.frame.
reptime_all <- bind_rows(reptime_SM001, reptime_SM002, reptime_SM004, reptime_SM006, reptime_SM008,
                         reptime_SM011, reptime_SM012, reptime_SM015, reptime_SM017, reptime_SM018, reptime_UC917)
## Strip the samples of unnecessary information.
reptime_all$sample <- gsub(".deduplicated", "",  reptime_all$sample)
reptime_all$sample <- gsub("_pe", "",  reptime_all$sample)

## I calculated/re-calculated the promoter PDR.
meta_data$`10X_id_short` <- as.character(meta_data$`10X_id_short`)
epiallele_info_rt = epiallele_info %>% 
  inner_join(reptime_all, by=c("sample_barcode"="sample")) %>% 
  mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) %>% 
  left_join(meta_data, by=c("case_normal_barcode"="case_barcode")) %>% 
  mutate(idh_codel_subtype = ifelse(case_normal_barcode=="Non-tumor", "Non-tumor", subtype)) %>% 
  select(sample_barcode, case_barcode = case_normal_barcode, promoter_fraction, promoter_old_PDR = promoter_PDR.x,  promoter_PDR = promoter_PDR.y, promoter_rt1_fraction:promoter_rt4_PDR, idh_codel_subtype)

## There seems to be a positive, but not a perfect correlation despite using the same input files.
cor.test(epiallele_info_rt$promoter_old_PDR, epiallele_info_rt$promoter_PDR, method="s")
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/epimutation-promoter-RepTime-issue.pdf", width = 7, height = 5)
ggplot(epiallele_info_rt, aes(x= promoter_old_PDR, y = promoter_PDR, color=idh_codel_subtype)) + 
  geom_point() +
  ylim(0.1, 0.35) + 
  xlim(0.1,0.35) +
  labs(x="promoter PDR (epimut table)", y="promoter PDR (RT)") +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)
dev.off()

## Visualizations of the epimutation-replication boxplots per sample/subtype.
prom_pdr_rt = epiallele_info_rt %>% 
  select(sample_barcode, case_barcode, idh_codel_subtype, promoter_rt1_PDR, promoter_rt2_PDR, promoter_rt3_PDR, promoter_rt4_PDR) %>% 
  gather(reptime, pdr, c(promoter_rt1_PDR, promoter_rt2_PDR, promoter_rt3_PDR, promoter_rt4_PDR)) %>% 
  mutate(reptime = recode(reptime, "promoter_rt1_PDR" = "very early", "promoter_rt2_PDR" = "early", "promoter_rt3_PDR" = "late", 
                          "promoter_rt4_PDR" = "very late")) %>% 
  mutate(idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  filter(case_barcode!="Non-tumor")
subtype_order <- c("Non-tumor", "IDHmut_codel", "IDHmut_noncodel", "IDHwt")
case_order <- c("Non-tumor","SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
prom_pdr_rt <- prom_pdr_rt %>% mutate(reptime = factor(reptime, levels = c("very early", "early", "late", "very late")))
prom_pdr_rt <- prom_pdr_rt %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
prom_pdr_rt <- prom_pdr_rt %>% mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = subtype_order))

## On a per sample level.
ggplot(prom_pdr_rt, aes(x= case_barcode, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme  +
  scale_fill_manual(name='Replication timing', values=c('very early'='#efedf5', 'early'='#9e9ac8', "late"= "#6a51a3", "very late"="#3f007d")) +
  labs(y = "Epimutation burden", x="") +
  facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
  
## Across the three subtypes.
ggplot(prom_pdr_rt, aes(x= idh_codel_subtype, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme +
  scale_fill_manual(name='Replication timing', values=c('very early'='#efedf5', 'early'='#9e9ac8', "late"= "#6a51a3", "very late"="#3f007d")) +
  labs(y = "Epimutation rate", x="") +
  stat_compare_means(method = "kruskal.test")

## Across the two larger categories of IDHmut vs. IDHwt.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/Epimutation-RepTime-Prom.pdf", width = 7, height = 5)
ggplot(prom_pdr_rt, aes(x= idh_status, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5)) +
  scale_fill_manual(name='Replication timing', values=c('very early'='#f1eef6', 'early'='#bdc9e1', "late"= "#74a9cf", "very late"="#0570b0")) +
  labs(y = "Epimutation burden (promoter)", x="") +
  stat_compare_means(method = "kruskal.test")
dev.off()
  
## Separate into IDHmut and IDHwt and perform correlation analysis.
prom_pdr_rt_wt = prom_pdr_rt %>% filter(idh_status=="IDHwt")
prom_pdr_rt_mut = prom_pdr_rt %>% filter(idh_status!="IDHwt")
cor.test(prom_pdr_rt$pdr, as.numeric(prom_pdr_rt$reptime), method="kendall")
cor.test(prom_pdr_rt_wt$pdr, as.numeric(prom_pdr_rt_wt$reptime), method="kendall")
cor.test(prom_pdr_rt_mut$pdr, as.numeric(prom_pdr_rt_mut$reptime), method="kendall")



####################
#### Gene body
###################
reptime_SM001 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-001-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM002 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-002-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM004 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-004-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM006 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-006-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM008 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-008-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM011 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-011-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM012 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-012-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM015 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-015-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM017 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-017-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_SM018 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-SM-018-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_UC917 <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/SCGP-UC-917-01D-S5M_cells.deduplicated-full_list-final_RT-genebody_PDR.txt", sep="\t", header=T, stringsAsFactors = F)
reptime_all <- bind_rows(reptime_SM001, reptime_SM002, reptime_SM004, reptime_SM006, reptime_SM008,
                         reptime_SM011, reptime_SM012, reptime_SM015, reptime_SM017, reptime_SM018, reptime_UC917)
reptime_all$sample <- gsub(".deduplicated", "",  reptime_all$sample)
reptime_all$sample <- gsub("_pe", "",  reptime_all$sample)

## Assess whether the re-calculated gene_body PDR equals the gold-standard.
epiallele_info_rt = epiallele_info %>% 
  inner_join(reptime_all, by=c("sample_barcode"="sample")) %>% 
  mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) %>% 
  left_join(meta_data, by=c("case_normal_barcode"="case_barcode")) %>% 
  mutate(idh_codel_subtype = ifelse(case_normal_barcode=="Non-tumor", "Non-tumor", subtype)) %>% 
  select(sample_barcode, case_barcode = case_normal_barcode, gene_body_fraction_reads, genebody_old_PDR = gene_body_PDR,  genebody_PDR, genebody_rt1_fraction:genebody_rt4_PDR, idh_codel_subtype)

## There seems to be a weaker association than expected.
cor.test(epiallele_info_rt$genebody_old_PDR, epiallele_info_rt$genebody_PDR, method="s")
ggplot(epiallele_info_rt, aes(x= genebody_old_PDR, y = genebody_PDR)) + geom_point()

## Prepare the gene_body replication timing and epimutation analysis for visualization.
genebody_pdr_rt = epiallele_info_rt %>% 
  select(sample_barcode, case_barcode, idh_codel_subtype, genebody_rt1_PDR, genebody_rt2_PDR, genebody_rt3_PDR, genebody_rt4_PDR) %>% 
  gather(reptime, pdr, c(genebody_rt1_PDR, genebody_rt2_PDR, genebody_rt3_PDR, genebody_rt4_PDR)) %>% 
  mutate(reptime = recode(reptime, "genebody_rt1_PDR" = "very early", "genebody_rt2_PDR" = "early", "genebody_rt3_PDR" = "late", 
                          "genebody_rt4_PDR" = "very late")) %>% 
  mutate(idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  filter(case_barcode!="Non-tumor")
subtype_order <- c("Non-tumor", "IDHmut_codel", "IDHmut_noncodel", "IDHwt")
case_order <- c("Non-tumor","SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
genebody_pdr_rt <- genebody_pdr_rt %>% mutate(reptime = factor(reptime, levels = c("very early", "early", "late", "very late")))
genebody_pdr_rt <- genebody_pdr_rt %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
genebody_pdr_rt <- genebody_pdr_rt %>% mutate(idh_codel_subtype = factor(idh_codel_subtype, levels = subtype_order))

## Separate by patient.
ggplot(genebody_pdr_rt, aes(x= case_barcode, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme  +
  scale_fill_manual(name='Replication timing', values=c('very early'='#efedf5', 'early'='#9e9ac8', "late"= "#6a51a3", "very late"="#3f007d")) +
  labs(y = "Epimutation rate", x="") +
  facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")

## Separate by idh_codel_subtype.
ggplot(genebody_pdr_rt, aes(x= idh_codel_subtype, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme +
  scale_fill_manual(name='Replication timing', values=c('very early'='#efedf5', 'early'='#9e9ac8', "late"= "#6a51a3", "very late"="#3f007d")) +
  stat_compare_means(method = "kruskal.test") +
  labs(y = "Epimutation rate", x="")

## Separate by IDHmut vs. IDHwt.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/Epimutation-RepTime-GeneBody.pdf", width = 7, height = 5)
ggplot(genebody_pdr_rt, aes(x= idh_status, y = pdr, fill=reptime)) + 
  geom_boxplot(outlier.shape = NA) +
  plot_theme +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5)) +
  scale_fill_manual(name='Replication timing', values=c('very early'='#f1eef6', 'early'='#bdc9e1', "late"= "#74a9cf", "very late"="#0570b0")) +
  stat_compare_means(method = "kruskal.test") +
  labs(y = "Epimutation burden (gene body)", x="")
dev.off()

## Separate into IDHmut and IDHwt and perform correlation analysis.
genebody_pdr_rt_wt = genebody_pdr_rt %>% filter(idh_status=="IDHwt")
genebody_pdr_rt_mut = genebody_pdr_rt %>% filter(idh_status!="IDHwt")
cor.test(genebody_pdr_rt$pdr, as.numeric(genebody_pdr_rt$reptime), method="kendall")
cor.test(genebody_pdr_rt_wt$pdr, as.numeric(genebody_pdr_rt_wt$reptime), method="kendall")
cor.test(genebody_pdr_rt_mut$pdr, as.numeric(genebody_pdr_rt_mut$reptime), method="kendall")

#### Combine the two metrics:
genebody_pdr_rt$region <- "Gene body"
prom_pdr_rt$region <- "Promoter"

comb_pdr_rt <- bind_rows(genebody_pdr_rt, prom_pdr_rt)
comb_pdr_rt$region <- factor(comb_pdr_rt$region, levels = c("Promoter", "Gene body"))


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig5/Fig5b-replication-timing.pdf", width = 7, height = 4)
ggplot(comb_pdr_rt, aes(x= idh_status, y = pdr, fill=reptime)) + 
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5)) +
  scale_fill_manual(name='Replication timing', values=c('very early'='#f1eef6', 'early'='#bdc9e1', "late"= "#74a9cf", "very late"="#0570b0")) +
  stat_compare_means(method = "kruskal.test") +
  labs(y = "DNAme disorder (PDR)", x="") +
  facet_grid(~ region, scales = "free_x", space = "free")
dev.off()

### END ###

