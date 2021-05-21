##################################
# Analyze single cell relationship between DNAme disorder and SCNA
# Updated: 2021.05.14
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(EnvStats)
library(RColorBrewer)
library(gplots)
library(grDevices)
library(GenomicRanges)
library(openxlsx)
library(ggpubr)
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

####################################################
## Combine DNAme disorder with SCNA
####################################################
# Load the SCGP subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
meta_data = meta_data %>%
  mutate(case_barcode = gsub("-", "", substr(subject_id, 6, 11)))

## Final epiallele information.
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)

# Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number. Do not include polyploid cells.
tumor_epimut <- epiallele_info %>% 
  filter(tumor_cnv == 1, fga < 0.5) %>% 
  #filter(tumor_cnv == 1, ploidy < 3) %>% 
  mutate(case_barcode = gsub("-", "", substr(sample, 6, 11))) %>% 
  left_join(meta_data, by="case_barcode") %>% 
  mutate(IDH_status = ifelse(subtype=="IDHwt", "IDHwt", "IDHmut"))
cor.test(tumor_epimut$fga, tumor_epimut$PDR, method = "s")
cor.test(tumor_epimut$fga, tumor_epimut$promoter_PDR, method = "s")
cor.test(tumor_epimut$fga, tumor_epimut$cgi_PDR, method = "s")


# Specify the order of cases.
#case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
#tumor_epimut <- tumor_epimut %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-SCNA.pdf", width = 8, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  ylim(0.25, 0.45) +
  labs(y = "Global Epimutation burden", x= "SCNA burden", color = "Subject")
dev.off()

pdf(file = "github/results/Fig5/Fig5a-promoter-disorder-SCNA.pdf", width = 7, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(tumor_epimut, aes(x=fga, y = promoter_PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  theme(panel.spacing.x = unit(1.5, "lines")) +
  geom_smooth(method='lm', se=FALSE) +
  labs(y = "Promoter DNAme disorder (PDR)", x= "SCNA burden", color = "Subject")
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-promoter-SCNA-corr.pdf", width = 7, height = 4)
ggplot(tumor_epimut, aes(x=fga, y = promoter_PDR)) +
  geom_point(alpha=0.8) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman") + 
  geom_smooth(method='lm', se=FALSE)
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-cgi-SCNA-corr.pdf", width = 8, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = cgi_PDR)) +
  geom_point(alpha=0.8) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman") + 
  geom_smooth(method='lm', se=FALSE)
dev.off()


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-SCNA-ind-corr.pdf", width = 8, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=FALSE) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  #ylim(0.25, 0.5) +
  labs(y = "Epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman", )
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-cgi-scna-patient-corr.pdf", width = 9, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = cgi_PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=FALSE) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "CGI epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman")
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-global-scna-patient-corr.pdf", width = 9, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=FALSE) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Global epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman")
dev.off()



pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-promoter-scna-patient-corr.pdf", width = 9, height = 5)
ggplot(tumor_epimut, aes(x=fga, y = promoter_PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=FALSE) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Promoter epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman")
dev.off()

### END ###


