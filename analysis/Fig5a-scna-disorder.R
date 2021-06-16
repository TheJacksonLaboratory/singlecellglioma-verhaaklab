##################################
# Analyze single cell relationship between DNAme disorder and SCNA
# Updated: 2021.05.14
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
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
### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Load in the DNAme disorder table. 
epiallele_info <- read.table(file="data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)

## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

tumor_epiallele_info <- epiallele_info %>% 
  inner_join(rrbs_qc, by=c("cell_barcode", "case_barcode")) %>% 
  left_join(meta_data, by="case_barcode") %>% 
  filter(tumor_status == 1, sc_somatic_copy_number_alt_burden < 0.5) %>%
  mutate(IDH_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))
cor.test(tumor_epiallele_info$sc_somatic_copy_number_alt_burden, tumor_epiallele_info$PDR, method = "s")
cor.test(tumor_epiallele_info$sc_somatic_copy_number_alt_burden, tumor_epiallele_info$promoter_PDR, method = "s")
cor.test(tumor_epiallele_info$sc_somatic_copy_number_alt_burden, tumor_epiallele_info$cgi_PDR, method = "s")

## Plot global PDR across individual SCNA values for each sample.
ggplot(tumor_epiallele_info, aes(x=sc_somatic_copy_number_alt_burden, y = PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "SM019" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  ylim(0.25, 0.45) +
  labs(y = "Global Epimutation burden", x= "SCNA burden", color = "Subject")

pdf(file = "results/Fig5/Fig5a-promoter-disorder-SCNA.pdf", width = 7, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(tumor_epiallele_info, aes(x=sc_somatic_copy_number_alt_burden, y = promoter_PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "SM019" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  theme(panel.spacing.x = unit(1.5, "lines")) +
  geom_smooth(method='lm', se=FALSE) +
  labs(y = "Promoter DNAme disorder (PDR)", x= "SCNA burden", color = "Subject")
dev.off()

pdf(file = "results/Fig5/subtype-level-promoter-DNAme-SCNA-corr.pdf", width = 7, height = 4)
ggplot(tumor_epiallele_info, aes(x=sc_somatic_copy_number_alt_burden, y = promoter_PDR)) +
  geom_point(alpha=0.8) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman") + 
  geom_smooth(method='lm', se=FALSE)
dev.off()


pdf(file = "results/Fig5/promoter-DNAme-disorder-SCNA-ind-corr.pdf", width = 8, height = 5)
ggplot(tumor_epiallele_info, aes(x=sc_somatic_copy_number_alt_burden, y = promoter_PDR, color = case_barcode)) +
  geom_point(alpha=0.8) +
  geom_smooth(method='lm', se=FALSE) +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "SM019" = "#FF63B6")) +
  facet_grid(~IDH_status, scales = "free", space="free") +
  plot_theme + 
  labs(y = "Promoter epimutation burden", x= "SCNA burden", color = "Subject") +
  stat_cor(method = "spearman")
dev.off()

### END ###