##################################
# Plot the CpG coverage and methylation across Alu and CpG island elements.
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

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Load the context-specific scRRBS DNAme disorder data.
disorder <- read.csv("data/analysis_scRRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)

## Load the context-specific epiallele DNA methylation data
meth <- read.csv("data/analysis_scRRBS_context_specific_methylation.csv", sep = ",", header = TRUE)
meth_wide <- meth %>% 
  dplyr::select(-c(case_barcode,num_cpgs)) %>% 
  pivot_wider(names_from = genomic_context, values_from = beta_value)

## Load the scRRBS qc data
rrbs_qc <- read.csv("data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

## Epiallele tumor cell data.
epiallele_tumor <- rrbs_qc %>%
  inner_join(meta_data, by="case_barcode") %>% 
  filter(tumor_status==1) %>%
  inner_join(disorder, by=c("cell_barcode", "case_barcode")) %>%
  inner_join(meth_wide, by="cell_barcode") %>% 
  mutate(idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut")) %>% 
  dplyr::select(cell_barcode, case_barcode, idh_status, cgi_PDR, cgi, alu_repeat_PDR, alu_repeat)

### alu *****
pdf("github/results/Fig1/Fig1f-alu-scatter.pdf", width = 5, height = 3.5)
ggplot(epiallele_tumor, aes(x=alu_repeat_PDR, y=alu_repeat, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "SM019" = "#FF63B6")) +
  guides(color=FALSE, fill=FALSE) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  stat_cor(method="spearman")
dev.off()

### cgi ****
pdf("results/Fig1//Fig1g-scatter.pdf", width = 5, height = 3.5)
ggplot(epiallele_tumor, aes(x=cgi_PDR, y=cgi, color=case_barcode)) + 
  geom_point(alpha=0.8) +
  labs(x="Mean DNAme disorder (PDR)", y="Mean DNA methylation", color="Subject") +
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "SM019" = "#FF63B6")) +
  geom_smooth(method = "lm", se=FALSE) + 
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid( ~ idh_status, scales = "free_y", space = "free") +
  stat_cor(method="spearman")
dev.off()

#### END ####