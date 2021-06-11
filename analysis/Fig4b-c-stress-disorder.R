##################################
# Assess the stress-associated changes in DNAme disorder/
# Updated: 2021.05.13
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github"
setwd(mybasedir)


###################################
# Necessary packages:
library(tidyverse)
library(parallel)
library(data.table)
library(ggpubr)
library(RColorBrewer)
library(EnvStats)
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

## Read in the DNAme disorder table and other sample data from Synapse download.
rrbs_20x_disorder <- read.csv("data/analysis_RRBS_context_specific_DNAme_disorder.csv", sep = ",", header = TRUE)
rrbs_qc <- read.csv("data/analysis_RRBS_sequencing_qc.csv", sep = ",", header = TRUE)

## Filter to the proportion of discordant read estimates that were calculated at locations with 20x coverage.
rrbs_20x_pdr_plot <- rrbs_qc %>% 
  inner_join(rrbs_20x_disorder, by="sample_barcode") %>% 
  dplyr::select(sample_barcode:rt_level, ends_with('PDR')) %>% 
  gather(feature_name, PDR, ends_with('PDR')) %>% 
  mutate(feature_name = gsub("_PDR", "", feature_name)) %>% 
  filter(feature_name%in%c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer"))

## Restrict to normoxia and hypoxia samples (exclude irradiated samples).
rrbs_20x_hypoxia_pdr_plot <- rrbs_20x_pdr_plot %>% 
  filter(!rt_level=="10Gy")

## Define factors.
oxygen_levels = c("Normoxia", "Hypoxia")
feature_levels <- c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer")
rrbs_20x_hypoxia_pdr_plot$feature_name <-  factor(rrbs_20x_hypoxia_pdr_plot$feature_name, levels = feature_levels)
rrbs_20x_hypoxia_pdr_plot$o2_level <-  factor(rrbs_20x_hypoxia_pdr_plot$o2_level, levels = oxygen_levels)

## Determine the relative DNAme disorder comparing hypoxia to normoxia replicateds.
rrbs_20x_hypoxia_pdr_plot_norm <- rrbs_20x_hypoxia_pdr_plot %>% 
  group_by(cell_line, timepoint, feature_name) %>%
  mutate_each(funs(./median(.[o2_level == "Normoxia"])), PDR)

## Filter to the genomic regions of interest that have been presented across the paper.
rrbs_20x_hypoxia_pdr_plot_norm_filt = rrbs_20x_hypoxia_pdr_plot_norm %>% 
  filter(feature_name%in%c("cgi", "promoter", "intergenic", "cgi_shore"))
feature_levels_filt <- c("cgi", "promoter", "intergenic", "cgi_shore")
rrbs_20x_hypoxia_pdr_plot_norm_filt$feature_name <-  factor(rrbs_20x_hypoxia_pdr_plot_norm_filt$feature_name, levels = feature_levels_filt)

pdf(file = "results/Fig4/Fig4b-relative-disorder-hypoxia.pdf", width = 6, height = 4, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt, aes(x=feature_name, y = PDR, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative local DNAme disorder", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(cell_line ~ timepoint) +
  stat_compare_means(label = "p.format", size=1.5)
dev.off()

pdf(file = "results/Fig4/Fig4b-relative-disorder-hypoxia-no-legend.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt, aes(x=feature_name, y = PDR, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative DNAme disorder\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  guides(fill=FALSE) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(cell_line ~ timepoint) +
  stat_compare_means(label = "p.format", size = 1.5)
dev.off()

## Filter to 3-day and present the oxygen dose concentrations for the same selected features.
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d <- rrbs_20x_hypoxia_pdr_plot_norm_filt %>% 
  filter(timepoint == "3 days") 

## Set the factor ordering.
o2_levels <- c("21%", "02%", "01%")
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc <-  factor(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc, levels = o2_levels)
my_comparisons <- list(c("21%", "01%"), c("21%", "02%"), c("01%", "02%"))

pdf(file = "results/Fig4/SuppFig8b-hypoxia-dosage.pdf", width = 6, height = 4, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d, aes(x=feature_name, y = PDR, fill=o2_conc)) +
  geom_boxplot() +
  labs(y = "Relative DNAme disorder\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                             "02%" = "#fdae61",  
                             "01%" = "#d7191c")) + 
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~.) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_compare_means(label = "p.format" , method = "kruskal.test", size = 1.5) 
dev.off()

###################
### Irradiation ###
###################
## Remove any 3-day timepoint data and restrict to normoxia (i.e., exclude 9-day hypoxia).
rrbs_20x_rt_pdr_plot <- rrbs_20x_pdr_plot %>% 
  filter(timepoint=="9 days", o2_conc =="21%")

feature_levels <- c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer")
rrbs_20x_rt_pdr_plot$feature_name <-  factor(rrbs_20x_rt_pdr_plot$feature_name, levels = feature_levels)
treatment_levels = c("0Gy", "10Gy")
rrbs_20x_rt_pdr_plot$rt_level <-  factor(rrbs_20x_rt_pdr_plot$rt_level, levels = treatment_levels)

## For each cell line determine the relative DNAme disorder (irradiation to normoxia).
rrbs_20x_rt_pdr_plot_norm <- rrbs_20x_rt_pdr_plot %>% 
  group_by(cell_line, timepoint, feature_name) %>%
  mutate_each(funs(./median(.[rt_level == "0Gy"])), PDR)

## Restrict to the features of interest that have been used elsewhere throughout the manuscript.
rrbs_20x_rt_pdr_plot_norm_filt = rrbs_20x_rt_pdr_plot_norm %>% 
  filter(feature_name%in%c("cgi", "promoter", "intergenic", "cgi_shore"))
feature_levels_filt <- c("cgi", "promoter", "intergenic",  "cgi_shore")
rrbs_20x_rt_pdr_plot_norm_filt$feature_name <-  factor(rrbs_20x_rt_pdr_plot_norm_filt$feature_name, levels = feature_levels_filt)

## Create the same plot for the relative PDR estimates for irradiation.
ggobj= ggplot(rrbs_20x_rt_pdr_plot_norm_filt, aes(x=feature_name, y = PDR, fill=rt_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative DNAme disorder\n(median 0 Gy normalized values)", x = "", fill = "Irradiation\n(dose)") +
  scale_fill_manual(values=c("0Gy" = "#abd9e9",
                             "10Gy" = "#0dba86")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  theme(legend.position="bottom") +
  facet_grid(cell_line ~ timepoint) +
  ylim(0.9, 1.25) +
  stat_compare_means(label = "p.format", size=1.5)
layer_scales(ggobj)$y$range$range

pdf(file = "results/Fig4/Fig4c-irradiation-stress.pdf", width = 3, height = 4, useDingbats = FALSE, bg="transparent")
ggobj
dev.off()

pdf(file = "results/Fig4/Fig4c-irradiation-stress-no-legend.pdf", width = 3, height = 5, useDingbats = FALSE, bg="transparent")
ggobj + guides(fill=FALSE)
dev.off()


### END ###