##################################
# Stacked barplots for stress scRNAseq with Neftel classification
# Updated: 2021.05.12
# Author: Kevin J.
###################################

## Project working directory
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)


###################################
library(tidyverse)
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



###############################################
## Create a stacked barplot for tumor cells
###############################################
## Import table with Neftel classification.
hf2354_umap_annot_neftel <- read.table("data/neftel-classification-hf2354-20210205.txt", sep="\t", header=T, stringsAsFactors = F)

hf2354_tumor_class = hf2354_umap_annot_neftel %>%
  group_by(condition_revalue, class) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(timepoint = ifelse(grepl("9d", condition_revalue), "9 days", "3 days"),
         exposure = sapply(strsplit(condition_revalue, "-"), "[[", 1),
         treatment = recode(exposure, `normoxia` = "Normoxia - 0Gy",
                            `hypoxia` = "Hypoxia",
                            `Irradiation` = "Normoxia - 10Gy"))

treatment_levels <- c("Normoxia - 0Gy", "Hypoxia", "Normoxia - 10Gy")
cell_state_order <- c("OPC", "NPC", "MES", "AC")

hf2354_tumor_class <- hf2354_tumor_class %>% mutate(treatment = factor(treatment, levels = treatment_levels))
hf2354_tumor_class <- hf2354_tumor_class %>% mutate(class = factor(class, levels = cell_state_order))

pdf("results/Fig4/SuppFig8i-stress-neftel.pdf", width = 7, height = 5, useDingbats = FALSE)
ggplot(hf2354_tumor_class, aes(x = treatment, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES" = "#d7191c", "AC" = "#fdae61",
                             "NPC" = "#abd9e9", "OPC" = "#2c7bb6")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()



###############################################
## Create a stacked bar plot for tumor cells
###############################################
## Import table with Neftel classification.
hf3016_umap_annot_neftel <- read.table("data/neftel-classification-hf3016-20210204.txt", sep="\t", header=T, stringsAsFactors = F)

hf3016_tumor_class = hf3016_umap_annot_neftel %>%
  group_by(condition_revalue, class) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(timepoint = ifelse(grepl("9d", condition_revalue), "9 days", "3 days"),
         exposure = sapply(strsplit(condition_revalue, "-"), "[[", 1),
         treatment = recode(exposure, `normoxia` = "Normoxia - 0Gy",
                            `hypoxia` = "Hypoxia",
                            `Irradiation` = "Normoxia - 10Gy"))

treatment_levels <- c("Normoxia - 0Gy", "Hypoxia", "Normoxia - 10Gy")
cell_state_order <- c("OPC", "NPC", "MES", "AC")

hf3016_tumor_class <- hf3016_tumor_class %>% mutate(treatment = factor(treatment, levels = treatment_levels))
hf3016_tumor_class <- hf3016_tumor_class %>% mutate(class = factor(class, levels = cell_state_order))

pdf("results/Fig4/SuppFig8j-stress-neftel.pdf", width = 7, height = 5, useDingbats = FALSE)
ggplot(hf3016_tumor_class, aes(x = treatment, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES" = "#d7191c", "AC" = "#fdae61",
                             "NPC" = "#abd9e9", "OPC" = "#2c7bb6")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()


### END ###