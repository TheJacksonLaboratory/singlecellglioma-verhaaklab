##################################
# Visualize specific genes for Supplementary Figure 6b - cluster gene expression violin plots.
# Updated: 2021.05.10
# Author: Kevin J.
###################################

mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(cowplot)
###################################

## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  axis.text.x = element_text(angle=45, hjust=1),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


###############################################
## Create a stacked barplot for non-tumor cells
###############################################
## 2D UMAP coordinates.
umap_metadata <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

non_tumor_clust_annot = umap_metadata %>%
  filter(!cell_state%in%c("Stem-like", "Diff.-like", "Prolif. stem-like")) %>% 
  group_by(case_barcode, cell_state) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(idh_status = ifelse(case_barcode%in%c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008"), "IDHmut", "IDHwt"))

case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
cell_state_order <- c("Oligodendrocyte", "Fibroblast", "Pericyte", "Endothelial", "B cell", "T cell", "Granulocyte", "Dendritic cell", "Myeloid")

non_tumor_clust_annot <- non_tumor_clust_annot %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
non_tumor_clust_annot <- non_tumor_clust_annot %>% mutate(cell_state = factor(cell_state, levels = cell_state_order))


pdf("results/Fig3/SuppFig6c-non-tumor-states.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(non_tumor_clust_annot, aes(x = case_barcode, y = freq, fill=cell_state)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion cell state", fill="Non-tumor cell state") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c")) +
  plot_theme +
  facet_grid(. ~ idh_status, scales = "free_x", space = "free")
dev.off()  


### END ####
