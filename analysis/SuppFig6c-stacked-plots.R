##################################
# Visualize specific genes for Supplementary Figure 6b - cluster gene expression violin plots.
# Updated: 2021.05.10
# Author: Kevin J.
###################################

mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
library(openxlsx)
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

## Load in pre-processed 10X data for all samples.
load("10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## 2D UMAP coordinates from cluster identification.
umap_coords_2d = read.csv("10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_umap_2d_embedding.csv")

## Change to HUGO gene name.
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are the 3D UMAP coordinates.
## Extract the numeric sample identifier.
tsne.data$sample_id = sapply(strsplit(tsne.data$cell_name, "-"), "[[", 3)

## Combine with the 2D results.
tsne.data = tsne.data %>% 
  left_join(umap_coords_2d, by=c("cell_name" = "barcode"))

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
clust_annot = tsne.data %>% 
  mutate(cell_type = recode(dbCluster, `1` = "Diff.-like",  `2` = "Myeloid", `3` = "Stem-like",
                            `4` = "Oligodendrocyte", `5` = "Prolif. stem-like", `6` = "Granulocyte", `7` = "Endothelial",
                            `8` = "T cell", `9` = "Pericyte", `10` = "Fibroblast", `11` = "B cell", `12` = "Dendritic cell")) 

# Create sample-specific labels for each patient to make it easier to refer back.
clust_annot$sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
clust_annot$sample_id <- gsub("^10$", "SM018", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^0$", "SM019", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^1$", "SM001", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^2$", "SM002", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^3$", "SM004", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^4$", "SM006", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^5$", "SM011", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^8$", "SM015", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^6$", "SM008", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^7$", "SM012", clust_annot$sample_id)
clust_annot$sample_id <- gsub("^9$", "SM017", clust_annot$sample_id)

###############################################
## Create a stacked barplot for non-tumor cells
###############################################
non_tumor_clust_annot = clust_annot %>%
  filter(!cell_type%in%c("Stem-like", "Diff.-like", "Prolif. stem-like")) %>% 
  group_by(sample_id, cell_type) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(idh_status = ifelse(sample_id%in%c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008"), "IDHmut", "IDHwt"))

case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
cell_state_order <- c("Oligodendrocyte", "Fibroblast", "Pericyte", "Endothelial", "B cell", "T cell", "Granulocyte", "Dendritic cell", "Myeloid")

non_tumor_clust_annot <- non_tumor_clust_annot %>% mutate(sample_id = factor(sample_id, levels = case_order))
non_tumor_clust_annot <- non_tumor_clust_annot %>% mutate(cell_type = factor(cell_type, levels = cell_state_order))


pdf("github/results/Fig3/SuppFig6c-non-tumor-states.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(non_tumor_clust_annot, aes(x = sample_id, y = freq, fill=cell_type)) +
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
