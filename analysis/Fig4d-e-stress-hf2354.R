##############################################
# Assess the differences in gene expression and cell states based on stress exposure for HF2354
# Updated: 2021.05.13
# Author: Kevin J.
###############################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###########################
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(Seurat)
library(EnvStats)
###########################

## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 10),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


## Load Scanpy filtered, processed, unnormalized count matrix.
hf2354_obj <- ReadH5AD("data/analysis_scRNAseq_stress_hf2354_expression.h5ad")

## Make a Seurat object with the standard pipeline through PCA.
## Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
## This is then natural-log transformed using log1p.
hf2354_obj <- hf2354_obj %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = hf2354_obj@var.genes, npcs = 20, verbose = FALSE)

## Visualize.
hf2354_obj <- hf2354_obj %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

## Quick Seurat visualization. Note that there will be a difference between Scanpy and Seurat clustering.
DimPlot(hf2354_obj) 

## Load the metadata including cell states and UMAP coordinates.
hf2354_metadata = read.csv("data/analysis_scRNAseq_stress_hf2354_metadata.csv", sep = ",", stringsAsFactors = F, header = T)

## Add metadata to Seurat object.
all(rownames(hf2354_obj@meta.data)==hf2354_metadata$cell_barcode)
hf2354_obj@meta.data$cell_barcode <- hf2354_metadata$cell_barcode
hf2354_obj@meta.data$condition <- hf2354_metadata$condition
hf2354_obj@meta.data$oxygen_levels <- hf2354_metadata$oxygen_levels
hf2354_obj@meta.data$condition_revalue <- hf2354_metadata$condition_revalue
hf2354_obj@meta.data$cell_type <- hf2354_metadata$cell_type
hf2354_obj@meta.data$exposure_state <- paste(hf2354_metadata$cell_type, hf2354_metadata$condition_revalue, sep="-")

## Change the UMAP coordinates to those from scanpy.
all(rownames(hf2354_obj@reductions$umap@cell.embeddings)==hf2354_metadata$cell_barcode)
hf2354_obj@reductions$umap@cell.embeddings[,"UMAP_1"] <- hf2354_metadata$UMAP_1
hf2354_obj@reductions$umap@cell.embeddings[,"UMAP_2"] <- hf2354_metadata$UMAP_2 

## First plot cell state:
Idents(hf2354_obj) <- hf2354_obj@meta.data$cell_type
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = state_levels)

pdf("results/Fig4/Fig4d-hf2354-scgp-umap.pdf", width = 4.5, height = 2.5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes()
dev.off()

png("results/Fig4/Fig4d-hf2354-scgp-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes() & NoLegend()
dev.off()


condition_levels <- c("normoxia-3d", "hypoxia-3d", "normoxia-9d", "hypoxia-9d", "Irradiation-9d")
Idents(hf2354_obj) <- hf2354_obj@meta.data$condition_revalue
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = condition_levels)

png("results/Fig4/Fig4d-hf2354-conditions-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes() & NoLegend()
dev.off()


pdf("results/Fig4/Fig4d-hf2354-conditions-umap-legend.pdf", width = 9, height = 5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes()
dev.off()


hf2354_neftel <- read.table("neftel-classification-hf2354-20210205.txt", header = T)
all(rownames(hf2354_obj@meta.data)==hf2354_neftel$cell_name)
hf2354_obj@meta.data$neftel_class <- paste(hf2354_neftel$class, "like", sep="-")
Idents(hf2354_obj) <- hf2354_obj@meta.data$neftel_class
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = c("OPC-like", "NPC-like", "MES-like", "AC-like"))

pdf("results/Fig4/Fig4d-hf2354-neftel-umap-legend.pdf", width = 4, height = 3, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes()
dev.off()

png("results/Fig4/Fig4d-hf2354-neftel-umap-no-legend.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes() & NoLegend()
dev.off()

## Create stacked barplots for Neftel cell classifications across the different experimental conditions.
hf2354_tumor_class = hf2354_metadata %>%
  inner_join(hf2354_neftel, by=c("cell_barcode" ="cell_name")) %>% 
  group_by(condition, class) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(timepoint = ifelse(grepl("-09", condition), "9 days", "3 days"),
         exposure = substr(condition, 10, 11),
         treatment = recode(exposure, `21` = "Normoxia - 0Gy",
                            `01` = "Hypoxia",
                            `02` = "Hypoxia",
                            `RT` = "Normoxia - 10Gy"), 
         class = paste(class, "like", sep="-"))

treatment_order <- c("Normoxia - 0Gy", "Hypoxia", "Normoxia - 10Gy")
hf2354_tumor_class <- hf2354_tumor_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
cell_state_order <- c("OPC-like", "NPC-like", "MES-like", "AC-like")
hf2354_tumor_class <- hf2354_tumor_class %>% mutate(class = factor(class, levels = cell_state_order))

pdf("results/Fig4/SuppFig8-stress-neftel-cellstates-hf2354.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf2354_tumor_class, aes(x = treatment, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                             "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()


#### Create stacked bar plots ####
hf2354_scgp_class = hf2354_metadata %>%
  group_by(condition, cell_type) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(timepoint = ifelse(grepl("-09", condition), "9 days", "3 days"),
         exposure = substr(condition, 10, 11),
         treatment = recode(exposure, `21` = "Normoxia - 0Gy",
                            `01` = "Hypoxia",
                            `02` = "Hypoxia",
                            `RT` = "Normoxia - 10Gy"))

treatment_order <- c("Normoxia - 0Gy", "Hypoxia", "Normoxia - 10Gy")
hf2354_scgp_class <- hf2354_scgp_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf2354_scgp_class <- hf2354_scgp_class %>% mutate(cell_type = factor(cell_type, levels = state_levels))

pdf("results/Fig4/Fig4e-stress-scgp-cellstates-hf2354.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf2354_scgp_class, aes(x = treatment, y = freq, fill=cell_type)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="SCGP\ncell state") +
  scale_fill_manual(values=c("Stem-like" = "#fb6a4a", 
                             "Diff.-like" = "#fcbba1", 
                             "Prolif. stem-like" = "#a50f15")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()

chisq.test(table(hf2354_metadata$cell_type[hf2354_metadata$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")], hf2354_metadata$exposure[hf2354_metadata$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")]))
chisq.test(table(hf2354_metadata$cell_type[hf2354_metadata$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")], hf2354_metadata$exposure[hf2354_metadata$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")]))
chisq.test(table(hf2354_metadata$cell_type[hf2354_metadata$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")], hf2354_metadata$exposure[hf2354_metadata$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")]))

### END ###