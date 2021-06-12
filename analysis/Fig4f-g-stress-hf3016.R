##############################################
# Assess the differences in gene expression and cell states based on stress exposure for HF3016
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
hf3016_obj <- ReadH5AD("data/analysis_scRNAseq_stress_hf3016_expression.h5ad")

## Make a Seurat object with the standard pipeline through PCA.
## Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
## This is then natural-log transformed using log1p.
hf3016_obj <- hf3016_obj %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = hf3016_obj@var.genes, npcs = 20, verbose = FALSE)

## Visualize.
hf3016_obj <- hf3016_obj %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

## Quick Seurat visualization. Note that there will be a difference between Scanpy and Seurat clustering.
DimPlot(hf3016_obj) 

## Load the metadata including cell states and UMAP coordinates.
hf3016_metadata = read.csv("/Users/johnsk/github/data/analysis_scRNAseq_stress_hf3016_metadata.csv", sep = ",", stringsAsFactors = F, header = T)

## Add metadata to Seurat object.
all(rownames(hf3016_obj@meta.data)==hf3016_metadata$cell_barcode)
hf3016_obj@meta.data$cell_barcode <- hf3016_metadata$cell_barcode
hf3016_obj@meta.data$condition <- hf3016_metadata$condition
hf3016_obj@meta.data$oxygen_levels <- hf3016_metadata$oxygen_levels
hf3016_obj@meta.data$condition_revalue <- hf3016_metadata$condition_revalue
hf3016_obj@meta.data$cell_type <- hf3016_metadata$cell_type
hf3016_obj@meta.data$exposure_state <- paste(hf3016_metadata$cell_type, hf3016_metadata$condition_revalue, sep="-")

## Change the UMAP coordinates to those from scanpy.
all(rownames(hf3016_obj@reductions$umap@cell.embeddings)==hf3016_metadata$cell_barcode)
hf3016_obj@reductions$umap@cell.embeddings[,"UMAP_1"] <- hf3016_metadata$UMAP_1
hf3016_obj@reductions$umap@cell.embeddings[,"UMAP_2"] <- hf3016_metadata$UMAP_2 

## First plot cell state:
Idents(hf3016_obj) <- hf3016_obj@meta.data$cell_type
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = state_levels)

pdf("results/Fig4/Fig4d-hf3016-scgp-umap.pdf", width = 4.5, height = 2.5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes()
dev.off()

png("results/Fig4/Fig4d-hf3016-scgp-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes() & NoLegend()
dev.off()

## Set the experimental conditions.
condition_levels <- c("normoxia-3d", "hypoxia-3d", "normoxia-9d", "hypoxia-9d", "Irradiation-9d")
Idents(hf3016_obj) <- hf3016_obj@meta.data$condition_revalue
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = condition_levels)

png("results/Fig4/Fig4d-hf3016-conditions-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes() & NoLegend()
dev.off()


pdf("results/Fig4/Fig4d-hf3016-conditions-umap-legend.pdf", width = 9, height = 5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf3016_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes()
dev.off()

## Assess the Neftel meta-module scores and classification. These can be calculated from gene expression tables on Synapse.
hf3016_neftel <- read.table("neftel-classification-hf3016-20210204.txt", header = T)
all(rownames(hf3016_obj@meta.data)==hf3016_neftel$cell_name)
## Add "-like" suffix to these cancer cell states.
hf3016_obj@meta.data$neftel_class <- paste(hf3016_neftel$class, "like", sep="-")
Idents(hf3016_obj) <- hf3016_obj@meta.data$neftel_class
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = c("OPC-like", "NPC-like", "MES-like", "AC-like"))

pdf("results/Fig4/Fig4d-hf3016-neftel-umap-legend.pdf", width = 4, height = 3, bg = "transparent", useDingbats = FALSE)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes()
dev.off()

png("results/Fig4/Fig4d-hf3016-neftel-umap-no-legend.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes() & NoLegend()
dev.off()

## Plot the proportions of the Neftel classification.
hf3016_tumor_class = hf3016_metadata %>%
  inner_join(hf3016_neftel, by=c("cell_barcode" ="cell_name.x", "condition")) %>% 
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
hf3016_tumor_class <- hf3016_tumor_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
cell_state_order <- c("OPC-like", "NPC-like", "MES-like", "AC-like")
hf3016_tumor_class <- hf3016_tumor_class %>% mutate(class = factor(class, levels = cell_state_order))

pdf("Fig4/SuppFig8-stress-neftel-cellstates-hf3016.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf3016_tumor_class, aes(x = treatment, y = freq, fill=class)) +
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
hf3016_scgp_class = hf3016_metadata %>%
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
hf3016_scgp_class <- hf3016_scgp_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf3016_scgp_class <- hf3016_scgp_class %>% mutate(cell_type = factor(cell_type, levels = state_levels))

pdf("results/Fig4/Fig4e-stress-scgp-cellstates-hf3016.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf3016_scgp_class, aes(x = treatment, y = freq, fill=cell_type)) +
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

chisq.test(table(hf3016_metadata$cell_type[hf3016_metadata$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")], hf3016_metadata$exposure[hf3016_metadata$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")]))
chisq.test(table(hf3016_metadata$cell_type[hf3016_metadata$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")], hf3016_metadata$exposure[hf3016_metadata$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")]))
chisq.test(table(hf3016_metadata$cell_type[hf3016_metadata$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")], hf3016_metadata$exposure[hf3016_metadata$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")]))

### END ###