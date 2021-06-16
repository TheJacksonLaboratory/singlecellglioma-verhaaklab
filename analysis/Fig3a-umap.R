##################################
# Visualize the dimensionality reduction on all 10X cells
# Updated: 2020.05.13
# Author: Kevin J.
###################################

mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(plotly)
library(openxlsx)
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


## Load the 10X data for all tumor samples.
load("/Users/johnsk/github/data/analysis_scRNAseq_tumor_gene_expression.Rds")

## 2D UMAP coordinates.
umap_coords_2d <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Use Seurat to create UMAP visualization.
raw_cpm = exp(expr_norm_data[c(1:24703), ])-1
scgp <- CreateSeuratObject(counts = raw_cpm, min.cells = 1, project = "scgp", names.field = 1, names.delim = "_")
scgp@meta.data$cell_barcode <- umap_coords_2d$cell_barcode
scgp@meta.data$cell_state <- umap_coords_2d$cell_state

## Make a Seurat object with the standard pipeline through PCA.
scgp <- scgp %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = scgp@var.genes, npcs = 20, verbose = FALSE)

## Visualize out of curiosity. Note that this is going to yield a different output from Scanpy, which uses different parameters.
scgp <- scgp %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

## Change the UMAP coordinates to those from scanpy.
all(rownames(scgp@reductions$umap@cell.embeddings)==umap_coords_2d$cell_barcode)
scgp@reductions$umap@cell.embeddings[,"UMAP_1"] <- umap_coords_2d$umap_1
scgp@reductions$umap@cell.embeddings[,"UMAP_2"] <- umap_coords_2d$umap_2

## Change to cell states already defined.
Idents(scgp) <- scgp@meta.data$cell_state

## Create desired levels of cell states defined with marker gene expression.
my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")
scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)

png("results/Fig3a-all-cells-umap-states.png", width = 8, height = 5, units = 'in', res = 300)
DimPlot(scgp) +
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        line = element_blank())
dev.off()

## or PDF for paper.
pdf("results/Fig3a-all-cells-umap-states.pdf", width = 8, height = 5, useDingbats = FALSE)
DimPlot(scgp) +
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        line = element_blank())
dev.off()

### END ###