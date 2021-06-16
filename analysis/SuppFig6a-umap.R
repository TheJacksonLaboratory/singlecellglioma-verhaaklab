##################################
# Visualize patient labels for UMAP in Supplementary Figure 6a.
# Updated: 2021.05.10
# Author: Kevin J.
###################################

mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
library(openxlsx)
###################################
## Load the 10X data for all tumor samples.
load("/Users/johnsk/github/data/analysis_scRNAseq_tumor_gene_expression.Rds")

## 2D UMAP coordinates.
umap_coords_2d <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Use Seurat to create UMAP visualization.
expr_unnorm_data = exp(expr_norm_data[c(1:24703), ])-1
scgp <- CreateSeuratObject(counts = expr_unnorm_data, min.cells = 1, project = "scgp", names.field = 1, names.delim = "_")
scgp@meta.data$cell_barcode <- umap_coords_2d$cell_barcode
scgp@meta.data$cell_state <- umap_coords_2d$cell_state
scgp@meta.data$case_barcode <- umap_coords_2d$case_barcode

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

png("github/results/all-cells-umap-samples.png", width = 7, height = 5, units = 'in', res = 300)
DimPlot(scgp, group.by = "case_barcode") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        line = element_blank())
dev.off()

### END ####