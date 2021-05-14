##################################
# Visualize patient labels for UMAP in Supplementary Figure 6a.
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

## Only raw counts.
raw_cpm = exp(log2cpm[c(1:24703), ])-1

## Create a Seurat object using raw counts.
scgp <- CreateSeuratObject(counts = raw_cpm, min.cells = 1, project = "scgp_all", names.field = 1, names.delim = "_")
scgp@meta.data$sample_id <- clust_annot$sample_id 
scgp@meta.data$cell_type <- clust_annot$cell_type

## Make a Seurat object with the standard pipeline through PCA to normalize the data.
scgp <- scgp %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = scgp@var.genes, npcs = 20, verbose = FALSE)

## Visualize.
scgp <- scgp %>% 
  RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 35, min.dist = 0.5) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20, ) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

## Change the UMAP coordinates to those from scanpy.
all(rownames(scgp@reductions$umap@cell.embeddings)==tsne.data$cell_name)
scgp@reductions$umap@cell.embeddings[,"UMAP_1"] <- tsne.data$x.coordinate
scgp@reductions$umap@cell.embeddings[,"UMAP_2"] <- tsne.data$y.coordinate

png("github/results/all-cells-umap-samples.png", width = 7, height = 5, units = 'in', res = 300)
DimPlot(scgp, group.by = "sample_id") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        line = element_blank())
dev.off()


### END ####