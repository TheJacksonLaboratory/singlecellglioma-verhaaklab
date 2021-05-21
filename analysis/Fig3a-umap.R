##################################
# Perform the dimensionality reduction on all 10X cells
# Updated: 2020.05.13
# Author: Kevin J.
###################################

mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
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


## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## 2D UMAP coordinates.
umap_coords_2d = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_umap_2d_embedding.csv")

### Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
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
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^10$", "SM018", sample_id)
sample_id <- gsub("^0$", "SM019", sample_id)
sample_id <- gsub("^1$", "SM001", sample_id)
sample_id <- gsub("^2$", "SM002", sample_id)
sample_id <- gsub("^3$", "SM004", sample_id)
sample_id <- gsub("^4$", "SM006", sample_id)
sample_id <- gsub("^5$", "SM011", sample_id)
sample_id <- gsub("^8$", "SM015", sample_id)
sample_id <- gsub("^6$", "SM008", sample_id)
sample_id <- gsub("^7$", "SM012", sample_id)
sample_id <- gsub("^9$", "SM017", sample_id)


## Only raw counts.
raw_cpm = exp(log2cpm[c(1:24703), ])

## Create a Seurat object using raw counts.
scgp <- CreateSeuratObject(counts = raw_cpm, min.cells = 1, project = "scgp", names.field = 1, names.delim = "_")
scgp@meta.data$sample_id <- sample_id
scgp@meta.data$cell_type <- clust_annot$cell_type

## Make a Seurat object with the standard pipeline through PCA.
scgp <- scgp %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = scgp@var.genes, npcs = 20, verbose = FALSE)

## Visualize.
scgp <- scgp %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()


## Change the UMAP coordinates to those from scanpy.
all(rownames(scgp@reductions$umap@cell.embeddings)==tsne.data$cell_name)
scgp@reductions$umap@cell.embeddings[,"UMAP_1"] <- tsne.data$x.coordinate
scgp@reductions$umap@cell.embeddings[,"UMAP_2"] <- tsne.data$y.coordinate

## Change to cell state.
Idents(scgp) <- scgp@meta.data$cell_type

my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")

scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig3a-all-cells-umap-states.png", width = 8, height = 5, units = 'in', res = 300)
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

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig3a-all-cells-umap-states.pdf", width = 8, height = 5, useDingbats = FALSE)
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