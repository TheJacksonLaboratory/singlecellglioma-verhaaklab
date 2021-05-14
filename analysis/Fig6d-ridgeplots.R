##################################
# Visualize specific features for tumor SM012.
# Updated: 2021.05.07
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

## Create a Seurat object using raw counts to leverage plotting functionality.
scgp <- CreateSeuratObject(counts = raw_cpm, min.cells = 1, project = "scgp_all", names.field = 1, names.delim = "_")
scgp@meta.data$sample_id <- clust_annot$sample_id
scgp@meta.data$cell_type <- clust_annot$cell_type

## Make a Seurat object with the standard pipeline through PCA.
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

DimPlot(scgp, group.by = "seurat_clusters") 
DimPlot(scgp, group.by = "cell_type") +
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) 
DimPlot(scgp, group.by = "sample_id") 
Idents(scgp) <- scgp@meta.data$cell_type

## Change the UMAP coordinates to those from scanpy.
all(rownames(scgp@reductions$umap@cell.embeddings)==tsne.data$cell_name)
scgp@reductions$umap@cell.embeddings[,"UMAP_1"] <- tsne.data$x.coordinate
scgp@reductions$umap@cell.embeddings[,"UMAP_2"] <- tsne.data$y.coordinate

my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")

scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)

#########################################
#### Identify subclonal events for SM012
#########################################
## Determined by inferCNV.
cell_annotation <- read.table("results/cnv/SM012-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Subset to only SM012 tumor cells.
sm012_scgp = subset(scgp, idents = c("Stem-like", "Prolif. stem-like", "Diff.-like"), subset = sample_id == "SM012")
DimPlot(object = sm012_scgp, reduction = "umap", label = T)
FeaturePlot(object = sm012_scgp, features = c("HES6", "EZH2", "EPAS1", "CD44")) + labs(x="UMAP 1",y = "UMAP 2")

## Need to add annotation of tumor subclone.
sm012_scgp@meta.data$cell_name <- rownames(sm012_scgp@meta.data)
## Add extra information (Neftel classification). 
all(rownames(sm012_scgp@meta.data)==cell_annotation$cell_name)
sm012_scgp@meta.data$class <- as.factor(cell_annotation$class)
DimPlot(sm012_scgp, group.by = "class", reduction = "umap", label = T) 

sm012_scgp@meta.data$clone <- as.factor(cell_annotation$clone)
DimPlot(sm012_scgp, group.by = "clone", reduction = "umap", label = T) 
## Examine the distribution of clones by Neftel class.
prop.table(table(cell_annotation$clone, cell_annotation$class), margin = 1)

# Find differentially expressed features between each clone and all other cells.
Idents(sm012_scgp) <- sm012_scgp@meta.data$clone
clone1_markers = FindMarkers(sm012_scgp, ident.1 = "subclone1", ident.2 = NULL, only.pos = TRUE)
clone2_markers = FindMarkers(sm012_scgp, ident.1 = "subclone2", ident.2 = NULL, only.pos = TRUE)
clone3_markers = FindMarkers(sm012_scgp, ident.1 = "subclone3", ident.2 = NULL, only.pos = TRUE)
clone4_markers = FindMarkers(sm012_scgp, ident.1 = "subclone4", ident.2 = NULL, only.pos = TRUE)
ecDNA_markers = FindMarkers(sm012_scgp, ident.1 = c("subclone3", "subclone4"), ident.2 = c("subclone1", "subclone2"), only.pos = TRUE)
ecDNA_neg_markers = FindMarkers(sm012_scgp, ident.1 = c("subclone1", "subclone2"), ident.2 = c("subclone3", "subclone4"),  only.pos = TRUE)

## Plot ecDNA gene (EGFR) as well as hypoxia-associated genes.
features = c("EGFR", "PDGFRA", "VEGFA", "PGK1")
pdf("github/results/Fig6d_SM012_ridgeplot.pdf", width = 7, height = 5)
RidgePlot(obj = sm012_scgp, features = features, ncol = 2, cols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  labs(x="", y="")
dev.off()

### END ###
