##################################
# Visualize specific features for tumor SM012.
# Updated: 2021.05.07
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
###################################

## Load the 10X data for all tumor samples.
load("data/analysis_scRNAseq_tumor_gene_expression.Rds")

## Change to HUGO gene name.
rownames(expr_norm_data)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## 2D UMAP coordinates.
clust_annot <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

## Untransformed the data.
expr_data = exp(expr_norm_data[c(1:24703), ])-1

## Create a Seurat object using raw counts to leverage plotting functionality.
scgp <- CreateSeuratObject(counts = expr_data, min.cells = 1, project = "scgp_all", names.field = 1, names.delim = "_")
scgp@meta.data$sample_id <- clust_annot$case_barcode
scgp@meta.data$cell_type <- clust_annot$cell_state

## Make a Seurat object with the standard pipeline through PCA.
scgp <- scgp %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = scgp@var.genes, npcs = 20, verbose = FALSE)

## Visualize with Seurat out of curiosity - clearly distinct from the Scanpy results.
scgp <- scgp %>% 
  RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 35, min.dist = 0.5) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20, ) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

DimPlot(scgp, group.by = "cell_type") +
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) 

Idents(scgp) <- scgp@meta.data$cell_type

## Change the UMAP coordinates to those from scanpy.
all(rownames(scgp@reductions$umap@cell.embeddings)==clust_annot$cell_barcode)
scgp@reductions$umap@cell.embeddings[,"UMAP_1"] <- clust_annot$umap_1
scgp@reductions$umap@cell.embeddings[,"UMAP_2"] <- clust_annot$umap_2

my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")

scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)



#########################################
#### Identify subclonal events for SM012
#########################################
## CNV subclonal groupings identified using inferCNV.
cell_annotation <- read.table("data/SM012-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Subset to only SM012 tumor cells.
sm012_scgp = subset(scgp, idents = c("Stem-like", "Prolif. stem-like", "Diff.-like"), subset = sample_id == "SM012")
DimPlot(object = sm012_scgp, reduction = "umap", label = T)
FeaturePlot(object = sm012_scgp, features = c("OLIG2", "EGFR", "PDGFRA", "EPAS1")) + labs(x="UMAP 1",y = "UMAP 2")

## Add annotation of tumor subclone.
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

## Plot ecDNA gene (EGFR/PGFRA) as well as hypoxia-associated genes (VEGFA + PGK1) upregulated in subclone 2.
features = c("EGFR", "PDGFRA", "VEGFA", "PGK1")
pdf("results/Fig6d_SM012_ridgeplot.pdf", width = 7, height = 5)
RidgePlot(obj = sm012_scgp, features = features, ncol = 2, cols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  labs(x="", y="")
dev.off()

### END ###
