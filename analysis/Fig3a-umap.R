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

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne.data$sample_id = sapply(strsplit(tsne.data$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
clust_annot = tsne.data %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 

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

DimPlot(scgp, group.by = "seurat_clusters") 
DimPlot(scgp, group.by = "cell_type") +
  scale_color_manual(values=c("b_cell" = "#eff3ff", "granulocyte" = "#bdd7e7", "t_cell" = "#6baed6", "dendritic_cell" = "#3182bd", "myeloid" = "#08519c",
                             "oligodendrocyte" = "#2ca25f",
                             "endothelial" = "#ffffd4", "pericyte" = "#fee391",
                             "fibroblast" = "#feb24c",
                             "stemcell_tumor" = "#fb6a4a", "differentiated_tumor" = "#fcbba1", "prolif_stemcell_tumor" = "#a50f15")) 
DimPlot(scgp, group.by = "sample_id") 


clust_annot = tsne.data %>% 
  mutate(cell_type = recode(dbCluster, `1` = "Diff-like (n = 16,537)",  `2` = "Myeloid (n = 15,667)", `3` = "Stem-like (n = 10,644)",
                            `4` = "Oligodendrocyte (n = 6,954)", `5` = "Prolif. stem-like (n = 3,650)", `6` = "Granulocyte (n = 787)", `7` = "Endothelial (n = 144)",
                            `8` = "T cell (n = 333)", `9` = "Pericyte (n = 103)", `10` = "Fibroblast (n = 78)", `11` = "B cell  (n = 54)", `12` = "Dendritic cell (n = 333)")) 

cell_state_order <- c("B cell  (n = 54)", "Granulocyte (n = 787)", "T cell (n = 333)", "Dendritic cell (n = 333)", "Myeloid (n = 15,667)", "Oligodendrocyte (n = 6,954)", "Endothelial (n = 144)",
                      "Pericyte (n = 103)", "Fibroblast (n = 78)", "Diff-like (n = 16,537)", "Stem-like (n = 10,644)", "Prolif. stem-like (n = 3,650)")
clust_annot <- clust_annot %>% mutate(cell_type = factor(cell_type, levels = cell_state_order))

ggplot(clust_annot, aes(x=V1, y=V2, color=cell_type)) + geom_point(size = 0.1) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                              "Oligodendrocyte" = "#2ca25f",
                              "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                              "Fibroblast" = "#feb24c",
                              "Stem-like" = "#fb6a4a", "Diff-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  plot_theme

cell_state_order <- c("B cell  (n = 54)", "Granulocyte (n = 787)", "T cell (n = 333)", "Dendritic cell (n = 333)", "Myeloid (n = 15,667)", "Oligodendrocyte (n = 6,954)", "Endothelial (n = 144)",
                      "Pericyte (n = 103)", "Fibroblast (n = 78)", "Diff-like (n = 16,537)", "Stem-like (n = 10,644)", "Prolif. stem-like (n = 3,650)")


## Filter out some cells in UMAP that impede understanding:
to_drop1 = which(clust_annot$cell_type == "Myeloid (n = 15,667)" & clust_annot$V1 < 1)
to_drop2 = which(clust_annot$cell_type == "Diff-like (n = 16,537)" & clust_annot$V1 > 1)
to_drop3 = which(clust_annot$cell_type == "Diff-like (n = 16,537)" & clust_annot$V2 > 5)
to_drop4 = which(clust_annot$cell_type == "Prolif. stem-like (n = 3,650)" & clust_annot$V1 > 1)
to_drop5 = which(clust_annot$cell_type == "Oligodendrocyte (n = 6,954)" & clust_annot$V2 < 5)
to_drop6 = which(clust_annot$cell_type == "Oligodendrocyte (n = 6,954)" & clust_annot$V2 < 5)
to_drop = c(to_drop1, to_drop2, to_drop3, to_drop4, to_drop5, to_drop6)
clust_annot_filt = clust_annot[-to_drop, ]

pal <- c("#eff3ff", "#bdd7e7", "#6baed6",  "#3182bd", "#08519c", "#2ca25f", "#ffffd4", "#fee391", "#feb24c", "#fb6a4a", "#fcbba1", "#a50f15")
pal <- setNames(pal, c("B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", "Oligodendrocyte", "Endothelial", "Pericyte", "Fibroblast",
                       "Stem-like", "Diff-like", "Prolif. stem-like"))

fig <- plot_ly(clust_annot_filt, x = ~V1, y = ~V2, z = ~V3, color = ~cell_type, colors = pal, size = 0.0001)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'UMAP 1'),
                                   yaxis = list(title = 'UMAP 2'),
                                   zaxis = list(title = 'UMAP 3'))) 
fig



png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/scgp-all-cells-umap.png", width = 7, height = 5, units = 'in', res = 300)
ggplot(clust_annot_filt, aes(x=V1, y=V2, color=cell_type)) + geom_point(size = 0.1) +
  labs(x = "UMAP 1", y = "UMAP 2", color ="Cell state\n(n = cell number)") +
  scale_color_manual(values=c("B cell  (n = 54)" = "#eff3ff","Granulocyte (n = 787)" = "#bdd7e7", "T cell (n = 333)" = "#6baed6", "Dendritic cell (n = 333)" = "#3182bd", "Myeloid (n = 15,667)" = "#08519c",
                              "Oligodendrocyte (n = 6,954)" = "#2ca25f",
                              "Endothelial (n = 144)" = "#ffffd4", "Pericyte (n = 103)" = "#fee391",
                              "Fibroblast (n = 78)" = "#feb24c",
                              "Stem-like (n = 10,644)" = "#fb6a4a", "Diff-like (n = 16,537)" = "#fcbba1", "Prolif. stem-like (n = 3,650)" = "#a50f15")) +
  plot_theme +
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

