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


## Credit to Tommy Tang: https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-2, -2, -2,-2), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-2, -2, -2, -2), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


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

## Create the different levels of cell states to be presented.
my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Pericyte", "Endothelial",
              "B cell", "T cell", "Granulocyte", "Dendritic cell", "Myeloid")

Idents(scgp) <- scgp@meta.data$cell_type
scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)

#### Provide a set of genes that discriminates cell states:
features_set1 <- c("SOX2", "OLIG2", "ASCL1", "TOP2A", "EGFR", "AQP4",
              "ID4", "MOG", "FBLN1")
features_set2 <- c("CLDN5", "DCN", "CD79A", "S100A9", 
              "CD3D", "HLA-DPB1", "CD14", "C1QA")



pdf("github/results/cell-states-violin_set1.pdf", width = 8, height = 11)
StackedVlnPlot(obj = scgp, features = features_set1, cols=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                                                       "Oligodendrocyte" = "#2ca25f",
                                                       "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                                                       "Fibroblast" = "#feb24c",
                                                       "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme(axis.text = element_text(size = 12, angle=45))
dev.off()

pdf("github/results/cell-states-violin_set2.pdf", width = 8, height = 11)
StackedVlnPlot(obj = scgp, features = features_set2, cols=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                                                            "Oligodendrocyte" = "#2ca25f",
                                                            "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                                                            "Fibroblast" = "#feb24c",
                                                            "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme(axis.text = element_text(size = 12, angle=45))
dev.off()

### END ####
