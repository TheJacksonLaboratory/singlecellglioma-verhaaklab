##################################
# Visualize specific genes for Supplementary Figure 6b - cluster gene expression violin plots.
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


## Load the 10X data for all tumor samples.
load("/Users/johnsk/github/data/analysis_scRNAseq_tumor_gene_expression.Rds")

## Change to HUGO gene name.
rownames(expr_norm_data)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

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

## Create the different levels of cell states to be presented.
my_levels = c("Stem-like", "Prolif. stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Pericyte", "Endothelial",
              "B cell", "T cell", "Granulocyte", "Dendritic cell", "Myeloid")

Idents(scgp) <- scgp@meta.data$cell_state
scgp@active.ident <- factor(x = scgp@active.ident, levels = my_levels)

#### Provide a set of genes that helps discriminate cell states:
features_set1 <- c("SOX2", "OLIG2", "ASCL1", "TOP2A", "EGFR", "AQP4",
              "ID4", "MOG", "HLA-DPB1")
features_set2 <- c("FBLN1", "CLDN5", "DCN", "CD79A", "S100A9", 
              "CD3D", "HLA-DPB1", "CD14", "C1QA")

pdf("results/Fig3/cell-states-violin_set1.pdf", width = 8, height = 11)
StackedVlnPlot(obj = scgp, features = features_set1, cols=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                                                       "Oligodendrocyte" = "#2ca25f",
                                                       "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                                                       "Fibroblast" = "#feb24c",
                                                       "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))
dev.off()

pdf("results/Fig3/cell-states-violin_set2.pdf", width = 8, height = 11)
StackedVlnPlot(obj = scgp, features = features_set2, cols=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                                                            "Oligodendrocyte" = "#2ca25f",
                                                            "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                                                            "Fibroblast" = "#feb24c",
                                                            "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15"))
dev.off()

### END ####