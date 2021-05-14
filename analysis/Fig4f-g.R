##################################
# Enumerate the cell states in HF3016 GSCs under different stress exposures
# Updated: 2021.02.06
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
library(tidyverse)
library(Seurat)
###################################

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



## Load filtered, processed, unnormalized count matrix.
hf3016_obj <- ReadH5AD("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20020-22-33/RV20020-22-33-qc_20210203.h5ad")

## Make a Seurat object with the standard pipeline through PCA.
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

## Test plot:
DimPlot(hf3016_obj) 

## Break the gene expression into 3 even groups. 
avg_expression_cluster <- AverageExpression(hf3016_obj)
avg_expression_mat <- as.matrix(avg_expression_cluster$RNA)
avg_expression_num <- apply(avg_expression_mat, 1, mean)
avg_express_bins = ntile(avg_expression_num, 3)
names(avg_express_bins) <- names(avg_expression_num)
## Output the file 
saveRDS(avg_express_bins, file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/hf3016_gene_expression_bins.rds")


## Load the 10X data for all the stressed and control samples.
## HF3016 ##
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20020-22-33/RV20020-22-33_20210203.Rds")
hf3016_featuredata <- featuredata
hf3016_log2cpm <- log2cpm

# Replace ENSEMBL gene names with Hugo gene names.
rownames(hf3016_log2cpm) <- hf3016_featuredata$Associated.Gene.Name
hf3016_tsne.data <- tsne.data
hf3016_tsne.data$library <- sapply(strsplit(rownames(hf3016_tsne.data), "-"), "[[", 3)
hf3016_tsne.data$cell_name <- rownames(hf3016_tsne.data)
hf3016_tsne.data$cell_barcode <- substr(rownames(hf3016_tsne.data), 1, 18)
hf3016_tsne.data$scbl_id[hf3016_tsne.data$library==0] <- "RV20020"
hf3016_tsne.data$scbl_id[hf3016_tsne.data$library==1] <- "RV20022"
hf3016_tsne.data$scbl_id[hf3016_tsne.data$library==2] <- "RV20033"


## 2D coordinates:
hf3016_2d_umap = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20020-22-33/RV20020-22-33_umap_2d_embedding_20210203.csv", sep = ",", stringsAsFactors = F, header = T)
hf3016_2d_umap$library <- sapply(strsplit(hf3016_2d_umap$barcode, "-"), "[[", 3)
hf3016_2d_umap$scbl_id[hf3016_2d_umap$library==0] <- "RV20020"
hf3016_2d_umap$scbl_id[hf3016_2d_umap$library==1] <- "RV20022"
hf3016_2d_umap$scbl_id[hf3016_2d_umap$library==2] <- "RV20033"
hf3016_2d_umap$cell_barcode <- substr(hf3016_2d_umap$barcode, 1, 18)

## Load in hashtag map:
hashtag_map <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20017_HTO_demux/stress-hashtag-10X-map.txt", sep="\t", header=T, stringsAsFactors = F)

## The HTO tags for 3 days.
hf3016_3d_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20017_HTO_demux/RV20020/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf3016_3d_tags <- hf3016_3d_tags %>% 
  mutate(cell_name = paste(X, "-0", sep=""),
         scbl_id = "RV20020",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)

## The HTO tags for 9 days hypoxia-normoxia.
hf3016_9d_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20017_HTO_demux/RV20022/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf3016_9d_tags <- hf3016_9d_tags %>% 
  mutate(cell_name = paste(X, "-1", sep=""),
         scbl_id = "RV20022",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)

## The HTO tags for 9 days irradiation.
hf3016_9d_rt_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20034-hto_demux/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf3016_9d_rt_tags <- hf3016_9d_rt_tags %>% 
  mutate(cell_name = paste(X, "-2", sep=""),
         scbl_id = "RV20033",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  filter(hash.ID == "HHTO2") %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)


hf3016_tags <- bind_rows(hf3016_3d_tags, hf3016_9d_tags, hf3016_9d_rt_tags)
hf3016_tags_annot <- hf3016_tags %>% 
  left_join(hashtag_map, by=c("hashtag", "scbl_id"))


## Loupe-defined clusters:
hf3016_graph_based_clusters = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20020-22-33/hf3016-graph-based-clustering.csv", sep = ",", stringsAsFactors = F, header = T)
hf3016_graph_based_clusters = hf3016_graph_based_clusters %>%
  mutate(library = sapply(strsplit(Barcode, "-"), "[[", 2),
         lib_adjust = recode(library, `1` = "0",  `2` = "1", `3` = "2"),
         barcode_merge = paste(substr(Barcode, 1, 16), "1", lib_adjust, sep="-"))

hf3016_tsne_annot <- hf3016_tsne.data %>% 
  inner_join(hf3016_tags_annot, by=c("cell_barcode", "scbl_id")) %>% 
  inner_join(hf3016_2d_umap, by=c("cell_barcode", "scbl_id")) %>% 
  inner_join(hf3016_graph_based_clusters, by=c("barcode"="barcode_merge")) %>% 
  mutate(cell_type_db = recode(dbCluster, `1` = "Prolif. stem-like",  `2` = "Diff.-like", `3` = "Prolif. stem-like",
                               `4` = "Prolif. stem-like", `5` = "Stem-like", `6` = "Diff.-like", `7` = "Diff.-like",
                               `8` = "Diff.-like", `9` = "Diff.-like"),
         cell_type_graph = recode(Graph.based, `Cluster 1` = "Diff.-like",  `Cluster 2` = "Diff.-like", `Cluster 3` = "Diff.-like",
                                  `Cluster 4` = "Stem-like", `Cluster 5` = "Diff.-like", `Cluster 6` = "Diff.-like", `Cluster 7` = "Prolif. stem-like",
                                  `Cluster 8` = "Diff.-like", `Cluster 9` = "Prolif. stem-like", `Cluster 10` = "Prolif. stem-like", `Cluster 11` = "Diff.-like",
                                  `Cluster 12` = "Prolif. stem-like", `Cluster 13` = "Stem-like", `Cluster 14` = "Prolif. stem-like"),
         condition_revalue = recode(condition, `HF3016-7221` = "normoxia-3d",
                                    `HF3016-7201` = "hypoxia-3d",
                                    `HF3016-0921` = "normoxia-9d",
                                    `HF3016-0901` = "hypoxia-9d",
                                    `HF3016-09RT` = "Irradiation-9d"))

all(rownames(hf3016_obj@meta.data)==hf3016_tsne_annot$cell_name)
hf3016_obj@meta.data$cell_name <- hf3016_tsne_annot$cell_name
hf3016_obj@meta.data$condition <- hf3016_tsne_annot$condition
hf3016_obj@meta.data$oxygen_levels <- hf3016_tsne_annot$oxygen_levels
hf3016_obj@meta.data$condition_revalue <- hf3016_tsne_annot$condition_revalue
hf3016_obj@meta.data$cell_type <- hf3016_tsne_annot$cell_type_graph
hf3016_obj@meta.data$exposure_state <- paste(hf3016_tsne_annot$cell_type_graph, hf3016_tsne_annot$condition_revalue, sep="-")

## Change the UMAP coordinates to those from scanpy.
all(rownames(hf3016_obj@reductions$umap@cell.embeddings)==hf3016_tsne_annot$cell_name)
hf3016_obj@reductions$umap@cell.embeddings[,"UMAP_1"] <- hf3016_tsne_annot$x.coordinates
hf3016_obj@reductions$umap@cell.embeddings[,"UMAP_2"] <- hf3016_tsne_annot$y.coordinates 


## First plot cell state:
Idents(hf3016_obj) <- hf3016_obj@meta.data$cell_type
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = state_levels)

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/stress-cellstates-umap-hf3016-scgp-nolegend.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes() & NoLegend()
dev.off()

selected_features <- c("ASCL1", "OLIG1", "APOE", "CD44", "TOP2A", "CA9")
selected_features <- c("APOE", "LRP1", "LDLR", "TOP2A")

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/stress-genes-umap-hf3016-corrected.png", width = 9, height = 5, units = 'in', res = 300)
FeaturePlot(object = hf3016_obj, features = selected_features, ncol=3) & NoLegend() & NoAxes()
dev.off()

condition_levels <- c("normoxia-3d", "hypoxia-3d", "normoxia-9d", "hypoxia-9d", "Irradiation-9d")
Idents(hf3016_obj) <- hf3016_obj@meta.data$condition_revalue
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = condition_levels)

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/stress-conditions-umap-hf3016-nolegend.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes() & NoLegend()
dev.off()




hf3016_neftel <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/neftel-classification-hf3016-20210204.txt", header = T)
all(rownames(hf3016_obj@meta.data)==hf3016_neftel$cell_name.x)
hf3016_obj@meta.data$neftel_class <- paste(hf3016_neftel$class, "like", sep="-")
Idents(hf3016_obj) <- hf3016_obj@meta.data$neftel_class
hf3016_obj@active.ident <- factor(x = hf3016_obj@active.ident, levels = c("OPC-like", "NPC-like", "MES-like", "AC-like"))

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/stress-cellstates-umap-hf3016-neftel.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf3016_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes() & NoLegend()
dev.off()


#### Create stacked barplots ####
hf3016_scgp_class = hf3016_tsne_annot %>%
  group_by(condition, cell_type_graph) %>%
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
hf3016_scgp_class <- hf3016_scgp_class %>% mutate(cell_type_graph = factor(cell_type_graph, levels = state_levels))

pdf("github/results/Fig4/Fig4e-stress-scgp-cellstates-hf3016.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf3016_scgp_class, aes(x = treatment, y = freq, fill=cell_type_graph)) +
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

### Neftel classification ###
hf3016_tumor_class = hf3016_neftel %>%
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

pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/stress-neftel-cellstates-hf3016.pdf", width = 9, height = 5)
ggplot(hf3016_tumor_class, aes(x = treatment, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                             "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme( axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()


chisq.test(table(hf3016_tsne_annot$cell_type_graph[hf3016_tsne_annot$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")], hf3016_tsne_annot$exposure[hf3016_tsne_annot$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")]))
chisq.test(table(hf3016_tsne_annot$cell_type_graph[hf3016_tsne_annot$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")], hf3016_tsne_annot$exposure[hf3016_tsne_annot$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")]))
chisq.test(table(hf3016_tsne_annot$cell_type_graph[hf3016_tsne_annot$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")], hf3016_tsne_annot$exposure[hf3016_tsne_annot$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")]))
chisq.test(table(hf3016_tsne_annot$cell_type_graph[hf3016_tsne_annot$condition_revalue%in%c("normoxia-3d", "normoxia-9d")], hf3016_tsne_annot$timepoint[hf3016_tsne_annot$condition_revalue%in%c("normoxia-3d", "normoxia-9d")]))

chisq.test(table(hf3016_neftel$class[hf3016_neftel$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")], hf3016_neftel$exposure[hf3016_neftel$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")]))
chisq.test(table(hf3016_neftel$class[hf3016_neftel$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")], hf3016_neftel$exposure[hf3016_neftel$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")]))
chisq.test(table(hf3016_neftel$class[hf3016_neftel$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")], hf3016_neftel$exposure[hf3016_neftel$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")]))
chisq.test(table(hf3016_neftel$class[hf3016_neftel$condition_revalue%in%c("normoxia-3d", "normoxia-9d")], hf3016_neftel$timepoint[hf3016_neftel$condition_revalue%in%c("normoxia-3d", "normoxia-9d")]))

###################################
### Differential expression     ###
###################################
## Interrogate conditions first.
Idents(hf3016_obj) <- hf3016_obj@meta.data$condition_revalue

## Find differentially expressed features between 9-day irradiation and normoxia.
hf3016_irradiation_de_markers <- FindMarkers(hf3016_obj, ident.1 = "Irradiation-9d", ident.2 = "normoxia-9d", 
                                      only.pos = TRUE,
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1)
hf3016_irradiation_de_markers_filt <- hf3016_irradiation_de_markers %>% 
  filter(p_val_adj < 0.05)
## Output table with irradiation results.
write.table(hf3016_irradiation_de_markers_filt, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/irradiation-deg-hf3016-20210217.txt", sep="\t", row.names = T, col.names = T, quote = F)


## Read in HF2354 differentially expressed genes for irradiation:
hf2354_irradiation_degs <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/irradiation-deg-hf3016-20210331.txt", sep="\t", header = TRUE, row.names = 1)



## What about hypoxia markers?
hf3016_hypoxia_de_markers <- FindMarkers(hf3016_obj, ident.1 = "hypoxia-9d", ident.2 = "normoxia-9d",
                                  only.pos = TRUE,
                                  logfc.threshold = 0.25,
                                  min.pct = 0.1)
hf3016_hypoxia_de_markers_filt <- hf3016_hypoxia_de_markers %>% 
  filter(p_val_adj < 0.05)
write.table(hf3016_hypoxia_de_markers_filt, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/hypoxia-deg-hf3016-20210217.txt", sep="\t", row.names = T, col.names = T, quote = F)


## Interrogate cell-state specific differences in conditions.
Idents(hf3016_obj) <- hf3016_obj@meta.data$exposure_state
irradiation_stem_de_markers <- FindMarkers(hf3016_obj, ident.1 = "Stem-like-Irradiation-9d", ident.2 = "Stem-like-normoxia-9d")
irradiation_diff_de_markers <- FindMarkers(hf3016_obj, ident.1 = "Diff.-like-Irradiation-9d", ident.2 = "Diff.-like-normoxia-9d")


###################################
### FindMarkers + Heatmap       ###
###################################
hf3016_obj@meta.data$condition_revalue  <- factor(x = hf3016_obj@meta.data$condition_revalue , levels = c("normoxia-3d", "hypoxia-3d", "normoxia-9d", "hypoxia-9d", "Irradiation-9d"))
Idents(hf3016_obj) <- hf3016_obj@meta.data$condition_revalue

stress_markers <- FindAllMarkers(hf3016_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- stress_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/heatmap_hf3016_timepoint_features.png", width = 9, height = 5, units = 'in', res = 300)
DoHeatmap(hf3016_obj, features = top10$gene,
          angle = 45, size = 3,
          group.colors = c("#bdc9e1", "#fd8d3c", "#045a8d", "#bd0026", "#0dba86"))+
  scale_fill_gradientn(colours = rev(mapal))
dev.off()



## Cell state:
Idents(hf3016_obj) <- hf3016_obj@meta.data$cell_type

state_markers <- FindAllMarkers(hf3016_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- state_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

png("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/heatmap_hf3016_cellstate_features.png", width = 9, height = 5, units = 'in', res = 300)
DoHeatmap(hf3016_obj, features = top10$gene,
          angle = 45, size = 3,
          group.colors = c("#fcbba1", "#a50f15", "#fb6a4a"))+
  scale_fill_gradientn(colours = rev(mapal))
dev.off()


## Write out the total denominator for this analysis.
write.table(featuredata, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/all-genes-measured-hf3016-20210331.txt", sep="\t", row.names = T, col.names = T, quote = F)


#### END ####
