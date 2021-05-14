##############################################
# Assess the differences in gene expression and cell states based on hypoxia exposure for HF2354
# Updated: 2021.05.13
# Author: Kevin J.
###############################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###########################
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(Seurat)
###########################

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
hf2354_obj <- ReadH5AD("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20016-18-33/RV20016-18-33-qc_20210203.h5ad")

## Make a Seurat object with the standard pipeline through PCA.
## Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
## This is then natural-log transformed using log1p.
hf2354_obj <- hf2354_obj %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = hf2354_obj@var.genes, npcs = 20, verbose = FALSE)

## Visualize.
hf2354_obj <- hf2354_obj %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

## Test plot:
DimPlot(hf2354_obj) 

## Break the gene expression into 3 even groups. 
avg_expression_cluster <- AverageExpression(hf2354_obj)
avg_expression_mat <- as.matrix(avg_expression_cluster$RNA)
avg_expression_num <- apply(avg_expression_mat, 1, mean)
avg_express_bins = ntile(avg_expression_num, 3)
names(avg_express_bins) <- names(avg_expression_num)
## Output the file 
saveRDS(avg_express_bins, file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/hf2354_gene_expression_bins.rds")

## Load the 10X data for all the stressed and control samples.
## HF2354 ##
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20016-18-33/RV20016-18-33_20210203.Rds")
hf2354_featuredata <- featuredata
hf2354_log2cpm <- log2cpm

# Replace ENSEMBL gene names with Hugo gene names.
rownames(hf2354_log2cpm) <- hf2354_featuredata$Associated.Gene.Name
hf2354_tsne.data <- tsne.data
hf2354_tsne.data$library <- sapply(strsplit(rownames(hf2354_tsne.data), "-"), "[[", 3)
hf2354_tsne.data$cell_name <- rownames(hf2354_tsne.data)
hf2354_tsne.data$cell_barcode <- substr(rownames(hf2354_tsne.data), 1, 18)
hf2354_tsne.data$scbl_id[hf2354_tsne.data$library==0] <- "RV20016"
hf2354_tsne.data$scbl_id[hf2354_tsne.data$library==1] <- "RV20018"
hf2354_tsne.data$scbl_id[hf2354_tsne.data$library==2] <- "RV20033"



## 2D coordinates:
hf2354_2d_umap = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20016-18-33/RV20016-18-33_umap_2d_embedding_20210203.csv", sep = ",", stringsAsFactors = F, header = T)
hf2354_2d_umap$library <- sapply(strsplit(hf2354_2d_umap$barcode, "-"), "[[", 3)
hf2354_2d_umap$scbl_id[hf2354_2d_umap$library==0] <- "RV20016"
hf2354_2d_umap$scbl_id[hf2354_2d_umap$library==1] <- "RV20018"
hf2354_2d_umap$scbl_id[hf2354_2d_umap$library==2] <- "RV20033"
hf2354_2d_umap$cell_barcode <- substr(hf2354_2d_umap$barcode, 1, 18)

## Load in hashtag map:
hashtag_map <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20017_HTO_demux/stress-hashtag-10X-map.txt", sep="\t", header=T, stringsAsFactors = F)

## The HTO tags for 3 days.
hf2354_3d_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20017_HTO_demux/RV20016/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf2354_3d_tags <- hf2354_3d_tags %>% 
  mutate(cell_name = paste(X, "-0", sep=""),
         scbl_id = "RV20016",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)

## The HTO tags for 9 days hypoxia-normoxia.
hf2354_9d_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20019-hto_demux/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf2354_9d_tags <- hf2354_9d_tags %>% 
  mutate(cell_name = paste(X, "-1", sep=""),
         scbl_id = "RV20018",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)

## The HTO tags for 9 days irradiation.
hf2354_9d_rt_tags = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/RV20034-hto_demux/hto_exome_2tags_matrix_hashid.csv", sep = ",", stringsAsFactors = F, header = T)
hf2354_9d_rt_tags <- hf2354_9d_rt_tags %>% 
  mutate(cell_name = paste(X, "-2", sep=""),
         scbl_id = "RV20033",
         cell_barcode = substr(cell_name, 1, 18)) %>% 
  filter(hash.ID == "HHTO1") %>% 
  select(cell_name, cell_barcode, hashtag = hash.ID, scbl_id)


hf2354_tags <- bind_rows(hf2354_3d_tags, hf2354_9d_tags, hf2354_9d_rt_tags)
hf2354_tags_annot <- hf2354_tags %>% 
  left_join(hashtag_map, by=c("hashtag", "scbl_id"))

hf2354_tsne_annot <- hf2354_tsne.data %>% 
  inner_join(hf2354_tags_annot, by=c("cell_barcode", "scbl_id")) %>% 
  inner_join(hf2354_2d_umap, by=c("cell_barcode", "scbl_id")) %>% 
  mutate(cell_type = recode(dbCluster, `1` = "Diff.-like",  `2` = "Diff.-like", `3` = "Prolif. stem-like",
                            `4` = "Stem-like", `5` = "Stem-like", `6` = "Prolif. stem-like", `7` = "Prolif. stem-like",
                            `8` = "Prolif. stem-like", `9` = "Stem-like"),
         condition_revalue = recode(condition, `HF2354-7221` = "normoxia-3d",
                                    `HF2354-7201` = "hypoxia-3d",
                                    `HF2354-0921` = "normoxia-9d",
                                    `HF2354-0901` = "hypoxia-9d",
                                    `HF2354-09RT` = "Irradiation-9d"))


all(rownames(hf2354_obj@meta.data)==hf2354_tsne_annot$cell_name)
hf2354_obj@meta.data$cell_name <- hf2354_tsne_annot$cell_name
hf2354_obj@meta.data$condition <- hf2354_tsne_annot$condition
hf2354_obj@meta.data$oxygen_levels <- hf2354_tsne_annot$oxygen_levels
hf2354_obj@meta.data$condition_revalue <- hf2354_tsne_annot$condition_revalue
hf2354_obj@meta.data$cell_type <- hf2354_tsne_annot$cell_type
hf2354_obj@meta.data$exposure_state <- paste(hf2354_tsne_annot$cell_type, hf2354_tsne_annot$condition_revalue, sep="-")

## Change the UMAP coordinates to those from scanpy.
all(rownames(hf2354_obj@reductions$umap@cell.embeddings)==hf2354_tsne_annot$cell_name)
hf2354_obj@reductions$umap@cell.embeddings[,"UMAP_1"] <- hf2354_tsne_annot$x.coordinates
hf2354_obj@reductions$umap@cell.embeddings[,"UMAP_2"] <- hf2354_tsne_annot$y.coordinates 

## First plot cell state:
Idents(hf2354_obj) <- hf2354_obj@meta.data$cell_type
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = state_levels)

pdf("github/results/Fig4/Fig4d-hf2354-scgp-umap.pdf", width = 4.5, height = 2.5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes()
dev.off()

png("github/results/Fig4/Fig4d-hf2354-scgp-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("Stem-like" = "#fb6a4a", 
                              "Diff.-like" = "#fcbba1", 
                              "Prolif. stem-like" = "#a50f15")) & NoAxes() & NoLegend()
dev.off()


condition_levels <- c("normoxia-3d", "hypoxia-3d", "normoxia-9d", "hypoxia-9d", "Irradiation-9d")
Idents(hf2354_obj) <- hf2354_obj@meta.data$condition_revalue
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = condition_levels)

png("github/results/Fig4/Fig4d-hf2354-conditions-umap.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes() & NoLegend()
dev.off()


pdf("github/results/Fig4/Fig4d-hf2354-conditions-umap-legend.pdf", width = 9, height = 5, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj)  +
  scale_color_manual(values=c("hypoxia-3d" = "#fd8d3c",
                              "hypoxia-9d" = "#bd0026",
                              "normoxia-3d" = "#bdc9e1",
                              "normoxia-9d" = "#045a8d",
                              "Irradiation-9d" = "#0dba86")) & NoAxes()
dev.off()


hf2354_neftel <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/neftel-classification-hf2354-20210205.txt", header = T)
all(rownames(hf2354_obj@meta.data)==hf2354_neftel$cell_name)
hf2354_obj@meta.data$neftel_class <- paste(hf2354_neftel$class, "like", sep="-")
Idents(hf2354_obj) <- hf2354_obj@meta.data$neftel_class
hf2354_obj@active.ident <- factor(x = hf2354_obj@active.ident, levels = c("OPC-like", "NPC-like", "MES-like", "AC-like"))

pdf("github/results/Fig4/Fig4d-hf2354-neftel-umap-legend.pdf", width = 4, height = 3, bg = "transparent", useDingbats = FALSE)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes()
dev.off()

png("github/results/Fig4/Fig4d-hf2354-neftel-umap-no-legend.png", width = 4, height = 3, units = 'in', res = 300)
DimPlot(hf2354_obj) +
  scale_color_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                              "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) & NoAxes() & NoLegend()
dev.off()

hf2354_tumor_class = hf2354_tsne_annot %>%
  inner_join(hf2354_neftel, by=c("cell_name.x" = "cell_name")) %>% 
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
hf2354_tumor_class <- hf2354_tumor_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
cell_state_order <- c("OPC-like", "NPC-like", "MES-like", "AC-like")
hf2354_tumor_class <- hf2354_tumor_class %>% mutate(class = factor(class, levels = cell_state_order))

pdf("github/results/Fig4/SuppFig8-stress-neftel-cellstates-hf2354.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf2354_tumor_class, aes(x = treatment, y = freq, fill=class)) +
  geom_bar(stat="identity") +
  labs(x="", y = "Proportion tumor\ncell state", fill="Neftel\ncell state") +
  scale_fill_manual(values=c("MES-like" = "#d7191c", "AC-like" = "#fdae61",
                             "NPC-like" = "#abd9e9", "OPC-like" = "#2c7bb6")) +
  facet_grid(.~timepoint, scales="free_x", space = "free_x") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(freq, 2)), size = 2, hjust = 0.5, vjust = 3, position ="stack") 
dev.off()


#### Create stacked barplots ####
hf2354_scgp_class = hf2354_tsne_annot %>%
  group_by(condition, cell_type) %>%
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
hf2354_scgp_class <- hf2354_scgp_class %>% mutate(treatment = factor(treatment, levels = treatment_order))
state_levels <- c("Diff.-like", "Stem-like", "Prolif. stem-like")
hf2354_scgp_class <- hf2354_scgp_class %>% mutate(cell_type = factor(cell_type, levels = state_levels))

pdf("github/results/Fig4/Fig4e-stress-scgp-cellstates-hf2354.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(hf2354_scgp_class, aes(x = treatment, y = freq, fill=cell_type)) +
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

chisq.test(table(hf2354_tsne_annot$cell_type[hf2354_tsne_annot$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")], hf2354_tsne_annot$exposure[hf2354_tsne_annot$condition_revalue%in%c("normoxia-3d", "hypoxia-3d")]))
chisq.test(table(hf2354_tsne_annot$cell_type[hf2354_tsne_annot$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")], hf2354_tsne_annot$exposure[hf2354_tsne_annot$condition_revalue%in%c("normoxia-9d", "hypoxia-9d")]))
chisq.test(table(hf2354_tsne_annot$cell_type[hf2354_tsne_annot$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")], hf2354_tsne_annot$exposure[hf2354_tsne_annot$condition_revalue%in%c("normoxia-9d", "Irradiation-9d")]))

### END ###