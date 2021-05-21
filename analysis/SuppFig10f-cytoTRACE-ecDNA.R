##################################
# Analyze transcriptional diversity differences amongst subclones.
# Updated: 2021.05.19
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
library(tidyverse)
library(Matrix)
library(openxlsx)
#library(devtools)
#devtools::install_local("/Users/johnsk/Documents/Single-Cell-DNAmethylation/CytoTRACE-master.zip")
library(CytoTRACE)
library(EnvStats)
###################################
## Theme #####
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




## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Change to HUGO gene name.
featuredata$ensembl_id <- rownames(featuredata)
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne_data_ind <- tsne.data %>% 
  filter(grepl("-4$", rownames(tsne.data)))
tsne_data_ind$sample_id = sapply(strsplit(tsne_data_ind$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
clust_annot = tsne_data_ind %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_name

## Create labels for scRNAseq.
tumor_clust <- as.factor(clust_annot$cell_type)
names(tumor_clust) <- clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^4$", "SM006", sample_id)

## Create labels for each cell.
clust_annot = clust_annot %>% 
  mutate(sample_id = recode(sample_id, `4` = "SM006"))
sample_clust <- as.factor(clust_annot$sample_id)
names(sample_clust) <- clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_filt <- log2cpm[ 1:24703, colnames(log2cpm)%in%cell_names_keep]
clust_annot_filt = clust_annot %>% 
  select(cell_name, cell_type)

## Load in the 10X CNV annotation for SM012.
cell_annotation <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM006-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)


## Run CytoTRACE.
SM006_results <- CytoTRACE(log2cpm_filt, ncores = 8)
pheno_vec = cell_annotation$ecDNA
names(pheno_vec) = cell_annotation$cell_name

## Create an output for the data.
cyto_df = as.data.frame(SM006_results$CytoTRACE)
colnames(cyto_df) <- "CytoTRACE"
cyto_df$cell_name = rownames(cyto_df)

## Combine with other annotation to plot.
cell_annotation = cell_annotation %>% 
  inner_join(cyto_df, by="cell_name")
  
##########################
### SM012
##########################
## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Change to HUGO gene name.
featuredata$ensembl_id <- rownames(featuredata)
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne_data_ind <- tsne.data %>% 
  filter(grepl("-7$", rownames(tsne.data)))
tsne_data_ind$sample_id = sapply(strsplit(tsne_data_ind$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
clust_annot = tsne_data_ind %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_name

## Create labels for scRNAseq.
tumor_clust <- as.factor(clust_annot$cell_type)
names(tumor_clust) <- clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^7$", "SM012", sample_id)

## Create labels for each cell.
clust_annot = clust_annot %>% 
  mutate(sample_id = recode(sample_id, `7` = "SM012"))
sample_clust <- as.factor(clust_annot$sample_id)
names(sample_clust) <- clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_filt <- log2cpm[ 1:24703, colnames(log2cpm)%in%cell_names_keep]
all(colnames(log2cpm_filt)==clust_annot$cell_name)

## Load in the 10X CNV annotation for SM012.
cell_annotation <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM012-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
SM012_results <- CytoTRACE(log2cpm_filt, ncores = 8)
pheno_vec = cell_annotation$ecDNA
names(pheno_vec) = cell_annotation$cell_name

## Output the CytoTRACE metric.
cyto_df = as.data.frame(SM012_results$CytoTRACE)
colnames(cyto_df) <- "CytoTRACE"
cyto_df$cell_name = rownames(cyto_df)

cell_annotation = cell_annotation %>% 
  inner_join(cyto_df, by="cell_name")

ggplot(cell_annotation, aes(x=cell_type, y=CytoTRACE, fill= ecDNA)) + 
  geom_boxplot() +
  labs(x="Cell state", y="CytoTRACE", fill = "Cell state") +
  scale_fill_manual(values=c("differentiated_tumor" = "#fcbba1", "prolif_stemcell_tumor" = "#a50f15", "stemcell_tumor" = "#fb6a4a")) 

#####################
### SM017
#####################
## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Change to HUGO gene name.
featuredata$ensembl_id <- rownames(featuredata)
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne_data_ind <- tsne.data %>% 
  filter(grepl("-9$", rownames(tsne.data)))
tsne_data_ind$sample_id = sapply(strsplit(tsne_data_ind$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
clust_annot = tsne_data_ind %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_name

## Create labels for scRNAseq.
tumor_clust <- as.factor(clust_annot$cell_type)
names(tumor_clust) <- clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^9$", "SM017", sample_id)

## Create labels for each cell.
clust_annot = clust_annot %>% 
  mutate(sample_id = recode(sample_id, `9` = "SM017"))
sample_clust <- as.factor(clust_annot$sample_id)
names(sample_clust) <- clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_filt <- log2cpm[ 1:24703, colnames(log2cpm)%in%cell_names_keep]
all(colnames(log2cpm_filt)==clust_annot$cell_name)

## Load in the 10X CNV annotation for SM017.
cell_annotation <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM017-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
SM017_results <- CytoTRACE(log2cpm_filt, ncores = 8)

pheno_vec = cell_annotation$ecDNA
names(pheno_vec) = cell_annotation$cell_name

cell_annotation = cell_annotation %>% 
  inner_join(cyto_df, by="cell_name")


#####################
### SM018
#####################
## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Change to HUGO gene name.
featuredata$ensembl_id <- rownames(featuredata)
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne_data_ind <- tsne.data %>% 
  filter(grepl("-10$", rownames(tsne.data)))
tsne_data_ind$sample_id = sapply(strsplit(tsne_data_ind$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
clust_annot = tsne_data_ind %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_name

## Create labels for scRNAseq.
tumor_clust <- as.factor(clust_annot$cell_type)
names(tumor_clust) <- clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^10$", "SM018", sample_id)

## Create labels for each cell.
clust_annot = clust_annot %>% 
  mutate(sample_id = recode(sample_id, `10` = "SM018"))
sample_clust <- as.factor(clust_annot$sample_id)
names(sample_clust) <- clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_filt <- log2cpm[ 1:24703, colnames(log2cpm)%in%cell_names_keep]
all(colnames(log2cpm_filt)==clust_annot$cell_name)

## Load in the 10X CNV annotation for SM018.
cell_annotation <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM018-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
SM018_results <- CytoTRACE(log2cpm_filt, ncores = 8)
pheno_vec = cell_annotation$ecDNA
names(pheno_vec) = cell_annotation$cell_name

cell_annotation_sm006 <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM006-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)
cell_annotation_sm012 <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM012-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)
cell_annotation_sm017 <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM017-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)
cell_annotation_sm018 <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/cnv/SM018-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Output the CytoTRACE metric.
cyto_df_sm006 = as.data.frame(SM006_results$CytoTRACE)
colnames(cyto_df_sm006) <- "CytoTRACE"
cyto_df_sm006$cell_name = rownames(cyto_df_sm006)
cyto_df_sm006 = cyto_df_sm006 %>% 
  inner_join(cell_annotation_sm006, by="cell_name")
cyto_df_sm006_clone1 = cyto_df_sm006 %>% 
  filter(clone=="subclone1")
cyto_df_sm006_clone2 = cyto_df_sm006 %>% 
  filter(clone=="subclone2")
cyto_df_sm006_clone3 = cyto_df_sm006 %>% 
  filter(clone=="subclone3")
wilcox.test(cyto_df_sm006_clone1$CytoTRACE~cyto_df_sm006_clone1$ecDNA)
wilcox.test(cyto_df_sm006_clone2$CytoTRACE~cyto_df_sm006_clone2$ecDNA)
wilcox.test(cyto_df_sm006_clone3$CytoTRACE~cyto_df_sm006_clone3$ecDNA)


cyto_df_sm012 = as.data.frame(SM012_results$CytoTRACE)
colnames(cyto_df_sm012) <- "CytoTRACE"
cyto_df_sm012$cell_name = rownames(cyto_df_sm012)
cyto_df_sm012 = cyto_df_sm012 %>% 
  inner_join(cell_annotation_sm012, by="cell_name")
cyto_df_sm012_clone1 = cyto_df_sm012 %>% 
  filter(clone=="subclone1")
cyto_df_sm012_clone3 = cyto_df_sm012 %>% 
  filter(clone=="subclone3")
cyto_df_sm012_clone4 = cyto_df_sm012 %>% 
  filter(clone=="subclone4")
wilcox.test(cyto_df_sm012_clone1$CytoTRACE~cyto_df_sm012_clone1$ecDNA)
wilcox.test(cyto_df_sm012_clone3$CytoTRACE~cyto_df_sm012_clone3$ecDNA)
wilcox.test(cyto_df_sm012_clone4$CytoTRACE~cyto_df_sm012_clone4$ecDNA)



cyto_df_sm017 = as.data.frame(SM017_results$CytoTRACE)
colnames(cyto_df_sm017) <- "CytoTRACE"
cyto_df_sm017$cell_name = rownames(cyto_df_sm017)
cyto_df_sm017 = cyto_df_sm017 %>% 
  inner_join(cell_annotation_sm017, by="cell_name")
cyto_df_sm017_clone1 = cyto_df_sm017 %>% 
  filter(clone=="subclone1")
cyto_df_sm017_clone2 = cyto_df_sm017 %>% 
  filter(clone=="subclone2")
wilcox.test(cyto_df_sm017_clone1$CytoTRACE~cyto_df_sm017_clone1$ecDNA)
wilcox.test(cyto_df_sm017_clone2$CytoTRACE~cyto_df_sm017_clone2$ecDNA)

cyto_df_sm018 = as.data.frame(SM018_results$CytoTRACE)
colnames(cyto_df_sm018) <- "CytoTRACE"
cyto_df_sm018$cell_name = rownames(cyto_df_sm018)
cyto_df_sm018 = cyto_df_sm018 %>% 
  inner_join(cell_annotation_sm018, by="cell_name")
cyto_df_sm018_clone1 = cyto_df_sm018 %>% 
  filter(clone=="subclone1")
wilcox.test(cyto_df_sm018_clone1$CytoTRACE~cyto_df_sm018_clone1$ecDNA)

cyto_df_ecDNA = bind_rows(cyto_df_sm006, cyto_df_sm012, cyto_df_sm017, cyto_df_sm018)
cyto_df_ecDNA = cyto_df_ecDNA %>% 
  mutate(cell_type = recode(cell_type, "differentiated_tumor" = "Diff-like", "prolif_stemcell_tumor" = "Prolif. stem-like",
                            "stemcell_tumor" = "Stem-like"))
cell_order = c("Prolif. stem-like", "Diff-like", "Stem-like") 
cyto_df_ecDNA <- cyto_df_ecDNA %>% mutate(cell_type = factor(cell_type, levels = cell_order))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig6/SuppFig10h-ecDNA-clone-transcriptional-diversity.pdf", width = 9, height = 4, useDingbats = FALSE)
ggplot(cyto_df_ecDNA, aes(x=clone,  y=CytoTRACE, fill=ecDNA)) + 
  geom_boxplot(outlier.shape = NA)  +
  ylab("Transcriptional diversity (CytoTRACE)") + xlab("") +
  scale_fill_manual(name='ecDNA status', values=c('ecDNA+'='#377eb8', "ecDNA-"= "gray")) +
  plot_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  stat_n_text() +
  facet_grid(~sample_id, scales = "free", space="free_x") 
dev.off()
table(cyto_df_ecDNA$clone, cyto_df_ecDNA$ecDNA, cyto_df_ecDNA$sample_id)

### END ###