##############################################
# Investigate DNAme disorder and gene expression levels
# Updated: 2021.05.18
# Author: Kevin J.
##################################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
library(tidyverse)
library(Matrix)
library(openxlsx)
library(Seurat)
library(EnvStats)
library(ggpubr)
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

######################
### IDHmut
######################
## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier for IDHmut samples.
tsne_data_mut <- tsne.data %>% 
  filter(grepl("-0$|-1$|-2$|-3$|-6$|-8$", rownames(tsne.data)))
tsne_data_mut$sample_id = sapply(strsplit(tsne_data_mut$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
mut_clust_annot = tsne_data_mut %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(mut_clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
mut_clust_annot <- mut_clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- mut_clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
sample_id <- sapply(strsplit(mut_clust_annot$cell_name, "-"), "[[", 3)
sample_id <- gsub("^10$", "SM018", sample_id)
sample_id <- gsub("^0$", "UC917", sample_id)
sample_id <- gsub("^1$", "SM001", sample_id)
sample_id <- gsub("^2$", "SM002", sample_id)
sample_id <- gsub("^3$", "SM004", sample_id)
sample_id <- gsub("^4$", "SM006", sample_id)
sample_id <- gsub("^5$", "SM011", sample_id)
sample_id <- gsub("^8$", "SM015", sample_id)
sample_id <- gsub("^6$", "SM008", sample_id)
sample_id <- gsub("^7$", "SM012", sample_id)
sample_id <- gsub("^9$", "SM017", sample_id)

## Restrict the RNAseq data to the tumor cells.
log2cpm_mut <- log2cpm[ , colnames(log2cpm)%in%cell_names_keep] 
all(colnames(log2cpm_mut)==mut_clust_annot$cell_name)

## Only raw counts.
raw_cpm_mut = exp(log2cpm_mut[c(1:24703), ])-1

scgp_mut <- CreateSeuratObject(counts = raw_cpm_mut, min.cells = 1, project = "IDHmut", names.field = 1, names.delim = "_")
scgp_mut@meta.data$sample_id <- sample_id
scgp_mut@meta.data$cell_type <- mut_clust_annot$cell_type
avail_genes = rownames(scgp_mut@assays$RNA@data)

## Make a Seurat object with normalization and calculating gene expression variability.
scgp_mut <- scgp_mut %>% 
  Seurat::NormalizeData(verbose = TRUE) %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

## Create a data.frame with the expression variability metrics.
rna_var_mut = scgp_mut@assays$RNA@meta.features
rna_var_mut$gene_id <- rownames(rna_var_mut)

## Read in the previously calculated DNAme disorder metrics.
epimut_gb_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_genebody_total_filt-20210207.txt", header = T, sep="\t")

## Combine with epimutation data for all .
rna_var_mut_epi = rna_var_mut %>% 
  inner_join(epimut_gb_total_filt, by="gene_id")

# Or just the IDHmut cells.
# rna_var_mut_epi = rna_var_mut %>% 
#  inner_join(epimut_IDHm_total, by="gene_id")

## Break into groups based on levels of epimutation (low, intermediate, and high). 
rna_var_mut_epi$epimut_groups <-cut(rna_var_mut_epi$tumor_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1.1))
## How do these groups look in terms of total number of genes.
table(rna_var_mut_epi$epimut_groups) # Fewer genes with high epimutation.

## Plots for both gene expression and variability that is controlled for expression levels.
## The dispersion method helps control for the relationship between variability and average expression.
#my_comparisons <- list( c("(-0.01,0.1]", "(0.1,0.4]"), c("(-0.01,0.1]", "(0.4,1.1]"), c("(0.1,0.4]", "(0.4,1.1]") )

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/IDHmut-genebody-epimutation-dispersion-loess.pdf", width = 7, height = 5)
ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.dispersion.scaled)) + geom_boxplot(outlier.shape = NA) +
  ylim(-2.25, 2.25) + 
  labs(x="Promoter epimutation groups", y="Expression dispersion (mean-scaled)", fill="Epimutation groups") +
#  scale_fill_manual(values=c("(-0.01,0.1]" = "#377eb8", 
#                              "(0.1,0.2]" = "gray90",
#                              "(0.2,0.3]" = "gray80",
#                              "(0.3,0.4]" = "gray70",
#                              "(0.4,0.5]" = "gray60",
#                              "(0.5,1.1]" = "#CD4F39")) +
  plot_theme + 
  stat_compare_means(method="kruskal.test")  + 
  geom_smooth(aes(x=as.integer(epimut_groups), y = mvp.dispersion.scaled), method="loess") +
  ggtitle("IDHmut scRNA tumor cells")
  #stat_compare_means(comparisons = my_comparisons) +
  #stat_n_text(y.pos = -2.25) +
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/IDHmut-genebody-epimutation-expression.pdf", width = 7, height = 5)
ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.mean)) + geom_boxplot() +
  stat_compare_means(method="kruskal.test") + 
  #ylim(-2.25, 2.25) + 
  #ylim(0, 1.5) + 
  labs(x="Promoter epimutation groups", y="Mean expression") +
  plot_theme + 
  stat_n_text(y.pos = 1)
dev.off()


######################
### IDHwt
######################
## scRNAseq ##
## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
## Extract the numeric sample identifier.
tsne_data_wt <- tsne.data %>% 
  filter(grepl("-4$|-5$|-7$|-9$|-10$", rownames(tsne.data)))
tsne_data_wt$sample_id = sapply(strsplit(tsne_data_wt$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
## Summarize the number of cells per sample as well as the proportion of each cell type per patient.
wt_clust_annot = tsne_data_wt %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(wt_clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
wt_clust_annot <- wt_clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- wt_clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_wt <- log2cpm[ , colnames(log2cpm)%in%cell_names_keep] 
all(colnames(log2cpm_wt)==wt_clust_annot$cell_name)

## Only raw counts.
raw_cpm_wt = exp(log2cpm_wt[c(1:24703), ])-1
scgp_wt <- CreateSeuratObject(counts = raw_cpm_wt, min.cells = 10, project = "IDHwt", names.field = 1, names.delim = "_")
scgp_wt@meta.data$cell_type <- wt_clust_annot$cell_type
avail_genes = rownames(scgp_wt@assays$RNA@data)

## Make a Seurat object with the standard pipeline through PCA.
scgp_wt <- scgp_wt %>% 
  Seurat::NormalizeData(verbose = TRUE) %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

rna_var_wt = scgp_wt@assays$RNA@meta.features
rna_var_wt$gene_id <- rownames(rna_var_wt)

## Specific to IDHwt
## epimut_IDHwt_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_genebody_total_filt-20210207.txt", header = T, sep="\t")

## Combine with epimutation data for all .
rna_var_wt_epi = rna_var_wt %>% 
  inner_join(epimut_gb_total_filt, by="gene_id")

## Filter to about 10% of the single-cells from each subtype.
# epimut_IDHwt_total_filt = epimut_IDHwt_total %>% 
#  filter(IDHwt_cells_cov > 50)
# rna_var_epi = rna_var %>% 
#  inner_join(epimut_IDHwt_total_filt, by="gene_id")

## Break into groups based on levels of epimutation. 
rna_var_wt_epi$epimut_groups <-cut(rna_var_wt_epi$tumor_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1.1))
table(rna_var_wt_epi$epimut_groups)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/IDHwt-genebody-epimutation-expression.pdf", width = 7, height = 5)
ggplot(rna_var_wt_epi, aes(x=epimut_groups, y=mvp.mean)) + geom_boxplot() +
  labs(x="Gene body epimutation groups", y="Mean expression") +
  #ylim(0, 1.5) + 
  plot_theme +
  #stat_compare_means(comparisons = my_comparisons)
  stat_compare_means(method="kruskal.test") + 
  stat_n_text(y.pos = 1.25)
dev.off()




###########################
### Combine plots (color w/subtype designation)
###########################
rna_var_all_epi = rna_var_wt_epi %>%
  inner_join(rna_var_mut_epi, by=c("gene_id", "tumor_pdr", "tumor_cov", "tumor_class", "epimut_groups"))
colnames(rna_var_all_epi) <- gsub("\\.x", ".IDHwt", colnames(rna_var_all_epi))
colnames(rna_var_all_epi) <- gsub("\\.y", ".IDHmut", colnames(rna_var_all_epi))


rna_var_all_epi_express = rna_var_all_epi %>% 
  gather(IDHstatus, mean_expression, c("mvp.mean.IDHwt", "mvp.mean.IDHmut")) %>% 
  mutate(IDHstatus = recode(IDHstatus, "mvp.mean.IDHwt" = "IDHwt", "mvp.mean.IDHmut" = "IDHmut"))


pdf(file = "github/results/Fig2/Fig5a-violin-gb=expression-disorder-violin-restricted.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
ggplot(rna_var_all_epi_express, aes(x=epimut_groups, y=mean_expression, fill=IDHstatus)) + 
  geom_violin(width=1.1) +
  geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
  ylim(0, 2) + 
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Gene body DNAme disorder groups", y="Mean expression log(cpm)", fill="Glioma subtype") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free")
dev.off()


### END ####





