##############################################
# Investigate DNAme disorder and gene expression levels.
# Updated: 2021.05.13
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

## Load the SCGP raw counts matrix that has already been filtered for quality metrics.
scgp_obj <- ReadH5AD("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9-qc_20190827.h5ad")
scgp_obj@meta.data$cell_name <- rownames(scgp_obj@meta.data)
avail_genes = rownames(scgp_obj@assays$RNA@data)

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Filter to only IDHmut tumor cells.
meta_10x <- readRDS("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_metadata.Rds")

## Extract the numeric sample identifier for IDHmut samples.
meta_10x_mut <- meta_10x %>% 
  filter(grepl("-0$|-1$|-2$|-3$|-6$|-8$", rownames(meta_10x)))
meta_10x_mut$sample_id = sapply(strsplit(meta_10x_mut$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
mut_clust_annot = meta_10x_mut %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(mut_clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
mut_clust_annot <- mut_clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- mut_clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
mut_clust_annot = mut_clust_annot %>% 
  mutate(sample_id = sapply(strsplit(cell_name, "-"), "[[", 3),
         case_barcode = recode(sample_id, `0` = "SM019",  
                               `1` = "SM001", 
                               `2` = "SM002",
                               `3` = "SM004", 
                               `4` = "SM006", 
                               `5` = "SM011", 
                               `6` = "SM008",
                               `7` = "SM012", 
                               `8` = "SM015", 
                               `9` = "SM017", 
                               `10` = "SM018")) 

## Subset to IDHmut tumor cells.
scgp_obj_idh_mut_tumor <- subset(x = scgp_obj, subset = cell_name %in% mut_clust_annot$cell_name)


## Make a Seurat object with the standard pipeline through PCA.
scgp_obj_idh_mut_tumor <- scgp_obj_idh_mut_tumor %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

## Create a data.frame with the expression variability metrics.
rna_var_mut = scgp_obj_idh_mut_tumor@assays$RNA@meta.features
rna_var_mut$gene_id <- rownames(rna_var_mut)

## Read in previously calculated promoter epimutation matrix.
epimut_prom_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_promoter_total_filt-20210207.txt", header = T, sep="\t")
epimut_pomoter_IDHm_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_promoter_idhmut_filt-20210207.txt", header = T, sep="\t")

## Combine RNA data with promoter epimutation data.
rna_var_mut_epi = rna_var_mut %>% 
  inner_join(epimut_prom_total_filt, by="gene_id")
## Alternatively, just IDHmut genes.
#rna_var_mut_epi = rna_var_mut %>% 
#  inner_join(epimut_pomoter_IDHm_total_filt, by="gene_id")

## Break into groups based on levels of epimutation (low, intermediate, and high). 
rna_var_mut_epi$epimut_groups <-cut(rna_var_mut_epi$tumor_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))
#rna_var_mut_epi$epimut_groups <-cut(rna_var_mut_epi$IDHmut_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))

## How do these groups look in terms of total number of genes.
table(rna_var_mut_epi$epimut_groups) # Fewer genes with high epimutation.

## Plots for both gene expression and variability that is controlled for expression levels.
## The dispersion method helps control for the relationship between variability and average expression.
my_comparisons <- list( c("(-0.01,0.1]", "(0.1,0.2]"), c("(-0.01,0.1]", "(0.3,0.4]"), c("(0.1,0.2]", "(0.4,1.1]") )

## Calculate the summary values across groups:
rna_var_mut_epi %>% 
  group_by(epimut_groups) %>% 
  summarise(avg_expression = median(mvp.mean, na.rm = T), 
            avg_dispersion = median(mvp.dispersion.scaled, na.rm = T))

ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.dispersion.scaled)) + 
  geom_boxplot(outlier.shape = NA) +
  #ylim(-2.25, 2.25) +  
  labs(x="Promoter epimutation groups", y="Dispersion (mean-scaled)") +
  plot_theme + 
  stat_compare_means(comparisons = my_comparisons)  + 
  stat_n_text(y.pos = -3.25) 

ggplot(rna_var_mut_epi, aes(x=epimut_groups, y=mvp.mean)) + 
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method="kruskal.test") + 
  ylim(0, 0.75) + 
  labs(x="Promoter epimutation groups", y="Mean expression") +
  plot_theme + 
  stat_n_text(y.pos = 1)


######################
### IDHwt
######################
## Re-load the SCGP raw counts matrix that has already been filtered for quality metrics.
scgp_obj <- ReadH5AD("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9-qc_20190827.h5ad")
scgp_obj@meta.data$cell_name <- rownames(scgp_obj@meta.data)
avail_genes = rownames(scgp_obj@assays$RNA@data)

## Extract the numeric sample identifier for IDHmut samples.
meta_10x_wt <- meta_10x %>% 
  filter(grepl("-4$|-5$|-7$|-9$|-10$", rownames(meta_10x)))
meta_10x_wt$sample_id = sapply(strsplit(meta_10x_wt$cell_name, "-"), "[[", 3)

## Cell populations were manually annotated using marker genes. Note: that these are relatively broad classifications.
wt_clust_annot = meta_10x_wt %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(wt_clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
wt_clust_annot <- wt_clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- wt_clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
wt_clust_annot = wt_clust_annot %>% 
  mutate(sample_id = sapply(strsplit(cell_name, "-"), "[[", 3),
         case_barcode = recode(sample_id, `0` = "SM019",  
                               `1` = "SM001", 
                               `2` = "SM002",
                               `3` = "SM004", 
                               `4` = "SM006", 
                               `5` = "SM011", 
                               `6` = "SM008",
                               `7` = "SM012", 
                               `8` = "SM015", 
                               `9` = "SM017", 
                               `10` = "SM018")) 

## Subset to IDHwt tumor cells.
scgp_obj_idh_wt_tumor <- subset(x = scgp_obj, subset = cell_name %in% wt_clust_annot$cell_name)

## Make a Seurat object with the standard pipeline through PCA.
scgp_obj_idh_wt_tumor <- scgp_obj_idh_wt_tumor %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "mvp", nfeatures = length(avail_genes)) 

## Create a data.frame with the expression variability metrics.
rna_var_wt = scgp_obj_idh_wt_tumor@assays$RNA@meta.features
rna_var_wt$gene_id <- rownames(rna_var_wt)

## Read in previously calculated promoter epimutation matrix.
epimut_prom_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_promoter_total_filt-20210207.txt", header = T, sep="\t")
epimut_pomoter_IDHwt_total_filt <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_promoter_idhwt_filt-20210207.txt", header = T, sep="\t")

## Combine RNA data with promoter epimutation data.
rna_var_wt_epi = rna_var_wt %>% 
  inner_join(epimut_prom_total_filt, by="gene_id")
## Alternatively, just IDHmut defined epimut classes.
#rna_var_wt_epi = rna_var_wt %>% 
#  inner_join(epimut_pomoter_IDHwt_total_filt, by="gene_id") %>% 
#  filter(!is.na(IDHwt_pdr))

## Break into groups based on levels of epimutation (low, intermediate, and high). 
rna_var_wt_epi$epimut_groups <-cut(rna_var_wt_epi$tumor_pdr,c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))
## Alternative: just IDHwt cells.
#rna_var_wt_epi$epimut_groups <-cut(rna_var_wt_epi$IDHwt_pdr, c(-0.01, 0.1, 0.2, 0.3, 0.4, 1.1))

## How do these groups look in terms of total number of genes.
table(rna_var_wt_epi$epimut_groups) # Fewer genes with high epimutation.

## Plots for both gene expression and variability that is controlled for expression levels.
## The dispersion method helps control for the relationship between variability and average expression.
my_comparisons <- list( c("(-0.01,0.1]", "(0.1,0.2]"), c("(-0.01,0.1]", "(0.3,0.4]"), c("(0.1,0.2]", "(0.4,1.1]") )

## Calculate the summary values across groups:
rna_var_wt_epi %>% 
  group_by(epimut_groups) %>% 
  summarise(avg_expression = median(mvp.mean, na.rm = T), 
            avg_dispersion = median(mvp.dispersion.scaled, na.rm = T))

ggplot(rna_var_wt_epi, aes(x=epimut_groups, y=mvp.dispersion.scaled)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x="Promoter epimutation groups", y="Dispersion (mean-scaled)") +
  plot_theme + 
  stat_compare_means(comparisons = my_comparisons)  + 
  stat_n_text(y.pos = -3.25) 

ggplot(rna_var_wt_epi, aes(x=epimut_groups, y=mvp.mean)) + 
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method="kruskal.test") + 
  ylim(0, 0.75) + 
  labs(x="Promoter epimutation groups", y="Mean expression") +
  plot_theme + 
  stat_n_text(y.pos = 1)


###########################
### Combine plots (color w/subtype designation)
###########################
rna_var_all_epi = rna_var_wt_epi %>%
  inner_join(rna_var_mut_epi, by=c("gene_id", "tumor_pdr", "epimut_groups"))
colnames(rna_var_all_epi) <- gsub("\\.x", ".IDHwt", colnames(rna_var_all_epi))
colnames(rna_var_all_epi) <- gsub("\\.y", ".IDHmut", colnames(rna_var_all_epi))

rna_var_all_epi_express = rna_var_all_epi %>% 
  gather(IDHstatus, mean_expression, c("mvp.mean.IDHwt", "mvp.mean.IDHmut")) %>% 
  mutate(IDHstatus = recode(IDHstatus, "mvp.mean.IDHwt" = "IDHwt", "mvp.mean.IDHmut" = "IDHmut"))


pdf(file = "github/results/Fig2/Fig2a-violin-expression-disorder-restricted-violin.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
ggplot(rna_var_all_epi_express, aes(x=epimut_groups, y=mean_expression, fill=IDHstatus)) + 
  geom_violin(width=1.1) +
  geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
  ylim(0, 0.9) + 
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Promoter DNAme disorder groups", y="Mean expression log(cpm)", fill="Glioma subtype") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free")
dev.off()

pdf(file = "github/results/Fig2/Fig2a-violin-expression-disorder-restricted.pdf", width = 6, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(rna_var_all_epi_express, aes(x=epimut_groups, y=mean_expression, fill=IDHstatus)) + 
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 0.9) + 
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Promoter DNAme disorder groups", y="Mean expression log(cpm)", fill="Glioma subtype") +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free")
dev.off()



rna_var_all_epi_dispersion = rna_var_all_epi %>% 
  gather(IDHstatus, dispersion_scaled, c("mvp.dispersion.scaled.IDHwt", "mvp.dispersion.scaled.IDHmut")) %>%
  mutate(IDHstatus = recode(IDHstatus, "mvp.dispersion.scaled.IDHwt" = "IDHwt", "mvp.dispersion.scaled.IDHmut" = "IDHmut"))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/RNA-dispersion-promoter-epimutation-outliers-20210207.pdf", width = 9, height = 5)
ggplot(rna_var_all_epi_dispersion, aes(x=epimut_groups, y=dispersion_scaled, fill=IDHstatus)) + geom_boxplot() +
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Gene epimutation groups", y="Expression dispersion (mean-scaled)",  fill="Glioma subtype") +
  plot_theme +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free") 
dev.off()


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/RNA-dispersion-promoter-epimutation-20210207.pdf", width = 9, height = 5)
ggplot(rna_var_all_epi_dispersion, aes(x=epimut_groups, y=dispersion_scaled, fill=IDHstatus)) + geom_boxplot(outlier.shape = NA) +
  ylim(-1.75, 1.75) +  
  scale_fill_manual(values=c("IDHwt" = "#7fbf7b", 
                             "IDHmut" = "#af8dc3")) +
  labs(x="Gene epimutation groups", y="Expression dispersion (mean-scaled)",  fill="Glioma subtype") +
  plot_theme +
  facet_grid(. ~ IDHstatus, scales = "free_x", space = "free") 
dev.off()


### END ####
