##################################
# Analyze transcriptional diversity differences amongst subclones.
# Updated: 2021.05.19
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
library(tidyverse)
library(Matrix)
library(openxlsx)
#library(devtools)
#devtools::install_local("data/CytoTRACE-master.zip")
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


## Load the 10X data for all tumor samples.
load("data/analysis_scRNAseq_tumor_gene_expression.Rds")

## Change to HUGO gene name.
rownames(expr_norm_data)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## 2D UMAP coordinates.
umap_metadata <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

# Limit to the sample of interest.
umap_data_sm006 <- umap_metadata %>% 
  filter(case_barcode=="SM006")

## Keep only tumor cells.
cells_to_keep = which(umap_data_sm006$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot <- umap_data_sm006[cells_to_keep, ]

# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_barcode

## Restrict the RNAseq data to the tumor cells.
expr_norm_data_filt <- expr_norm_data[ 1:24703, colnames(expr_norm_data)%in%cell_names_keep]


## Load in the 10X CNV annotation for SM012.
cell_annotation_sm006 <- read.table("data/SM006-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)


## Run CytoTRACE.
all(colnames(expr_norm_data_filt)==cell_annotation_sm006$cell_name)
SM006_results <- CytoTRACE(expr_norm_data_filt, ncores = 8)
pheno_vec = cell_annotation_sm006$ecDNA
names(pheno_vec) = cell_annotation_sm006$cell_name

## Create an output for the data.
cyto_df = as.data.frame(SM006_results$CytoTRACE)
colnames(cyto_df) <- "CytoTRACE"
cyto_df$cell_name = rownames(cyto_df)



###################
###### SM012 ######
###################
# Limit to the sample of interest.
umap_data_sm012 <- umap_metadata %>% 
  filter(case_barcode=="SM012")

## Keep only tumor cells.
cells_to_keep_sm012 = which(umap_data_sm012$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot_sm012 <- umap_data_sm012[cells_to_keep_sm012, ]

# Extract the exact name of cells.
cell_names_keep_sm012 <- clust_annot_sm012$cell_barcode

## Restrict the RNAseq data to the tumor cells.
expr_norm_data_sm012 <- expr_norm_data[ 1:24703, colnames(expr_norm_data)%in%cell_names_keep_sm012]


## Load in the 10X CNV annotation for SM012.
cell_annotation_sm012 <- read.table("data/SM012-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
all(colnames(expr_norm_data_sm012)==cell_annotation_sm012$cell_name)
SM012_results <- CytoTRACE(expr_norm_data_sm012, ncores = 8)
pheno_vec = cell_annotation_sm012$ecDNA
names(pheno_vec) = cell_annotation_sm012$cell_name

## Output the CytoTRACE metric.
cyto_df_sm012 = as.data.frame(SM012_results$CytoTRACE)
colnames(cyto_df_sm012) <- "CytoTRACE"
cyto_df_sm012$cell_name = rownames(cyto_df_sm012)


###################
###### SM017 ######
###################
# Limit to the sample of interest.
umap_data_sm017 <- umap_metadata %>% 
  filter(case_barcode=="SM017")

## Keep only tumor cells.
cells_to_keep_sm017 = which(umap_data_sm017$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot_sm017 <- umap_data_sm017[cells_to_keep_sm017, ]

# Extract the exact name of cells.
cell_names_keep_sm017 <- clust_annot_sm017$cell_barcode

## Restrict the RNAseq data to the tumor cells.
expr_norm_data_sm017 <- expr_norm_data[ 1:24703, colnames(expr_norm_data)%in%cell_names_keep_sm017]


## Load in the 10X CNV annotation for SM017.
cell_annotation_sm017 <- read.table("data/SM017-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
all(colnames(expr_norm_data_sm017)==cell_annotation_sm017$cell_name)
SM017_results <- CytoTRACE(expr_norm_data_sm017, ncores = 8)
pheno_vec = cell_annotation_sm017$ecDNA
names(pheno_vec) = cell_annotation_sm017$cell_name

## Output the CytoTRACE metric.
cyto_df_sm017 = as.data.frame(SM017_results$CytoTRACE)
colnames(cyto_df_sm017) <- "CytoTRACE"
cyto_df_sm017$cell_name = rownames(cyto_df_sm017)


###################
###### SM018 ######
###################
# Limit to the sample of interest.
umap_data_sm018 <- umap_metadata %>% 
  filter(case_barcode=="SM018")

## Keep only tumor cells.
cells_to_keep_sm018 = which(umap_data_sm018$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot_sm018 <- umap_data_sm018[cells_to_keep_sm018, ]

# Extract the exact name of cells.
cell_names_keep_sm018 <- clust_annot_sm018$cell_barcode

## Restrict the RNAseq data to the tumor cells.
expr_norm_data_sm018 <- expr_norm_data[ 1:24703, colnames(expr_norm_data)%in%cell_names_keep_sm018]


## Load in the 10X CNV annotation for SM018.
cell_annotation_sm018 <- read.table("data/SM018-10X-cnv-classification.txt", sep="\t", header=T, stringsAsFactors = F)

## Run CytoTRACE.
all(colnames(expr_norm_data_sm018)==cell_annotation_sm018$cell_name)
SM018_results <- CytoTRACE(expr_norm_data_sm018, ncores = 8)
pheno_vec = cell_annotation_sm018$ecDNA
names(pheno_vec) = cell_annotation_sm018$cell_name

## Output the CytoTRACE metric.
cyto_df_sm018 = as.data.frame(SM018_results$CytoTRACE)
colnames(cyto_df_sm018) <- "CytoTRACE"
cyto_df_sm018$cell_name = rownames(cyto_df_sm018)


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

pdf(file = "results/Fig6/SuppFig10h-ecDNA-clone-transcriptional-diversity.pdf", width = 9, height = 4, useDingbats = FALSE)
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