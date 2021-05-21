##############################################
# Assess single cell stress scores in 10X data.
# Updated: 2021.05.17
# Author: Kevin J.
###############################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###########################
library(tidyverse)
library(openxlsx)
library(AUCell)
library(GSEABase)
library(ggpubr)
###########################
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
## Martine provided the 2D coordinates from that run.
umap_coords_2d = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10X/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_umap_2d_embedding.csv")

## Test out the plotting coordinates.
ggplot(umap_coords_2d, aes(x= x.coordinate, y= y.coordinate)) + geom_point()

## Change to HUGO gene name.
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are original UMAP coordinates.
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
sample_id <- sapply(strsplit(clust_annot$cell_name, "-"), "[[", 3)
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
clust_annot$sample_id <- sample_id

## Keep only the tumor cells.
cells_to_keep = which(clust_annot$cell_type%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot_filt <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot_filt$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm_filt <- log2cpm[1:24703, colnames(log2cpm)%in%cell_names_keep]
all(colnames(log2cpm_filt)==clust_annot_filt$cell_name)

## Downsample for input.
set.seed(34)
down_sample = sample(ncol(log2cpm_filt), 2500)
log2cpm_filt_downsample <- log2cpm_filt[ , down_sample]
clust_annot_filt_downsample <- clust_annot_filt[down_sample, ]

## Load some example gene sets.
gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
geneSets <- getGmt(gmtFile)


## Define gene sets
## Load in the referent gene sets for stress, oxidative stress, hypoxia..
stress_genes <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/GO_REGULATION_OF_RESPONSE_TO_STRESS_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
stress_genes_list = list(stress_genes$V1)
ox_stress_genes <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/GO_RESPONSE_TO_OXIDATIVE_STRESS_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
ox_stress_genes_list = list(ox_stress_genes$V1)
harris_hypoxia_genes <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/HARRIS_HYPOXIA_geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
harris_hypoxia_list = list(harris_hypoxia_genes$V1)
ELVIDGE_hypoxia_genes <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/ELVIDGE_HYPOA_UP-geneset.txt", sep="\t", header=F, skip=2, stringsAsFactors = F)
ELVIDGE_hypoxia_list = list(ELVIDGE_hypoxia_genes$V1)

## Define the stress gene sets that were used in the bulk analyses.
stressGeneSets <- c(GSEABase::GeneSet(stress_genes$V1, setName="stress_set"),
  GSEABase::GeneSet(ox_stress_genes$V1, setName="ox_stress_set"),
  GSEABase::GeneSet(harris_hypoxia_genes$V1, setName="hypoxia_harris_set"),
  GSEABase::GeneSet(ELVIDGE_hypoxia_genes$V1, setName="hypoxia_elvidge_set"))

## How many genes are in each set?
#stressGeneSets <- GSEABase::subsetGeneSets(stressGeneSets[[1]], rownames(log2cpm_filt_downsample))

# Random gene sets.
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(log2cpm_filt), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(log2cpm_filt), 500), setName="Random (500g)"))

countsPerGene <- apply(log2cpm_filt, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(stressGeneSets,extraGeneSets))
names(geneSets)

## Step 1: Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(data.matrix(log2cpm_filt), nCores=1, plotStats=TRUE)
save(cells_rankings, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cells_stress_rankings.RData")

cells_rankings


## Step 2: Calculate enrichment for the gene signatures (AUC). 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/cells_stress_AUC.RData")

## Since the AUC represents the proportion of expressed genes in the signature, we can use the relative AUCs 
## across the cells to explore the population of cells that are present in the dataset according to the expression of the gene-set


set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

## Determine the stress scores for each sample.

## Create a data.frame with the gene sets/TFs and cells.
cells_AUC_df = as.data.frame(t(getAUC(cells_AUC)))
cells_AUC_df$cell_name <- rownames(cells_AUC_df)

clust_annot_filt_comb <- clust_annot_filt %>% 
  inner_join(cells_AUC_df, by="cell_name")

clust_annot_filt_comb$idh_status <- ifelse(clust_annot_filt_comb$sample_id%in%c("SM004", "SM001", "SM002", "SM015", "SM008", "UC917"), "IDHmut", "IDHwt")
case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
clust_annot_filt_comb$sample_id <- factor(clust_annot_filt_comb$sample_id, levels = case_order)
state_levels <- c("Prolif. stem-like", "Stem-like", "Diff.-like")
clust_annot_filt_comb$cell_type <- factor(clust_annot_filt_comb$cell_type, levels = state_levels)

clust_annot_filt_comb_non_prolif <- clust_annot_filt_comb %>% 
  filter(cell_type!="Prolif. stem-like")

clust_annot_filt_comb_long = clust_annot_filt_comb %>% 
  dplyr::select(sample_id, idh_status, cell_type, stress_set, hypoxia_elvidge_set, `Random (500g)`) %>% 
  pivot_longer(stress_set:`Random (500g)`,
             names_to = "gene_set", 
             values_to =  "auc") %>% 
  mutate(gene_set = recode(gene_set, stress_set = "Response to stress", hypoxia_elvidge_set = "Hypoxia (Elvidge)", `Random (500g)` = "Random"))
clust_annot_filt_comb_long$gene_set <- factor(clust_annot_filt_comb_long$gene_set, levels = c("Response to stress", "Hypoxia (Elvidge)", "Random"))

pdf("github/results/10X/auc_cell_stress_cellstates.pdf", width = 9, height = 5, bg = "transparent", useDingbats = FALSE)
ggplot(clust_annot_filt_comb_long, aes(x= sample_id, y=auc, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("Prolif. stem-like" = "#a50f15", 
                             "Stem-like" = "#fb6a4a",
                             "Diff.-like" = "#fcbba1")) +
  facet_grid(gene_set ~ idh_status, scales = "free") +
  stat_compare_means(method="kruskal", label = "p.signif", size = 2) +
  labs(x="", y="AUC - gene set (scRNAseq)", fill = "Cell state") +
  plot_theme
dev.off()

## Normalize to the median value to be able to combine across different patients.
clust_annot_filt_comb_long_norm <- clust_annot_filt_comb_long %>% 
  group_by(sample_id, gene_set) %>%
  mutate_each(funs(./median(., na.rm = TRUE)), auc) %>% 
  dplyr::select(sample_id:gene_set, auc_normalized = auc)
  
my_comparisons <- list( c("Diff.-like" , "Stem-like"), c("Diff.-like", "Prolif. stem-like")) 
clust_annot_filt_comb_long_norm$idh_status <- factor(clust_annot_filt_comb_long_norm$idh_status, levels = c("IDHmut", "IDHwt"))

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank())


pdf("github/results/Fig3/SuppFig6f-cell-stress.pdf", width = 6, height = 4, bg = "transparent", useDingbats = FALSE)
ggplot(clust_annot_filt_comb_long_norm, aes(x= cell_type, y=auc_normalized, fill = cell_type)) +
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  geom_hline(yintercept = 1, alpha=0.8, linetype=2 ) +
  scale_fill_manual(values=c("Prolif. stem-like" = "#a50f15", 
                             "Stem-like" = "#fb6a4a",
                             "Diff.-like" = "#fcbba1")) +
  facet_grid(idh_status~gene_set, scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.format", size = 3) +
  labs(x="", y="AUC gene set (scRNAseq)\nsubject median normalized", fill = "Pan-glioma\ncell state") +
  plot_theme +
  theme(legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust=1)) +
  null_x
dev.off()


## Cutoff the outliers for plotting purposes, but retain the p-values above.
clust_annot_filt_comb_long_norm_cutoff = clust_annot_filt_comb_long_norm %>% 
  filter(auc_normalized < 1.8)

pdf("github/results/Fig3/SuppFig6f-cell-stress-cutoff.pdf", width = 6, height = 4, bg = "transparent", useDingbats = FALSE)
ggplot(clust_annot_filt_comb_long_norm_cutoff, aes(x= cell_type, y=auc_normalized, fill = cell_type)) +
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  geom_hline(yintercept = 1, alpha=0.8, linetype=2 ) +
  scale_fill_manual(values=c("Prolif. stem-like" = "#a50f15", 
                             "Stem-like" = "#fb6a4a",
                             "Diff.-like" = "#fcbba1")) +
  facet_grid(idh_status~gene_set, scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.format", size = 3) +
  labs(x="", y="AUC gene set (scRNAseq)\nsubject median normalized", fill = "Pan-glioma\ncell state") +
  plot_theme +
  theme(legend.position="bottom", 
        axis.text.x = element_text(angle=45, hjust=1)) +
  null_x
dev.off()

### END ###