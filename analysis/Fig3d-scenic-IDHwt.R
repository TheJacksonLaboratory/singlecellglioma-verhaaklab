##################################
# Determine + visualize TFs that regulate cell states (IDHwt samples).
# Updated: 2021.05.16
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(NMF)
library(psycho)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(matrixTests)
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

## Load the 10X data for all tumor samples.
load("data/analysis_scRNAseq_tumor_gene_expression.Rds")

## Change to HUGO gene name.
rownames(expr_norm_data)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## 2D UMAP coordinates.
umap_coords_2d <- read.csv("data/analysis_scRNAseq_tumor_metadata.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

# Limit to IDHwt samples.
umap_data_wt <- umap_coords_2d %>% 
  filter(case_barcode%in%c("SM006", "SM012", "SM017", "SM018", "SM011"))

## Keep only tumor cells.
cells_to_keep = which(umap_data_wt$cell_state%in%c("Diff.-like", "Stem-like", "Prolif. stem-like"))
clust_annot <- umap_data_wt[cells_to_keep, ]

# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_barcode

## Restrict the RNAseq data to the tumor cells.
expr_norm_data <- expr_norm_data[ ,colnames(expr_norm_data)%in%cell_names_keep] 
all(colnames(expr_norm_data)==clust_annot$cell_barcode)

## Downsample for input into SCENIC to make the cell numbers consistent across analyses + run time with R SCENIC version.
set.seed(34)
down_sample = sample(ncol(expr_norm_data), 5000)
expr_norm_data_sample <- expr_norm_data[ , down_sample]
clust_annot_sample <- clust_annot[down_sample, ]
# Make sure that the same cells are being subsetted.
all(colnames(expr_norm_data_sample)==clust_annot_sample$cell_barcode)

## Untransformed data.
raw_cpm = exp(expr_norm_data_sample[c(1:24703), ])-1

## Create a Seurat object using raw counts.
scgp <- CreateSeuratObject(counts = raw_cpm, min.cells = 1, project = "scgp_wt", names.field = 1, names.delim = "_")

## Use broad tumor cell classification.
scgp@meta.data$case_barcode <- clust_annot_sample$case_barcode
scgp@meta.data$cell_state <- clust_annot_sample$cell_state


##########################
### Begin SCENIC approach
##########################
## Building the **gene regulatory network (GRN)**: 
## 1. Identify potential targets for each TF based on co-expression.
# - Filtering the expression matrix and running GENIE3/GRNBoost. 
# - Formatting the targets from GENIE3/GRNBoost into co-expression modules. 

## Initialize SCENIC settings:
org="hgnc" 
dbDir="reference/cisTarget_databases" 
myDatasetTitle="SCENIC IDHwt" 
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Save to use at a later time.
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## Load expression matrix.
exprMat <- data.matrix(scgp@assays$RNA@counts)
cellInfo <- data.frame(scgp@meta.data$case_barcode)
rownames(cellInfo) <- rownames(scgp@meta.data)
colnames(cellInfo) <- "CellType"

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("SM006"="#64B200", 
                           "SM012"="#00BADE", 
                           "SM017"="#B385FF", 
                           "SM018"="#EF67EB", 
                           "SM011"="#00C1A7"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

## Save outputs for cell annotation and color code. 
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="int/colVars.Rds")


## Examine how many genes have greater than 0.
cellInfo$nGene <- colSums(exprMat>0)

## Filter based on the number of genes.
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)


## Filter the expression matrix only to keep these genes.
exprMat_filtered <- exprMat[genesKept, ]


## Split the targets into positive- and negative-correlated targets 
## (i.e. Spearman correlation between the TF and the potential target).
runCorrelation(exprMat_filtered, scenicOptions)

## Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

### Run GENIE3 (this is computationally intensive). 
runGenie3(exprMat_filtered, scenicOptions)


## 2.  Select potential direct-binding targets (regulons) based on DNA-motif analysis (*RcisTarget*: TF motif analysis) 
## Build and score the GRN.
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

## 1. Get co-expression modules.
runSCENIC_1_coexNetwork2modules(scenicOptions)

## 2. Get regulons (with RcisTarget): TF motif analysis).
runSCENIC_2_createRegulons(scenicOptions)

## 3. Score GRN (regulons) in the cells (with AUCell).
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

## 4. Determine the binarized activities.
runSCENIC_4_aucell_binarize(scenicOptions, skipBoxplot = FALSE, skipHeatmaps = FALSE,
                            skipTsne = FALSE, exprMat = exprMat_filtered)



################################
### Re-load SCENIC results - that can be derived from scripts above
################################
auc_rankings <- readRDS("data/SCENIC/IDHwt/3.3_aucellRankings.Rds")
regulonAUC <- readRDS("data/SCENIC/IDHwt/3.4_regulonAUC.Rds")


## Create a data.frame with the gene sets/TFs and cells.
regulonAUC_df = as.data.frame(getAUC(regulonAUC))

## The annotation files we have match the regulonAUC data.
all(clust_annot_sample$cell_barcode==colnames(regulonAUC_df))

## generate z-scores for variable A using the scale() function
## scale(A, center = TRUE, scale = TRUE). These are the defaults. 
regulonAUC_scaled = t(apply(as.matrix(regulonAUC_df), 1, scale))

## Provide the case_barcode, cell_state, and subtype annotations for each cell.
annot_df = data.frame(clust_annot_sample$cell_barcode, rep("IDHwt", dim(clust_annot_sample)[1]), clust_annot_sample$case_barcode, clust_annot_sample$cell_state)
colnames(annot_df) <- c("barcode", "subtype", "case_barcode", "cell_state")
epimut_cols <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1.0), c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))

## Define the annotation colors:
ha = HeatmapAnnotation(df = annot_df,
                       col = list(sample_id = c("SM006" = "#64B200",
                                                "SM011" = "#00C1A7",
                                                "SM012" = "#00BADE",
                                                "SM017" = "#B385FF",
                                                "SM018" = "#EF67EB"),
                                  cell_state = c("differentiated" = "#fcbba1",
                                                 "stemcell" = "#fb6a4a",
                                                 "prolif_stemcell" = "#a50f15"),
                                  subtype = c("IDHwt" = "#7FBF7B")))


## What are the TFs considered in this analysis?
scenic_tfs = data.frame(tf = sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1)) 

## Quick visualization.
set.seed(43)
Heatmap(regulonAUC_scaled, name = "TF activity\nZ-score",row_km = 5, column_km = 3, col = colorRamp2(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), magma(17)),
        top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = 8))

## Apply kruskal.wallis across three cell types. Select X number of TFs
kw_test_results = row_kruskalwallis(regulonAUC_scaled, clust_annot_sample$cell_state)
kw_test_results$adj_pvalue = p.adjust(kw_test_results$pvalue, method = "fdr", n = length(kw_test_results$pvalue))
kw_test_results_diff = kw_test_results[which(kw_test_results$pvalue< 1e-127), ]
tfs_to_keep = rownames(kw_test_results_diff)

## Remove any "_extended" TFs that are also represented.
tmp = sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1) %>% as.data.frame() 
colnames(tmp) <- "tf"
duplicated_tf = tmp %>% 
  group_by(tf) %>% 
  summarise(tf_counts = n()) %>% 
  filter(tf_counts >= 2)
dup_tfs_remove = rownames(regulonAUC_scaled)[which(sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1)%in%duplicated_tf$tf & grepl("_extended", rownames(regulonAUC_scaled)))]
tfs_to_keep_uniq = rownames(regulonAUC_scaled)[!rownames(regulonAUC_scaled)%in%dup_tfs_remove]

## Take the 15 most enriched TFs per cell state. Some may overlap leading there to be fewer than 45 across the two cell states
diff_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, clust_annot_sample$cell_state=="Diff.-like"], 1, median), decreasing = TRUE)
stem_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, clust_annot_sample$cell_state=="Stem-like"], 1, median), decreasing = TRUE)
prol_stem_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, clust_annot_sample$cell_state=="Prolif. stem-like"], 1, median), decreasing = TRUE)
ranked_tfs_to_keep = unique(c(names(prol_stem_rank[1:15]), names(stem_rank[1:15]), names(diff_rank[1:15])))

## Filter the regulonAUC plot.
regulonAUC_scaled_filt = regulonAUC_scaled[ranked_tfs_to_keep, ]
rownames(regulonAUC_scaled_filt) <- gsub("_extended", "", rownames(regulonAUC_scaled_filt))

## Gather activity by cell state.
prolif_stem_df <- as.data.frame(prol_stem_rank)
prolif_stem_df$tf <- sapply(strsplit(rownames(prolif_stem_df), " "), "[[", 1)
prolif_stem_df$tf  <- reorder(prolif_stem_df$tf , prolif_stem_df$prol_stem_rank)

diff_rank_df <- as.data.frame(diff_rank)
diff_rank_df$tf <- sapply(strsplit(rownames(diff_rank_df), " "), "[[", 1)
diff_rank_df$tf  <- reorder(diff_rank_df$tf , diff_rank_df$diff_rank)

stem_rank_df <- as.data.frame(stem_rank)
stem_rank_df$tf <- sapply(strsplit(rownames(stem_rank_df), " "), "[[", 1)
stem_rank_df$tf  <- reorder(stem_rank_df$tf , stem_rank_df$stem_rank)

### Create a matrix of summary values:
median_activity_df <- stem_rank_df %>% 
  inner_join(prolif_stem_df, by="tf") %>% 
  inner_join(diff_rank_df, by="tf") %>% 
  select(tf, stem_rank, prol_stem_rank, diff_rank)

median_activity_df_filt <- median_activity_df %>% 
  pivot_longer(cols= c(stem_rank:diff_rank),
               names_to = "cell_state",
               values_to = "activity") %>% 
  mutate(cell_state = recode(cell_state, `stem_rank` = "Stem-like",
                             `prol_stem_rank` = "Prolif. stem-like",
                             `diff_rank` = "Diff.-like")) %>% 
  mutate(tf = gsub("_extended", "", tf)) %>% 
  filter(tf%in%c(sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)))


## Set the TF order to be ranked for each cell state.
tf_order <- sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)

tf_levels = tf_order
state_levels = c("Prolif. stem-like", "Stem-like", "Diff.-like")
median_activity_df_filt$tf <-  factor(median_activity_df_filt$tf, levels = tf_levels)
median_activity_df_filt$cell_state <-  factor(median_activity_df_filt$cell_state, levels = rev(state_levels))

null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

pdf(file = "results/Fig3/Fig3d-IDHwt-TF-activity.pdf", height = 4, width = 7.5, bg = "transparent", useDingbats = FALSE)
ggplot(median_activity_df_filt, aes(tf, cell_state)) +
  geom_tile(aes(fill=activity)) +
  scale_fill_viridis() +
  labs(x="", y= "", fill="Relative TF activity\n(Z-score)") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10),
        legend.position="bottom") +
  null_y
dev.off()

### END ###