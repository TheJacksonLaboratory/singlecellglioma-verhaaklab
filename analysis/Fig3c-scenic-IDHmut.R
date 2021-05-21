##################################
# Plot TFs that regulate gene activity (IDHmut samples).
# Updated: 2021.05.16
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/mnt/verhaak-lab/scgp/results/scenic/IDHmut/"
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
library(ggrepel)
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


## Load the 10X data for all samples.
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds")

## Change to HUGO gene name.
rownames(log2cpm)[1:24703] <- featuredata[1:24703, "Associated.Gene.Name"]

## Define the cell clusters of interest (i.e., tumor cells inferred from marker genes using CellView and confirmed by CNVs).
tsne.data$cell_name = rownames(tsne.data)

## The object says, "tSNE" but really these are the original UMAP coordinates.
# Limit to IDHmut samples.
tsne_data_mut <- tsne.data %>% 
  filter(grepl("-0$|-1$|-2$|-3$|-6$|-8$", rownames(tsne.data)))
tsne_data_mut$sample_id = sapply(strsplit(tsne_data_mut$cell_name, "-"), "[[", 3)

## Keep only tumor cells.
clust_annot = tsne_data_mut %>% 
  mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                            `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                            `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) 
cells_to_keep = which(clust_annot$cell_type%in%c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor"))
clust_annot <- clust_annot[cells_to_keep, ]
# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_name

# Create sample-specific labels for each patient to make it easier to refer back.
clust_annot = clust_annot %>% 
  mutate(sample_id = recode(sample_id, `0` = "UC917", `1` = "SM001", `2` = "SM002", `3` = "SM004",
                            `6` = "SM008", `8` = "SM015"))
sample_clust <- as.factor(clust_annot$sample_id)
names(sample_clust) <- clust_annot$cell_name

## Restrict the RNAseq data to the tumor cells.
log2cpm <- log2cpm[ , colnames(log2cpm)%in%cell_names_keep] 
all(colnames(log2cpm)==clust_annot$cell_name)

## Downsample for input into SCENIC.
set.seed(43)
down_sample = sample(ncol(log2cpm), 5000)
log2cpm <- log2cpm[ , down_sample]
clust_annot <- clust_annot[down_sample, ]
# Make sure that the same cells are being subsetted.
all(colnames(log2cpm)==clust_annot$cell_name)

## Enumerate the mitochondrial reads.
clust_annot$mito = as.numeric(log2cpm[24706, ])


################################
### Re-load SCENIC results
################################
#scenicOptions <- readRDS("int/scenicOptions.Rds")
#cellInfo <- readRDS("int/cellInfo.Rds")
#colVars <- readRDS("int/colVars.Rds")
auc_rankings <- readRDS("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/data/SCENIC/IDHmut/3.3_aucellRankings.Rds")
regulonAUC <- readRDS("/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/data/SCENIC/IDHmut/3.4_regulonAUC.Rds")

## Create a data.frame with the gene sets/TFs and cells.
regulonAUC_df = as.data.frame(getAUC(regulonAUC))

## The annotation files we have match the regulonAUC data.
all(clust_annot$cell_name==colnames(regulonAUC_df))

## generate z-scores for variable A using the scale() function
## scale(A, center = TRUE, scale = TRUE). These are the defaults. 
regulonAUC_scaled = t(apply(as.matrix(regulonAUC_df), 1, scale))

## Load the IDHmut pseudotime data:
pseudotime_IDHmut = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/10X/monocle/monocle3_IDHmut_pseudotime.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)

## Add the pseudotime data as a covariate.
clust_annot = clust_annot %>% 
  left_join(pseudotime_IDHmut, by="cell_name")
all(clust_annot$cell_name==colnames(regulonAUC_scaled))

## Provide the "sample_id", "cell_type", and mitochondrial percent annotations for each cell.
cell_state = gsub("_tumor", "", clust_annot$cell_type)
sample_id = clust_annot$sample_id
mito_pct = clust_annot$mito
pseudotime = clust_annot$pseudotime
subtype = rep("IDHmut", dim(clust_annot)[1])
annot_df = data.frame(subtype, sample_id, cell_state)
#mito_cols <- colorRamp2(c(0, 5, 10, 15, 20), c("#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe"))
## viridis(7)
#pt_cols <- colorRamp2(c(0, 2, 4, 6, 8, 10, 12), c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF",
#                                                  "#8FD744FF", "#FDE725FF"))
epimut_cols <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1.0), c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))


## Define the annotation colors:
ha = HeatmapAnnotation(df = annot_df,
                       col = list(sample_id = c("SM001" = "#F8766D",
                                    "SM002" = "#DB8E00",
                                    "SM004" = "#AEA200",
                                    "SM008" = "#00BD5C",
                                    "SM015" = "#00A6FF",
                                    "UC917" = "#FF63B6"),
                                  cell_state = c("differentiated" = "#fcbba1",
                                                "stemcell" = "#fb6a4a",
                                                "prolif_stemcell" = "#a50f15"),
                                  subtype = c("IDHmut" = "#AF8DC3")))


## What are the high epimutation binding sites for TFs?
tfbs_epimut = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/tfbs_epimutation_subtype.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
tfbs_epimut = tfbs_epimut %>% 
  select(tf, IDHmut) %>% 
  distinct()
scenic_tfs = data.frame(tf = sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1)) 

scenic_tfs_epimut = scenic_tfs %>% 
  left_join(tfbs_epimut, by="tf") %>% 
  mutate(high_tf = ifelse(IDHmut > 0.399, "high", NA)) %>% 
  select(tf, epimut = IDHmut, high_tf)
row_ha = rowAnnotation(epimut = scenic_tfs_epimut$epimut)

## Set z-score limits:
#regulonAUC_scaled[regulonAUC_scaled < -4] = -4
#regulonAUC_scaled[regulonAUC_scaled > 4] = 4
set.seed(43)
Heatmap(regulonAUC_scaled, name = "TF activity\nZ-score",row_km = 5, column_km = 3, col = colorRamp2(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), magma(17)),
        top_annotation = ha,
        right_annotation = row_ha,
        show_row_names = FALSE, show_column_names = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = 8))

## Apply kruskal.wallis across three cell types. Select X number of TFs
kw_test_results = row_kruskalwallis(regulonAUC_scaled, cell_state)
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

## Or take the 10 most enriched TFs per cell state
diff_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="differentiated"], 1, median), decreasing = TRUE)
stem_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="stemcell"], 1, median), decreasing = TRUE)
prol_stem_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="prolif_stemcell"], 1, median), decreasing = TRUE)
ranked_tfs_to_keep = unique(c(names(prol_stem_rank[1:15]), names(stem_rank[1:15]), names(diff_rank[1:15])))

## Filter the regulonAUC plot.
regulonAUC_scaled_filt = regulonAUC_scaled[ranked_tfs_to_keep, ]
rownames(regulonAUC_scaled_filt) <- gsub("_extended", "", rownames(regulonAUC_scaled_filt))

## Include epimutation as side panel. Keep TFs in the same order as above.
scenic_tfs_epimut_sub = scenic_tfs_epimut[rownames(regulonAUC_scaled)%in%ranked_tfs_to_keep, ]
## Set the TFs so that they are in the same order:
sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)==scenic_tfs_epimut_sub$tf[match(sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1), scenic_tfs_epimut_sub$tf)]
scenic_tfs_epimut_sub_ord = scenic_tfs_epimut_sub[match(sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1), scenic_tfs_epimut_sub$tf), ]
## The annotation should now be correct:
row_ha = rowAnnotation(epimut = scenic_tfs_epimut_sub_ord$epimut,
                       col = list(epimut = epimut_cols))


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig3/Fig3c-IDHmut-TF-activity.pdf", height = 6.545454, width = 9, bg = "transparent", useDingbats = FALSE)
set.seed(43)
Heatmap(regulonAUC_scaled_filt, name = "TF activity\nZ-score", row_km = 1, column_km = 1, col = colorRamp2(c(-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), viridis(17)),
        top_annotation = ha, 
        right_annotation = row_ha,
        show_row_names = TRUE, show_column_names = FALSE,
        show_column_dend = FALSE, show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 8))
dev.off()

## Objective: Create boxplots for IDHmut TFBS epimutation data at active TFs per cell state:
diff_tfs = data.frame(tf = unique(sapply(strsplit(names(diff_rank), "_| "), "[[", 1))) 
diff_tfs$diff_activity_rank = seq_len(nrow(diff_tfs))
diff_tfs_epimut = diff_tfs %>% 
  left_join(tfbs_epimut, by="tf") %>% 
  filter(!is.na(IDHmut)) %>% 
  mutate(state = "diff") %>% 
  select(tf, activity_rank = diff_activity_rank, state, epimut = IDHmut)

stem_tfs = data.frame(tf = unique(sapply(strsplit(names(stem_rank), "_| "), "[[", 1))) 
stem_tfs$stem_activity_rank = seq_len(nrow(stem_tfs))
stem_tfs_epimut = stem_tfs %>% 
  left_join(tfbs_epimut, by="tf") %>% 
  filter(!is.na(IDHmut)) %>% 
  mutate(state = "stem") %>% 
  select(tf = tf, activity_rank =stem_activity_rank, state, epimut = IDHmut)

prol_stem_tfs = data.frame(tf = unique(sapply(strsplit(names(prol_stem_rank), "_| "), "[[", 1))) 
prol_stem_tfs$prol_stem_activity_rank = seq_len(nrow(prol_stem_tfs))
prol_stem_tfs_epimut = prol_stem_tfs %>% 
  left_join(tfbs_epimut, by="tf") %>% 
  filter(!is.na(IDHmut)) %>%
  mutate(state = "prolif_stem") %>% 
  select(tf, activity_rank = prol_stem_activity_rank, state, epimut = IDHmut)


top_tf_epimut <- bind_rows(diff_tfs_epimut[1:15, ],
                           stem_tfs_epimut[1:15, ],
                           prol_stem_tfs_epimut[1:15, ])

## Test whether the general epimutation rate differs between these three cell states.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/IDHmut-states-TF-epimutation.pdf", height = 4, width = 6, bg = "transparent", useDingbats = FALSE)
ggplot(top_tf_epimut, aes(x=state, y=epimut)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(alpha = 0.4), width = 0.1) +
  stat_compare_means(method = "kruskal") +
  labs(y="High TF activity per cell state (RNA)\nTFBS epimutation burden (DNAm)", x="Cell state") +
  plot_theme +
  guides(alpha=FALSE)
dev.off()


## Recode plot so that it's comparbale.
top_tf_epimut = top_tf_epimut %>% 
  mutate(state = recode(state, "diff" = "Diff.-like",
                        "stem" = "Stem-like",
                        "prolif_stem" = "Prolif. stem-like"))

mu = top_tf_epimut %>% 
  group_by(state) %>% 
  summarise(grp_mu = mean(epimut))

## ks.test for diff.-like versus others.
ks.test(top_tf_epimut$epimut[top_tf_epimut$state=="Diff.-like"], top_tf_epimut$epimut[top_tf_epimut$state=="Stem-like"])
ks.test(top_tf_epimut$epimut[top_tf_epimut$state=="Diff.-like"], top_tf_epimut$epimut[top_tf_epimut$state=="Prolif. stem-like"])

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig3/SuppFig6-IDHmut-states-TF-disorder-density.pdf", height = 4, width = 5, bg = "transparent", useDingbats = FALSE)
ggplot(top_tf_epimut, aes(x= epimut, fill=state)) +
  geom_vline(data=mu, aes(xintercept=grp_mu, color=state),
             linetype="dashed") +
  scale_fill_manual(values = c("Diff.-like" = "#fcbba1",
                               "Stem-like" = "#fb6a4a",
                               "Prolif. stem-like" = "#a50f15")) +
  scale_color_manual(values = c("Diff.-like" = "#fcbba1",
                                "Stem-like" = "#fb6a4a",
                                "Prolif. stem-like" = "#a50f15")) +
  geom_density(alpha=0.75) +
  labs(y="Density", x ="TFBS motif DNAme disorder (DNAm)\nHigh TF activity per cell state (RNA)", fill="Cell state") +
  guides(color=FALSE) +
  plot_theme +
  theme(legend.position="bottom") +
  annotate(geom="text", x=.5, y=18, label="Kolmogorov-Smirnov",
           color="black") +
  annotate(geom="text", x=.5, y=16, label="Diff vs. Stem p=0.66",
           color="black") +
  annotate(geom="text", x=.5, y=14, label="Diff vs. Polif. stem p=0.92",
           color="black")
dev.off()


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

ggplot(median_activity_df_filt, aes(tf, cell_state)) +
  geom_tile(aes(fill=activity)) +
  scale_fill_gradientn(colours=c("#0571b0", "#92c5de", "#f4a582","#ca0020"), values=c(0, 0.5, 1, 1.5), na.value="white") +
  labs(x="", y= "", fill="TF activity") +
  plot_theme

## Set the TF order to be ranked for each cell state.
tf_order <- sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)

tf_levels = tf_order
state_levels = c("Prolif. stem-like", "Stem-like", "Diff.-like")
median_activity_df_filt$tf <-  factor(median_activity_df_filt$tf, levels = tf_levels)
median_activity_df_filt$cell_state <-  factor(median_activity_df_filt$cell_state, levels = rev(state_levels))

null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig3/Fig3c-IDHmut-TF-activity.pdf", height = 4, width = 7.5, bg = "transparent", useDingbats = FALSE)
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