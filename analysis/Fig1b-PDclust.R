##################################
# Unsupervised clustering of tumor scRRBS data.
# Updated: 2020.03.21
# Author: Kevin J.
###################################

## Based on tutorial: https://github.com/hui-tony-zk/PDclust

# HPC working directory for this analysis on Helix.
mybasedir = "/Users/johnsk/mnt/verhaak-lab/scgp/"
datadir  = "/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/rerun-dedup_bed_graph"
pattern   = '_pe.deduplicated.bismark.cov.gz$'
setwd(mybasedir)

###################################
library(tidyverse)
library(devtools)
library(openxlsx)
library(PDclust)
library(parallel)
library(gridExtra)
library(pheatmap)
###################################
## Plot theme to be used for consistency.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                 axis.text = element_text(size = 12),
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(),
                                                 panel.background = element_rect(fill = "transparent"),
                                                 axis.line = element_blank(),
                                                 strip.background =element_rect(fill="white"))

## Data are available on Synapse for single cells and 50-cell populations under the bismark coverage files.
## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)
rrbs_50_qc <- read.table(file="data/analysis_RRBS_50cell_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  bind_rows(rrbs_50_qc) %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(meta_data, by="case_barcode")  %>% 
  mutate(file_path = paste("/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/rerun-dedup_bed_graph/", cell_barcode, "_pe.deduplicated.bismark.cov.gz", sep=""))

## Note: For Synapse, all of these files have been condensed into a single coverage file.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)
files <- gsub("\\//", "\\/", files)
# Restrict to the single and 50-cell pops that passed our QC metrics.
files = files[files%in%rrbs_qc_pass$file_path]

# If it is desirable to include the sample names.
samples = data.frame(sample_id=gsub("_pe.deduplicated.bismark.cov.gz", "", basename(files)), library_type = substring(basename(files), 23, 23))

### Step 1: Set-up the data for PDclust.
scgp_files = mclapply(files, function(f){
  dat = tryCatch(data.table::fread(cmd = sprintf("zcat < %s", f),
                                   verbose = FALSE, showProgress = FALSE), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Create the `dat` to fit the PDclust designation. Remove scaffolds and sex chromosomes.
  dat <- dat[ , 1:4] 
  colnames(dat) <- c("chr", "start", "end", "meth")
  dat <- dat %>%
    filter(!grepl("GL|X|Y|MT", chr)) %>% 
    mutate(chr = paste("chr", chr, sep=""))
  
  return(dat)
  
}, mc.cores=20)

# Assign cell names to individual data.frames.
names(scgp_files) <- samples$sample_id

### Step 2: Perform pairwise comparisons:
system.time(
  scgp_files_pairwise <- create_pairwise_master(scgp_files, cores_to_use = 4, digital = FALSE)
)

# Save output as pairwise comparisons were a computationally intensive process.
#saveRDS(scgp_files_pairwise, "/Users/johnsk/mnt//verhaak-lab/scgp/results/scRRBS/human/PDclust/all-samples-PDclust-pairwise-20200320.rds")

### Step 3: Convert to a symmetrical matrix.
scgp_files_pairwise_matrix <- convert_to_dissimilarity_matrix(scgp_files_pairwise)
#saveRDS(scgp_files_pairwise_matrix, "/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/PDclust/all-samples-PDclust-pairwise-mat-20200320.rds")

## Read in results from 1) single-cells plus small populations:
pd_clust_pairwise = readRDS("/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/PDclust/all-samples-PDclust-pairwise-20200320.rds")
pd_clust_pairwise_mat = readRDS("/Users/johnsk/mnt/verhaak-lab/scgp/results/scRRBS/human/PDclust/all-samples-PDclust-pairwise-mat-20200320.rds")

## Plot the average number of CpGs overlapping across all comparisons.
summary(pd_clust_pairwise$num_cpgs) # average: 8,346

## Filter out the 50-cell samples to determine the average overlap amongst single cells.
pd_clust_pairwise_filt = pd_clust_pairwise %>% 
  filter(!grepl("-50P[0-9]", x)) %>% 
  filter(!grepl("-50P[0-9]", y))
summary(pd_clust_pairwise_filt$num_cpgs)

## Generate clusters with PDclust.
scgp_cluster_results <- cluster_dissimilarity(pd_clust_pairwise_mat, num_clusters = 6)

## Generate clusters for single-cells only
pd_clust_sc <- pd_clust_pairwise_mat[!grepl("-50P[0-9]", rownames(pd_clust_pairwise_mat)), !grepl("-50P[0-9]", colnames(pd_clust_pairwise_mat))]
scgp_cluster_sc_results <- cluster_dissimilarity(pd_clust_sc, num_clusters = 6)
hc <- scgp_cluster_sc_results$hclust_obj
hc$labels[hc$order]

## Determine cluster assignments.
scgp_clusters <- as.data.frame(scgp_cluster_results$cluster_assignments)
scgp_clusters$sample_barcode <- rownames(scgp_clusters)

## Plot a heatmap with annotations:
heatmap_pallete <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21)
pheatmap(pd_clust_pairwise_mat,
         cluster_rows = scgp_cluster_results$hclust_obj,
         cluster_cols = scgp_cluster_results$hclust_obj,
         treeheight_row = 0,
         border_color = NA,
         color = heatmap_pallete,
         show_colnames = F,
         show_rownames = F, 
         annotation_col = scgp_cluster_results$cluster_assignments)

## Visualize the clusters from a dimensionality reduction point of view:
viz_df <- visualize_clusters(pd_clust_pairwise_mat, cluster_labels = scgp_cluster_results$cluster_assignments)
viz_df_revised <- as.data.frame(viz_df)

## This will provide assessment regarding other quality metrics (including: bisulfite conversion, read count, epimutation, etc).
scgp_cluster_results_annot <- viz_df_revised %>% 
  inner_join(rrbs_qc_pass, by = c("Row.names"="cell_barcode")) %>% 
  mutate(cell_num = ifelse(grepl("50P", Row.names), 50, 1))

## Plot clusters by a few different variables, including the Molecular NeuroPathology tool.
ggplot(scgp_cluster_results_annot, aes(V1, V2, color = mnp_classification_epic, shape = factor(cell_num), size = factor(cell_num))) +
  geom_point(alpha=0.8) + 
  labs(color = "MNP classification", shape = "Cell num", x="MDS Dimension 1", y="MDS Dimension 2") +
  guides(size=FALSE) +
  plot_theme +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

pdf(file = "/Users/johnsk/github/results/Fig1/Fig1b-PDclust.pdf", width = 7, height = 5)
ggplot(scgp_cluster_results_annot, aes(V1, V2, color = case_barcode, shape = factor(cell_num), size = factor(cell_num))) +
  geom_point(alpha = 0.8) + 
  scale_color_manual(values=c("SM001" = "#F8766D", "SM002" = "#DB8E00", "SM004" = "#AEA200", "SM006" = "#64B200",
                              "SM008" = "#00BD5C" , "SM011" = "#00C1A7", "SM012" =  "#00BADE", "SM015" = "#00A6FF",
                              "SM017" = "#B385FF", "SM018" = "#EF67EB", "UC917" = "#FF63B6")) +
  labs(x= "MDS Dimension 1", y = "MDS Dimension 2", color="Subject", shape = "Tumor cell\nnumber") +
  guides(size=FALSE) +
  plot_theme +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "/Users/johnsk/github/results/Fig1/bisulfite-conversion-rates.pdf", width = 9, height = 5)
ggplot(scgp_cluster_results_annot, aes(V1, V2, color = bisulfite_conversion_rate, shape = factor(cell_num), size = factor(cell_num))) +
  geom_point() + 
  xlab("MDS Dimension 1") + 
  ylab("MDS Dimension 2") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title =  element_text(size = 14)) +
  labs(color="Bisulfite Conv.\nEfficiency") + 
  labs(shape = "Cell number") +
  guides(size=FALSE) 
dev.off()


### END ###