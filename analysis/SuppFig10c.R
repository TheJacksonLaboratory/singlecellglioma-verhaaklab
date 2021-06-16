##################################
# Clustering of CNV and DNA methylation from scRRBS.
# Updated: 2021.05.19
# Author: Kevin J.
###################################

## Refresher on hclust:
# https://uc-r.github.io/hc_clustering
## Hints to annotate dendrogram
# https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend.html

# Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)


###################################
library(tidyverse)  # data organization
library(RColorBrewer) # color options
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(openxlsx) # for opening xlsx files.
library(ggpubr)
#################################################
## Plotting theme:
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




## Load the SCGP subject-level metadata.
full_meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Need to extract subtype and IDHmut status.
meta_data = full_meta_data %>% 
  select(case_barcode = subject_id, subtype) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(case_barcode_short = gsub("-", "", substr(case_barcode, 6, 11)))

### Epimutation
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
## Create a new variable to indicate "non_tumor" or shortened case barcode.
gg_epiallele = epiallele_info %>% 
  mutate(case_normal_barcode = ifelse(tumor_cnv==0, "Non-tumor", gsub("-","", substr(sample, 6, 11)))) %>% 
  left_join(meta_data, by=c("sample"="case_barcode")) %>% 
  mutate(idh_codel_subtype = ifelse(case_normal_barcode=="Non-tumor", "Non-tumor", subtype)) %>% 
  select(sample_barcode, case_barcode = case_normal_barcode, PDR, idh_codel_subtype, library_id)

## Summarized DNA methylation data across the 914 cells passing QC.
tiles_10kb = read.table("data/Samples_annotation_hg19_tiling10kb_methylation-filtered_passedQC.txt", sep="\t", header=T, stringsAsFactors = F)
colnames(tiles_10kb) <- gsub("\\.", "-",  colnames(tiles_10kb))

## TILED DNA methylation.
## Create a mean methylation value across the 10-kb tiled region.
bin_counts_10kb <- colSums(!is.na(tiles_10kb[, 5:918]))
bin_meth_10kb <- colMeans(tiles_10kb[, 5:918], na.rm = TRUE)
broad_meth <- as.data.frame(t(bind_rows(bin_counts_10kb, bin_meth_10kb)))
colnames(broad_meth) <- c("bins_covered_10kb", "mean_methylation_10kb")
broad_meth$sample_barcode <- rownames(broad_meth)
broad_meth$case_barcode <- gsub("-", "", substr(broad_meth$sample_barcode, 6, 11))

## LIGER classes.
sm012_classes = read.table("data/SM012-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm012_classes = sm012_classes %>% 
  filter(dataset=="dna") %>% 
  mutate(binary_liger_state = ifelse(liger_state%in%c("stemcell", "prolif_stemcell"), "stemcell", "differentiated"))

## Combine all meta-data
cnv_plot_meta = broad_meth %>% 
  inner_join(sm012_classes, by=c("sample_barcode"="cell_name")) %>% 
  inner_join(gg_epiallele, by="sample_barcode")

## Final object to be used in downstream analyses for metadata.
rrbs_SM012_meta = cnv_plot_meta %>% 
  arrange(sample_barcode)

### 1Mb variable bins as processed by Ginkgo.
SM012 <- read.table("data/SM012-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
SM012 <- SM012[-which(SM012$CHR%in%c("chrX", "chrY")), ]
colnames(SM012) <- gsub("\\.", "-", colnames(SM012))

# Subset the SM012 CNV data to the cells with passing QC.
index <- colnames(SM012)%in%rrbs_SM012_meta$sample_barcode
SM012_filtered <- SM012[ ,index]
colnames(SM012_filtered) <- gsub("SCGP-SM-012-01D-S5M-", "SM012-", colnames(SM012_filtered))

## Create a new color for library and for cell state.
all(rownames(df)==rrbs_SM012_meta$match_id)
rrbs_SM012_meta = rrbs_SM012_meta %>% 
  mutate(lib_col = recode(library_id, "SCGP-SM-012-1-1" = "#2FB3CA",
                          "SCGP-SM-012-1-3" = "#CA2F66",
                          "SCGP-SM-012-1-4" = "#CA932F"))
rrbs_SM012_meta = rrbs_SM012_meta %>% 
  mutate(state_col = recode(liger_state_dna, "differentiated" = "#fcbba1",
                          "stemcell" = "#fb6a4a"))

rrbs_SM012_meta$meth_col = colorRampPalette(c('blue', 'yellow'))(length(rrbs_SM012_meta$mean_methylation_10kb))[rank(rrbs_SM012_meta$mean_methylation_10kb)]
rrbs_SM012_meta$epimut_col = colorRampPalette(c('blue', 'red'))(length(rrbs_SM012_meta$PDR))[rank(rrbs_SM012_meta$PDR)]

annot_bars <- cbind(rrbs_SM012_meta$epimut_col, rrbs_SM012_meta$meth_col, rrbs_SM012_meta$state_col)


####################################
#### CNV hierarchal clustering
####################################
df2 <- na.omit(SM012_filtered)
df2 <- scale(t(df2))

# Dissimilarity matrix. 
d2 <- dist(df2, method = "euclidean")

# Hierarchical clustering using Ward's.
hc_cnv <- hclust(d2, method = "ward.D2" )

# Plot the obtained dendrogram.
plot(hc_cnv, cex = 0.6, hang = -1)

# Color the annotation bar with library color or LIGER cell state.
dend_cnv = hc_cnv %>% 
  as.dendrogram()

## Multiple annotations:
dend_cnv %>% set("labels_cex", 0.01) %>%  plot(main="SM012 CNV clustering")
dend_cnv %>%  rect.dendrogram(k=5, border = 4, lty = 5, lwd = 2)    
colored_bars(colors = annot_bars, dend = dend_cnv, rowLabels = c("DNAme disorder", "DNA methylation", "Cell state"), y_shift = -25)

## Cut tree at k = 5 for this sample:
cutree(dend_cnv, k =5)
all(rrbs_SM012_meta$match_id==names(cutree(dend_cnv, k =5)))

## Define cluster based on cuts.
rrbs_SM012_meta$scna_k_clust <- cutree(dend_cnv, k =5)
rrbs_SM012_meta$scna_subclone_clust <- NA
rrbs_SM012_meta$scna_subclone_clust[rrbs_SM012_meta$scna_k_clust==1] <- "cluster1"
rrbs_SM012_meta$scna_subclone_clust[rrbs_SM012_meta$scna_k_clust==2] <- "cluster2"
rrbs_SM012_meta$scna_subclone_clust[rrbs_SM012_meta$scna_k_clust==3] <- "cluster3"
rrbs_SM012_meta$scna_subclone_clust[rrbs_SM012_meta$scna_k_clust==4] <- "cluster4"
rrbs_SM012_meta$scna_subclone_clust[rrbs_SM012_meta$scna_k_clust==5] <- "cluster1"

## Create a visualization for the differences in DNAme disorder.
pdf(file = "results/Fig6/SuppFig10c-SM012-DNAme-disorder-clusters.pdf", height = 2, width = 3, bg = "transparent", useDingbats = FALSE)
ggplot(rrbs_SM012_meta, aes(x=scna_subclone_clust, y=PDR)) +
  geom_boxplot() +
  labs(x = "SCNA cluster", y="DNAme disorder") +
  plot_theme +
  stat_compare_means(method="kruskal")
dev.off()
## Create a visualization for the differences in DNA methylation.
pdf(file = "results/Fig6/SuppFig10c-SM012-DNAme-clusters.pdf", height = 2, width = 3, bg = "transparent", useDingbats = FALSE)
ggplot(rrbs_SM012_meta, aes(x=scna_subclone_clust, y=mean_methylation_10kb)) +
  geom_boxplot() +
  labs(x = "SCNA cluster", y="DNAme (10-kb)") +
  plot_theme +
  stat_compare_means(method="kruskal")
dev.off()


### END ####


