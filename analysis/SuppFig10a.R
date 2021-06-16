##################################
# Clustering of CNV and DNA methylation from scRRBS.
# Updated: 2020.05.12
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
sm001_classes = read.table("data/SM001-liger-classification.txt", sep="\t", header=T, stringsAsFactors = F)
sm001_classes = sm001_classes %>% 
  filter(dataset=="dna") %>% 
  mutate(binary_liger_state = ifelse(liger_state%in%c("stemcell", "prolif_stemcell"), "stemcell", "differentiated"))

## Combine all meta-data
cnv_plot_meta = broad_meth %>% 
  inner_join(sm001_classes, by=c("sample_barcode"="cell_name")) %>% 
  inner_join(gg_epiallele, by="sample_barcode")

## Final object to be used in downstream analyses for metadata.
rrbs_sm001_meta = cnv_plot_meta %>% 
  arrange(sample_barcode)

### 1Mb variable bins as processed by Ginkgo.
sm001 <- read.table("data/SM001-SegCopy-1Mb.txt", sep="\t", header=T, stringsAsFactors = F)
sm001 <- sm001[-which(sm001$CHR%in%c("chrX", "chrY")), ]
colnames(sm001) <- gsub("\\.", "-", colnames(sm001))

# Subset the sm001 CNV data to the cells with passing QC.
index <- colnames(sm001)%in%rrbs_sm001_meta$sample_barcode
sm001_filtered <- sm001[ ,index]
colnames(sm001_filtered) <- gsub("SCGP-SM-001-01D-S5M-", "SM001-", colnames(sm001_filtered))


## Create a new color for library and for cell state.
all(rownames(df)==rrbs_sm001_meta$match_id)
rrbs_sm001_meta = rrbs_sm001_meta %>% 
  mutate(lib_col = recode(library_id, "SCGP-SM-001-1-3" = "#ff7f00",
                          "SCGP-SM-001-2-1" = "#ffff33",
                          "SCGP-SM-001-5-1" = "#b2e2e2",
                          "SCGP-SM-001-5-3" = "#2ca25f",
                          "SCGP-SM-001-5-4" = "#006d2c",
                          "SCGP-SM-001-6-1" = "#fdcc8a",
                          "SCGP-SM-001-6-3" = "#e34a33",
                          "SCGP-SM-001-6-4" = "#b30000"))
rrbs_sm001_meta = rrbs_sm001_meta %>% 
  mutate(state_col = recode(liger_state_dna, "differentiated" = "#fcbba1",
                          "stemcell" = "#fb6a4a",
                          "prolif_stemcell" = "#fb6a4a"))

rrbs_sm001_meta$meth_col = colorRampPalette(c('blue', 'yellow'))(length(rrbs_sm001_meta$mean_methylation_10kb))[rank(rrbs_sm001_meta$mean_methylation_10kb)]
rrbs_sm001_meta$epimut_col = colorRampPalette(c('blue', 'red'))(length(rrbs_sm001_meta$PDR))[rank(rrbs_sm001_meta$PDR)]

annot_bars <- cbind(rrbs_sm001_meta$epimut_col, rrbs_sm001_meta$meth_col, rrbs_sm001_meta$state_col)


####################################
#### CNV hierarchal clustering
####################################
df2 <- na.omit(sm001_filtered)
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

## Library ID:
plot(dend_cnv)
colored_bars(colors = rrbs_sm001_meta$lib_col, dend = dend_cnv, rowLabels = "lib")

## Multiple annotations:
pdf(file = "results/Fig6/SM001-CNV-clus-annotation.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
dend_cnv %>% set("labels_cex", 0.01) %>% 
  plot(main="SM001 CNV clustering") 
dend_cnv %>%  rect.dendrogram(k=7, border = 4, lty = 5, lwd = 2)    
colored_bars(colors = annot_bars, dend = dend_cnv, rowLabels = c( "DNAme disorder", "DNA methylation", "Cell state"), y_shift = -25)
dev.off()

## Cell state ID -
plot(dend_cnv)
colored_bars(colors = rrbs_sm001_meta$state_col, dend = dend_cnv, rowLabels = "lib")

## Cut tree at k = 7 for this sample:
cutree(dend_cnv, k =7)
all(rrbs_sm001_meta$match_id==names(cutree(dend_cnv, k =7)))

## Define cluster based on cuts.
rrbs_sm001_meta$scna_k_clust <- cutree(dend_cnv, k =7)
rrbs_sm001_meta$scna_subclone_clust <- ifelse(rrbs_sm001_meta$scna_k_clust==1, "cluster1", "cluster2")

## Create a visualization for the differences in DNAme disorder.
pdf(file = "results/Fig6/SuppFig10a-SM001-DNAme-disorder-clusters.pdf", height = 2, width = 3, bg = "transparent", useDingbats = FALSE)
ggplot(rrbs_sm001_meta, aes(x=scna_subclone_clust, y=PDR)) +
  geom_boxplot() +
  labs(x = "SCNA cluster", y="DNAme disorder") +
  plot_theme +
  stat_compare_means(method="wilcox")
dev.off()

## Create a visualization for the differences in DNA methylation.
pdf(file = "results/Fig6/SuppFig10a-SM001-DNAme-clusters.pdf", height = 2, width = 3, bg = "transparent", useDingbats = FALSE)
ggplot(rrbs_sm001_meta, aes(x=scna_subclone_clust, y=mean_methylation_10kb)) +
  geom_boxplot() +
  labs(x = "SCNA cluster", y="DNAme (10-kb)") +
  plot_theme +
  stat_compare_means(method="wilcox")
dev.off()

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=0.25, height=1)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
color.bar(colorRampPalette(c("blue", "yellow"))(100), -1)
color.bar(colorRampPalette(c("blue", "red"))(100), -1)

### END ####