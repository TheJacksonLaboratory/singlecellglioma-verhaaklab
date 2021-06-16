##################################
# Examine the DNAme disorder summarized at the gene-level
# Updated: 2021.05.13
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(topGO)
# Method to determine "high" levels based on ranked data.
source("singlecellglioma-verhaaklab/analysis//ROSE-high-delineation.R")
###################################

## Generate plot theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 10),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())

## Determine promoter-level PDR levels based on categories (i.e., high vs. low):
prom_epimut <- read.table("data/Samples-passQC_single_cells_individual_promoter-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)
## Revise the sample names to match between the promoter methylation data.
colnames(prom_epimut) <- gsub("\\.", "-",  colnames(prom_epimut))


## Additional information about single-cells passing QC.
rrbs_qc <- read.table(file="data/analysis_scRRBS_sequencing_qc.csv", sep = ",", header = TRUE)

### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and tumor status (inferred CNVs).
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>% 
  left_join(meta_data, by="case_barcode") %>% 
  mutate(IDH_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))

## Define subgroups based on IDHmut status.
rrbs_qc_pass_wt = rrbs_qc_pass %>% filter(IDH_status=="IDHwt")
rrbs_qc_pass_mut = rrbs_qc_pass %>% filter(IDH_status=="IDHmut")

## Restrict to only tumor cells.
prom_epimut_tumor = prom_epimut[, colnames(prom_epimut)%in%rrbs_qc_pass$cell_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(prom_epimut_tumor)/dim(prom_epimut_tumor)[2])
## Identify the number of gene promoters measured in at least 100 cells.
sum(missingness<(1-100/844))

## Restrict to genes covered in at least 100 tumor cells samples:
epimut_prom_cov <- as.data.frame(rowSums(!is.na(prom_epimut_tumor)))
epimut_prom_avg <- as.data.frame(rowMeans(prom_epimut_tumor, na.rm = TRUE))
epimut_prom_total = cbind(epimut_prom_avg, epimut_prom_cov)
colnames(epimut_prom_total) <- c("tumor_pdr", "tumor_cov")
epimut_prom_total$gene_id <- rownames(epimut_prom_total)
epimut_prom_total_filt = epimut_prom_total %>% 
  filter(tumor_cov>100)

## The distribution is heavily tilted toward low epimutation rate.
hist(epimut_prom_total_filt$tumor_pdr)
ggplot(epimut_prom_total_filt, aes(x=tumor_cov, y=tumor_pdr)) + geom_point(alpha=0.2) +
  xlim(0, 844) + stat_cor(method="spearman") +
  plot_theme

## Create sort-able object for plotting by increasing gene promoter epimutation.
sort_df <- epimut_prom_total_filt %>%
  arrange(desc(tumor_pdr)) 
gene_order <- rev(unique(sort_df$gene_id))
epimut_prom_total_filt <- epimut_prom_total_filt %>% mutate(gene_id = factor(gene_id, levels = gene_order))

## Quick visualization.
ggplot(epimut_prom_total_filt, aes(x=gene_id)) +
  geom_point(aes(y=tumor_pdr)) +
  plot_theme +
  null_x

## Determine the CUTOFF for high gene promoter epimutation.
epimut_vector = epimut_prom_total_filt$tumor_pdr
calculate_cutoff(epimut_vector) # 0.399 is the cut-off for "high" DNAme disorder (aka epimutation).

## Examine the "high" DNAme disorder gene set.
epimut_high <- epimut_prom_total_filt %>% 
  filter(tumor_pdr >= 0.399)

## What about promoter DNAme disorder rates less than 0.1, manually setting a lower bound.
epimut_low <- epimut_prom_total_filt %>% 
  filter(tumor_pdr < 0.1)

####################
### topGO analysis - all tumor cells
####################
GO_prep = epimut_prom_total_filt 
gene_list <- GO_prep$tumor_pdr
names(gene_list) <- GO_prep$gene_id

## Function for enrichment of "high" DNAme disorder genes against covered background.
selFun = function(x) {
  ifelse(x>=0.399, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'symbol', # We're using gene symbols.
              nodeSize = 10
)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")

## Extract GOs
allGO = usedGO(object = GOdata) 

# Make the table
resultTable_all_high <- GenTable(GOdata, 
                        classicFisher = resultFisher, 
                        classicKS = resultKS, 
                        topNodes = 30, 
                        ranksOf = 'classicFisher'
)

# Use extracted GO terms as follows:
resultTable_high_epimut <- GenTable(GOdata, 
                                 classicFisher = resultFisher, 
                                 classicKS = resultKS, 
                                 topNodes = length(allGO), 
                                 ranksOf = 'classicFisher'
)
hist(as.numeric(resultTable_high_epimut$classicFisher))
resultTable_high_epimut$adj_pvalue <- p.adjust(resultTable_high_epimut$classicFisher, method = "fdr", n = length(resultTable_high_epimut$classicFisher))
hist(resultTable_high_epimut$adj_pvalue)
sum(resultTable_high_epimut$adj_pvalue<0.05)

#############################
### Low DNAme disorder promoters - all
#############################
## Function for enrichment of "low" DNAme disorder genes against covered background.
selFun = function(x) {
  ifelse(x<0.1, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'symbol', # We're using gene symbols.
              nodeSize = 10
)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
## Extract GOs
allGO = usedGO(object = GOdata) 

# Make the table containing the results for low-DNAme disorder.
resultTable_all_low <- GenTable(GOdata, 
                             classicFisher = resultFisher, 
                             classicKS = resultKS, 
                             topNodes = 30, 
                             ranksOf = 'classicFisher'
)

# Use extracted GO terms as follows:
resultTable_low_epimut <- GenTable(GOdata, 
                                    classicFisher = resultFisher, 
                                    classicKS = resultKS, 
                                    topNodes = length(allGO), 
                                    ranksOf = 'classicFisher'
)

hist(as.numeric(resultTable_low_epimut$classicFisher))
resultTable_low_epimut$adj_pvalue <- p.adjust(resultTable_low_epimut$classicFisher, method = "fdr", n = length(resultTable_low_epimut$classicFisher))
hist(resultTable_low_epimut$adj_pvalue)
sum(resultTable_low_epimut$adj_pvalue<0.05)


##################################
### Generate enrichment visualizations for two sets
##################################
resultTable_high_epimut$GO.ID <- factor(resultTable_high_epimut$GO.ID)
resultTable_high_epimut$adj_pvalue <- as.numeric(resultTable_high_epimut$adj_pvalue)
resultTable_filtered <- resultTable_high_epimut[1:10, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

## Just based on Gene Ontology IDs.
pdf(file = "results/Fig2/Fig2d-high-DNAme disorder-promoter-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 4) +
  coord_flip()
dev.off()

pdf(file = "results/Fig2/Fig2d-high-DNAme disorder-promoter-GOterms.pdf", height = 6, width = 4, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 4) +
  coord_flip()
dev.off()


## Enrichment for low-DNAme disorder genes.
resultTable_low_epimut$GO.ID <- factor(resultTable_low_epimut$GO.ID)
resultTable_low_epimut$adj_pvalue <- as.numeric(resultTable_low_epimut$adj_pvalue)
resultTable_filtered <- resultTable_low_epimut[1:10, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

pdf(file = "results/Fig2/Fig2d-low-DNAme disorder-promoter-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 3) +
  coord_flip()
dev.off()

pdf(file = "results/Fig2/Fig2d-low-DNAme disorder-promoter-GOterms.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 3) +
  coord_flip()
dev.off()

### END ###