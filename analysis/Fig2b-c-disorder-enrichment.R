##################################
# Examine the DNAme disorder summarized at the gene-level
# Updated: 2021.05.13
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(topGO)
# Method to determine "high" levels based on ranked data.
source("/Users/johnsk/Documents/Single-Cell-DNAmethylation/scripts/misc/ROSE-high-delineation.R")
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


## Promoter DNAme disorder rates.
prom_epimut <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-passQC_single_cells_individual_promoter-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)

## Load in the quality control data for these samples.
rrbs_qc <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/scgp_cnv_status.txt", sep="\t", header=T, stringsAsFactors = F)

## Load the subject-level metadata.
meta_data = readWorkbook("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Restrict to the samples that pass the general QC guidelines of CpG count, bisulfite conversion, and cell number.
rrbs_qc_pass <- rrbs_qc %>% 
  filter(cell_num == 1, cpg_unique > 40000, conversion_rate > 95, tumor_cnv == 1) %>% 
  left_join(meta_data, by=c("sample"="subject_id")) %>% 
  mutate(IDH_status = ifelse(subtype=="IDHwt", "IDHwt", "IDHmut"))

## Define subgroups based on IDHmut status.
rrbs_qc_pass_wt = rrbs_qc_pass %>% filter(IDH_status=="IDHwt")
rrbs_qc_pass_mut = rrbs_qc_pass %>% filter(IDH_status=="IDHmut")

## Revise the sample names to match between the promoter methylation data.
colnames(prom_epimut) <- gsub("\\.", "-",  colnames(prom_epimut))

## Restrict to only tumor cells.
prom_epimut_tumor = prom_epimut[, colnames(prom_epimut)%in%rrbs_qc_pass$sample_barcode]

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

## The distribution is heavily tilted toward low DNAme disorder rate.
hist(epimut_prom_total_filt$tumor_pdr)
ggplot(epimut_prom_total_filt, aes(x=tumor_cov, y=tumor_pdr)) + geom_point(alpha=0.2) +
  xlim(0, 844) + stat_cor(method="spearman")

## Create sortable object for plotting by increasing gene promoter DNAme disorder.
sort_df <- epimut_prom_total_filt %>%
  arrange(desc(tumor_pdr)) 
gene_order <- rev(unique(sort_df$gene_id))
epimut_prom_total_filt <- epimut_prom_total_filt %>% mutate(gene_id = factor(gene_id, levels = gene_order))

## Quick visualization.
ggplot(epimut_prom_total_filt, aes(x=gene_id)) +
  geom_point(aes(y=tumor_pdr)) +
  plot_theme +
  null_x

## Determine the CUTOFF for high gene promoter DNAme disorder.
epimut_vector = epimut_prom_total_filt$tumor_pdr
calculate_cutoff(epimut_vector)

## Add the new variable for classification:
epimut_prom_total_filt = epimut_prom_total_filt %>% 
  mutate(tumor_class = ifelse(tumor_pdr > 0.399, "high", "low"))

## Examine the "high" DNAme disorder gene set.
epimut_high <- epimut_prom_total_filt %>% 
  filter(tumor_pdr >= 0.399)
## When splitting on the median.
#epimut_high <- epimut_prom_total_filt %>% 
#  filter(tumor_pdr >= 0.2)
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


###################################
## Restrict to only IDHwt tumors
###################################
prom_epimut_wt = prom_epimut[rownames(prom_epimut)%in%epimut_prom_total_filt$gene_id, colnames(prom_epimut)%in%rrbs_qc_pass_wt$sample_barcode]

## Calculate the same metrics for IDHwt tumors as all tumor samples.
epimut_IDHwt_cov <- as.data.frame(rowSums(!is.na(prom_epimut_wt)))
epimut_IDHwt_avg <- as.data.frame(rowMeans(prom_epimut_wt, na.rm = TRUE))
epimut_IDHwt_total = cbind(epimut_IDHwt_avg, epimut_IDHwt_cov)
colnames(epimut_IDHwt_total) <- c("IDHwt_pdr", "IDHwt_cells_cov")
epimut_IDHwt_total$gene_id <- rownames(epimut_IDHwt_total)
## Most genes in IDHwt cells are covered at ~10%.
sum(epimut_IDHwt_total$IDHwt_cells_cov>34)
## Apply the same classification labels to "high" gene promoters as all.
epimut_IDHwt_total$IDHwt_class <- ifelse(epimut_IDHwt_total$IDHwt_pdr>=0.399, "high", "low")

## Determine the CUTOFF for IDHmut high gene promoter DNAme disorder.
IDHwt_pdr_vector = epimut_IDHwt_total$IDHwt_pdr
calculate_cutoff(IDHwt_pdr_vector) # It's close at 0.35 with more promoters that have "very low" DNAme disorder rates.

## Output table with IDHwt promoter DNAme disorder results.
write.table(epimut_IDHwt_total, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/DNAme disorder/epimut_promoter_idhwt_filt-20210207.txt", sep="\t", row.names = F, col.names = T, quote = F)


###################################
## Restrict to only IDHmut tumors
###################################
## Restrict to only tumor cells.
prom_epimut_mut = prom_epimut[rownames(prom_epimut)%in%epimut_prom_total_filt$gene_id, colnames(prom_epimut)%in%rrbs_qc_pass_mut$sample_barcode]

epimut_IDHm_cov <- as.data.frame(rowSums(!is.na(prom_epimut_mut)))
epimut_IDHm_avg <- as.data.frame(rowMeans(prom_epimut_mut, na.rm = TRUE))
epimut_IDHm_total = cbind(epimut_IDHm_avg, epimut_IDHm_cov)
colnames(epimut_IDHm_total) <- c("IDHmut_pdr", "IDHmut_cells_cov")
epimut_IDHm_total$gene_id <- rownames(epimut_IDHm_total)
epimut_IDHm_total$IDHmut_class <- ifelse(epimut_IDHm_total$IDHmut_pdr>=0.399, "high", "low")

## Determine the CUTOFF for IDHmut high gene promoter DNAme disorder.
IDHmut_pdr_vector = epimut_IDHm_total$IDHmut_pdr
calculate_cutoff(IDHmut_pdr_vector) # It's very close at 0.379.

## Output table with IDHmut promoter DNAme disorder results.
write.table(epimut_IDHm_total, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/DNAme disorder/epimut_promoter_idhmut_filt-20210207.txt", sep="\t", row.names = F, col.names = T, quote = F)


#################################
### Concordance between All, IDHwt-only, and IDHmut-only.
#################################
all_analyses = epimut_prom_total_filt %>% 
  inner_join(epimut_IDHm_total, by ="gene_id") %>% 
  inner_join(epimut_IDHwt_total, by = "gene_id")

## What's the overlap between the "high" promoter DNAme disorder classification.
table(all_analyses$tumor_class, all_analyses$IDHmut_class)
table(all_analyses$tumor_class, all_analyses$IDHwt_class)
table(all_analyses$IDHwt_class, all_analyses$IDHmut_class)



####################
### topGO analysis - HIGH - all IDHwt tumor cells
####################
GO_prep_IDHwt = epimut_IDHwt_total %>% 
  filter(!is.na(IDHwt_pdr))
gene_list <- GO_prep_IDHwt$IDHwt_pdr
names(gene_list) <- GO_prep_IDHwt$gene_id
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

# Make the table
resultTable_IDHwt_high <- GenTable(GOdata, 
                                 classicFisher = resultFisher, 
                                 classicKS = resultKS, 
                                 topNodes = 30, 
                                 ranksOf = 'classicFisher'
)

####################
### topGO analysis - LOW - all IDHwt tumor cells
####################
GO_prep_IDHwt = epimut_IDHwt_total %>% 
  filter(!is.na(IDHwt_pdr))
gene_list <- GO_prep_IDHwt$IDHwt_pdr
names(gene_list) <- GO_prep_IDHwt$gene_id
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

# Make the table
resultTable_IDHwt_low <- GenTable(GOdata, 
                                   classicFisher = resultFisher, 
                                   classicKS = resultKS, 
                                   topNodes = 30, 
                                   ranksOf = 'classicFisher'
)



####################
### topGO analysis - HIGH - all IDHmut tumor cells
####################
GO_prep_IDHmut = epimut_IDHm_total 
gene_list <- GO_prep_IDHmut$IDHmut_pdr
names(gene_list) <- GO_prep_IDHmut$gene_id
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

# Make the table for the results. 
resultTable_IDHmut_high <- GenTable(GOdata, 
                                   classicFisher = resultFisher, 
                                   classicKS = resultKS, 
                                   topNodes = 30, 
                                   ranksOf = 'classicFisher'
)

####################
### topGO analysis - LOW - all IDHmut tumor cells
####################
GO_prep_IDHmut = epimut_IDHm_total 
gene_list <- GO_prep_IDHmut$IDHmut_pdr
names(gene_list) <- GO_prep_IDHmut$gene_id
## Function for enrichment of "high" DNAme disorder genes against covered background.
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

# Make the table for the results. 
resultTable_IDHmut_low <- GenTable(GOdata, 
                                    classicFisher = resultFisher, 
                                    classicKS = resultKS, 
                                    topNodes = 30, 
                                    ranksOf = 'classicFisher'
)

##################################
### Generate enrichment analyses for all genes
##################################

resultTable_high_epimut$GO.ID <- factor(resultTable_high_epimut$GO.ID)
resultTable_high_epimut$adj_pvalue <- as.numeric(resultTable_high_epimut$adj_pvalue)
resultTable_filtered <- resultTable_high_epimut[1:10, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

## Just based on Gene Ontology IDs.
pdf(file = "github/results/Fig2/Fig2d-high-DNAme disorder-promoter-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 4) +
  coord_flip()
dev.off()

pdf(file = "github/results/Fig2/Fig2d-high-DNAme disorder-promoter-GOterms.pdf", height = 6, width = 4, useDingbats = FALSE, bg="transparent")
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

pdf(file = "github/results/Fig2/Fig2d-low-DNAme disorder-promoter-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 3) +
  coord_flip()
dev.off()

pdf(file = "github/results/Fig2/Fig2d-low-DNAme disorder-promoter-GOterms.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  ylim(0, 3) +
  coord_flip()
dev.off()

### END ###

