##################################
# Examine the epimutation rates summarized at the gene-level
# Updated: 2021.02.07
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

null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


## Gene body epimutation rates.
gb_epimut <- read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Samples-passQC_single_cells_individual_gene_body-specific_PDR-filtered.txt", sep="\t", row.names=1, header=T, stringsAsFactors = F)

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
colnames(gb_epimut) <- gsub("\\.", "-",  colnames(gb_epimut))

## Restrict to only tumor cells.
gb_epimut_tumor = gb_epimut[, colnames(gb_epimut)%in%rrbs_qc_pass$sample_barcode]

## Determine the missingness by colSums.
missingness <- rowSums(is.na(gb_epimut_tumor)/dim(gb_epimut_tumor)[2])
## Identify the number of gene promoters measured in at least 100 cells.
sum(missingness<(100/844))

## Restrict to genes covered in at least 100 tumor cells samples:
epimut_gb_cov <- as.data.frame(rowSums(!is.na(gb_epimut_tumor)))
epimut_gb_avg <- as.data.frame(rowMeans(gb_epimut_tumor, na.rm = TRUE))
epimut_gb_total = cbind(epimut_gb_avg, epimut_gb_cov)
colnames(epimut_gb_total) <- c("tumor_pdr", "tumor_cov")
epimut_gb_total$gene_id <- rownames(epimut_gb_total)
epimut_gb_total_filt = epimut_gb_total %>% 
  filter(tumor_cov>100)

## The distribution is heavily tilted toward low epimutation rate.
hist(epimut_gb_total_filt$tumor_pdr)
ggplot(epimut_gb_total_filt, aes(x=tumor_cov, y=tumor_pdr)) + geom_point(alpha=0.2) +
  xlim(0, 844) + stat_cor(method="spearman")

## Create sortable object for plotting by increasing gene gboter epimutation.
sort_df <- epimut_gb_total_filt %>%
  arrange(desc(tumor_pdr)) 
gene_order <- rev(unique(sort_df$gene_id))
epimut_gb_total_filt <- epimut_gb_total_filt %>% mutate(gene_id = factor(gene_id, levels = gene_order))


## Determine the CUTOFF for high gene gboter epimutation.
epimut_vector = epimut_gb_total_filt$tumor_pdr
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/high-epimutation-vector.pdf", height = 5, width = 7)
calculate_cutoff(epimut_vector)
dev.off()

## Be able to intersect with the HUGO IDs:
featuredata = readRDS(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/ensembl-hugo.Rds")
featuredata$gene_id <- rownames(featuredata)

## Add the new variable for classification:
epimut_gb_total_filt = epimut_gb_total_filt %>% 
  mutate(tumor_class = ifelse(tumor_pdr > 0.492, "high", "low")) %>% 
  inner_join(featuredata, by="gene_id")

## Output table with gene body epimutation results.
write.table(epimut_gb_total_filt, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_genebody_total_filt-20210207.txt", sep="\t", row.names = F, col.names = T, quote = F)


## Quick visualization.
ggplot(epimut_gb_total_filt, aes(x=gene_id)) +
  geom_point(aes(y=tumor_pdr, col = tumor_class)) +
  scale_color_manual(values = c("high" = "#CD4F39",
                                "low" ="black")) +
  plot_theme +
  null_x

## Examine the "high" epimutation gene set.
epimut_high <- epimut_gb_total_filt %>% 
  filter(tumor_pdr >= 0.492) %>% 
  inner_join(featuredata, by="gene_id")
## Output gene-body high epimutation across all single cell RRBS samples with the gene being covered by at least 100 samples.
saveRDS(epimut_high, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/high-epimutation-genebody-all-scRRBS.Rds")

## What about gb epimutation rates less than 0.1.
epimut_low <- epimut_gb_total_filt %>% 
  filter(tumor_pdr < 0.1) %>% 
  inner_join(featuredata, by="gene_id")
saveRDS(epimut_low, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/low-epimutation-genebody-all-scRRBS.Rds")

epimut_nothigh <- epimut_gb_total_filt %>% 
  filter(tumor_pdr <= 0.492) %>% 
  inner_join(featuredata, by="gene_id")
saveRDS(epimut_nothigh, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/nothigh-epimutation-genebody-all-scRRBS.Rds")

all_epimut_genes <- epimut_gb_total_filt %>% 
  inner_join(featuredata, by="gene_id")
saveRDS(all_epimut_genes, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/all-genebody-all-scRRBS.Rds")

####################
### topGO analysis - all tumor cells
####################
GO_prep = epimut_gb_total_filt 
gene_list <- GO_prep$tumor_pdr
names(gene_list) <- GO_prep$gene_id

## Function for enrichment of "high" epimutation genes against covered background.
selFun = function(x) {
  ifelse(x>=0.492, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'ensembl', # We're using gene symbols.
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
### Low epimutation gene body - all
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
              ID = 'ensembl', # We're using gene symbols.
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
### Generate enrichment analyses for all genes
##################################
resultTable_high_epimut$GO.ID <- factor(resultTable_high_epimut$GO.ID)
resultTable_high_epimut$adj_pvalue <- as.numeric(resultTable_high_epimut$adj_pvalue)
resultTable_filtered <- resultTable_high_epimut[1:10, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

## Just based on Gene Ontology IDs.
pdf(file = "github/results/Fig2/SuppFig5d-high-DNAme disorder-genebody-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  coord_flip()
dev.off()

pdf(file = "github/results/Fig2/SuppFig5d-high-DNAme disorder-genebody-GOterms.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +

  coord_flip()
dev.off()


## Enrichment for low-DNAme disorder genes.
resultTable_low_epimut$GO.ID <- factor(resultTable_low_epimut$GO.ID)
resultTable_low_epimut$adj_pvalue <- as.numeric(resultTable_low_epimut$adj_pvalue)
resultTable_filtered <- resultTable_low_epimut[1:10, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

pdf(file = "github/results/Fig2/SuppFig5e-low-DNAme disorder-genebody-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  coord_flip()
dev.off()

pdf(file = "github/results/Fig2/SuppFig5e-low-DNAme disorder-genebody-GOterms.pdf", height = 5, width = 3.5, useDingbats = FALSE, bg="transparent")
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  coord_flip()
dev.off()

### END ###



###################################
## Restrict to only IDHwt tumors
###################################
gb_epimut_wt = gb_epimut[rownames(gb_epimut)%in%epimut_gb_total_filt$gene_id, colnames(gb_epimut)%in%rrbs_qc_pass_wt$sample_barcode]

## Calculate the same metrics for IDHwt tumors as all tumor samples.
epimut_IDHwt_cov <- as.data.frame(rowSums(!is.na(gb_epimut_wt)))
epimut_IDHwt_avg <- as.data.frame(rowMeans(gb_epimut_wt, na.rm = TRUE))
epimut_IDHwt_total = cbind(epimut_IDHwt_avg, epimut_IDHwt_cov)
colnames(epimut_IDHwt_total) <- c("IDHwt_pdr", "IDHwt_cells_cov")
epimut_IDHwt_total$gene_id <- rownames(epimut_IDHwt_total)
epimut_IDHwt_total = epimut_IDHwt_total %>% 
  inner_join(featuredata, by="gene_id")
## Most genes in IDHwt cells are covered at ~10%.
sum(epimut_IDHwt_total$IDHwt_cells_cov>34)

## Determine the CUTOFF for IDHmut high gene promoter epimutation.
IDHwt_pdr_vector = epimut_IDHwt_total$IDHwt_pdr
calculate_cutoff(IDHwt_pdr_vector) # It's close at 0.35 with a more promoters that have "very low" epimutation rates.

## Apply the same classification labels to "high" gene promoters as all.
epimut_IDHwt_total$IDHwt_class <- ifelse(epimut_IDHwt_total$IDHwt_pdr>=0.536, "high", "low")

## Output table with IDHwt gene body epimutation results.
write.table(epimut_IDHwt_total, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_genebody_idhwt_filt-20210207.txt", sep="\t", row.names = F, col.names = T, quote = F)


###################################
## Restrict to only IDHmut tumors
###################################
## Restrict to only tumor cells.
gb_epimut_mut = gb_epimut[rownames(gb_epimut)%in%epimut_gb_total_filt$gene_id, colnames(gb_epimut)%in%rrbs_qc_pass_mut$sample_barcode]

epimut_IDHm_cov <- as.data.frame(rowSums(!is.na(gb_epimut_mut)))
epimut_IDHm_avg <- as.data.frame(rowMeans(gb_epimut_mut, na.rm = TRUE))
epimut_IDHm_total = cbind(epimut_IDHm_avg, epimut_IDHm_cov)
colnames(epimut_IDHm_total) <- c("IDHmut_pdr", "IDHmut_cells_cov")
epimut_IDHm_total$gene_id <- rownames(epimut_IDHm_total)
epimut_IDHm_total = epimut_IDHm_total %>% 
  inner_join(featuredata, by="gene_id")
## Determine the CUTOFF for IDHmut high gene promoter epimutation.
IDHmut_pdr_vector = epimut_IDHm_total$IDHmut_pdr
calculate_cutoff(IDHmut_pdr_vector) # It's very close.

epimut_IDHm_total$IDHmut_class <- ifelse(epimut_IDHm_total$IDHmut_pdr>=0.467, "high", "low")
## Output table with IDHwt gene body epimutation results.
write.table(epimut_IDHm_total, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimut_genebody_idhmut_filt-20210207.txt", sep="\t", row.names = F, col.names = T, quote = F)


#################################
### Concordance between All, IDHwt-only, and IDHmut-only.
#################################
all_analyses = epimut_gb_total_filt %>% 
  inner_join(epimut_IDHm_total, by ="gene_id") %>% 
  inner_join(epimut_IDHwt_total, by = "gene_id")

## What's the overlap between the "high" promoter epimutation classification.
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
## Function for enrichment of "high" epimutation genes against covered background.
selFun = function(x) {
  ifelse(x>=0.536, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'ensembl', # We're using gene symbols.
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
## Function for enrichment of "low" epimutation genes against covered background.
selFun = function(x) {
  ifelse(x<0.1, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'ensembl', # We're using gene symbols.
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
## Function for enrichment of "high" epimutation genes against covered background.
selFun = function(x) {
  ifelse(x>= 0.467, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'ensembl', # We're using gene symbols.
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
## Function for enrichment of "high" epimutation genes against covered background.
selFun = function(x) {
  ifelse(x<0.1, TRUE, FALSE)
}

GOdata <- new('topGOdata',
              ontology='BP', # We'll use Biological Process.
              allGenes=gene_list,
              geneSel = selFun,
              annot=annFUN.org,
              mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
              ID = 'ensembl', # We're using gene symbols.
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
### Generate enrichment analyses plots across ALL cells
##################################
resultTable_all_high$GO.ID <- factor(resultTable_all_high$GO.ID)
resultTable_all_high$classicFisher <- as.numeric(resultTable_all_high$classicFisher)
resultTable_filtered <- resultTable_all_high[1:15, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[order(resultTable_filtered$classicFisher)])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[order(resultTable_filtered$classicFisher)])

## Just based on Gene Ontology IDs.
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 12) +
  coord_flip()


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/high-epimutation-genebody-GO.pdf", height = 5, width = 7)
ggplot(resultTable_filtered, aes(x=Term, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#CD4F39") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 12) +
  coord_flip()
dev.off()


## Enrichment for low-epimutation genes.
resultTable_all_low$GO.ID <- factor(resultTable_all_low$GO.ID)
resultTable_all_low$classicFisher <- as.numeric(resultTable_all_low$classicFisher)
resultTable_filtered <- resultTable_all_low[1:15, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[order(resultTable_filtered$classicFisher)])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[order(resultTable_filtered$classicFisher)])

ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 12) +
  coord_flip()


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/low-epimutation-genebody-GO.pdf", height = 5, width = 7)
ggplot(resultTable_filtered, aes(x=Term, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#377eb8") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 12) +
  coord_flip()
dev.off()



### END ###

