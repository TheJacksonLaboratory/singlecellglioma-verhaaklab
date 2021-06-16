##############################################
# Results from linear mixed effects model.
# Updated: 2020.06.02
# Author: Kevin J.
##################################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
## Load the necessary packages.
library(tidyverse)
library(minfi)
library(openxlsx)
library(LOLA)
library(GenomicRanges)
library(devtools)
library(ggrepel)
library(topGO)
##########################################
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


#########################################
### Pre-processing and QC
#########################################
## These are publicly available datasets each with their own metadata.
## These include De Souza et al, TCGA, Mazor et al, and GLASS.
## The `all_data` object represents the list of IDAT files to be processed.

## Start off by creating a RG channel set.
all_450k = read.metharray(all_data, verbose=FALSE)

## Determine some QC measures.
il450k_detp = detectionP(all_450k)

## Perform functional normalization.
glass_mset <- preprocessFunnorm(all_450k, nPCs=2) 

## Ensure probes are in the same order in the `glass_mset` and `il450k_detp` objects.
il450k_detp_filt = il450k_detp[match(featureNames(glass_mset), rownames(il450k_detp)), ]

## Cross reactive probes identified by Chen et al. 
xReactiveProbes_450 = openxlsx::read.xlsx("/48639-non-specific-probes-Illumina450k.xlsx")

## Drop samples with high average detection P-values.
#il450k_g_filt = glass_mset[,which(colMeans(il450k_detp_filt) < 0.01)]
il450k_g_filt = glass_mset[rowSums(il450k_detp_filt < 0.01) == ncol(glass_mset) & !(featureNames(glass_mset) %in% xReactiveProbes_450$TargetID) & !(featureNames(glass_mset) %in% getAnnotation(glass_mset)$Name[getAnnotation(glass_mset)$chr %in% c("chrX", "chrY")]),]

## Add known SNP information to the MethylSet.
il450k_g_filt <- addSnpInfo(il450k_g_filt)
il450k_g_filt <- dropLociWithSnps(il450k_g_filt, snps=c("SBE","CpG"), maf=0)

############################
## Re-load this processed data. 
il450k_g_filt = readRDS(file = "data/glass-QC-filtered-normalized.rds")
annot_450K <- getAnnotation(il450k_g_filt)

## Load in the results from the linear mixed effects model below applied to each individual CpG treated as a M-value.
## fit = lme(CpG probe ~ timepoint + idh_codel_subtype + Cancer + Immune, random = ~1|case_barcode)
Results = readRDS(file = "data/differential-methylation.rds")

# Reduce the results to remove probes associated with: SNP, cross-hybridizing, low detection P-values, and sex chromosomes.
Results_filt = Results[rownames(Results)%in%rownames(annot_450K), ]
Results_filt$timepoint_pval_adj = p.adjust(Results_filt$timepoint_pval, method = "fdr", n = length(Results_filt$timepoint_pval))
hist(Results_filt$timepoint_pval_adj)

## We are interested in the CpGs that change most consistently across patients while controlling for differences in glioma subtype, purity, different immune cell proportions, and patient.
## Given the lme (adjusting for patient) and testing across thousands of probes, we want to balance correction for multiple tests with discovery.
hist(Results_filt$timepoint_pval_adj)
## Setting the FDR to 0.3 as an threshold that strikes a balance between enough regions for enrichment (3,419 CpGs) and trying to limit the false discoveries.
sig_cpgs = rownames(Results_filt)[which(Results_filt$timepoint_pval_adj < 0.3)]
progression_cpgs = annot_450K[rownames(annot_450K)%in%sig_cpgs, ]

## Investigations into gene-level epimutation estimates.
n_distinct(progression_cpgs$UCSC_RefGene_Name)

## Investigate the overlap with high DNAme disorder genes.
## Load the high DNAme disorder genes (calculated from gene body epimutation).
epimut_high = readRDS(file="data/high-epimutation-genebody-all-scRRBS.Rds")
epimut_low = readRDS(file="data/low-epimutation-genebody-all-scRRBS.Rds")
all_epimut_genes = readRDS(file="data/all-genebody-all-scRRBS.Rds")

## Retrieve the gene_id from the Illumina array, but a single CpG can track to multiple genes. Extract the first gene.
annot_450K$gene_id <- sapply(strsplit(annot_450K$UCSC_RefGene_Name, ";"), "[", 1)


annotated_df <- merge(Results_filt, annot_450K, by='row.names')
rownames(annotated_df) <- annotated_df$Row.names
annotated_df$Row.names <- NULL
annotated_df <- data.frame(annotated_df)

# Select the results that meet the threshold for significance
max(annotated_df$timepoint_pval[annotated_df$timepoint_pval_adj<0.3], na.rm = T)
annotated_df$sig <- ifelse(annotated_df$timepoint_pval_adj<0.3, "P-value < 0.003", "Not Sig")

# Genes with labeled regions
sum(annotated_df$timepoint_coeff>0 & annotated_df$timepoint_pval_adj<0.3, na.rm = T)
sum(annotated_df$timepoint_coeff<0 & annotated_df$timepoint_pval_adj<0.3, na.rm = T)

## Volcano plot visualization of the results as they pertain to CpG DNAme changes at recurrence.
ggplot(annotated_df, aes(timepoint_coeff, -log10(timepoint_pval))) +
  geom_point(aes(col=sig)) + #xlim(-.075, 0.075) + 
  xlab("Recurrence coefficient") + 
  ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red")) +
  plot_theme

###############################
### Perform enrichment analyses
###############################
# Create start and stop genomic location indicators
progression_cpgs$hg19_start = as.numeric(progression_cpgs$pos)
# Create another new variable with pos for CpG end
progression_cpgs$hg19_end = as.numeric(progression_cpgs$pos)
# Create a genomic ranges object.
progression_cpgs <- progression_cpgs[, c('chr', 'hg19_start', 'hg19_end', "Name", 'UCSC_RefGene_Name', 'Regulatory_Feature_Group')]
progression_cpgs_gr <- makeGRangesFromDataFrame(progression_cpgs, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# Create start and stop genomic location indicators
annot_450K$hg19_start = as.numeric(annot_450K$pos)
# Create another new variable with pos for CpG end
annot_450K$hg19_end = as.numeric(annot_450K$pos)
# Create a genomic ranges object.
annot_450K <- annot_450K[, c('chr', 'hg19_start', 'hg19_end', "Name", 'UCSC_RefGene_Name')]
annot_450K_gr <- makeGRangesFromDataFrame(annot_450K, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# Set Universe as CpGs considered in GLASS longitudinal analysis.
universe <- GRanges(annot_450K_gr)
geneset_progress <- GRanges(progression_cpgs_gr)

## What are the high epimutation binding sites for TFs?
tfbs_epimut = read.table("data/tfbs_epimutation_subtype.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
tfbs_epimut = tfbs_epimut %>% 
  dplyr::select(tf, IDHwt) %>% 
  distinct()

# Compute gene enrichment set (see document "Using LOLA Core")
# Load the LOLA core (cached version; hg19/38, ENCODE TFBS, UCSC CGIs, Citrome epigenome)
#install_github("sheffien/simpleCache")
lolaDB <- loadRegionDB("data/LOLAJaspar/hg19/")
locResults_progress <- runLOLA(geneset_progress, universe, lolaDB)

## convert -log10(p-value) to p-value and correct for FDR. 
locResults_progress$pvalue <- 10^-locResults_progress$pValueLog
## Perform multiple-hypotheses testing correction.
locResults_progress$pvalue_adj = p.adjust(locResults_progress$pvalue , method = "fdr", n = length(locResults_progress$pvalue))
locResults_progress_filt <- locResults_progress[which(locResults_progress$pvalue_adj<0.05) ,]

locResults_progress_filt$TF <- gsub(".bed", "", locResults_progress_filt$filename)
locResults_progress_filt$TF <- toupper(locResults_progress_filt$TF)
locResults_progress_filt$TF <- sapply(strsplit(locResults_progress_filt$TF, "__"), "[[", 1)
locResults_progress_filt$TF <- sapply(strsplit(locResults_progress_filt$TF, "_"), "[[", 1)
threshold <- -log10(0.05)

locResults_progress_filt_annot = locResults_progress_filt %>% 
  left_join(tfbs_epimut, by=c("TF"="tf")) %>% 
  filter(filename!="Arnt.bed") %>% 
  mutate(epimut = ifelse(IDHwt >= 0.39, "yes", ifelse(is.na(IDHwt), NA, "no")))


TFBS <- ggplot(locResults_progress_filt_annot, aes(x = reorder(TF, -pvalue_adj), y = -log10(pvalue_adj), label = TF)) +
  geom_point(aes(size = 2, color = IDHwt, alpha=.5)) + 
  geom_text_repel(size = 3) + ylim(0,12) +
  #scale_size(range = c(1,15))  +
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed", size = 1.1) +
  labs(color="scRRBS TFBS motif\nDNAme disorder")

pdf(file = "results/Fig7/Fig7f-tf-enrichment.pdf", height = 4, width = 6.5)
TFBS + 
  plot_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_alpha(guide = 'none') + 
  xlab("Transcription Factor (JASPAR)") + 
  ylab("-log10(adj. p-value)") + 
  guides(size=FALSE)
dev.off()

#########################################
### Gene Ontology enrichment analysis
### for significant TFs
#########################################
## Not reported in the paper, but helpful to understand the biological processes impacted by methylation changes at specific TFs. 
locResults_progress$TF <- gsub(".bed", "", locResults_progress$filename)
locResults_progress$TF <- toupper(locResults_progress$TF)
locResults_progress$TF <- sapply(strsplit(locResults_progress$TF, "__"), "[[", 1)
locResults_progress$TF <- sapply(strsplit(locResults_progress$TF, "_"), "[[", 1)

## Enrichment of TFs.
gene_list <- locResults_progress$pvalue_adj
names(gene_list) <- locResults_progress$TF
## Function for enrichment of significant TFs against covered background.
selFun = function(x) {
  ifelse(x<0.05, TRUE, FALSE)
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
resultTable_GLASS_LOLA <- GenTable(GOdata, 
                                   classicFisher = resultFisher, 
                                   classicKS = resultKS, 
                                   topNodes = 30, 
                                   ranksOf = 'classicFisher'
)

resultTable_GLASS_LOLA$GO.ID <- factor(resultTable_GLASS_LOLA$GO.ID)
resultTable_GLASS_LOLA$classicFisher <- as.numeric(resultTable_GLASS_LOLA$classicFisher)
resultTable_filtered <- resultTable_GLASS_LOLA[1:15, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[order(resultTable_filtered$classicFisher)])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[order(resultTable_filtered$classicFisher)])

## Just based on Gene Ontology IDs.
ggplot(resultTable_filtered, aes(x=Term, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#CA2F66") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 3) +
  coord_flip()

### END ###