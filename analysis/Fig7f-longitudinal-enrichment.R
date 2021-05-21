##############################################
# Results from linear mixed effects model.
# Updated: 2020.06.02
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/"
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

## Skip to loading the data:

##########################################
# Step 1: Identify files to process
##########################################
## de Souza data. Make sure that matching TCGA pairs are retrieved (n = 3).
desouza_metadata = readWorkbook("data/desouza/Metadata_Mendeley.xlsx", sheet = 1, startRow = 4, colNames = TRUE)
desouza_metadata = desouza_metadata %>% 
  mutate(Slide = substr(desouza_metadata$Files, 1, 12),
         Array = substr(desouza_metadata$Files, 14, 19),
         idat = paste(Slide, Array, sep = "_"),
         file_path = paste("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/desouza/", idat, sep = ""))

# TCGA data available (450K).
tcga_metadata = read.csv("data/tcga/tcga_idats_sample_sheet.csv", sep = ",", stringsAsFactors = F, header = T)
# Location of IDAT files for study (Place the IDAT files in a designated file location) for minfi

# Isolate primary tumors (only three categories. Other is "Solid Tissue Normal").
p1 = tcga_metadata %>% 
  dplyr::select(Sample_Name, Status) %>%
  filter(Status %in% "Primary solid Tumor") %>% 
  mutate(pair_id = substring(Sample_Name, 1, 12)) 

# Isolate recurrent tumors (only three categories. Other is "Solid Tissue Normal").
p2 = tcga_metadata %>% 
  dplyr::select(Sample_Name, Status) %>%
  filter(Status %in% "Recurrent Solid Tumor") %>%
  mutate(pair_id = substring(Sample_Name, 1, 12)) 

# NOTE: There are three additional TCGA samples available in Houtan's Cell Reports paper.
# 26 unique pairs contained within TCGA.
pairs = inner_join(p1,p2, by="pair_id")
pair_names <- c(pairs$Sample_Name.x, pairs$Sample_Name.y)

## Matched 450K pairs (n = 26 subjects, 56 samples).
tcga_matched_450k = tcga_metadata[tcga_metadata$Sample_Name%in%pair_names, ]
tcga_matched_450k$idat = paste(tcga_matched_450k$Slide, tcga_matched_450k$Array, sep="_")

# Find the primary TCGA tumors (450K) that match the de Souza recurrences.
tcga_matched_primary = desouza_metadata %>% 
  inner_join(p1, by = c("Patient.ID" = "pair_id")) %>% 
  dplyr::select(pair_id = Patient.ID, Sample_Name, Status)

## Matched 450K pairs (n = 3 subjects, 3 samples).
desouza_matched_450K = tcga_metadata[tcga_metadata$Sample_Name%in%tcga_matched_primary$Sample_Name, ]
desouza_matched_450K$idat = paste(desouza_matched_450K$Slide, desouza_matched_450K$Array, sep = "_")

## The mazor dataset profiled both multisector and longitudinal glioma samples.
mazor_metadata = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/mazor-cancer-cell/mazor_idat/mazor-metadata.csv", sep = ",", stringsAsFactors = F, header = T)
mazor_metadata = mazor_metadata %>% 
  mutate(idat = paste(Slide, Array, sep="_"),
         file_path = paste("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/mazor-cancer-cell/mazor_idat/", idat, sep=""))

## The dkfz dataset that accompanied the WXS data.
dkfz_metadata = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/dkfz/methylation_array_samplesheet.csv", sep = ",", stringsAsFactors = F, header = T)
dkfz_metadata = dkfz_metadata %>% 
  mutate(idat = paste(sentrix_id, sentrix_position, sep="_"),
         file_path = paste("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/dkfz/", idat, sep=""))

## The Leeds dataset was provided by Lucy Stead.
leeds_metadata = read.csv("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/leeds/methylation-leeds-20190401.csv", sep = ",", stringsAsFactors = F, header = T)
leeds_metadata = leeds_metadata %>% 
  mutate(idat = Idat.filename.stem,
         sample_id = substr(Glass.barcode, 1, 15),
         file_path = paste("/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/leeds/", idat, sep=""))


# Combine all files to apply the classifier.
all_data = c(desouza_metadata$file_path, desouza_matched_450K$Basename, tcga_matched_450k$Basename, mazor_metadata$file_path, 
             dkfz_metadata$file_path, leeds_metadata$file_path)

#########################################
### Step 2: QC and pre-process
#########################################
## Start off by creating a RG channel set.
all_450k = read.metharray(all_data, verbose=FALSE)

## Determine some QC measures.
il450k_detp = detectionP(all_450k)

## Perform functional normalization.
glass_mset <- preprocessFunnorm(all_450k, nPCs=2) 

## Ensure probes are in the same order in the `glass_mset` and `il450k_detp` objects.
il450k_detp_filt = il450k_detp[match(featureNames(glass_mset), rownames(il450k_detp)), ]

## Filter probes.
xReactiveProbes_450 = openxlsx::read.xlsx("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/ref/48639-non-specific-probes-Illumina450k.xlsx")

## Drop samples with high average detection P-values.
#il450k_g_filt = glass_mset[,which(colMeans(il450k_detp_filt) < 0.01)]
il450k_g_filt = glass_mset[rowSums(il450k_detp_filt < 0.01) == ncol(glass_mset) & !(featureNames(glass_mset) %in% xReactiveProbes_450$TargetID) & !(featureNames(glass_mset) %in% getAnnotation(glass_mset)$Name[getAnnotation(glass_mset)$chr %in% c("chrX", "chrY")]),]


## Add known SNP information to the MethylSet.
il450k_g_filt <- addSnpInfo(il450k_g_filt)
il450k_g_filt <- dropLociWithSnps(il450k_g_filt, snps=c("SBE","CpG"), maf=0)
## Save output for the final R data file.
saveRDS(il450k_g_filt, file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/glass-QC-filtered-normalized.rds")

annot_450K <- getAnnotation(il450k_g_filt)
## What's the final set of CpGs to be included in an analysis?
dim(annot_450K)


######### *****************************************************************
## Linear mixed effect model results including all CpGs.
il450k_g_filt = readRDS(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/public/glass-QC-filtered-normalized.rds")
annot_450K <- getAnnotation(il450k_g_filt)

## fit = lme(y ~ timepoint + idh_codel_subtype + Cancer + Immune, random = ~1|case_barcode)
Results = readRDS(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/public/differential-methylation.rds")

# Reduce the results to remove probes associated with: SNP, cross-hybridizing, low detection P-values, and sex chromosomes.
Results_filt = Results[rownames(Results)%in%rownames(annot_450K), ]
Results_filt$timepoint_pval_adj = p.adjust(Results_filt$timepoint_pval, method = "fdr", n = length(Results_filt$timepoint_pval))
hist(Results_filt$timepoint_pval)
sig_cpgs = rownames(Results_filt)[which(Results_filt$timepoint_pval_adj < 0.3)]
progression_cpgs = annot_450K[rownames(annot_450K)%in%sig_cpgs, ]

## Investigations into gene-level epimutation estimates.
n_distinct(progression_cpgs$UCSC_RefGene_Name)

## Investigate the overlap with high epimutation genes.
## Load the high epimutation genes (calculated from gene body epimutation).
epimut_high = readRDS(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/high-epimutation-genebody-all-scRRBS.Rds")
epimut_low = readRDS(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/low-epimutation-genebody-all-scRRBS.Rds")
all_epimut_genes = readRDS(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/all-genebody-all-scRRBS.Rds")

annot_450K$gene_id <- sapply(strsplit(annot_450K$UCSC_RefGene_Name, ";"), "[[", 1)

sum(unique(annot_450K$UCSC_RefGene_Name)%in%epimut_high$Associated.Gene.Name)
sum(unique(annot_450K$UCSC_RefGene_Name)%in%epimut_low$Associated.Gene.Name)
sum(unique(progression_cpgs$UCSC_RefGene_Name)%in%epimut_high$Associated.Gene.Name)
sum(unique(progression_cpgs$UCSC_RefGene_Name)%in%epimut_low$Associated.Gene.Name)
unique(progression_cpgs$UCSC_RefGene_Name)[unique(progression_cpgs$UCSC_RefGene_Name)%in%epimut_high$Associated.Gene.Name]
unique(progression_cpgs$UCSC_RefGene_Name)[unique(progression_cpgs$UCSC_RefGene_Name)%in%epimut_low$Associated.Gene.Name]


annotated_df <- merge(Results_filt, annot_450K, by='row.names')
rownames(annotated_df) <- annotated_df$Row.names
annotated_df$Row.names <- NULL
annotated_df <- data.frame(annotated_df)

# Select the results that meet the threshold for significance
annotated_df$sig <- ifelse(annotated_df$timepoint_pval_adj<0.3, "P-value < 0.002", "Not Sig")
max(annotated_df$timepoint_pval[annotated_df$timepoint_pval_adj<0.3], na.rm = T)

# Figure 1. Genes with labeled regions
sum(annotated_df$timepoint_coeff>0 & annotated_df$timepoint_pval_adj<0.3, na.rm = T)
sum(annotated_df$timepoint_coeff<0 & annotated_df$timepoint_pval_adj<0.3, na.rm = T)

#p = ggplot(annotated_df, aes(timepoint_coeff, -log10(timepoint_pval))) +
#  geom_point(aes(col=sig)) + #xlim(-.075, 0.075) + 
#  xlab("Recurrence coefficient") + ylab("-log10(P-value)") +
#  scale_color_manual(values=c("black", "red"))
#pdf(file = "/Users/johnsk/Documents/glass-diff-meth-volcano.pdf", height = 5, width = 7)
#p + plot_theme
#dev.off()

######
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
tfbs_epimut = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/tfbs_epimutation_subtype.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
tfbs_epimut = tfbs_epimut %>% 
  dplyr::select(tf, IDHwt) %>% 
  distinct()

# Compute gene enrichment set (see document "Using LOLA Core")
# Load the LOLA core (cached version; hg19/38, ENCODE TFBS, UCSC CGIs, Citrome epigenome)
#install_github("sheffien/simpleCache")
lolaDB <- loadRegionDB("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/LOLAJaspar/hg19/")
locResults_progress <- runLOLA(geneset_progress, universe, lolaDB)
ocResults_progress_filt <- locResults_progress[which(locResults_progress$pValueLog>1.30103) ,]
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
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig7/Fig7f-tf-enrichment.pdf", height = 4, width = 6.5)
TFBS + 
  plot_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_alpha(guide = 'none') + 
  xlab("Transcription Factor (JASPAR)") + 
  ylab("-log10(adj. p-value)") + 
  guides(size=FALSE)
dev.off()

### Gene Ontology enrichment analysis:
locResults_progress$TF <- gsub(".bed", "", locResults_progress$filename)
locResults_progress$TF <- toupper(locResults_progress$TF)
locResults_progress$TF <- sapply(strsplit(locResults_progress$TF, "__"), "[[", 1)
locResults_progress$TF <- sapply(strsplit(locResults_progress$TF, "_"), "[[", 1)

## Enrichment of TFs.
gene_list <- locResults_progress$qValue
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

library(ALL)
data(ALL)
affyLib <- paste(annotation(ALL), "db", sep = ".")

goID <- resultTable_GLASS_LOLA[1, "GO.ID"]


locResults_progress

geneGO_list = locResults_progress$TF
tmp = annFUN.gene2GO(whichOnto = "BP", feasibleGenes = NULL, gene2GO = geneGO_list)

resultTable_GLASS_LOLA$GO.ID <- factor(resultTable_GLASS_LOLA$GO.ID)
resultTable_GLASS_LOLA$classicFisher <- as.numeric(resultTable_GLASS_LOLA$classicFisher)
resultTable_filtered <- resultTable_GLASS_LOLA[1:15, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[order(resultTable_filtered$classicFisher)])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[order(resultTable_filtered$classicFisher)])

## Just based on Gene Ontology IDs.
pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/public/glass-lme-GO.pdf", height = 5, width = 7)
ggplot(resultTable_filtered, aes(x=Term, y=-log10(classicFisher))) + 
  geom_bar(stat="identity", fill="#CA2F66") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed", col="black") +
  plot_theme +
  labs(x="", y="-log10(p-value)") +
  ylim(0, 3) +
  coord_flip()
dev.off()

###### END ########

## Examine the
glass_betas <- getBeta(il450k_g_filt)

epimut_genes = which(annot_450K$UCSC_RefGene_Name%in%epimut_prom_total_filt$gene_id)
high_epimut_genes = which(annot_450K$UCSC_RefGene_Name%in%epimut_high$gene_id)
glass_epimut_betas = glass_betas[high_epimut_genes, ]
densityPlot(glass_epimut_betas)


