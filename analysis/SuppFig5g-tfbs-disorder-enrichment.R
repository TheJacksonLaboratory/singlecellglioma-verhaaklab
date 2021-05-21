##################################
# Test for enrichment of TFs across high TFBS DNAme disorder sites.
# Updated: 2021.05.18
# Author: Kevin J.
###################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)

###################################
## Load the essential packages.
library(openxlsx)
library(ggrepel)
library(corrplot)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(topGO)
library(ggpubr)
library(tidyverse)
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


## Load the CpG density file for the scRRBS data.
cpg_density <- read.table(file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/Samples-passQC_single_tumor_cells_epiallele_CpG_density_summary.txt", sep="\t", header=T, stringsAsFactors = F)
cpg_density$TF <- gsub("\\(", "\\.", cpg_density$TF)
cpg_density$TF <- gsub("\\)", "\\.", cpg_density$TF)

## Load the SCGP subject-level metadata.
full_meta_data = readWorkbook("/Users/johnsk/mnt/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

## Need to extract subtype and IDHmut status.
meta_data = full_meta_data %>% 
  select(case_barcode = subject_id, subtype) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  mutate(case_barcode_short = gsub("-", "", substr(case_barcode, 6, 11)))

## The epimutation rate was calculated for TFs in the JASPAR database, including normal samples.
epimut_TFBS <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC-single_tumor_cells_TFBS_epimutation_rate.txt", sep="\t", header=T, stringsAsFactors = F)
epimut_TFBS_pdr = epimut_TFBS %>% 
  select(sample_barcode = patient, Ar_fraction:ZSCAN29_PDR)

## Load in the data with the total number of epialleles.
epiallele_info <- read.table(file="/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/rerun-reformatted_deduplicated-final_recalculated/Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", sep="\t", header=T, stringsAsFactors = F)
epiallele_sub = epiallele_info %>% 
  select(sample_barcode, epiallele_total_reads) 

## Restrict to only tumor cells (n = 844) because the normal cells were relabelled to have sample_barcode of "normal".
epimut_total = epimut_TFBS_pdr %>% 
  inner_join(epiallele_sub, by="sample_barcode")

## Rearrange the column variable names to make it easier to pivot to long format.
colnames(epimut_total) <- sub("^(.*)_(.*)$", "\\2_\\1", colnames(epimut_total))

## Split into reads per epiallele AND TF binding site specific pdr/epimutation rate.
epimut_epiallele_tf = epimut_total %>% 
  select(sample_barcode = barcode_sample, epiallele_total_reads = reads_epiallele_total, starts_with("fraction_")) %>% 
  pivot_longer(
    cols = starts_with("fraction_"),
    names_to = "tf",
    names_prefix = "fraction_",
    values_to = "fraction") %>% 
  mutate(epiallele_tf = fraction*epiallele_total_reads)
# PDR.
epimut_pdr_tf = epimut_total %>% 
  select(sample_barcode = barcode_sample, starts_with("PDR_")) %>% 
  pivot_longer(
    cols = starts_with("PDR_"),
    names_to = "tf",
    names_prefix = "PDR_",
    values_to = "pdr")

## Combine into single dataframe. 
all(epimut_epiallele_tf$sample_barcode==epimut_pdr_tf$sample_barcode)
epimut_comb_tf = epimut_epiallele_tf %>% 
  inner_join(epimut_pdr_tf, by=c("sample_barcode", "tf")) %>% 
  mutate(case_barcode = substr(sample_barcode, 1, 11)) %>% 
  select(sample_barcode, case_barcode, tf, epiallele_tf, pdr) %>% 
  left_join(meta_data, by="case_barcode") %>% 
  mutate(case_barcode = gsub("-", "", substr(case_barcode, 6, 11)))

## Set the case_order as in all other figures.
case_order <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
epimut_comb_tf <- epimut_comb_tf %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

## Generate some visuals of the number of supporting reads per TFBS.
## What's the correlation?
cor.test(epimut_comb_tf$epiallele_tf, epimut_comb_tf$pdr)
## There seems to be a good number of outliers, but a relatively weak association.
ggplot(epimut_comb_tf, aes(x=epiallele_tf, y=pdr)) + geom_point()

## Generate boxplots for each sample.
ggplot(epimut_comb_tf, aes(x=case_barcode, y=epiallele_tf)) + geom_boxplot(outlier.shape = NA) +
  ylim(0, 500) + 
  plot_theme + 
  facet_grid(~subtype, scales = "free_x", space = "free")

## Overall, the first quartile seems to be listed at ~88 supporting an observation.
# For an accurate estimation, let's select at least 88 reads per TF.
summary(epimut_comb_tf$epiallele_tf)
## This only removes cell+TF-specific pairings of PDR estimates setting a minimum of 100 reads supporting TFBS.
## This should provide a more stable estimate of TFBS epimutation rate.
## Combine with K2's list of TFs to drop.
tfs_to_drop <-  c("HLTF", "EN1", "FOXD1", "FOXG1", "FOXN3",
                  "FOXO3", "MEIS3", "MSX1", "PRRX1", "Shox2", "VAX2")
unique(epimut_comb_tf$tf)[unique(epimut_comb_tf$tf)%in%tfs_to_drop]
epimut_comb_tf_filt = epimut_comb_tf %>% 
  filter(epiallele_tf > 100, !tf%in%tfs_to_drop)
n_distinct(epimut_comb_tf$sample_barcode)
n_distinct(epimut_comb_tf_filt$sample_barcode)
n_distinct(epimut_comb_tf_filt$tf)

## Save output for downstream analyses.
saveRDS(epimut_comb_tf_filt, file="/Users/johnsk/Documents/Single-Cell-DNAmethylation/epimut_comb_tf_filt.rds")

## This preliminary filtering also ensures that estimates across cells would be more likely to be stable
# with all cell-TFBS pairings having AT LEAST 100 reads. 
epimut_comb_tf_tab = epimut_comb_tf_filt %>% 
  group_by(tf, case_barcode) %>% 
  summarise(avg_supporting_reads = mean(epiallele_tf))
## Sanity check.
sum(epimut_comb_tf_tab$avg_supporting_reads<100)

## Calculate the average PDR per TFBS and per subject.
epimut_comb_tf_filt_pdr = epimut_comb_tf_filt %>% 
  group_by(tf, case_barcode) %>% 
  summarise(avg_epimutation = mean(pdr)) %>% 
  ungroup()

## Create a sort for each TF based on PDR.
sort_df <- epimut_comb_tf_filt_pdr %>%
  group_by(tf) %>% 
  summarise(epimutation = median(avg_epimutation)) %>% 
  arrange(desc(epimutation))  %>% 
  ungroup()
tf_order <- unique(sort_df$tf)
epimut_comb_tf_filt_pdr <- epimut_comb_tf_filt_pdr %>% 
  mutate(tf = factor(tf, levels = tf_order))


## What about creating a dot-plot with lines connecting two subtypes?
epimut_comb_tf_filt_merge = epimut_comb_tf_filt_pdr %>% 
  left_join(meta_data, by=c("case_barcode"="case_barcode_short")) %>% 
  select(-case_barcode.y) %>% 
  group_by(tf, idh_status) %>% 
  summarise(tf_subtype_pdr = median(avg_epimutation)) %>% 
  ungroup() %>% 
  spread(idh_status, tf_subtype_pdr) %>% 
  mutate(pdr_diff = IDHwt-IDHmut,
         pdr_change = ifelse(pdr_diff>=0, "+", "-")) %>% 
  ## Combine with CpG density
  inner_join(cpg_density, by=c("tf"="TF"))



## Remove "var2|var3" and "." from TFs.
n_distinct(epimut_comb_tf_filt_merge$tf)
epimut_comb_tf_filt_merge$tf = gsub("\\.", "", epimut_comb_tf_filt_merge$tf)
epimut_comb_tf_filt_merge$tf = gsub("var[0-9]", "", epimut_comb_tf_filt_merge$tf)

## Remove any non-human TFs (inclusion of lowercase letters).
epimut_comb_tf_filt_subtype = epimut_comb_tf_filt_merge %>% 
  filter(!grepl("[a-z]", tf))

## Are the IDHwt average TFBS epimutation burden estimates confounded by CG density?
epimut_comb_tf_filt_subtype_plot = epimut_comb_tf_filt_subtype %>% 
  gather(key = "density_metric", value = "cg_density", c(nearby_motif.25bp_AUC, nearby_motif.50bp_AUC, nearby_motif.2r_AUC, nearby_motif.3r_AUC))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/cg_density_tfbs_mean_epimutation_corr.pdf", width = 10, height = 7)
ggplot(epimut_comb_tf_filt_subtype_plot, aes(x=cg_density,  y=IDHwt)) +
  geom_point() +
  plot_theme +
  facet_grid(~density_metric, scales = "free") +
  stat_cor(method="s") +
  labs(x="TFBS CpG density metric", y="IDHwt subtype mean TFBS epimutation")
dev.off()


#######################################
### Annotate high TFBS epimutation TFs
#######################################
## Define the number of distinct TFs.
n_distinct(epimut_comb_tf_filt_subtype$tf)
epimut_comb_tf_filt_subtype = epimut_comb_tf_filt_subtype %>% 
  distinct()

## Create a vector of the IDHwt TFBS values (n = 204). There are some missing values for IDHwt + IDHmut, respectively.
epimut_vector <- epimut_comb_tf_filt_subtype$IDHwt


#######################################
## ROSE - super enhancer detection code
#######################################
#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    plot(1:length(inputVector), inputVector,type="l",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(b,slope),col=2)
    #title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    title(paste("High TFBS epimutation value=", signif(y_cutoff,3), sep=""), cex.main=0.8)
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="black",col="black") #Number of regions with zero signal
  }
  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

## Create visual for determining high TFBS epimutation via ROSE.
pdf("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/epimutation-TFBS-rose.pdf", width=7, height=5, useDingbats = FALSE)
calculate_cutoff(epimut_vector, ylab="TFBS epimutation", xlab="Ranked TFBS epimutation")
dev.off()

####################
### topGO analysis
####################
GO_prep = epimut_comb_tf_filt_subtype %>% 
  filter(!is.na(IDHwt))
gene_list <- GO_prep$IDHwt
names(gene_list) <- GO_prep$tf
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
resultTable <- GenTable(GOdata, 
                        classicFisher = resultFisher, 
                        classicKS = resultKS, 
                        topNodes = 30, 
                        ranksOf = 'classicFisher'
)

# Look at the first few lines
hist(as.numeric(resultTable$classicFisher))
resultTable$adj_pvalue <- p.adjust(resultTable$classicFisher, method = "fdr", n = length(resultTable$classicFisher))
head(resultTable)



resultTable$GO.ID <- factor(resultTable$GO.ID)
resultTable$classicFisher <- as.numeric(resultTable$classicFisher)
resultTable_filtered <- resultTable[1:15, ]
resultTable_filtered$GO.ID <- factor(resultTable_filtered$GO.ID, levels = resultTable_filtered$GO.ID[rev(order(resultTable_filtered$adj_pvalue))])
resultTable_filtered$Term <- factor(resultTable_filtered$Term, levels = resultTable_filtered$Term[rev(order(resultTable_filtered$adj_pvalue))])

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5g-TFBS-disorder-GOid.pdf", height = 5, width = 3.5, useDingbats = FALSE)
ggplot(resultTable_filtered, aes(x=GO.ID, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#7fbf7b") + 
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  coord_flip()
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig2/SuppFig5g-TFBS-disorder-GOterm.pdf", height = 5, width = 3.5, useDingbats = FALSE)
ggplot(resultTable_filtered, aes(x=Term, y=-log10(adj_pvalue))) + 
  geom_bar(stat="identity", fill="#7fbf7b") + 
  plot_theme +
  labs(x="", y="-log10(adj. p-value)") +
  coord_flip()
dev.off()


#### END #####

