##################################
# TFs with DNAme disorder + TF activity
# Updated: 2021.05.20
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(Seurat)
library(SCENIC)
library(AUCell)
library(ggpubr)
###################################

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

############################
#### TFBS DNAme disorder ###
############################
## Additional information about single-cells passing QC.
tfbs <- read.table(file="data/analysis_RRBS_individual_TFBS_motif_DNAme_disorder.csv", sep = ",", header = TRUE)
hypoxia_tfbs <- tfbs %>% 
  filter(treatment!="RT")
rt_tfbs <- tfbs %>% 
  filter(treatment=="RT")
  

## Split into Unique CpGs per TF AND TFBS PDR.
epimut_cpg_tf_hypoxia = hypoxia_tfbs %>% 
  filter(timepoint=="9 days") %>% 
  dplyr::select(sample, cell_line, timepoint, treatment, ends_with("_num_unique_CpGs")) %>% 
  pivot_longer(
    cols = ends_with("_num_unique_CpGs"),
    names_to = "tf",
    values_to = "num_unique") %>% 
  mutate(tf = gsub("_num_unique_CpGs", "", tf))

epimut_cpg_tf_rt = rt_tfbs %>% 
  dplyr::select(sample, cell_line, timepoint, treatment, ends_with("_num_unique_CpGs")) %>% 
  pivot_longer(
    cols = ends_with("_num_unique_CpGs"),
    names_to = "tf",
    values_to = "num_unique") %>% 
  mutate(tf = gsub("_num_unique_CpGs", "", tf))

## Combine two data.frames.
all(colnames(epimut_cpg_tf_hypoxia)==colnames(epimut_cpg_tf_rt))
epimut_cpg_tf <- bind_rows(epimut_cpg_tf_hypoxia, epimut_cpg_tf_rt)

# Do the same for PDR.
epimut_pdr_tf_hypoxia = hypoxia_tfbs %>% 
  filter(timepoint=="9 days") %>% 
  dplyr::select(sample, cell_line, timepoint, treatment,  ends_with("_PDR")) %>% 
  pivot_longer(
    cols = ends_with("_PDR"),
    names_to = "tf",
    values_to = "pdr") %>% 
  mutate(tf = gsub("_PDR", "", tf))

epimut_pdr_tf_rt = rt_tfbs %>% 
  dplyr::select(sample, cell_line, timepoint, treatment,  ends_with("_PDR")) %>% 
  pivot_longer(
    cols = ends_with("_PDR"),
    names_to = "tf",
    values_to = "pdr") %>% 
  mutate(tf = gsub("_PDR", "", tf))
## Combine two data.frames.
all(colnames(epimut_pdr_tf_hypoxia)==colnames(epimut_pdr_tf_rt))
epimut_pdr_tf <- bind_rows(epimut_pdr_tf_hypoxia, epimut_pdr_tf_rt)


## Combine into single dataframe. 
all(epimut_cpg_tf$sample==epimut_pdr_tf$sample)
all(epimut_cpg_tf$tf==epimut_pdr_tf$tf)
epimut_comb_tf = epimut_cpg_tf %>% 
  inner_join(epimut_pdr_tf, by=c("sample", "cell_line", "timepoint", "treatment", "tf")) 

## Tabulate the median number of CpGs per TF across cell lines.
epimut_comb_tf_tab = epimut_comb_tf %>% 
  group_by(tf) %>% 
  summarise(median_cpgs = median(num_unique))

## Set a filter to retain specific TFs.
epimut_comb_tf_tab_filt = epimut_comb_tf_tab %>% 
  filter(median_cpgs > 30)

## Identify the TFs to keep based on minimum 30 unique CpGs.
tf_keep <- unique(epimut_comb_tf_tab_filt$tf)

## Filter the TFs to keep.
epimut_comb_tf_filt = epimut_comb_tf %>% 
  filter(tf%in%tf_keep) 

## Examine the median TF values between normoxia and rt/hypoxia for each cell line.
epimut_group_med_compare = epimut_comb_tf_filt %>% 
  group_by(tf, cell_line, treatment) %>% 
  summarise(median_pdr = median(pdr)) %>% 
  pivot_wider(names_from = treatment, values_from = median_pdr) %>% 
  mutate(rt_delta_median_pdr = RT-Normoxia,
         hypoxia_delta_median_pdr = Hypoxia-Normoxia) 


## Perform analyses using the relative epimutation burden (Hypoxia/Normoxia).
epimut_group_med =  epimut_comb_tf_filt %>% 
  mutate_each(funs(./median(.[treatment == "Normoxia"])), pdr) %>% 
  group_by(cell_line, tf, treatment) %>%
  summarise(median_relative_epimutation = median(pdr)) %>% 
  pivot_wider(names_from = treatment, values_from = median_relative_epimutation) %>% 
  select(cell_line, tf, Normoxia_rel_epimut = Normoxia, Hypoxia_rel_epimut = Hypoxia, RT_rel_epimut = RT)

## Calculate a wilcoxon rank sum p-value for each RT/hypoxia vs. normoxia comparison.
epimut_group_wilcox_hypoxia = epimut_comb_tf_filt %>% 
  filter(treatment!="RT") %>% 
  group_by(cell_line, tf) %>% 
  do(w = wilcox.test(pdr~treatment, data=., paired=FALSE)) %>% 
  summarise(cell_line, tf, Wilcox = w$p.value) %>% 
  select(cell_line, tf, hypoxia_wilcox = Wilcox)

epimut_group_wilcox_rt = epimut_comb_tf_filt %>% 
  filter(treatment!="Hypoxia") %>% 
  group_by(cell_line, tf) %>% 
  do(w = wilcox.test(pdr~treatment, data=., paired=FALSE)) %>% 
  summarise(cell_line, tf, Wilcox = w$p.value) %>% 
  select(cell_line, tf, rt_wilcox = Wilcox)

## Combine two analyses together.
epimut_group_wilcox = epimut_group_wilcox_hypoxia %>% 
  inner_join(epimut_group_wilcox_rt, by=c("cell_line", "tf"))

## Combine data to identify common modes of shifts in epimutation burden.
epimut_group_med_wilcox <- epimut_group_med %>% 
  inner_join(epimut_group_wilcox, by=c("cell_line", "tf")) %>% 
  ungroup() 


## Which TFBS motifs are statistically significant in all conditions?
tfbs_pvalues <- epimut_group_med_wilcox %>% 
  dplyr::select(cell_line, tf, hypoxia_wilcox, rt_wilcox) %>% 
  pivot_wider(names_from = cell_line, values_from = c(hypoxia_wilcox, rt_wilcox)) %>% 
  filter(hypoxia_wilcox_HF2354 < 0.1, hypoxia_wilcox_HF3016 < 0.1, rt_wilcox_HF2354 < 0.1, rt_wilcox_HF3016 < 0.1)


#############################
##### TF activity      ######
#############################
## Load the SCENIC results for each sample:.
hf2354_umap <- read.table(file="data/analysis_scRNAseq_stress_hf2354_metadata.csv", sep = ",", header = TRUE)

## Load filtered, processed, unnormalized count matrix.
hf2354_obj <- ReadH5AD("data/analysis_scRNAseq_stress_hf2354_expression.h5ad")

## Make a Seurat object with the standard pipeline.
## Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
## This is then natural-log transformed using log1p.
hf2354_obj <- hf2354_obj %>% 
  Seurat::NormalizeData() 

## Use the experimental conditions.
hf2354_obj@meta.data$cell_name = rownames(hf2354_obj@meta.data)

## Downsample for input into SCENIC.
cells_to_sample <- 5000
set.seed(34)
cells_to_subset <- sample(x = hf2354_obj@meta.data$cell_name, size = cells_to_sample, replace = F)

# Subset Seurat object.
hf2354_obj_sub <- SubsetData(object = hf2354_obj, cells = cells_to_subset)

## Restrict the metadata data to that which was included in the SCENIC run.
hf2354_umap_filt <- hf2354_umap %>% 
  filter(cell_barcode%in%hf2354_obj_sub@meta.data$cell_name)

##########################
### Begin SCENIC approach
##########################
## Building the **gene regulatory network (GRN)**: 
## 1. Identify potential targets for each TF based on co-expression.
# - Filtering the expression matrix and running GENIE3/GRNBoost. 
# - Formatting the targets from GENIE3/GRNBoost into co-expression modules. 

## Running SCENIC on a subsample of the data for run time considerations.
## Initialize SCENIC settings:
org="hgnc" 
dbDir="/projects/verhaak-lab/scgp/reference/cisTarget_databases" 
myDatasetTitle="SCENIC HF2354" 
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Save to use at a later time.
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## Load expression matrix.
exprMat <- data.matrix(GetAssayData(hf2354_obj_sub))
cellInfo <- data.frame(hf2354_obj_sub@meta.data$library)
rownames(cellInfo) <- rownames(hf2354_obj_sub@meta.data)
colnames(cellInfo) <- "CellType"

# Color to assign to the variables (same format as for NMF::heatmap)
colVars <- list(CellType=c("RV20018" = "#fcbba1",
                           "RV20033" = "#00FF00"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]


## Save outputs for cell annotation and color code. 
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="int/colVars.Rds")


## Examine how many genes have greater than 0.
cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)

## Filter based on the number of genes.
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)


## Filter the expression matrix only to keep these genes.
## The data is already logged / normalized.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)


## Split the targets into positive- and negative-correlated targets 
## (i.e. Spearman correlation between the TF and the potential target).
runCorrelation(exprMat_filtered, scenicOptions)

### Run GENIE3 (this is computationally intensive). 
runGenie3(exprMat_filtered, scenicOptions)

## 2.  Select potential direct-binding targets (regulons) based on DNA-motif analysis (*RcisTarget*: TF motif analysis) 
## Build and score the GRN.
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

## 1. Get co-expression modules.
runSCENIC_1_coexNetwork2modules(scenicOptions)

## 2. Get regulons (with RcisTarget): TF motif analysis).
runSCENIC_2_createRegulons(scenicOptions)

## 3. Score GRN (regulons) in the cells (with AUCell).
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

## 4. Determine the binarized activities.
runSCENIC_4_aucell_binarize(scenicOptions, skipBoxplot = FALSE, skipHeatmaps = FALSE,
                            skipTsne = FALSE, exprMat = exprMat_filtered)

##### The outputs from the SCENIC runs.
## HF2354 single cells exposed to different stressprs.
auc_rankings_hf2354 <- readRDS("data/SCENIC/HF2354/3.3_aucellRankings.Rds")
regulonAUC_hf2354 <- readRDS("data/SCENIC/HF2354/3.4_regulonAUC.Rds")

## Create a data.frame with the gene sets/TFs and cells.
regulonAUC_hf2354_df = as.data.frame(getAUC(regulonAUC_hf2354))


## generate z-scores for variable A using the scale() function
## scale(A, center = TRUE, scale = TRUE). These are the defaults. 
regulonAUC_hf2354_scaled = t(apply(as.matrix(regulonAUC_hf2354_df), 1, scale))


## Extract the TFs of interest.
regulonAUC_hf2354_scaled_filt <- as.data.frame(t(regulonAUC_hf2354_scaled[grep("RELA|TFDP1|ELK4|KLF16", rownames(regulonAUC_hf2354_scaled)), ]))
regulonAUC_hf2354_scaled_filt$cell_name <- colnames(regulonAUC_hf2354_df)

regulon_filt_9d_annot_hf2354 <- hf2354_umap_filt %>% 
  inner_join(regulonAUC_hf2354_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`TFDP1 (3124g)`, `RELA (406g)`, `ELK4_extended (913g)`, `KLF16_extended (242g)`), 
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days")  %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1))


regulon_filt_9d_annot_hf2354$condition_revalue[regulon_filt_9d_annot_hf2354$condition_revalue=="Irradiation-9d"] <- "Irradiation"
regulon_filt_9d_annot_hf2354$condition_revalue[regulon_filt_9d_annot_hf2354$condition_revalue=="hypoxia-9d"] <- "Hypoxia"
regulon_filt_9d_annot_hf2354$condition_revalue[regulon_filt_9d_annot_hf2354$condition_revalue=="normoxia-9d"] <- "Normoxia"
regulon_filt_9d_annot_hf2354$condition_revalue <- factor(x = regulon_filt_9d_annot_hf2354$condition_revalue, levels = c("Normoxia", "Hypoxia", "Irradiation"))
my_comparisons <- list( c("Normoxia", "Irradiation"), c("Normoxia", "Hypoxia"))

ggplot(regulon_filt_9d_annot_hf2354, aes(x=condition_revalue, y=activity, fill=condition_revalue)) +
  geom_boxplot(outlier.shape = NA) + 
  #ylim(-3, 4) +
  facet_wrap(~cell_line, scales = "free_y") +
  labs(y = "TF activity (Z-score)", x = "Oxygen\nconc.", fill = "Oxygen\nconc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  facet_grid( ~ tf) +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format") 

###################################
##### HF3016 TF activity      #####
###################################
## Load the SCENIC results for each sample:.
hf3016_umap <- read.table(file="data/analysis_scRNAseq_stress_hf3016_metadata.csv", sep = ",", header = TRUE)

## Load filtered, processed, unnormalized count matrix.
hf3016_obj <- ReadH5AD("data/analysis_scRNAseq_stress_hf3016_expression.h5ad")

# Make a Seurat object with the standard pipeline.
## Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
## This is then natural-log transformed using log1p.
hf3016_obj <- hf3016_obj %>% 
  Seurat::NormalizeData() 

## Use the experimental conditions.
hf3016_obj@meta.data$cell_name = rownames(hf3016_obj@meta.data)

## Downsample for input into SCENIC.
cells_to_sample <- 5000
set.seed(34)
cells_to_subset <- sample(x = hf3016_obj@meta.data$cell_name, size = cells_to_sample, replace = F)

# Subset Seurat object.
hf3016_obj_sub <- SubsetData(object = hf3016_obj, cells = cells_to_subset)

## Restrict the metadata data to that which was included in the SCENIC run.
hf3016_umap_filt <- hf3016_umap %>% 
  filter(cell_barcode%in%hf3016_obj_sub@meta.data$cell_name)

##########################
### Begin SCENIC approach
##########################
## Building the **gene regulatory network (GRN)**: 
## 1. Identify potential targets for each TF based on co-expression.
# - Filtering the expression matrix and running GENIE3/GRNBoost. 
# - Formatting the targets from GENIE3/GRNBoost into co-expression modules. 

## Initialize SCENIC settings:
org="hgnc" 
dbDir="/projects/verhaak-lab/scgp/reference/cisTarget_databases" 
myDatasetTitle="SCENIC HF3016" 
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Save to use at a later time.
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## Load expression matrix.
exprMat <- data.matrix(GetAssayData(hf3016_obj_sub))
cellInfo <- data.frame(hf3016_obj_sub@meta.data$library)
rownames(cellInfo) <- rownames(hf3016_obj_sub@meta.data)
colnames(cellInfo) <- "CellType"

# Color to assign to the variables (same format as for NMF::heatmap)
colVars <- list(CellType=c("RV20020" = "#fb6a4a", 
                           "RV20022" = "#fcbba1",
                           "RV20033" = "#00FF00"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]


## Save outputs for cell annotation and color code. 
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="int/colVars.Rds")


## Examine how many genes have greater than 0.
cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)

## Filter based on the number of genes.
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)


## Filter the expression matrix only to keep these genes.
## The data is already logged / normalized.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)


## Split the targets into positive- and negative-correlated targets 
## (i.e. Spearman correlation between the TF and the potential target).
runCorrelation(exprMat_filtered, scenicOptions)

### Run GENIE3 (this is computationally intensive). 
runGenie3(exprMat_filtered, scenicOptions)

## 2.  Select potential direct-binding targets (regulons) based on DNA-motif analysis (*RcisTarget*: TF motif analysis) 
## Build and score the GRN.
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

## 1. Get co-expression modules.
runSCENIC_1_coexNetwork2modules(scenicOptions)

## 2. Get regulons (with RcisTarget): TF motif analysis).
runSCENIC_2_createRegulons(scenicOptions)

## 3. Score GRN (regulons) in the cells (with AUCell).
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

## 4. Determine the binarized activities.
runSCENIC_4_aucell_binarize(scenicOptions, skipBoxplot = FALSE, skipHeatmaps = FALSE,
                            skipTsne = FALSE, exprMat = exprMat_filtered)


#########################
### Load SCENIC results
#########################
auc_rankings_hf3016 <- readRDS("data/SCENIC/HF3016/3.3_aucellRankings.Rds")
regulonAUC_hf3016 <- readRDS("data/SCENIC/HF3016/3.4_regulonAUC.Rds")

## Create a data.frame with the gene sets/TFs and cells.
regulonAUC_hf3016_df = as.data.frame(getAUC(regulonAUC_hf3016))

## generate z-scores for variable A using the scale() function
## scale(A, center = TRUE, scale = TRUE). These are the defaults. 
regulonAUC_hf3016_scaled = t(apply(as.matrix(regulonAUC_hf3016_df), 1, scale))

## Extract the TFs of interest.
regulonAUC_hf3016_scaled_filt <- as.data.frame(t(regulonAUC_hf3016_scaled[grep("RELA|TFDP1|ELK4|KLF16", rownames(regulonAUC_hf3016_scaled)), ]))
regulonAUC_hf3016_scaled_filt$cell_name <- colnames(regulonAUC_hf3016_df)

regulon_filt_9d_annot_hf3016 <- hf3016_umap_filt %>% 
  inner_join(regulonAUC_hf3016_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`TFDP1 (2910g)`, `RELA (650g)`, `ELK4_extended (1222g)`, `KLF16_extended (40g)`), 
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days")  %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1))


regulon_filt_9d_annot_hf3016$condition_revalue[regulon_filt_9d_annot_hf3016$condition_revalue=="Irradiation-9d"] <- "Irradiation"
regulon_filt_9d_annot_hf3016$condition_revalue[regulon_filt_9d_annot_hf3016$condition_revalue=="hypoxia-9d"] <- "Hypoxia"
regulon_filt_9d_annot_hf3016$condition_revalue[regulon_filt_9d_annot_hf3016$condition_revalue=="normoxia-9d"] <- "Normoxia"
regulon_filt_9d_annot_hf3016$condition_revalue <- factor(x = regulon_filt_9d_annot_hf3016$condition_revalue, levels = c("Normoxia", "Hypoxia", "Irradiation"))
my_comparisons <- list( c("Normoxia", "Irradiation"), c("Normoxia", "Hypoxia"))

ggplot(regulon_filt_9d_annot_hf3016, aes(x=condition_revalue, y=activity, fill=condition_revalue)) +
  geom_boxplot(outlier.shape = NA) + 
  #ylim(-3, 4) +
  facet_wrap(~cell_line, scales = "free_y") +
  labs(y = "TF activity (Z-score)", x = "Oxygen\nconc.", fill = "Oxygen\nconc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  facet_grid( ~ tf) +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format") 


##############################
#### Combine all data    #####
##############################
## Print the TFs with consistently stress-altered TFBS motif DNAme disorder
tfbs_pvalues$tf

## How many are covered in the SCENIC analysis?
rownames(regulonAUC_hf2354_scaled[grep("ELK4|FOSL2|KLF10|KLF16|MEF2D|NFYB|NR2C1|REL_|RELA|RFX3|TFDP1|VEZF1|ZNF263|ZNF384", rownames(regulonAUC_hf2354_scaled)), ])
rownames(regulonAUC_hf3016_scaled[grep("ELK4|FOSL2|KLF10|KLF16|MEF2D|NFYB|NR2C1|REL_|RELA|RFX3|TFDP1|VEZF1|ZNF263|ZNF384", rownames(regulonAUC_hf3016_scaled)), ])


##### HF2354 ############
regulonAUC_hf2354_scaled_filt <- as.data.frame(t(regulonAUC_hf2354_scaled[grep("ELK4|FOSL2|KLF10|KLF16|MEF2D|NFYB|NR2C1|REL_|RELA|RFX3|TFDP1|VEZF1|ZNF263|ZNF384", rownames(regulonAUC_hf2354_scaled)), ]))
regulonAUC_hf2354_scaled_filt$cell_name <- colnames(regulonAUC_hf2354_df)

tf_activity_wilcox_hypoxia_hf2354 <- hf2354_umap_filt %>% 
  inner_join(regulonAUC_hf2354_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`NFYB_extended (1712g)`:`ELK4_extended (913g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days", condition_revalue!="Irradiation-9d") %>% 
  group_by(tf, cell_line) %>% 
  do(w = wilcox.test(activity~condition_revalue, data=., paired=FALSE)) %>% 
  summarise(tf, cell_line, Wilcox = w$p.value) %>% 
  dplyr::select(cell_line, tf, hypoxia_wilcox = Wilcox) 

tf_activity_wilcox_rt_hf2354 <- hf2354_umap_filt %>% 
  inner_join(regulonAUC_hf2354_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`NFYB_extended (1712g)`:`ELK4_extended (913g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days", condition_revalue!="hypoxia-9d") %>% 
  group_by(tf, cell_line) %>% 
  do(w = wilcox.test(activity~condition_revalue, data=., paired=FALSE)) %>% 
  summarise(tf, cell_line, Wilcox = w$p.value) %>% 
  dplyr::select(cell_line, tf, rt_wilcox = Wilcox) 

tf_activity_median_wilcox_hf2354 <- hf2354_umap_filt %>% 
  inner_join(regulonAUC_hf2354_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`NFYB_extended (1712g)`:`ELK4_extended (913g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days") %>% 
  group_by(tf, cell_line, condition_revalue) %>%
  summarise(avg_activity = median(activity)) %>% 
  pivot_wider(names_from = condition_revalue, values_from = avg_activity) %>% 
  mutate(rt_delta_median_activity = `Irradiation-9d`-`normoxia-9d`,
         hypoxia_delta_median_activity = `hypoxia-9d`-`normoxia-9d`) %>% 
  inner_join(tf_activity_wilcox_hypoxia_hf2354, by=c("tf", "cell_line")) %>% 
  inner_join(tf_activity_wilcox_rt_hf2354, by=c("tf", "cell_line")) %>% 
  mutate(rt_group = ifelse(rt_delta_median_activity>0 & rt_wilcox < 0.05, "upregulation",  ifelse(rt_delta_median_activity<0 & rt_wilcox < 0.05, "downregulation", "no change")),
         hypoxia_group = ifelse(hypoxia_delta_median_activity>0 & hypoxia_wilcox < 0.05, "upregulation",  ifelse(hypoxia_delta_median_activity<0 & hypoxia_wilcox < 0.05, "downregulation", "no change"))) %>% 
  dplyr::select(cell_line, tf, hypoxia_group, rt_group) %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1)) %>% 
  distinct()

tf_activity_median_wilcox_hf2354 <- hf2354_umap_filt %>% 
  inner_join(regulonAUC_hf2354_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`NFYB (1192g)`,`TFDP1 (3124g)`, `FOSL2 (62g)`, `KLF16_extended (242g)`,
                 `ZNF384 (116g)`, `REL_extended (56g)`,  `RELA (406g)`,
                 `RFX3 (16g)`, `VEZF1_extended (38g)`, `ELK4_extended (913g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days") %>% 
  group_by(tf, cell_line, condition_revalue) %>%
  summarise(avg_activity = median(activity)) %>% 
  pivot_wider(names_from = condition_revalue, values_from = avg_activity) %>% 
  mutate(rt_delta_median_activity = `Irradiation-9d`-`normoxia-9d`,
         hypoxia_delta_median_activity = `hypoxia-9d`-`normoxia-9d`) %>% 
  inner_join(tf_activity_wilcox_hypoxia_hf2354, by=c("tf", "cell_line")) %>% 
  inner_join(tf_activity_wilcox_rt_hf2354, by=c("tf", "cell_line")) %>% 
  mutate(rt_group = ifelse(rt_delta_median_activity>0, "upregulation", "downregulation"),
         hypoxia_group = ifelse(hypoxia_delta_median_activity>0, "upregulation", "downregulation")) %>% 
  dplyr::select(cell_line, tf, hypoxia_group, hypoxia_wilcox, rt_group, rt_wilcox) %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1)) %>% 
  distinct()


###### HF3016 ########
regulonAUC_hf3016_scaled_filt <- as.data.frame(t(regulonAUC_hf3016_scaled[grep("ELK4|FOSL2|KLF10|KLF16|MEF2D|NFYB|NR2C1|REL_|RELA|RFX3|TFDP1|VEZF1|ZNF263|ZNF384", rownames(regulonAUC_hf3016_scaled)), ]))
regulonAUC_hf3016_scaled_filt$cell_name <- colnames(regulonAUC_hf3016_df)


tf_activity_wilcox_hypoxia_hf3016 <- hf3016_umap_filt %>% 
    inner_join(regulonAUC_hf3016_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
    pivot_longer(c(`FOSL2_extended (649g)`:`RELA (650g)`),             
                 names_to = "tf", 
                 values_to =  "activity") %>% 
    filter(timepoint=="9days", condition_revalue!="Irradiation-9d") %>% 
  group_by(tf, cell_line) %>% 
  do(w = wilcox.test(activity~condition_revalue, data=., paired=FALSE)) %>% 
  summarise(tf, cell_line, Wilcox = w$p.value) %>% 
  dplyr::select(cell_line, tf, hypoxia_wilcox = Wilcox) 

tf_activity_wilcox_rt_hf3016 <- hf3016_umap_filt %>% 
  inner_join(regulonAUC_hf3016_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`FOSL2_extended (649g)`:`RELA (650g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days", condition_revalue!="hypoxia-9d") %>% 
  group_by(tf, cell_line) %>% 
  do(w = wilcox.test(activity~condition_revalue, data=., paired=FALSE)) %>% 
  summarise(tf, cell_line, Wilcox = w$p.value) %>% 
  dplyr::select(cell_line, tf, rt_wilcox = Wilcox) 

tf_activity_median_wilcox_hf3016 <- hf3016_umap_filt %>% 
  inner_join(regulonAUC_hf3016_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`FOSL2_extended (649g)`:`RELA (650g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days") %>% 
  group_by(tf, cell_line, condition_revalue) %>%
  summarise(avg_activity = median(activity)) %>% 
  pivot_wider(names_from = condition_revalue, values_from = avg_activity) %>% 
  mutate(rt_delta_median_activity = `Irradiation-9d`-`normoxia-9d`,
         hypoxia_delta_median_activity = `hypoxia-9d`-`normoxia-9d`) %>% 
  inner_join(tf_activity_wilcox_hypoxia_hf3016, by=c("tf", "cell_line")) %>% 
  inner_join(tf_activity_wilcox_rt_hf3016, by=c("tf", "cell_line")) %>% 
  mutate(rt_group = ifelse(rt_delta_median_activity>0 & rt_wilcox < 0.05, "upregulation",  ifelse(rt_delta_median_activity<0 & rt_wilcox < 0.05, "downregulation", "no change")),
         hypoxia_group = ifelse(hypoxia_delta_median_activity>0 & hypoxia_wilcox < 0.05, "upregulation",  ifelse(hypoxia_delta_median_activity<0 & hypoxia_wilcox < 0.05, "downregulation", "no change"))) %>% 
  dplyr::select(cell_line, tf, hypoxia_group, rt_group) %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1)) %>% 
  distinct()

tf_activity_median_wilcox_hf3016 <- hf3016_umap_filt %>% 
  inner_join(regulonAUC_hf3016_scaled_filt, by=c("cell_barcode" ="cell_name")) %>% 
  pivot_longer(c(`FOSL2 (490g)`,`TFDP1 (2910g)`, `KLF16_extended (40g)`, `RFX3_extended (16g)`,
                 `ELK4_extended (1222g)`, `NFYB (13g)`, `MEF2D (11g)`, `ZNF263_extended (166g)`,
                 `RELA (650g)`),             
               names_to = "tf", 
               values_to =  "activity") %>% 
  filter(timepoint=="9days") %>% 
  group_by(tf, cell_line, condition_revalue) %>%
  summarise(avg_activity = median(activity)) %>% 
  pivot_wider(names_from = condition_revalue, values_from = avg_activity) %>% 
  mutate(rt_delta_median_activity = `Irradiation-9d`-`normoxia-9d`,
         hypoxia_delta_median_activity = `hypoxia-9d`-`normoxia-9d`) %>% 
  inner_join(tf_activity_wilcox_hypoxia_hf3016, by=c("tf", "cell_line")) %>% 
  inner_join(tf_activity_wilcox_rt_hf3016, by=c("tf", "cell_line")) %>% 
  mutate(rt_group = ifelse(rt_delta_median_activity>0, "upregulation", "downregulation"),
         hypoxia_group = ifelse(hypoxia_delta_median_activity>0, "upregulation", "downregulation")) %>% 
  dplyr::select(cell_line, tf, hypoxia_group, hypoxia_wilcox, rt_group, rt_wilcox) %>% 
  mutate(tf = gsub("_extended", "", tf),
         tf = sapply(strsplit(tf, " "), "[[", 1)) %>% 
  distinct()

comb_tf_activity <- tf_activity_median_wilcox_hf3016 %>% 
  bind_rows(tf_activity_median_wilcox_hf2354) %>% 
  dplyr::select(cell_line, tf, hypoxia_group, rt_group) %>% 
  pivot_longer(
    cols = c(hypoxia_group, rt_group),
    names_to = "stress",
    values_to = "activity") %>%
  mutate(stress = recode(stress, "hypoxia_group" = "Hypoxia", "rt_group" = "Irradiation")) %>% 
  filter(tf%in%c("TFDP1", "RFX3", "RELA", "NFYB", "KLF16", "FOSL2", "ELK4"))

comb_tf_pvalue <- tf_activity_median_wilcox_hf3016 %>% 
  bind_rows(tf_activity_median_wilcox_hf2354) %>% 
  dplyr::select(cell_line, tf, hypoxia_wilcox, rt_wilcox) %>% 
  pivot_longer(
    cols = c(hypoxia_wilcox, rt_wilcox),
    names_to = "stress",
    values_to = "pvalue") %>% 
  mutate(stress = recode(stress, "hypoxia_wilcox" = "Hypoxia", "rt_wilcox" = "Irradiation")) %>% 
  filter(tf%in%c("TFDP1", "RFX3", "RELA", "NFYB", "KLF16", "FOSL2", "ELK4"))

comb_tf_all <- comb_tf_pvalue %>% 
  inner_join(comb_tf_activity, by=c("cell_line", "tf", "stress")) %>% 
  mutate(pvalue = ifelse(pvalue <2.2e-16, 2.2e-16, pvalue))

comb_tf_all$activity <- factor(x = comb_tf_all$activity, levels = c("upregulation","downregulation"))

pdf(file = "results/Fig4/SuppFig8d-tf-activity.pdf", height = 6, width = 4, bg = "transparent", useDingbats = FALSE)
ggplot() +
  geom_tile(data = comb_tf_all, aes(x = stress, y = tf, fill = activity, alpha = -log10(pvalue)), color = "black") +
  labs(y = "", fill = "9-day TF activity \nchange from normoxia") +
  scale_fill_manual(values = c("upregulation" = "#CD4F39", "downregulation" = "#27408B", "no change" = "#d3d3d3",
                               "NA" = "#FFFFFF")) +
  plot_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(. ~ cell_line, scales="free") +
  ggtitle("Selected TFs with significant TFBS motif\nDNAme disorder changes")
dev.off()

####################################################################

regulon_filt_9d_annot_hf3016_sub <- regulon_filt_9d_annot_hf3016 %>% 
  dplyr::select(cell_line, stress = condition_revalue, tf, value = activity) %>% 
  mutate(type = "tf_activity")
regulon_filt_9d_annot_hf2354_sub <- regulon_filt_9d_annot_hf2354 %>% 
  dplyr::select(cell_line, stress = condition_revalue, tf, value = activity) %>% 
  mutate(type = "tf_activity")
tf_activity <- regulon_filt_9d_annot_hf2354_sub %>% 
  bind_rows(regulon_filt_9d_annot_hf3016_sub)


tfbs_disorder <- epimut_comb_tf_filt %>% 
  group_by(cell_line, tf) %>%
  mutate_each(funs(./median(.[treatment == "Normoxia"])), pdr) %>% 
  filter(tf%in%c("ELK4", "TFDP1")) %>% 
  mutate(type = "tf_disorder",
         treatment = recode(treatment, `RT` = "Irradiation")) %>% 
  select(cell_line, stress = treatment, tf, value = pdr, type)

tfbs_disorder$stress <- factor(x = tfbs_disorder$stress, levels = c("Normoxia", "Hypoxia", "Irradiation"))
my_comparisons <- list( c("Normoxia", "Irradiation"), c("Normoxia", "Hypoxia"))

pdf(file = "results/Fig4/SuppFig8e-f-tf-disorder.pdf", height = 3, width = 6, bg = "transparent", useDingbats = FALSE)
ggplot(tfbs_disorder, aes(x=tf, y=value, fill= stress)) +
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y = "Relative TFBS motif\nDNAme disorder", x = "", fill = "9-day\nstress exposure") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ cell_line, scales ="free") 
dev.off()

## Extract p-values
ggplot(tfbs_disorder, aes(x=stress, y=value, fill= stress)) +
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y = "Relative TFBS motif\nDNAme disorder", x = "", fill = "9-day\nstress exposure") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(tf ~ cell_line, scales ="free") +
  stat_compare_means(comparisons = my_comparisons, label = "p.format")


tf_activity <- regulon_filt_9d_annot_hf2354_sub %>% 
  bind_rows(regulon_filt_9d_annot_hf3016_sub) %>% 
  filter(tf%in%c("ELK4", "TFDP1"))

pdf(file = "results/Fig4/SuppFig8e-tf-activity.pdf", height = 3, width = 6, bg = "transparent", useDingbats = FALSE)
ggplot(tf_activity, aes(x=tf, y=value, fill=stress)) +
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y = "SCENIC TF activity\n(Z-score)", x = "", fill = "9-day\nstress exposure") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(. ~ cell_line, scales="free") 
dev.off()

## Extract p-values.
ggplot(tf_activity, aes(x=stress, y=value, fill=stress)) +
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  labs(y = "SCENIC TF activity (Z-score)", x = "", fill = "9-day\nstress exposure") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c",
                             "Irradiation" = "#0dba86")) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(tf ~ cell_line, scales="free") +
  stat_compare_means(comparisons = my_comparisons, label = "p.format") 

#### END #####
