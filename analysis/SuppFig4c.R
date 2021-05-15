# mutRate_vs_epimutRate.R - Correlates cell DNAme disorder with bulk mutation rates

require(VariantAnnotation)
require(genomation)
require(GenomicRanges)
require(plyr)
require(tidyverse)
require(RColorBrewer)
require(ggplot2)


### Input arguments ###

# Input directory for variant calls
variantDir

# Input directory for reference annotations
annotationDir

# Input directory for WgsMetrics coverage statistics output for sample BAM files
wgsMetricsDir

# Input directory for genomic context-specific coverage statistics
# pre-calculated by read_context_WgsMetrics.R
context_coverageDir 

# Input directory for calculated DNAme disorder
PDRDir 

# Output directory
outDir 
###


# Identify vcf files for each sample
variantFiles <- list.files(path = variantDir, pattern = "*WGS.filtered.vcf$")

# Extract first portion of file names to determine sample names
bulkSamples <- lapply(strsplit(variantFiles,"-"),'[',c(1:3))
bulkSamples <- unlist(lapply(bulkSamples, function(x) paste(x,collapse = "-")))

# Load reference annotations
# NOTE: due to low feature length, removing TSS from analysis
annos <- list(
  exon=readBed(paste0(annotationDir,"/","Ensembl_b37_exon-sorted.bed")),
  gene_body=readBed(paste0(annotationDir,"/","Ensembl_b37_gene_body-sorted.bed")),
  intergenic=readBed(paste0(annotationDir,"/","Ensembl_b37_intergenic-sorted.bed")),
  intron=readBed(paste0(annotationDir,"/","Ensembl_b37_intron-sorted.bed")),
  enhancer=readBed(paste0(annotationDir,"/","FANTOM5_b37_human_permissive_enhancers_phase_1_and_2-sorted.bed")),
  promoter=readBed(paste0(annotationDir,"/","FANTOM5_b37_gene_matched_TSS.padded_1500u_500d-sorted.bed")),
#  tss=readBed(paste0(annotationDir,"/","FANTOM5_b37_gene_matched_TSS-sorted.bed")),
  cgi=readBed(paste0(annotationDir,"/","UCSC_b37_CGI-sorted.bed")),
  dnaseI=readBed(paste0(annotationDir,"/","UCSC_b37_DNaseI_hypersensitive_sites-sorted.bed"))
)

# Example variants@metadata$header@samples
# [1] "SCGP-SM-001-NB-01D-WGS-BHJXI2" "SCGP-SM-001-R1-01D-WGS-FC04O9"


##### For each sample, load global and context-specific output from WgsMetrics #####
##### and use to calculate global and context-specific read coverage           #####

# List sample wgsMetrics files
wgsMetricsFiles <- list.files(path = wgsMetricsDir, pattern = "*.WgsMetrics.txt", full.names = TRUE)
# Exclude normal samples
if (any(grep("NB",wgsMetricsFiles))) {wgsMetricsFiles <- wgsMetricsFiles[-grep("NB",wgsMetricsFiles)]}

# Extract sample name from wgsMetrics filename
wgsMetricsSamples <- lapply(strsplit(wgsMetricsFiles,"-"),'[',c(1:3))
wgsMetricsSamples <- unlist(lapply(wgsMetricsSamples, function(x) paste(x,sep = "",collapse = "-")))

# Initialize output sample coverage table
sample_coverage <- as.data.frame(matrix(nrow = length(bulkSamples), ncol = length(annos) + 1))
colnames(sample_coverage) <- c("global",names(annos))
rownames(sample_coverage) <- wgsMetricsSamples

# Load global coverage metrics; ignore the first 10 lines of input 
# and only retain output for coverage split by base depth
wgs_coverage.global <- lapply(wgsMetricsFiles, function(x) read.delim(x,header = TRUE, skip = 10))
names(wgs_coverage.global) <- wgsMetricsSamples

# Calculate sample coverage by summing high quality coverage count 
# for all bases between 14x and 250x
sample_coverage$global <- unlist(lapply(wgs_coverage.global,function(x) sum(x[,2][c(15:251)])))

### Load context specific coverage outputs and add to sample coverage table
# NOTE: due to low feature length, removing TSS from analysis

# Load context-specific coverage metrics
wgsContextCoverageFiles <- list.files(path = context_coverageDir, pattern = "*context_coverage.txt", full.names = TRUE)
wgs_coverage.context <- lapply(wgsContextCoverageFiles, function(x) read.delim(x,header = FALSE, stringsAsFactors = FALSE))

# Adjust context-specific coverage metrics table context names to match reference annotations
# and reorder to match reference annotation order before adding to sample coverage table
for (i in 1:length(wgs_coverage.context)) {
  if (any(grep("tss",tolower(wgs_coverage.context[[i]][,1])))) {
    wgs_coverage.context[[i]] <- wgs_coverage.context[[i]][-grep("tss",tolower(wgs_coverage.context[[i]][,1])),]
  }
  wgs_coverage.context[[i]]$V1[which(wgs_coverage.context[[i]]$V1 == "CGI")] <- "cgi"
  wgs_coverage.context[[i]]$V1[which(wgs_coverage.context[[i]]$V1 == "DNaseI")] <- "dnaseI"
  wgs_coverage.context[[i]]$V1[which(wgs_coverage.context[[i]]$V1 == "gene")] <- "gene_body"
  wgs_coverage.context[[i]]$V1 <- factor(wgs_coverage.context[[i]]$V1, levels = names(annos))
  wgs_coverage.context[[i]] <- wgs_coverage.context[[i]][order(wgs_coverage.context[[i]]$V1),]
  sample_coverage[i,-1] <- t(wgs_coverage.context[[i]]$V2)
}

write.table(sample_coverage, file = paste0(outDir,"/","Samples_WGS_coverage.txt"), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
#####


##### For each sample, calculate global and context-specific mutation counts #####

# Initialize output sample mutation counts
sample_varCount <- as.data.frame(matrix(nrow = length(bulkSamples), ncol = length(annos) + 1))
rownames(sample_varCount) <- bulkSamples
colnames(sample_varCount) <- c("global",names(annos))

# Count mutations for each patient
for (i in 1:length(bulkSamples)) {
  
  # Load sample vcf file
  variants <- readVcf(paste0(variantDir,"/",variantFiles[grep(bulkSamples[i],variantFiles)]))
  
  # Apply quality filters to vcfs
  median_mapping_quality <- unlist(lapply(variants@info$MMQ,'[',2))
  median_mapping_quality.pass <- which(median_mapping_quality >= 20)
  median_base_quality <- unlist(lapply(variants@info$MBQ,'[',2))
  median_base_quality.pass <- which(median_base_quality >= 20)
  read_depth <- geno(variants)$DP[,2]
  read_depth.pass <- which(read_depth >= 14 & read_depth <= 250)
  
  ind_variants.pass <- Reduce(intersect,list(median_mapping_quality.pass,median_base_quality.pass,read_depth.pass))
  
  # Calculate global variant counts (all "PASS" variants that pass filters)
  varCount.global <- length(which(rowRanges(variants)$FILTER[ind_variants.pass] == "PASS"))
  
  # Calculate variant count per genomic context
  varCount.context <- vector(mode = "numeric", length = length(annos))
  names(varCount.context) <- names(annos)
  
  for (j in 1:length(annos)) {
    context_overlap <- findOverlaps(rowRanges(variants)[ind_variants.pass],annos[[j]])
    
    varCount.context[j] <- length(which(rowRanges(variants)$FILTER[ind_variants.pass][context_overlap@from] == "PASS"))
  }
  
  # Combine global and context-specific variant counts
  varCount <- c(global=varCount.global,varCount.context)
  sample_varCount[i,] <- varCount
}

# Write global and context-specific variant counts to file
write.table(sample_varCount, file = paste0(outDir,"/","Samples_WGS_mutation_counts.txt"), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Calculate global and context-specific mutation rates as mutations per megabase
sample_mutRate <- (sample_varCount/sample_coverage) * 10^6

# Write mutation rates to file
write.table(sample_mutRate, file = paste0(outDir,"/","Samples_WGS_mutation_rates.txt"), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
#####


##### Compare global and context-specific mutation rates with cell DNAme disorder #####

# Load cell DNAme disorder
cells_epr <- read.delim("Samples-passQC_single_cells_global_and_context-specific_pdr_score_table.txt", 
                        stringsAsFactors = FALSE)

# Load patient metadata
meta_data <- openxlsx::readWorkbook("scgp-subject-metadata.xlsx", 
                                    sheet = 1, startRow = 1, colNames = TRUE)

### Join DNAme disorder and metadata table
cells_epr <- cells_epr %>% 
  # Join merged table with metadata table
  left_join(meta_data, by=c("sample"="subject_id")) %>% 
  mutate(idh_status = ifelse(subtype == "IDHwt", "IDHwt", "IDHmut")) 

# Rename the second "sample" column to "Patient"
names(cells_epr)[grep("sample",names(cells_epr))[2]] <- "patient"

# Reformat "patient" column to shorthand form
cells_epr$patient <- gsub("-","",gsub("SCGP-","",cells_epr$patient))
                 
# Join DNAme disorder and mutation rate tables
sample_mutRate <- cbind(Patient=gsub("-","",gsub("SCGP-","",rownames(sample_mutRate))),sample_mutRate)
cells_epr_w_mutRate <- cells_epr %>% left_join(sample_mutRate, by=c("patient"="Patient")) 

# Separate hypermutants from nonhypermutants
cells_epr_w_mutRate.hypermutants <- cells_epr_w_mutRate %>% filter(global >= 10)
cells_epr_w_mutRate.nonhypermutants <- cells_epr_w_mutRate %>% filter(global < 10)

# Reshape hypermutant table for plotting
cells_epr_w_mutRate.hypermutants.long_form <- reshape2::melt(cells_epr_w_mutRate.hypermutants,
                                                              id.vars = c("sample_barcode","patient","PDR","library_id","idh_status"),
                                                              measure.vars = colnames(sample_mutRate)[-1])

# Plot mutation rate vs DNAme disorder for nonhypermutants
for (i in 1:(dim(sample_mutRate)[2] - 1L)) {
  figure <- ggplot(cells_epr_w_mutRate.nonhypermutants) + 
    geom_point(aes_string(x = names(cells_epr_w_mutRate.nonhypermutants)[i+58], 
                          y = "PDR", 
                          shape = "idh_status", color = "patient", alpha = 0.3, size = 8)) +
    guides(size = FALSE, alpha = FALSE) + 
    labs(x = "Mutations/Mb", y = "DNAme disorder", 
         title = paste0(names(cells_epr_w_mutRate.nonhypermutants)[i+58]," regions"), shape = "Subtype") +
    scale_color_manual(name='patient', values=c('SM001'='#F8766D',
                                                 'SM002'='#DB8E00',
                                                 'SM004'='#AEA200',
                                                 'SM006'='#64B200',
                                                 'SM008'='#00BD5C',
                                                 'SM011'='#00C1A7',
                                                 'SM012'='#00BADE',
                                                 'SM015'='#00A6FF',
                                                 'SM017'='#B385FF',
                                                 'SM018'='#EF67EB',
                                                 'SM019'='#FF63B6'))
  ggsave(paste0(outDir,"/",names(sample_mutRate)[i+1],"_region_mutation_vs_DNAme_disorder.pdf"), 
         figure, device = "pdf", units = "in", width = 12, height = 6)
}

# Plot mutation rate vs DNAme disorder for hypermutants
hypmut_figure <- ggplot(cells_epr_w_mutRate.hypermutants.long_form) + 
  geom_point(aes(x = value, y = PDR, color = variable, alpha = 0.3, size = 8),shape = 15) +
  guides(size = FALSE, alpha = FALSE) + 
  labs(x = "Mutations/Mb", y = "DNAme disorder", title = "SM011", shape = "Subtype", color = "Genomic\nContext") +
  scale_colour_brewer(palette = "Set1", direction = -1)

### Aggregate DNAme disorder by patient and compare to patient mutation rates ###

# Calculate patient average global DNAme disorder
avg_PDR <- ddply(cells_epr,.(patient,initial_recurrence,idh_status),summarise,PDR = median(PDR))
names(avg_PDR)[grep("PDR",names(avg_PDR))] <- "avg_PDR"

# Calculate standard deviation of patient global DNAme disorder
sd_PDR <- aggregate(cells_epr$PDR,by=list(cells_epr$patient),sd)
names(sd_PDR) <- c("patient","sd_PDR")

# Calculate patient average promoter DNAme disorder
avg_promoter_PDR <- aggregate(cells_epr$promoter_PDR,by=list(cells_epr$patient),median)
names(avg_promoter_PDR) <- c("patient","avg_promoter_PDR")

# Calculate standard deviation of patient average promoter DNAme disorder
sd_promoter_PDR <- aggregate(cells_epr$promoter_PDR,by=list(cells_epr$patient),sd)
names(sd_promoter_PDR) <- c("patient","sd_promoter_PDR")

# Calculate patient average gene body DNAme disorder
avg_gene_body_PDR <- aggregate(cells_epr$gene_body_PDR,by=list(cells_epr$patient),median)
names(avg_gene_body_PDR) <- c("patient","avg_gene_body_PDR")

# Calculate standard deviation of patient average gene body DNAme disorder
sd_gene_body_PDR <- aggregate(cells_epr$gene_body_PDR,by=list(cells_epr$patient),sd)
names(sd_gene_body_PDR) <- c("patient","sd_gene_body_PDR")

# Calculate patient average intergenic DNAme disorder
avg_intergenic_PDR <- aggregate(cells_epr$intergenic_PDR,by=list(cells_epr$patient),median)
names(avg_intergenic_PDR) <- c("patient","avg_intergenic_PDR")

# Calculate standard deviation of patient average intergenic DNAme disorder
sd_intergenic_PDR <- aggregate(cells_epr$intergenic_PDR,by=list(cells_epr$patient),sd)
names(sd_intergenic_PDR) <- c("patient","sd_intergenic_PDR")

# Combine aggregated mutation rate statistics with epimutaiton rate statistics
patient_epr_w_mutRate <- avg_PDR %>% left_join(sd_PDR, by=c("patient"="patient")) %>%
  left_join(avg_promoter_PDR, by=c("patient"="patient")) %>%
  left_join(sd_promoter_PDR, by=c("patient"="patient")) %>%
  left_join(avg_gene_body_PDR, by=c("patient"="patient")) %>%
  left_join(sd_gene_body_PDR, by=c("patient"="patient")) %>%
  left_join(avg_intergenic_PDR, by=c("patient"="patient")) %>%
  left_join(sd_intergenic_PDR, by=c("patient"="patient")) %>%
  left_join(sample_mutRate, by=c("patient"="Patient"))

# Write combined table to file
write.table(patient_epr_w_mutRate, file = paste0(outDir,"/","Samples_WGS_context-specific_DNAme_disorder_and_mutation_rates.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Separate hypermutants from non-hypermutants
patient_epr_w_mutRate.hypermutants <- patient_epr_w_mutRate %>% filter(global >= 10)
patient_epr_w_mutRate.nonhypermutants <- patient_epr_w_mutRate %>% filter(global < 10)

# Reshape aggregated mutation/DNAme disorder tables for plotting
patient_epr_w_mutRate.long_form <- reshape2::melt(patient_epr_w_mutRate,
                                                  id.vars = c("patient","initial_recurrence","idh_status","avg_PDR","sd_PDR",
                                                              "avg_promoter_PDR","sd_promoter_PDR","avg_gene_body_PDR",
                                                              "sd_gene_body_PDR","avg_intergenic_PDR","sd_intergenic_PDR"),
                                                  measure.vars = colnames(sample_mutRate)[-1])

patient_epr_w_mutRate.hypermutants.long_form <- reshape2::melt(patient_epr_w_mutRate.hypermutants,
                                                               id.vars = c("patient","initial_recurrence","idh_status","avg_PDR","sd_PDR",
                                                                           "avg_promoter_PDR","sd_promoter_PDR","avg_gene_body_PDR",
                                                                           "sd_gene_body_PDR","avg_intergenic_PDR","sd_intergenic_PDR"),
                                                               measure.vars = colnames(sample_mutRate)[-1])

patient_epr_w_mutRate.nonhypermutants.long_form <- reshape2::melt(patient_epr_w_mutRate.nonhypermutants,
                                                                  id.vars = c("patient","initial_recurrence","idh_status","avg_PDR","sd_PDR",
                                                                              "avg_promoter_PDR","sd_promoter_PDR","avg_gene_body_PDR",
                                                                              "sd_gene_body_PDR","avg_intergenic_PDR","sd_intergenic_PDR"),
                                                                  measure.vars = colnames(sample_mutRate)[-1])

# Reorder sample levels so that subtypes are grouped together,
# and primary samples are ordered before recurrent
patient_order <- c("SM004","SM001","SM015","SM019","SM002","SM008","SM006","SM012","SM017","SM018","SM011")

patient_epr_w_mutRate.long_form$patient <- factor(patient_epr_w_mutRate.long_form$patient, levels = patient_order)
patient_epr_w_mutRate.hypermutants.long_form$patient <- factor(patient_epr_w_mutRate.hypermutants.long_form$patient, levels = patient_order)
patient_epr_w_mutRate.nonhypermutants.long_form$patient <- factor(patient_epr_w_mutRate.nonhypermutants.long_form$patient, levels = patient_order)

# Labels for coloring sample names by primary-recurrent status (PRS)
patient_PRS_color <- c("#CA2F66","#2FB3CA","#CA2F66","#CA2F66","#2FB3CA","#2FB3CA","#CA2F66","#CA2F66","#CA2F66","#CA2F66","#2FB3CA")

# Fix order for subtype as IDHmut_codel -> IDHmut_noncodel -> IDHwt
subtype_order <- c("IDHmut","IDHwt")

patient_epr_w_mutRate.long_form$idh_status <- factor(patient_epr_w_mutRate.long_form$idh_status, levels = subtype_order)
patient_epr_w_mutRate.hypermutants.long_form$idh_status <- factor(patient_epr_w_mutRate.hypermutants.long_form$idh_status, levels = subtype_order)
patient_epr_w_mutRate.nonhypermutants.long_form$idh_status <- factor(patient_epr_w_mutRate.nonhypermutants.long_form$idh_status, levels = subtype_order)


# Plot mutation rates for non-hypermutants
# ggplot(patient_epr_w_mutRate.long_form) + geom(point(aes))
mutRate.nonhypermutants_figure <- ggplot(patient_epr_w_mutRate.long_form) + 
  geom_point(aes(x = patient, 
                 y = value, 
                 shape = variable, color = variable, alpha = 0.75, size = avg_PDR)) +
  #facet_wrap(~subtype, nrow = 1) +
  facet_grid(~idh_status, scales = "free", space="free_x") +
  #geom_tile(aes(fill="Primary"), shape = NA) +
  #geom_tile(aes(fill="Recurrent"), shape = NA) +
  ylim(0,10) +
  scale_size(range = c(3, 12)) +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  scale_shape_manual(name="Genomic\nContext", values=c(18,16,16,16,16,16,16,16,16),
                     labels=c("Global","Exon","Gene Body","Intergenic","Intron","Enhancer","Promoter","CGI","DNaseI")) +
  guides(shape = guide_legend("Genomic\nContext",override.aes=list(color = rev(brewer.pal(9,"Set1")),size = 8))) +
  #scale_fill_manual(name='Time Point', values=c('Primary'='#CA2F66', 'Recurrent'='#2FB3CA')) +
  #guides(fill = guide_legend("Time Point",override.aes=list(shape=15, size = 6, color = c("#CA2F66","#2FB3CA")))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.title = element_text(face = "bold"), strip.background =element_rect(fill="white")) +
  guides(alpha = FALSE, color = FALSE) + 
  labs(x = "Patient", y = "Mutations/Mb", color = "Genomic\nContext",
       size = "Average\nDNAme\nDisorder", shape = "Subtype", title = NULL) 

ggsave(paste0(outDir,"/","Samples_bulk.nonhypermutants-region_mutation_rate.pdf"), 
       mutRate.nonhypermutants_figure, device = "pdf", units = "in", width = 12, height = 8)

# Plot mutation rates for hypermutants
# ggplot(patient_epr_w_mutRate.long_form) + geom(point(aes))
mutRate.hypermutants_figure <- ggplot(patient_epr_w_mutRate.long_form) + 
  geom_point(aes(x = Patient, 
                 y = value, 
                 shape = variable, color = variable, alpha = 0.75, size = avg_PDR)) +
  #facet_wrap(~subtype, nrow = 1) +
  facet_grid(~idh_status, scales = "free", space="free_x") +
  #geom_tile(aes(fill="Primary"), shape = NA) +
  #geom_tile(aes(fill="Recurrent"), shape = NA) +
  ylim(115,410) +
  scale_size(range = c(3, 12)) +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  scale_shape_manual(name="Genomic\nContext", values=c(18,16,16,16,16,16,16,16,16),
                     labels=c("Global","Exon","Gene Body","Intergenic","Intron","Enhancer","Promoter","CGI","DNaseI")) +
  guides(shape = guide_legend("Genomic\nContext",override.aes=list(color = rev(brewer.pal(9,"Set1")),size = 8))) +
  #scale_fill_manual(name='Time Point', values=c('Primary'='#CA2F66', 'Recurrent'='#2FB3CA')) +
  #guides(fill = guide_legend("Time Point",override.aes=list(shape=15, size = 6, color = c("#CA2F66","#2FB3CA")))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.title = element_text(face = "bold"), strip.background =element_rect(fill="white")) +
  guides(alpha = FALSE, color = FALSE) + 
  labs(x = "Patient", y = "Mutations/Mb", color = "Genomic\nContext",
       size = "Average\nDNAme\nDisorder", shape = "Subtype", title = NULL) 

ggsave(paste0(outDir,"/","Samples_bulk.hypermutants-region_mutation_rate.pdf"), 
       mutRate.hypermutants_figure, device = "pdf", units = "in", width = 12, height = 8)

# # Calculate spearman's correlation coefficient for DNAme disorder vs. mutation rates
# epr_mutRate_corr <- cor.test(patient_epr_w_mutRate$avg_PDR,patient_epr_w_mutRate$global, method = "spearman")

# Create labels for facet wrap plots
context_names <- list(
  'global'="Global",
  'exon'="Exon",
  'gene_body'="Gene body",
  'intergenic'="Intergenic",
  'intron'="Intron",
  'enhancer'="Enhancer",
  'promoter'="Promoter",
  'cgi'="CGI",
  'dnaseI'="DNaseI"
)

context_labeller <- function(variable,value){
  return(context_names[value])
}

### Plot DNAme disorder vs mutation rates for nonhypermutants in select genomic contexts
# Subset epr_w_mutRate table to select genomic contexts
patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts <- 
  patient_epr_w_mutRate.nonhypermutants.long_form[which(patient_epr_w_mutRate.nonhypermutants.long_form$variable %in% c("global","promoter","gene_body","intergenic")),]

# Reformat table
patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable <- gsub("global","Global",patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable)
patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable <- gsub("promoter","Promoter",patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable)
patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable <- gsub("gene_body","Gene body",patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable)
patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable <- gsub("intergenic","Intergenic",patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable)

patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable <- factor(patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable, 
                                                                                   levels = c("Global","Promoter","Gene body","Intergenic"))

# epr_mutRate.nonhypermutants_figures <- vector(mode = "list", length = length(unique(patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts$variable)))

epr_mutRate.nonhypermutants_figures <- ggplot(patient_epr_w_mutRate.nonhypermutants.long_form.select_contexts, 
                                              aes(x = value, y = avg_PDR)) + 
  geom_point() +
  geom_smooth(aes(x = value, y = avg_PDR), method = "lm") +
  ggpubr::stat_cor(method="spearman",label.sep = "\n") + 
  facet_grid(~variable, scales = "free", space="free_y") +
  guides(alpha = FALSE) +
  ylim(0.2,0.5) +
  theme_bw(base_size = 12) + 
  theme(axis.line = element_line(color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.5,"lines")) +
  labs(y = "Context-specific DNAme disorder", x = "Context-specific mutations / Mb")

pdf("Samples_bulk.nonhypermutants-region_mutation_vs_DNAme_disorder.pdf", 
    width = 6, height = 4, useDingbats = FALSE, bg="transparent")
epr_mutRate.nonhypermutants_figures
dev.off()
###
#####


