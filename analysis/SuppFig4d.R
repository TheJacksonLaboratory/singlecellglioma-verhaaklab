# - Generate mutational signatures from vcf data

require(VariantAnnotation)
require(MutationalPatterns)



### Input arguments ###
# Input directory for variant calls
variantDir <- "~/Documents/scgp/mutect2/post_processing"

# Reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Patient metadata 
meta_data <- openxlsx::readWorkbook("/Users/anderk/Documents/SCGP/metadata/scgp-subject-metadata.xlsx", 
                                    sheet = 1, startRow = 1, colNames = TRUE)

# COSMIC mutational signature version (either "v2" or "v3")
cosmic_ver <- "v3"

# Output directory
outDir <- setwd("~/Documents/scgp/mutational_patterns")
#######################


# # Identify vcf files for each sample
# variantFiles <- list.files(path = variantDir, 
#                            pattern = "*-WGS.filtered.normalized.vep.vcf.gz$", full.names = TRUE)
# 
# # Extract first portion of file names to determine sample names
# bulkSamples <- lapply(strsplit(basename(variantFiles),"-"),'[',c(1:3))
# bulkSamples <- unlist(lapply(bulkSamples, function(x) paste(x,collapse = "-")))
# 
# # Reorder metadata to match order of bulkSamples
# meta_data$subject_id <- factor(meta_data$subject_id, levels = bulkSamples)
# meta_data <- meta_data[order(meta_data$subject_id),]
# 
# # Load vcf files as a GRanges list
# vcfs <- read_vcfs_as_granges(variantFiles, bulkSamples, ref_genome)
# 
# 
# ### Identify mutational profile characteristics ###
# # Retrieve the base substitution types and contexts for all positions in the VCF GRanges object
# type_context <- lapply(vcfs, function(x) type_context(x, ref_genome))
# 
# # Count mutation type occurrences for all vcfs
# type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
# # Write to file
# write.table(type_occurrences, file = paste0(outDir,"/","Samples_bulk-mutation_type_occurrences.txt"),
#             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# 
# # Plot the mutation spectrum across samples with distinction between C>T at CpG sites and other sites
# fig.mut_spectrum <- plot_spectrum(type_occurrences, CT = TRUE)
# # Write to file
# ggsave(paste0(outDir,"/","Samples_bulk-consensus_mutation_spectrum.pdf"), 
#        fig.mut_spectrum, device = "pdf", units = "in", width = 12, height = 8)
# 
# # Plot the mutation spectrum across samples with distinction between C>T at CpG sites and other sites
# fig.patient_mut_spectrum <- plot_spectrum(type_occurrences, by = meta_data$subject_id, CT = TRUE, legend = TRUE)
# # Write to file
# ggsave(paste0(outDir,"/","Samples_bulk-patient_mutation_spectrum.pdf"), 
#        fig.patient_mut_spectrum, device = "pdf", units = "in", width = 12, height = 8)
# 
# # Plot the mutation spectrum across subtypes with distinction between C>T at CpG sites and other sites
# fig.subtype_mut_spectrum <- plot_spectrum(type_occurrences, by = meta_data$subtype, CT = TRUE, legend = TRUE)
# # Write to file
# ggsave(paste0(outDir,"/","Samples_bulk-subtype_mutation_spectrum.pdf"), 
#        fig.subtype_mut_spectrum, device = "pdf", units = "in", width = 12, height = 8)
# 
# # Make a 96 trinucleodide mutation count matrix
# mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# 
# # Plot 96 trinucleotide mutation profile
# fig.96_profile <- plot_96_profile(mut_mat[,c(1,7)], condensed = TRUE)
# ###################################################
# 
# 
# ### Record loci with C>T mutations at CpG loci ###
# # Initialize output list of target mutation indices
# transition_loci <- vector(mode = "list", length = length(vcfs))
# names(transition_loci) <- names(vcfs)
# 
# # For each sample, Identify C>T mutations within the context NCG
# for (i in 1:length(transition_loci)) {
#   # Identify C>T mutations
#   CtoTs <- which(type_context[[i]]$types == "C>T")
#   
#   # Identify mutations at CpGs
#   CpGs <- which(type_context[[i]]$context %in% c("ACG","CCG","GCG","TCG"))
#   
#   # Subset patient VCFs to C>T transitions at CpG loci
#   transition_loci[[i]] <- vcfs[[i]][intersect(CtoTs,CpGs)]
#   
#   # Convert subsetted GRanges to BED formatted df
#   transition_df <- data.frame(chr=as.character(seqnames(transition_loci[[i]])),
#                               start=start(transition_loci[[i]]) - 1L,
#                               end=end(transition_loci[[i]]),
#                               name=names(transition_loci[[i]]),
#                               score=rep(0L,length(names(transition_loci[[i]]))),
#                               strand=rep(".",length(names(transition_loci[[i]]))))
#   
#   # Write formatted df to file
#   write.table(transition_df, file = paste0(outDir,"/",names(transition_loci[i]),"-C>T_at_CpG_loci.bed"),
#               sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# }
# ##################################################
# 
# 
# ##### Find optimal contribution of known signatures #####
# 
# ### Obtain COSMIC mutational signatures
# # Load mutational signatures; previously downloaded from the COSMIC website
# if (cosmic_ver == "v3") {
#   cancer_signatures <- as.data.frame(readRDS("~/Documents/scgp/reference/COSMIC_v3_cancer_signatures.rds"))
# } else {
#   cancer_signatures <- as.data.frame(readRDS("~/Documents/scgp/reference/COSMIC_v2_cancer_signatures.rds"))
# }
# 
# 
# # Match the order of the mutation types to MutationalPatterns standard
# new_order <- match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# 
# # Reorder cancer signatures dataframe
# cancer_signatures <- cancer_signatures[as.vector(new_order),]
# 
# # Add trinucletiode changes names as row.names
# row.names(cancer_signatures) <- cancer_signatures$Somatic.Mutation.Type
# 
# # Keep only 96 contributions of the signatures in matrix
# cancer_signatures <- as.matrix(cancer_signatures[,-c(1:3)])
# 
# # Hierarchically cluster the COSMIC signatures based on their similarity with average linkage
# hclust_cosmic <- cluster_signatures(cancer_signatures, method = "average")
# 
# # Store signatures in new order
# cosmic_order <- colnames(cancer_signatures)[hclust_cosmic$order]
# 
# # Plot clustering results
# pdf(file = paste0(outDir,"/","Samples_bulk-pairwise_COSMIC_",cosmic_ver,"_signature_clustering.pdf"),
#     width = 12, height = 8)
# plot(hclust_cosmic)
# dev.off()
# ###
# 
# # Calculate similarity between mutational profiles and COSMIC signatures
# cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, cancer_signatures)
# 
# # Plot heatmap with specified signature order
# pdf(file = paste0(outDir,"/","Samples_bulk-pairwise_COSMIC_",cosmic_ver,"_signature_cosine_similarity.pdf"),
#     width = 12, height = 8)
# plot_cosine_heatmap(cos_sim_samples_signatures,
#                     col_order = cosmic_order,
#                     cluster_rows = TRUE)
# dev.off()
# 
# ### Find optimal contribution of COSMIC signatures to reconstruct 96 mutational profiles
# # Fit mutation matrix to the COSMIC mutational signatures
# fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
# 
# # Plot the optimal contribution of the COSMIC signatures in each sample as a stacked barplot
# # Select signatures with some contribution
# select <- which(rowSums(fit_res$contribution) > 10)
# # Plot contribution barplot
# plot_contribution(fit_res$contribution[select,],
#                   cancer_signatures[,select],
#                   coord_flip = FALSE,
#                   mode = "relative")
# ### 
# 
# ### Extract the contribution of the top n signatures in each sample, 
# ### and collapse the remaining signatures into an other category
# # Number of signatures to extract
# num_sigs <- 5 
# 
# # Top num_sigs signatures in each sample
# top_sigs <- apply(fit_res$contribution,2,function(x) names(sort(x, decreasing = TRUE))[c(1:num_sigs)])
# 
# # Convert mutational signatures output from absolute to relative value
# mut_sigs.relative <- apply(fit_res$contribution,2,function(x) x/sum(x))
# 
# # Extracting top n signatures and collapsing remaining signatures to "Other"
# for (i in 1:ncol(mut_sigs.relative)) {
#   # Identify top signature contributions
#   target_sigs <- which(rownames(mut_sigs.relative) %in% top_sigs[,i])
#   
#   # Create subset of signatures
#   temp <- data.frame(Signature=rownames(mut_sigs.relative)[target_sigs],
#                      Value=as.numeric(mut_sigs.relative[target_sigs,i]),
#                      stringsAsFactors = FALSE)
#   # Add "Other" category to signature subset, with a contribution value 
#   # equal to the sum of remaining signatures
#   temp[num_sigs + 1,1] <- "Other"
#   temp[num_sigs + 1,2] <- sum(mut_sigs.relative[-target_sigs,i])
#   # Add category for patient ID
#   temp <- data.frame(Patient=rep(colnames(mut_sigs.relative)[i],nrow(temp)),Signature=temp$Signature,Value=temp$Value)
#   if (i == 1) {mut_sigs.top_relative <- temp} else {mut_sigs.top_relative <- rbind(mut_sigs.top_relative,temp)}
# }
# 
# # Join mutational signature table with metadata
# mut_sigs.top_relative.annotated <- mut_sigs.top_relative %>% left_join(meta_data, by = c("Patient"="subject_id"))
# 
# # Write table to file
# write.table(mut_sigs.top_relative.annotated, 
#             file = paste0(outDir,"/","Samples_bulk-top_",num_sigs,"_COSMIC_",cosmic_ver,"_signatures.txt"),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

mut_sigs.top_relative.annotated <- read.delim("~/Documents/scgp/mutational_patterns/Samples_bulk-top_5_COSMIC_v3_signatures.txt", check.names = FALSE, stringsAsFactors = FALSE)

# Reformat table for plotting
mut_sigs.top_relative.annotated$Patient <- gsub("-","",gsub("SCGP-","",as.character(mut_sigs.top_relative.annotated$Patient)))

mut_sigs.top_relative.annotated$Patient <- gsub("UC917","SM019",mut_sigs.top_relative.annotated$Patient)

mut_sigs.top_relative.annotated$Patient <- factor(mut_sigs.top_relative.annotated$Patient, 
                                                  levels = c("SM004","SM001","SM015","SM019","SM002","SM008",
                                                             "SM006","SM012","SM017","SM018","SM011"))

mut_sigs.top_relative.annotated$Signature <- gsub("\\."," ",as.character(mut_sigs.top_relative.annotated$Signature))

# Manually order/rename Signature factor levels. Current parameters are for num_sigs=5
if (cosmic_ver == "v3") {
  mut_sigs.top_relative.annotated$Signature <- mapvalues(mut_sigs.top_relative.annotated$Signature,
                                                         from = c("Signature 1","Signature 7d","Signature 8","Signature 9",
                                                                  "Signature 16","Signature 20","Signature 24","Signature 25",
                                                                  "Signature 32","Signature 37","Signature 39","Signature 57",
                                                                  "Other"),
                                                         to = c("Signature 1 (Clock-like)",
                                                                "Signature 7d (Translesion DNA synthesis)",
                                                                "Signature 8 (Unknown)",
                                                                "Signature 9 (Polymerase eta errors)",
                                                                "Signature 16 (Unknown)",
                                                                "Signature 20 (Defective mismatch repair)",
                                                                "Signature 24 (Aflatoxin exposure)",
                                                                "Signature 25 (Unknown)",
                                                                "Signature 32 (Azathioprine exposure)",
                                                                "Signature 37 (Unknown)",
                                                                "Signature 39 (Unknown)",
                                                                "Signature 57 (Transcriptional strand bias)",
                                                                "Other"))
  mut_sigs.top_relative.annotated$Signature <- factor(mut_sigs.top_relative.annotated$Signature,
                                                      levels = c("Signature 1 (Clock-like)",
                                                                 "Signature 7d (Translesion DNA synthesis)",
                                                                 "Signature 8 (Unknown)",
                                                                 "Signature 9 (Polymerase eta errors)",
                                                                 "Signature 16 (Unknown)",
                                                                 "Signature 20 (Defective mismatch repair)",
                                                                 "Signature 24 (Aflatoxin exposure)",
                                                                 "Signature 25 (Unknown)",
                                                                 "Signature 32 (Azathioprine exposure)",
                                                                 "Signature 37 (Unknown)",
                                                                 "Signature 39 (Unknown)",
                                                                 "Signature 57 (Transcriptional strand bias)",
                                                                 "Other"))
} else {
  mut_sigs.top_relative.annotated$Signature <- factor(mut_sigs.top_relative.annotated$Signature,
                                                      levels = c("Signature 1","Signature 2","Signature 3","Signature 5",
                                                                 "Signature 8","Signature 9","Signature 11","Signature 12",
                                                                 "Signature 13","Signature 14","Signature 15","Signature 16",
                                                                 "Signature 19","Signature 21","Signature 22","Signature 25",
                                                                 "Signature 26","Other"))
}


# Set steps for signature color ramp
if (cosmic_ver == "v3") {color_steps = 3} else {color_steps = 4}

# Plot signature contributions
mut_sigs.top_relative.annotated$subtype <- ifelse(mut_sigs.top_relative.annotated$subtype == "IDHwt","IDHwt","IDHmut")

fig.mut_sigs.top_relative <- ggplot(mut_sigs.top_relative.annotated) + 
  geom_bar(aes(x = Patient, y = Value, fill = Signature),position="fill", stat="identity") +
  facet_grid(~subtype, scales = "free", space = "free_x") +
  scale_fill_manual(values=c(colorRamps::primary.colors(n = length(unique(mut_sigs.top_relative.annotated$Signature)), 
                                                        steps = color_steps)[-1],"#999999")) +
  #scale_fill_manual(name='Signature', values=c(c(sig_color[c(1,2,3,5,8,9,11:16,19,21,22,25,26)],"#999999"))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.5, "lines")) +
  labs(x = NULL, y = "Relative Signature Contribution", fill = "COSMIC V3 Mutational Signature\n(potential etiology)")

pdf(paste0(outDir,"/",
           "Samples_bulk-top_",num_sigs,"_COSMIC_",cosmic_ver,"_signatures.pdf"), 
    width = 6, height = 4, useDingbats = FALSE, bg="transparent")
fig.mut_sigs.top_relative
dev.off()



