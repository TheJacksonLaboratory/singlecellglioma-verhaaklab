---
title: "Figure 1h"
---
  
  ### Preprocessing
  
  #### Load data
  
  ```{r}
  library(genomation)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)

  # Load sample binned feature DNAme disorder and methylation tables
  # NOTE: current object format is a list of data frames, one per patient, 
  # containing binned CGI feature (upstream_CGI_shore, CGI, downstream_CGI_shore) values per cell
  samples_pdr <- readRDS("~/Documents/scgp/synapse_tables/recoded_id/files/analysis_scRRBS_binned_CGI_and_shores_DNAme_disorder.Rds")
  samples_meth <- readRDS("~/Documents/scgp/synapse_tables/recoded_id/files/analysis_scRRBS_binned_CGI_and_shores_methylation.Rds")
  
  # Load clinical metadata
  clinical_meta <- read.csv("~/Documents/scgp/synapse_tables/recoded_id/tables/clinical_metadata.csv")
  
  # Combine IDHmut subtypes
  clinical_meta <- clinical_meta %>% mutate(idh_status = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut"))
  ```
  
  #### Average metrics across cells to compare sample level aggregation of epimutation burden and methylation status across CGI sets
  
  ```{r}
  # Initialize output lists
  samples_pdr.cell_avg <- vector(mode = "list", length = length(samples_pdr))
  samples_meth.cell_avg <- vector(mode = "list", length = length(samples_pdr))
  
  for (i in 1:length(samples_pdr)) {
    # For each feature, average epimutation burden across cells
    samples_pdr.cell_avg[[i]] <- data.frame(samples_pdr[[i]][,c(1:3)],
                                            names(samples_pdr)[i],
                                            apply(samples_pdr[[i]][,-c(1:3)],1,mean, na.rm = TRUE),
                                            stringsAsFactors = FALSE,check.names = FALSE)
    names(samples_pdr.cell_avg[[i]]) <- c("anno_name","rel_pos","scaled_rel_pos","sample","mean_pdr")
    
    # For convenience, remove the end of rel_pos and set as a numeric vector
    samples_pdr.cell_avg[[i]]$rel_pos <- as.numeric(unlist(lapply(strsplit(samples_pdr.cell_avg[[i]]$rel_pos,":"),'[',1)))
    
    # For convenience, remove the end of scaled_rel_pos and set as a numeric vector
    samples_pdr.cell_avg[[i]]$scaled_rel_pos <- as.numeric(unlist(lapply(strsplit(samples_pdr.cell_avg[[i]]$scaled_rel_pos,":"),'[',1)))
    
    
    # For each feature, average methylation beta value across cells
    samples_meth.cell_avg[[i]] <- data.frame(samples_meth[[i]][,c(1:3)],
                                             names(samples_meth)[i],
                                             apply(samples_meth[[i]][,-c(1:3)],1,mean, na.rm = TRUE),
                                             stringsAsFactors = FALSE,check.names = FALSE)
    names(samples_meth.cell_avg[[i]]) <- c("anno_name","rel_pos","scaled_rel_pos","sample","mean_beta_value")
    
    # For convenience, remove the end of rel_pos and set as a numeric vector
    samples_meth.cell_avg[[i]]$rel_pos <- as.numeric(unlist(lapply(strsplit(samples_meth.cell_avg[[i]]$rel_pos,":"),'[',1)))
    
    # For convenience, remove the end of scaled_res_pos and set as a numeric vector
    samples_meth.cell_avg[[i]]$scaled_rel_pos <- as.numeric(unlist(lapply(strsplit(samples_meth.cell_avg[[i]]$scaled_rel_pos,":"),'[',1)))
  }
  
  # Join and reformat individual sample tables
  samples_pdr.cell_avg.df <- bind_rows(samples_pdr.cell_avg)
  samples_pdr.cell_avg.df$anno_name <- factor(samples_pdr.cell_avg.df$anno_name, levels = c("upstream_CGI_shore","CGI","downstream_CGI_shore"))
  samples_pdr.cell_avg.df$sample <- gsub("-","",substr(samples_pdr.cell_avg.df$sample,6,11))
  samples_pdr.cell_avg.df <- samples_pdr.cell_avg.df %>% left_join(clinical_meta[,c(1,16)], by = c("sample"="case_barcode"))
  
  samples_meth.cell_avg.df <- bind_rows(samples_meth.cell_avg)
  samples_meth.cell_avg.df$anno_name <- factor(samples_meth.cell_avg.df$anno_name, levels = c("upstream_CGI_shore","CGI","downstream_CGI_shore"))
  samples_meth.cell_avg.df$sample <- gsub("-","",substr(samples_meth.cell_avg.df$sample,6,11))
  samples_meth.cell_avg.df <- samples_meth.cell_avg.df %>% left_join(clinical_meta[,c(1,16)], by = c("sample"="case_barcode"))
  ```
  
  ### Figure generation
  
  ```{r}
  plot_theme  <- theme_bw(base_size = 12) + theme(axis.title = element_text(12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background =element_rect(fill="white"),
                                                  panel.grid = element_blank())
  
  # Plot a fitted curve of average methylation beta value for each subtype
  merged_feature_sample_meth_fig.subtype <- ggplot(samples_meth.cell_avg.df) +
    geom_smooth(aes(x=rel_pos, y=mean_beta_value, group = sample, color = idh_status, alpha = 0.5), 
                formula = "y ~ x", method = "loess", se = FALSE) +
    geom_vline(xintercept = -1000, linetype = "dashed") +
    geom_vline(xintercept = 1000, linetype = "dashed") +
    coord_cartesian(xlim = c(-2300,2300)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.line = element_line(color = "black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_color_manual(name='Subject averaged\nsingle cell profiles', values=c('IDHmut'='#AF8DC3','IDHwt'='#7FBF7B'),
                       labels = c("IDHmut (n = 6 subjects)","IDHwt (n = 5 subjects)")) +
    guides(alpha = FALSE) +
    labs(x = NULL, y = "Mean cell DNA methylation") 
  
  # Plot a fitted curve of average epimutation burden for each subtype
  merged_feature_sample_pdr_fig.subtype <- ggplot(samples_pdr.cell_avg.df) +
    geom_smooth(aes(x=rel_pos, y=mean_pdr, group = sample, color = idh_status, alpha = 0.5), 
                formula = "y ~ x", method = "loess", se = FALSE) +
    geom_vline(xintercept = -1000, linetype = "dashed") +
    geom_vline(xintercept = 1000, linetype = "dashed") +
    coord_cartesian(xlim = c(-2300,2300)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.line = element_line(color = "black"),
          strip.background = element_rect(fill="white"),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_color_manual(name='Glioma subtype', values=c('IDHmut'='#AF8DC3','IDHwt'='#7FBF7B')) +
    guides(color = FALSE, alpha = FALSE) +
    labs(x = NULL, y = "Mean cell DNA\nmethylation disorder (PDR)") 
  
  # Combine plots
  combined_subtype_fig <- egg::ggarrange(merged_feature_sample_meth_fig.subtype, merged_feature_sample_pdr_fig.subtype, nrow = 2, ncol = 1, heights = c(1,1))
  
  pdf(paste0("~/Documents/scgp/epimutation/results/rerun-reformatted_deduplicated-final_recalculated","/",
             "Samples-passQC_single_tumor_cells_CGI_and_shores_merged_feature_meth_and_pdr_curves-subtype_annotation.pdf"), 
      width = 5.5, height = 5.6, useDingbats = FALSE, bg="transparent")
  combined_subtype_fig
  dev.off()
  ```