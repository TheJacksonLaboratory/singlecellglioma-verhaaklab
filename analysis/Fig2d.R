---
title: "Figure 1h"
---
  
  ### Preprocessing
  
  #### Load data
  
  ```{r}
  library(topGO)
  library(tidyverse)
  library(openxlsx)
  library(ggrepel)
  library(corrplot)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(egg)
  library(reshape2)
  library(ggpubr)


  # Load clinical metadata (available on Synapse)
  clinical_meta <- read.csv("clinical_metadata.csv")
  
  # Combine IDHmut subtypes
  clinical_meta <- clinical_meta %>% mutate(idh_status = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut"))
  
  # Load scRRBS sequencing QC (available on Synapse)
  scRRBS_qc <- read.csv("analysis_scRRBS_sequencing_qc.csv")
  
  # Load metrics for CpG density overlapping TFBS motifs
  cpg_density_info <- read.delim("scRRBS_epiallele_CpG_density_summary.txt")
  
  # Load TFBS motif DNAme disorder (available on Synapse)
  tf_pdr <- read.csv("analysis_scRRBS_individual_TFBS_motif_DNAme_disorder.csv")
  
  # Add IDH status to DNAme disorder table
  tf_pdr <- tf_pdr %>% left_join(clinical_meta[c("case_barcode","idh_status")])
  
  # Remove (var.X) suffix from TF names
  tf_pdr$TF <- gsub("\\(var.2\\)|\\(var.3\\)","",tf_pdr$TF)
  
  # Set the case_barcode order as in manuscript figures.
  case_order <- c("SM004", "SM001", "SM015", "SM019", "SM002", "SM008", "SM006", "SM012", "SM017", "SM018", "SM011")
  tf_pdr$case_barcode <- factor(tf_pdr$case_barcode, levels = case_order)
  ```
  
  #### Process data for figure generation
  
  ```{r}
  # Filter data to remove TF observations for cells with < 100 reads for that TF,
  # TFs identified as lacking CpGs at binding motifs, and any TFs lacking observations in both subtypes
  tfs_to_drop <-  c("HLTF", "EN1", "FOXD1", "FOXG1", "FOXN3", "FOXO3", "MEIS3", "MSX1", "PRRX1", "Shox2", "VAX2")
  
  tf_pdr_filt <- tf_pdr %>% 
    filter(epiallele_count >= 100, !TF%in%tfs_to_drop)
  
  # Calculate the average PDR per TF and per subject
  tf_pdr_filt_pdr <- tf_pdr_filt %>% 
    group_by(TF, case_barcode) %>% 
    summarise(avg_PDR = mean(PDR)) %>% 
    ungroup()
  
  # Create a sort for each TF based on PDR
  sort_df <- tf_pdr_filt_pdr %>%
    group_by(TF) %>% 
    summarise(DNAme_disorder = median(avg_PDR)) %>% 
    arrange(desc(DNAme_disorder))  %>% 
    ungroup()
  
  tf_order <- unique(sort_df$TF)
  
  tf_pdr_filt_pdr <- tf_pdr_filt_pdr %>% mutate(TF = factor(TF, levels = tf_order))
  
  # Create a rank based on average DNAme disorder:
  tf_pdr_filt_ranked <- tf_pdr_filt_pdr %>% 
    group_by(case_barcode) %>% 
    mutate(rank = dense_rank(desc(avg_PDR))) %>% 
    ungroup() %>% 
    left_join(clinical_meta[c("case_barcode","idh_status")])
  
  # Create an average rank per IDHmut status
  tf_pdr_filt_ranked_avg <- tf_pdr_filt_ranked %>% 
    group_by(idh_status, TF) %>% 
    summarise(avg_subtype_rank = mean(rank, na.rm = T)) %>% 
    spread(idh_status, avg_subtype_rank) %>% 
    mutate(rank_diff = IDHwt-IDHmut)
  
  # Calculate difference between subtype-averaged TFBS motif DNAme disorder
  tf_pdr_filt_subtype <- tf_pdr_filt_pdr %>% 
    left_join(clinical_meta[c("case_barcode","idh_status")]) %>%
    group_by(TF, idh_status) %>% 
    summarise(TF_subtype_PDR = median(avg_PDR)) %>% 
    ungroup() %>% 
    spread(idh_status, TF_subtype_PDR) %>% 
    mutate(PDR_diff = IDHwt-IDHmut,
           PDR_change = ifelse(PDR_diff>=0, "+", "-"))
  
  # Remove any non-human TFs (inclusion of lowercase letters).
  tf_pdr_filt_subtype <- tf_pdr_filt_subtype %>% 
    filter(!grepl("[a-z]", TF))
  
  # Create sortable object for plotting by increasing TFBS epimutation.
  sort_df <- tf_pdr_filt_subtype %>%
    arrange(desc(IDHwt)) 
  
  tf_order <- rev(unique(sort_df$TF))
  
  tf_pdr_filt_subtype <- tf_pdr_filt_subtype %>% mutate(TF = factor(TF, levels = tf_order))
  
  tf_pdr_filt_subtype$placeholder <- c("IDHmut", rep("IDHwt", dim(tf_pdr_filt_subtype)[1]-1))
  
  # Write out results to be used downstream
  tf_pdr_filt_subtype_out = tf_pdr_filt_subtype %>% 
    distinct()
  
  # Filter cpg_density_info to TFBS motifs considered for analysis
  cpg_density_info_temp <- cpg_density_info
  cpg_density_info_temp$TF <- gsub("\\(var.2\\)|\\(var.3\\)","",cpg_density_info$TF)
  
  cpg_density_info_filt <- left_join(data.frame(TF=tf_pdr_filt_subtype_out$TF),
                                     cpg_density_info_temp, by = c("TF"="TF"))
  
  # Combine DNAme disorder with CpG density
  tf_pdr_filt_subtype_density <- tf_pdr_filt_subtype_out %>% 
    left_join(cpg_density_info_temp, by = c("TF"="TF"))
  rm(cpg_density_info_temp)
  
  # Reformat combined table for plotting 
  tf_pdr_filt_subtype_density.long <- melt(tf_pdr_filt_subtype_density[,c(1:3,7,10:11)],
                                           id.vars=names(tf_pdr_filt_subtype_density)[c(1,7,10:11)],
                                           melt.vars=c("IDHmut","IDHwt"))
  names(tf_pdr_filt_subtype_density.long)[c(5:6)] <- c("subtype","DNAme_disorder")
  
  # # Filter out the extra TF row that was created due to multiple SELEX information:
  # tf_pdr_filt_subtype <- tf_pdr_filt_subtype %>% 
  #   distinct()
  # dup_to_drop <- which(tf_pdr_filt_subtype$TF =="YY1" & tf_pdr_filt_subtype$`methyl-SELEX.call`=="MethylMinus and MethylPlus")
  # tf_pdr_filt_plot <- tf_pdr_filt_subtype[-dup_to_drop,]
  
  ### Add GSC essential + fitness genes
  # Load the CRISPR screen with Bayes Factor value (BAGEL algorithm) per cell line and per gene.
  stab1 <- readWorkbook("~/Documents/Papers Supplemental Information/Genome-Wide CRISPR-Cas9 Screens Expose Genetic Vulnerabilities and Mechanisms of Temozolomide Sensitivity in Glioblastoma Stem Cells/maccleod-supptable1.xlsx", sheet = 1, startRow = 1, rowNames=TRUE, colNames = TRUE)
  stab2 <- readWorkbook("~/Documents/Papers Supplemental Information/Genome-Wide CRISPR-Cas9 Screens Expose Genetic Vulnerabilities and Mechanisms of Temozolomide Sensitivity in Glioblastoma Stem Cells/maccleod-supptable2.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
  
  # Include all profiled cells.
  stab1_gsc <- stab1[,1:12]
  stab1_gsc[] <- sapply(stab1_gsc, function(x) ifelse(x>=3, 1, 0))
  
  # Find which genes are significant in more than 6 GSCs.
  ess_genes <- names(which(rowSums(stab1_gsc)>6))
  
  # Determine the genes described as GSC fitness genes.
  gsc_fitness_gene_df <- stab2 %>% 
    filter(`Z-score.GBM.vs..other.ex-NSC` > 1.65)
  
  gsc_fitness_genes <- gsc_fitness_gene_df$X1
  
  # Annotate DNAme disorder table with GSC fitness genes
  tf_pdr_filt_plot <- tf_pdr_filt_subtype %>% 
    mutate(match_TF = gsub("\\.var\\.[0-9]\\.", "", TF),
           crispr_screened = ifelse(match_TF%in%rownames(stab1), 1, 0),
           crispr_essential = ifelse(!match_TF%in%rownames(stab1), "Not Avail.", ifelse(match_TF%in%ess_genes, "Essential", "Not sig.")),
           crispr_fitness = ifelse(!match_TF%in%rownames(stab1), "Not Avail.", ifelse(match_TF%in%gsc_fitness_genes, "Fitness", "Not sig.")))
  ###
  
  # Annotate DNAme disorder table with CpG density
  tf_pdr_filt_plot.update <- left_join(tf_pdr_filt_plot,cpg_density_info_filt[,-c(3:5)], by = c("match_TF"="TF"))
  
  # Extract CpG density metrics and reformat for plotting as an annotation track
  tf_pdr_filt_plot.update.density_long <- tf_pdr_filt_plot.update[,c(1,12)]
  
  names(tf_pdr_filt_plot.update.density_long) <- gsub("nearby_motif.","\u00B1",names(tf_pdr_filt_plot.update.density_long))
  
  names(tf_pdr_filt_plot.update.density_long) <- gsub("_AUC","",names(tf_pdr_filt_plot.update.density_long))
  
  tf_pdr_filt_plot.update.density_long <- melt(tf_pdr_filt_plot.update.density_long,
                                               id.vars = "TF",
                                               melt.vars = names(tf_pdr_filt_plot.update.density_long)[-1])
  
  names(tf_pdr_filt_plot.update.density_long)[c(2:3)] <- c("distance_from_motif_center","AUC")
  
  
  # Create and reformat a table to plot correlation of CpG density AUC vs. TFBS motif DNAme disorder
  tf_pdr_filt_subtype_density.corr <- melt(tf_pdr_filt_plot.update[,c(1:3,12)],
                                           id.vars = names(tf_pdr_filt_plot.update)[c(1,12)],
                                           melt.vars = names(tf_pdr_filt_plot.update)[c(2:3)])
  names(tf_pdr_filt_subtype_density.corr)[c(3:4)] <- c("subtype","DNAme_disorder")
  
  tf_pdr_filt_subtype_density.corr <- melt(tf_pdr_filt_subtype_density.corr,
                                           id.vars = names(tf_pdr_filt_subtype_density.corr)[-2],
                                           melt.vars = names(tf_pdr_filt_subtype_density.corr)[2])
  names(tf_pdr_filt_subtype_density.corr)[c(4,5)] <- c("distance_from_motif_center","CpG_density_AUC")
  
  tf_pdr_filt_subtype_density.corr$distance_from_motif_center <- 
    gsub("nearby_motif.","",tf_pdr_filt_subtype_density.corr$distance_from_motif_center)
  
  tf_pdr_filt_subtype_density.corr$distance_from_motif_center <- 
    gsub("_AUC","",tf_pdr_filt_subtype_density.corr$distance_from_motif_center)
  
  tf_pdr_filt_subtype_density.corr$distance_from_motif_center <- 
    paste0("Distance from motif center: \u00B1",tf_pdr_filt_subtype_density.corr$distance_from_motif_center)
  
  tf_pdr_filt_subtype_density.corr$distance_from_motif_center <- 
    gsub("2r","2*motif radius",tf_pdr_filt_subtype_density.corr$distance_from_motif_center)
  
  tf_pdr_filt_subtype_density.corr$distance_from_motif_center <- 
    factor(tf_pdr_filt_subtype_density.corr$distance_from_motif_center,
           levels = unique(tf_pdr_filt_subtype_density.corr$distance_from_motif_center))
  ```
  
  ### Common plotting elements
  
  ```{r}
  # Set plot themes
  plot_theme <- theme_bw(base_size = 12) + 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size=12, angle = 90),
          axis.line = element_line(color = "black"),
          strip.background = element_rect(fill="white"),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  null_x <- theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())
  
  null_y        <- theme(axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank())
  
  bottom_x      <- theme(axis.text.x=element_blank())
  
  null_facet    <- theme(strip.background = element_blank(),
                         strip.text.x = element_blank())
  
  top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines"))
  
  middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
  
  bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))
  
  null_legend   <- theme(legend.position = 'none')
  
  # Extract legend function
  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  
  # Plot function
  testPlot <- function(gg, grid = TRUE) {
    if(grid)
      gg + plot_theme  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    else
      gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  ```
  
  
  
  ### Figure generation
  
  ```{r}
  # Generate plot of subtype mean TFBS motif DNAme disorder, ranked by DNAme disorder
  # with a dotted line connecting IDHmut and IDHwt, and TFs associated with GSC fitness genes annotated
  gg_epimut_tf_fitness <- ggplot(tf_pdr_filt_plot, aes(x=TF)) +
    geom_hline(yintercept = 0.4, alpha=0.8, linetype=2, color="gray70") +
    geom_point(aes(y=IDHwt), color = "#7fbf7b") +
    geom_linerange(aes(ymin = IDHwt, ymax = IDHmut, color = PDR_change), linetype = 2) +
    geom_point(aes(y=IDHmut), color = "#af8dc3") +
    scale_color_manual(values = c("+" = "gray50", "-" = "#fb8072")) +
    labs(y = "Subtype mean TFBS motif DNAme disorder", color="IDHwt-IDHmut diff.") + 
    plot_theme + 
    null_x +
    geom_text_repel(data=filter(tf_pdr_filt_plot, crispr_fitness=="Fitness"), aes(label=TF, y=IDHwt),
                    nudge_y = 0.05,
                    direction = "x",
                    angle = 90,
                    vjust = 0,
                    segment.size = 0.3, 
                    color="#2D2926FF")
  
  # Generate plot comparing DNAme disorders 
  epimut_density_corr_plot <- ggplot(tf_pdr_filt_subtype_density.corr,
                                     aes(x=CpG_density_AUC, y=DNAme_disorder, color = subtype, group = subtype)) +
    facet_wrap(~distance_from_motif_center, nrow = 2, ncol = 2) +
    geom_point() +
    stat_cor(method = "spearman") +
    scale_color_manual(values = c('IDHwt'='#7fbf7b','IDHmut'='#af8dc3')) +
    theme_bw() + 
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background =element_rect(fill="white")) +
    guides(color = guide_legend("Glioma\nSubtype",override.aes=list(shape=15))) +
    labs(x = "Area under CpG density curve", y = "DNAme disorder", color = "Glioma subtype")
  
  # Generate annotation for TFBS motif length
  gg_motif_length <- tf_pdr_filt_plot.update %>% 
    ggplot(aes(x=TF)) +
    geom_tile(aes(fill = motif_width, y = 1)) +
    scale_fill_gradient2(low = "#edf8b1",
                         mid = "#7fcdbb",
                         high = "#2c7fb8",
                         midpoint = median(tf_pdr_filt_plot.update$motif_width)) +
    labs(y="", fill = "TFBS motif\nlength") +
    plot_theme
  
  # Generate annotation for CpG density for +-2r motif distance
  gg_motif_density <-
    tf_pdr_filt_plot.update.density_long %>% 
    ggplot(aes(x=TF)) +
    geom_tile(aes(fill = AUC, y = distance_from_motif_center)) +
    scale_fill_gradient(low = "#ffffd4",
                        high = "#993404") +
    labs(y="", fill = "Area under CpG\ndensity curve") +
    plot_theme
  
  # Generate figure legends to plot separately
  gleg1 = g_legend(gg_epimut_tf_fitness)
  gleg2 = g_legend(epimut_density_corr_plot + plot_theme)
  gleg3 = g_legend(gg_motif_length)
  gleg4 = g_legend(gg_motif_density)
  
  # Generate plots with plot themes configured for combining plots
  g1 = ggplotGrob(gg_epimut_tf_fitness + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
  g2 = ggplotGrob(gg_motif_length + plot_theme + null_legend + null_y + null_x + null_facet + theme(plot.margin = unit(c(0,1,0,1), "lines")))  %>% gtable_frame()
  g3 = ggplotGrob(gg_motif_density + plot_theme + null_legend + null_y + null_x + null_facet + bottom_margin)  %>% gtable_frame()
  
  g = gtable_rbind(g1, g2, g3)
  gleg = gtable_rbind(gleg2, gleg1, gleg3, gleg4)
  
  # Adjust relative height of figure panels
  panels = g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- unit(c(4,0.2,0.2), "null")
  
  # Plot combined figure panels
  grid.newpage()
  grid.draw(g)
  
  # Plot combined legend panels
  plot(gleg)
  ```
  
  
  
