##################################
# Assess the stress-associated changes in DNAme disorder
# Updated: 2021.02.18
# Author: Kevin J.
###################################

## Working directory for this analysis.
mybasedir = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/"
setwd(mybasedir)


###################################
# Necessary packages:
library(tidyverse)
library(parallel)
library(data.table)
library(ggpubr)
library(RColorBrewer)
library(EnvStats)
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

## Define factors.
oxygen_levels = c("Normoxia", "Hypoxia")

### Download the 3-day hypoxia data.
files_20x = list.files("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/hypoxia_epimutation", include.dirs = TRUE, full.names = T, pattern = "_pe_sorted_ge20x_global_and_context_set1and2-specific_weighted_mean_PDR.txt", recursive=T)

# Each data frame is a single row.
rrbs_20x = mclapply(files_20x, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub("_pe_sorted_ge20x_global_and_context_set1and2-specific_weighted_mean_PDR.txt", "", basename(f)))
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples.
rrbs_20x_hypoxia = data.table::rbindlist(rrbs_20x)

rrbs_20x_hypoxia = rrbs_20x_hypoxia %>% 
  mutate(cell_line = gsub("-", "", substr(sample, 6, 12)),
         replicate = substr(sample, 14, 16),
         o2_conc = substr(sample, 24, 25),
         o2_level = ifelse(o2_conc==21, "Normoxia", "Hypoxia"),
         time = substr(sample, 22, 23),
         timepoint = ifelse(time==72, "3 days", "9 days"))



### Plot the 20X data by feature  + timepoint
rrbs_20x_hypoxia_pdr_plot <- rrbs_20x_hypoxia %>% 
  dplyr::select(sample_id:timepoint, global_PDR = global, ends_with('PDR')) %>% 
  gather(feature_name, epimutation_burden, ends_with('PDR')) %>% 
  mutate(feature_name = gsub("_PDR", "", feature_name)) %>% 
  filter(feature_name%in%c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer"))

feature_levels <- c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer")
rrbs_20x_hypoxia_pdr_plot$feature_name <-  factor(rrbs_20x_hypoxia_pdr_plot$feature_name, levels = feature_levels)
rrbs_20x_hypoxia_pdr_plot$o2_level <-  factor(rrbs_20x_hypoxia_pdr_plot$o2_level, levels = oxygen_levels)

ggplot(rrbs_20x_hypoxia_pdr_plot, aes(x=feature_name, y = epimutation_burden, fill=o2_level)) +
  geom_boxplot() +
  labs(y = "Epimutation Burden", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format") 


rrbs_20x_hypoxia_pdr_plot_norm <- rrbs_20x_hypoxia_pdr_plot %>% 
  group_by(cell_line, timepoint, feature_name) %>%
  mutate_each(funs(./median(.[o2_level == "Normoxia"])), epimutation_burden)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_hypoxia_ge20x_median_normoxia_normalize.pdf", width = 11, height = 6, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm, aes(x=feature_name, y = epimutation_burden, fill=o2_level)) +
  geom_boxplot() +
  labs(y = "Relative local DNAme disorder", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format", size = 1.5) 
  #stat_compare_means(label = "p.signif") 
dev.off()

rrbs_20x_hypoxia_pdr_plot_norm_filt = rrbs_20x_hypoxia_pdr_plot_norm %>% 
  filter(feature_name%in%c("cgi", "promoter", "intergenic", "cgi_shore"))
feature_levels_filt <- c("cgi", "promoter", "intergenic", "cgi_shore")
rrbs_20x_hypoxia_pdr_plot_norm_filt$feature_name <-  factor(rrbs_20x_hypoxia_pdr_plot_norm_filt$feature_name, levels = feature_levels_filt)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/presentations/DFCI/relative_DNAme_hypoxia_ge20x_median_normoxia_normalize-filtered.pdf", width = 8, height = 6, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt, aes(x=feature_name, y = epimutation_burden, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative local DNAme disorder", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  stat_compare_means(label = "p.format", size = 1.5)
  #stat_compare_means(label = "p.signif") 
dev.off()

## Filter to 3-day and present the oxygen concentrations for select features
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d <- rrbs_20x_hypoxia_pdr_plot_norm_filt %>% 
  filter(timepoint == "3 days") %>% 
  mutate(o2_conc = paste(o2_conc, "%", sep=""))

o2_levels <- c("21%", "02%", "01%")
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc <-  factor(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc, levels = o2_levels)
my_comparisons <- list(c("21%", "01%"), c("21%", "02%"), c("01%", "02%"))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig4/SuppFig8b-hypoxia-dosage.pdf", width = 6, height = 4, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d, aes(x=feature_name, y = epimutation_burden, fill=o2_conc)) +
  geom_boxplot() +
  labs(y = "Relative DNAme disorder\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                              "02%" = "#fdae61",  
                              "01%" = "#d7191c")) + 
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~.) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_compare_means(label = "p.format" , method = "kruskal.test") 
dev.off()

my_comparisons <- list( c("21%", "01%"), c("21%", "02%"), c("02%", "01%"))
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d, aes(x=o2_conc, y = epimutation_burden, fill=o2_conc)) +
  geom_boxplot() +
  labs(y = "Relative epimutation burden\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                             "02%" = "#fdae61",  
                             "01%" = "#d7191c")) + 
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~ feature_name) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format", size = 2 , comparisons = my_comparisons) +
  ggtitle("Hypoxia dosage effect - 3 day")

### END ###