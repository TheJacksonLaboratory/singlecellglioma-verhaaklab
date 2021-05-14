##################################
# Assess the irradiation-associated changes in DNAme disorder
# Updated: 2021.05.13
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

## Define factors.
oxygen_levels = c("Normoxia", "Hypoxia")

## Hypoxia and control normoxia data.
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

rrbs_20x_hypoxia_full =read.delim("/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/hypoxia_experiment_RRBS-ge20x_coverage/hypoxia_experiment_RRBS-9day-ge20x_coverage-global_and_context_set1and2-specific_weighted_mean_PDR.txt", as.is=T, header=T)


rrbs_20x_hypoxia = rrbs_20x_hypoxia %>% 
  mutate(cell_line = gsub("-", "", substr(sample, 6, 12)),
         replicate = substr(sample, 14, 16),
         o2_conc = substr(sample, 24, 25),
         o2_level = ifelse(o2_conc==21, "Normoxia", "Hypoxia"),
         time = substr(sample, 22, 23),
         timepoint = ifelse(time==72, "3 days", "9 days"))



## RT (10 Gy) exposed cells.
rt_files_20x = list.files("/Users/johnsk/mnt/verhaak-lab/scgp/results/epimutation/RT_experiment_RRBS-ge20x_coverage", include.dirs = TRUE, full.names = T, pattern = "_pe_sorted_ge20x_global_and_context_set1and2-specific_weighted_mean_PDR.txt", recursive=T)

# Each data frame is a single row.
rrbs_rt_20x = mclapply(rt_files_20x, function(f){
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
rrbs_20x_rt = data.table::rbindlist(rrbs_rt_20x)

rrbs_20x_rt = rrbs_20x_rt %>% 
  mutate(cell_line = gsub("-", "", substr(sample, 6, 12)),
         replicate = substr(sample, 14, 16),
         o2_conc = substr(sample, 24, 25),
         o2_level = ifelse(o2_conc=="RT", "Normoxia", NA),
         time = substr(sample, 22, 23),
         timepoint = ifelse(time==72, "3 days", "9 days"))

## Combine the normoxia 9-day data with RT data.
rt_9d_rrbs_compare <- rrbs_20x_hypoxia %>% 
  filter(o2_level=="Normoxia", timepoint=="9 days") %>% 
  bind_rows(rrbs_20x_rt) %>% 
  mutate(treatment = ifelse(o2_conc=="RT", "10Gy", "0Gy"))

### Plot the 20X data by feature  + timepoint
rrbs_20x_rt_pdr_plot <- rt_9d_rrbs_compare %>% 
  dplyr::select(sample_id:treatment, global_PDR = global, ends_with('PDR')) %>% 
  gather(feature_name, epimutation_burden, ends_with('PDR')) %>% 
  mutate(feature_name = gsub("_PDR", "", feature_name)) %>% 
  filter(feature_name%in%c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer"))

feature_levels <- c("global", "alu_repeat", "intergenic", "cgi_shore", "h1hesc_ezh2", "dnaseI", "h1hesc_ctcf_1", "cgi", "promoter", "enhancer")
rrbs_20x_rt_pdr_plot$feature_name <-  factor(rrbs_20x_rt_pdr_plot$feature_name, levels = feature_levels)
treatment_levels = c("0Gy", "10Gy")
rrbs_20x_rt_pdr_plot$treatment <-  factor(rrbs_20x_rt_pdr_plot$treatment, levels = treatment_levels)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_rt_ge20x_complete.pdf", width = 11, height = 6, useDingbats = FALSE)
ggplot(rrbs_20x_rt_pdr_plot, aes(x=feature_name, y = epimutation_burden, fill=treatment)) +
  geom_boxplot() +
  labs(y = "Epimutation Burden", x = "", fill = "Treatment\n(dose)") +
  scale_fill_manual(values=c("0Gy" = "#abd9e9",
                             "10Gy" = "#0dba86")) +
  plot_theme +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) #+
  #stat_compare_means(label = "p.format") 
dev.off()

rrbs_20x_rt_pdr_plot_norm <- rrbs_20x_rt_pdr_plot %>% 
  group_by(cell_line, timepoint, feature_name) %>%
  mutate_each(funs(./median(.[treatment == "0Gy"])), epimutation_burden)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_rt_ge20x_median_normoxia_normalize.pdf", width = 11, height = 6, useDingbats = FALSE)
ggplot(rrbs_20x_rt_pdr_plot_norm, aes(x=feature_name, y = epimutation_burden, fill=treatment)) +
  geom_boxplot() +
  labs(y = "Relative epimutation burden\n(median 0Gy normalized values)", x = "", fill = "Irradiation\n(dose)") +
  scale_fill_manual(values=c("0Gy" = "#abd9e9",
                             "10Gy" = "#0dba86")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format", size = 2) 
dev.off()

rrbs_20x_rt_pdr_plot_norm_filt = rrbs_20x_rt_pdr_plot_norm %>% 
  filter(feature_name%in%c("cgi", "promoter", "intergenic", "cgi_shore"))
feature_levels_filt <- c("cgi", "promoter", "intergenic",  "cgi_shore")
rrbs_20x_rt_pdr_plot_norm_filt$feature_name <-  factor(rrbs_20x_rt_pdr_plot_norm_filt$feature_name, levels = feature_levels_filt)


ggobj= ggplot(rrbs_20x_rt_pdr_plot_norm_filt, aes(x=feature_name, y = epimutation_burden, fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative DNAme disorder\n(median 0 Gy normalized values)", x = "", fill = "Irradiation\n(dose)") +
  scale_fill_manual(values=c("0Gy" = "#abd9e9",
                             "10Gy" = "#0dba86")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  theme(legend.position="bottom") +
  facet_grid(cell_line ~ timepoint) +
  ylim(0.9, 1.25) +
  stat_compare_means(label = "p.format", size = 1.5)

layer_scales(ggobj)$y$range$range

pdf(file = "github/results/Fig4/Fig4c-irradiation-disorder.pdf", width = 3, height = 4, useDingbats = FALSE, bg="transparent")
ggobj
dev.off()

pdf(file = "github/results/Fig4/Fig4c-irradiation-disorder-no-legend.pdf", width = 3, height = 5, useDingbats = FALSE, bg="transparent")
ggobj + guides(fill=FALSE)
dev.off()




### END ###

