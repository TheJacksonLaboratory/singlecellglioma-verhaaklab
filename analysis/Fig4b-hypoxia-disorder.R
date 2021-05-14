##################################
# Assess the stress-associated changes in DNAme disorder
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

pdf(file = "github/results/Fig4/Fig4b-relative-disorder-hypoxia.pdf", width = 6, height = 4, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt, aes(x=feature_name, y = epimutation_burden, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative local DNAme disorder", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(cell_line ~ timepoint) +
  stat_compare_means(label = "p.format", size = 1.5)
dev.off()

pdf(file = "github/results/Fig4/Fig4b-relative-disorder-hypoxia-no-legend.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt, aes(x=feature_name, y = epimutation_burden, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative DNAme disorder\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  guides(fill=FALSE) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
  facet_grid(cell_line ~ timepoint) +
  stat_compare_means(label = "p.format", size = 1.5)
dev.off()

## Filter to 3-day and present the oxygen concentrations for select features
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d <- rrbs_20x_hypoxia_pdr_plot_norm_filt %>% 
  filter(timepoint == "3 days") %>% 
  mutate(o2_conc = paste(o2_conc, "%", sep=""))

o2_levels <- c("21%", "02%", "01%")
rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc <-  factor(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d$o2_conc, levels = o2_levels)
my_comparisons <- list(c("21%", "01%"), c("21%", "02%"), c("01%", "02%"))

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_hypoxia_ge20x_median_normoxia_normalize_filtered_3d_dosage.pdf", width = 7, height = 5, useDingbats = FALSE)
ggplot(rrbs_20x_hypoxia_pdr_plot_norm_filt_3d, aes(x=feature_name, y = epimutation_burden, fill=o2_conc)) +
  geom_boxplot() +
  labs(y = "Relative epimutation burden\n(median normoxia normalized values)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                              "02%" = "#fdae61",  
                              "01%" = "#d7191c")) + 
  geom_hline(yintercept = 1, linetype = 2) +
  plot_theme +
  facet_grid(cell_line ~ .) +
  theme(panel.spacing.x = unit(1.5, "lines"),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format", size = 2 , method = "kruskal.test") +
  ggtitle("Hypoxia dosage effect - 3 day")
dev.off()


#### Calculate the median normoxia value
hypoxia_3d_comb <- rrbs_10x_hypoxia %>% 
  filter(timepoint=="3 days") %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) 

hypoxia_9d_comb <- rrbs_10x_hypoxia %>% 
  filter(timepoint=="9 days") %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) 


## Normalize by median normoxia:
hypoxia_3d_hf2354 <- rrbs_10x_hypoxia %>% 
  filter(cell_line == "HF2354"  & timepoint=="3 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.320) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)

hypoxia_9d_hf2354 <- rrbs_10x_hypoxia %>% 
  filter(cell_line == "HF2354" & timepoint=="9 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.323) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)


hypoxia_3d_hf3016 <- rrbs_10x_hypoxia %>% 
  filter(cell_line == "HF3016" & timepoint=="3 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.307) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)

hypoxia_9d_hf3016 <- rrbs_10x_hypoxia %>% 
  filter(cell_line == "HF3016" & timepoint=="9 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.299) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)



hypoxia_cgi <- hypoxia_3d_hf2354 %>% 
  bind_rows(hypoxia_9d_hf2354, hypoxia_3d_hf3016, hypoxia_9d_hf3016)

hypoxia_cgi$o2_level <-  factor(hypoxia_cgi$o2_level, levels = oxygen_levels)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_cgi_hypoxia_ge10x_median_normoxia_normalize.pdf", width =7, height = 5, useDingbats = FALSE)
ggplot(hypoxia_cgi, aes(x=o2_level, y = relative_cgi_epimutation, fill=o2_level)) +
  geom_boxplot() +
  labs(y = "Relative CGI Epimutation Burden", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  stat_compare_means(method="wilcox.test") +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  stat_n_text()
dev.off()  



##########################################
#### 20X
##########################################
#### Calculate the median normoxia value
hypoxia_3d_comb <- rrbs_20x_hypoxia %>% 
  filter(timepoint=="3 days") %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) 

hypoxia_9d_comb <- rrbs_20x_hypoxia %>% 
  filter(timepoint=="9 days") %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) 


## Normalize by median normoxia:
hypoxia_3d_hf2354 <- rrbs_12x_hypoxia %>% 
  filter(cell_line == "HF2354"  & timepoint=="3 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.320) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)

hypoxia_9d_hf2354 <- rrbs_20x_hypoxia %>% 
  filter(cell_line == "HF2354" & timepoint=="9 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.323) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)


hypoxia_3d_hf3016 <- rrbs_20x_hypoxia %>% 
  filter(cell_line == "HF3016" & timepoint=="3 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.307) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)

hypoxia_9d_hf3016 <- rrbs_20x_hypoxia %>% 
  filter(cell_line == "HF3016" & timepoint=="9 days") %>% 
  mutate(relative_cgi_epimutation = cgi_PDR/0.299) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_cgi_epimutation)



hypoxia_cgi <- hypoxia_3d_hf2354 %>% 
  bind_rows(hypoxia_9d_hf2354, hypoxia_3d_hf3016, hypoxia_9d_hf3016)

hypoxia_cgi$o2_level <-  factor(hypoxia_cgi$o2_level, levels = oxygen_levels)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_cgi_hypoxia_ge10x_median_normoxia_normalize.pdf", width =7, height = 5, useDingbats = FALSE)
ggplot(hypoxia_cgi, aes(x=o2_level, y = relative_cgi_epimutation, fill=o2_level)) +
  geom_boxplot() +
  labs(y = "Relative CGI Epimutation Burden", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  stat_compare_means(method="wilcox.test") +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  stat_n_text()
dev.off()  









###########################################
### Normalize with the specific experiment
###########################################
## Normalize by replicate-specific normoxia by subtraction:
hypoxia_3d_replicate_sum <- rrbs_3d_hypoxia %>% 
  group_by(cell_line, replicate, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) %>% 
  filter(o2_level=="Normoxia") 

hypoxia_3d_replicate <- rrbs_3d_hypoxia %>% 
  mutate(relative_cgi_epimutation = NA) %>% 
  dplyr::select(sample, cell_line, replicate, timepoint, o2_conc, o2_level, cgi_PDR, relative_cgi_epimutation)

normalized_values <- c(hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "01D"]/0.2844893,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "02D"]/0.3074093,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "03D"]/0.3397975,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "04D"]/0.3378776,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "01D"]/0.3093080,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "02D"]/0.3064193,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "03D"]/0.3089814,
                       hypoxia_3d_replicate$cgi_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "04D"]/0.3070904)

hypoxia_3d_replicate$relative_cgi_epimutation <- as.numeric(normalized_values)


## Normalize by replicate-specific normoxia by subtraction:
hypoxia_9d_replicate_sum <- rrbs_9d_hypoxia %>% 
  group_by(cell_line, replicate, o2_level) %>% 
  summarise(med_norm_cgi = median(cgi_PDR)) %>% 
  filter(o2_level=="Normoxia") 

hypoxia_9d_replicate <- rrbs_9d_hypoxia %>% 
  mutate(relative_cgi_epimutation = NA) %>% 
  select(sample, cell_line, replicate, timepoint, o2_conc, o2_level, cgi_PDR, relative_cgi_epimutation) 

normalized_values_9d <- c(hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "05D"]/0.3162809,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "06D"]/0.3208067,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "07D"]/0.3218552,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "08D"]/0.3277412,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "09D"]/0.3196707,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "10D"]/0.3271328,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "05D"]/0.2965759,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "06D"]/0.3016780,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "07D"]/0.2960708,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "08D"]/0.3079938,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "09D"]/0.2953324,
                          hypoxia_9d_replicate$cgi_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "10D"]/0.2933232)

hypoxia_9d_replicate$relative_cgi_epimutation <- as.numeric(normalized_values_9d)

hypoxia_batch_normalize <- bind_rows(hypoxia_3d_replicate, hypoxia_9d_replicate)

oxygen_conc= c("21%", "02%", "01%")
hypoxia_batch_normalize$o2_level <-  factor(hypoxia_batch_normalize$o2_level, levels = oxygen_levels)
hypoxia_batch_normalize$o2_conc <-  factor(hypoxia_batch_normalize$o2_conc, levels = oxygen_conc)


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_cgi_hypoxia_ge10x_replicate_normoxia_normalize.pdf", width =7, height = 5, useDingbats = FALSE)
ggplot(hypoxia_batch_normalize, aes(x=o2_level, y = relative_cgi_epimutation, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative CGI Epimutation Burden (3 day)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  stat_compare_means(method="wilcox.test") +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"), 
        axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_n_text()
dev.off()

#######################################
###### Promoter-specific calculations
#######################################
#### Calculate the median normoxia value
hypoxia_3d_comb <- rrbs_3d_hypoxia %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_promoter = median(promoter_PDR)) 

hypoxia_9d_comb <- rrbs_9d_hypoxia %>% 
  group_by(cell_line, o2_level) %>% 
  summarise(med_norm_promoter = median(promoter_PDR)) 


## Normalize by median normoxia:
hypoxia_3d_hf2354 <- rrbs_3d_hypoxia %>% 
  filter(cell_line == "HF2354") %>% 
  mutate(relative_promoter_epimutation = promoter_PDR/0.293) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_promoter_epimutation)

hypoxia_9d_hf2354 <- rrbs_9d_hypoxia %>% 
  filter(cell_line == "HF2354") %>% 
  mutate(relative_promoter_epimutation = promoter_PDR/0.295) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_promoter_epimutation)


hypoxia_3d_hf3016 <- rrbs_3d_hypoxia %>% 
  filter(cell_line == "HF3016") %>% 
  mutate(relative_promoter_epimutation = promoter_PDR/0.286) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_promoter_epimutation)

hypoxia_9d_hf3016 <- rrbs_9d_hypoxia %>% 
  filter(cell_line == "HF3016") %>% 
  mutate(relative_promoter_epimutation = promoter_PDR/0.277) %>% 
  dplyr::select(sample, cell_line, timepoint, o2_level, relative_promoter_epimutation)



hypoxia_promoter <- hypoxia_3d_hf2354 %>% 
  bind_rows(hypoxia_9d_hf2354, hypoxia_3d_hf3016, hypoxia_9d_hf3016)

hypoxia_promoter$o2_level <-  factor(hypoxia_promoter$o2_level, levels = oxygen_levels)

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_promoter_hypoxia_ge10x_median_normoxia_normalize.pdf", width =7, height = 5, useDingbats = FALSE)
ggplot(hypoxia_promoter, aes(x=o2_level, y = relative_promoter_epimutation, fill=o2_level)) +
  geom_boxplot() +
  labs(y = "Relative Promoter Epimutation Burden", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  stat_compare_means(method="wilcox.test") +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  stat_n_text()
dev.off()  



###########################################
### Normalize with the specific experiment
###########################################
## Normalize by replicate-specific normoxia by subtraction:
hypoxia_3d_replicate_sum <- rrbs_3d_hypoxia %>% 
  group_by(cell_line, replicate, o2_level) %>% 
  summarise(med_norm_promoter = median(promoter_PDR)) %>% 
  filter(o2_level=="Normoxia") 

hypoxia_3d_replicate <- rrbs_3d_hypoxia %>% 
  mutate(relative_promoter_epimutation = NA) %>% 
  dplyr::select(sample, cell_line, replicate, timepoint, o2_conc, o2_level, promoter_PDR, relative_promoter_epimutation)

normalized_values <- c(hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "01D"]/0.2614826,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "02D"]/0.2788018,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "03D"]/0.3073845,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF2354" & hypoxia_3d_replicate$replicate == "04D"]/0.3066636,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "01D"]/0.2865903,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "02D"]/0.2827244,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "03D"]/0.2871650,
                       hypoxia_3d_replicate$promoter_PDR[hypoxia_3d_replicate$cell_line=="HF3016" & hypoxia_3d_replicate$replicate == "04D"]/0.2863562)

hypoxia_3d_replicate$relative_promoter_epimutation <- as.numeric(normalized_values)


## Normalize by replicate-specific normoxia by subtraction:
hypoxia_9d_replicate_sum <- rrbs_9d_hypoxia %>% 
  group_by(cell_line, replicate, o2_level) %>% 
  summarise(med_norm_cgi = median(promoter_PDR)) %>% 
  filter(o2_level=="Normoxia") 

hypoxia_9d_replicate <- rrbs_9d_hypoxia %>% 
  mutate(relative_promoter_epimutation = NA) %>% 
  select(sample, cell_line, replicate, timepoint, o2_conc, o2_level, promoter_PDR, relative_promoter_epimutation) 

normalized_values_9d <- c(hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "05D"]/0.2979073,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "06D"]/0.2901325,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "07D"]/0.2923010,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "08D"]/0.3001898,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "09D"]/0.2902329,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF2354" & hypoxia_9d_replicate$replicate == "10D"]/0.2994855,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "05D"]/0.2784737,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "06D"]/0.2816359,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "07D"]/0.2752614,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "08D"]/0.2878591,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "09D"]/0.2760957,
                          hypoxia_9d_replicate$promoter_PDR[hypoxia_9d_replicate$cell_line=="HF3016" & hypoxia_9d_replicate$replicate == "10D"]/0.2748070)

hypoxia_9d_replicate$relative_promoter_epimutation <- as.numeric(normalized_values_9d)

hypoxia_batch_normalize <- bind_rows(hypoxia_3d_replicate, hypoxia_9d_replicate)

oxygen_conc= c("21%", "02%", "01%")
hypoxia_batch_normalize$o2_level <-  factor(hypoxia_batch_normalize$o2_level, levels = oxygen_levels)
hypoxia_batch_normalize$o2_conc <-  factor(hypoxia_batch_normalize$o2_conc, levels = oxygen_conc)


pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/gsc/relative_epimutation_promoter_hypoxia_ge10x_replicate_normoxia_normalize.pdf", width =7, height = 5, useDingbats = FALSE)
ggplot(hypoxia_batch_normalize, aes(x=o2_level, y = relative_promoter_epimutation, fill=o2_level)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Relative Promoter Epimutation Burden (3 day)", x = "", fill = "Oxygen conc.") +
  scale_fill_manual(values=c("Normoxia" = "#abd9e9",
                             "Hypoxia" = "#d7191c")) +
  plot_theme +
  stat_compare_means(method="wilcox.test") +
  facet_grid(cell_line ~ timepoint) +
  theme(panel.spacing.x = unit(1.5, "lines"), 
        axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_n_text()
dev.off()


