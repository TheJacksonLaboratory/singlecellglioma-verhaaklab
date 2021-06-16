##############################################
# Analysis of epimutation at different replication times
# Updated: 2020.06.01
# Author: Kevin J.
###############################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###############################################
## Load the necessary packages.
library(tidyverse)
library(ggpubr)
library(openxlsx)
###############################################
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

## Analysis of DNAme disorder at different replication times
## First, generate the context-specific PDR plot for reference.
### Load the subject-level metadata.
meta_data = read.csv("data/clinical_metadata.csv", sep = ",", header = TRUE)

## Load in the PDR metric calculated across different replication timings.
reptime_disorder = read.csv("data/analysis_scRRBS_replication_timing_DNAme_disorder.csv", sep = ",", header = TRUE)

reptime_disorder_annot <- reptime_disorder %>% 
  left_join(meta_data, by="case_barcode") %>% 
  mutate(idh_status = ifelse(idh_codel_subtype=="IDHwt", "IDHwt", "IDHmut"))

prom_pdr_rt = reptime_disorder_annot %>% 
  dplyr::select(cell_barcode, case_barcode, idh_status, promoter_rt1_PDR, promoter_rt2_PDR, promoter_rt3_PDR, promoter_rt4_PDR) %>% 
  gather(reptime, pdr, c(promoter_rt1_PDR, promoter_rt2_PDR, promoter_rt3_PDR, promoter_rt4_PDR)) %>% 
  mutate(reptime = recode(reptime, "promoter_rt1_PDR" = "very early", "promoter_rt2_PDR" = "early", "promoter_rt3_PDR" = "late", 
                          "promoter_rt4_PDR" = "very late"))

subtype_order <- c("IDHmut", "IDHwt")
prom_pdr_rt <- prom_pdr_rt %>% mutate(reptime = factor(reptime, levels = c("very early", "early", "late", "very late")))
prom_pdr_rt <- prom_pdr_rt %>% mutate(idh_status = factor(idh_status, levels = subtype_order))


## Separate into IDHmut and IDHwt and perform correlation analysis.
prom_pdr_rt_wt = prom_pdr_rt %>% filter(idh_status=="IDHwt")
prom_pdr_rt_mut = prom_pdr_rt %>% filter(idh_status!="IDHwt")
cor.test(prom_pdr_rt$pdr, as.numeric(prom_pdr_rt$reptime), method="kendall")
cor.test(prom_pdr_rt_wt$pdr, as.numeric(prom_pdr_rt_wt$reptime), method="kendall")
cor.test(prom_pdr_rt_mut$pdr, as.numeric(prom_pdr_rt_mut$reptime), method="kendall")



## Perform same analyses for gene body.
genebody_pdr_rt = reptime_disorder_annot %>% 
  select(cell_barcode, case_barcode, idh_status, genebody_rt1_PDR, genebody_rt2_PDR, genebody_rt3_PDR, genebody_rt4_PDR) %>% 
  gather(reptime, pdr, c(genebody_rt1_PDR, genebody_rt2_PDR, genebody_rt3_PDR, genebody_rt4_PDR)) %>% 
  mutate(reptime = recode(reptime, "genebody_rt1_PDR" = "very early", "genebody_rt2_PDR" = "early", "genebody_rt3_PDR" = "late", 
                          "genebody_rt4_PDR" = "very late")) 
genebody_pdr_rt <- genebody_pdr_rt %>% mutate(reptime = factor(reptime, levels = c("very early", "early", "late", "very late")))
genebody_pdr_rt <- genebody_pdr_rt %>% mutate(idh_status = factor(idh_status, levels = subtype_order))


## Separate into IDHmut and IDHwt and perform correlation analysis.
genebody_pdr_rt_wt = genebody_pdr_rt %>% filter(idh_status=="IDHwt")
genebody_pdr_rt_mut = genebody_pdr_rt %>% filter(idh_status!="IDHwt")
cor.test(genebody_pdr_rt$pdr, as.numeric(genebody_pdr_rt$reptime), method="kendall")
cor.test(genebody_pdr_rt_wt$pdr, as.numeric(genebody_pdr_rt_wt$reptime), method="kendall")
cor.test(genebody_pdr_rt_mut$pdr, as.numeric(genebody_pdr_rt_mut$reptime), method="kendall")

#### Combine the two metrics:
genebody_pdr_rt$region <- "Gene body"
prom_pdr_rt$region <- "Promoter"

comb_pdr_rt <- bind_rows(genebody_pdr_rt, prom_pdr_rt)
comb_pdr_rt$region <- factor(comb_pdr_rt$region, levels = c("Promoter", "Gene body"))


pdf(file = "results/Fig5/Fig5b-replication-timing.pdf", width = 7, height = 4)
ggplot(comb_pdr_rt, aes(x= idh_status, y = pdr, fill=reptime)) + 
  geom_violin() +
  geom_boxplot(width=0.9, color="black", alpha=0.1, outlier.shape = NA) +
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5)) +
  scale_fill_manual(name='Replication timing', values=c('very early'='#f1eef6', 'early'='#bdc9e1', "late"= "#74a9cf", "very late"="#0570b0")) +
  stat_compare_means(method = "kruskal.test") +
  labs(y = "DNAme disorder (PDR)", x="") +
  facet_grid(~ region, scales = "free_x", space = "free")
dev.off()

### END ###