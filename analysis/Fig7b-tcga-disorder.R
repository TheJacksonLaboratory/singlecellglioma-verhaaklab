##############################################
# Investigate the DNAme disorder metric in TCGA bulk samples.
# Updated: 2020.05.27
# Author: Kevin J.
##################################################

# Working directory for this analysis in the SCGP project. 
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)


########################################
# Necessary packages:
library(tidyverse)
library(RColorBrewer)
library(openxlsx)
library(DBI)
library(survminer)
library(survival)
library(ggpubr)
library(EnvStats)
########################################
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

## The TCGA beta-values were processed using minfi and adjusted with BMIQ.
tcga_betas_epimut = readRDS("data/tcga_betas_bmiq_epimut_genes.rds")

## DNA methylation disorder metric was calculated as the sum of the intermediate value
## over the the "high" epimutation genes.
eITH = readRDS("data/tcga_dna_methylation_instability.rds")

## In the TCGA data, there are some samples from multiple cases. Filter just the primary tumors.
eITH = eITH %>% 
  mutate(case_barcode = substr(idat, 1, 12)) %>% 
  filter(grepl("-01A-", idat))
## Is there only one sample per case? Yes.
dim(eITH)[1]==n_distinct(eITH$case_barcode)

## Load in the clinical data from the Supplement of Ceccarelli et al Cell 2016.
tcga_clinical = readWorkbook("data/tcga-glioma-2016-clinical.xlsx", sheet = 1, startRow = 2, colNames = TRUE)

# Join together with the eITH metric with the clinical data.
tcga_clinical_filt = tcga_clinical %>% 
  inner_join(eITH, by=c("Case"="case_barcode")) 

## How does the DNA methylation disorder index vary by IDHmut status and DNAm subtype.
ggplot(tcga_clinical_filt, aes(x=IDH.status, y=eITH)) + geom_boxplot()
ggplot(tcga_clinical_filt, aes(x=`IDH-specific.DNA.Methylation.Cluster`, y=eITH)) + geom_boxplot()

## Filter to those with a bulk DNAme classification.
tcga_clinical_filt_plot = tcga_clinical_filt %>% 
  filter(!is.na(`Supervised.DNA.Methylation.Cluster`))

## Define subtype order:
subtype_order <- c("Codel", "G-CIMP-high", "G-CIMP-low", "Classic-like", "Mesenchymal-like", "LGm6-GBM", "PA-like")
tcga_clinical_filt_plot <- tcga_clinical_filt_plot %>% mutate(Supervised.DNA.Methylation.Cluster = factor(Supervised.DNA.Methylation.Cluster, levels = subtype_order))

ggplot(tcga_clinical_filt_plot, aes(x=`Supervised.DNA.Methylation.Cluster`, y=eITH, fill=Supervised.DNA.Methylation.Cluster)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Codel" = "#af8dc3",
                           "G-CIMP-high" = "#af8dc3",
                           "G-CIMP-low" = "#af8dc3",
                           "Classic-like" = "#7fbf7b",
                           "Mesenchymal-like" = "#7fbf7b",
                           "LGm6-GBM" = "#7fbf7b",
                           "PA-like" = "#7fbf7b")) +
  labs(x="TCGA Pan-glioma DNA methylation classification", y="DNA methylation instability (high epimutation genes)", fill="Glioma subtype") +
  plot_theme + 
  stat_compare_means(method="kruskal.test") + 
  stat_n_text()

## Adding in ENCODE data types to represent pure normal cell states per reviewer request.
encode_eITH <- read.table(file="data/encode-eITH-metrics.txt", sep="\t", header=T, stringsAsFactors = F)
subtype_order_norm <- c("Normal (ENCODE)","Codel", "G-CIMP-high", "G-CIMP-low", "Classic-like", "Mesenchymal-like", "LGm6-GBM", "PA-like")

## Add this information to the object to be plotted.
tcga_clinical_filt_plot_encode <- tcga_clinical_filt_plot %>% 
  select(Case, Supervised.DNA.Methylation.Cluster, eITH) %>% 
  bind_rows(encode_eITH) %>% 
  mutate(Supervised.DNA.Methylation.Cluster = factor(Supervised.DNA.Methylation.Cluster, levels = subtype_order_norm))


pdf(file = "results/Fig7/Fig7b-tcga-subtype-encode-high-disorder-no-legend.pdf", height = 4, width = 6, useDingbats = FALSE, bg='transparent')
ggplot(tcga_clinical_filt_plot_encode, aes(x=`Supervised.DNA.Methylation.Cluster`, y=eITH, fill=Supervised.DNA.Methylation.Cluster)) + 
  geom_violin() +
  geom_boxplot(width=0.25, color="black", alpha=0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("Normal (ENCODE)" = "gray70",
                               "Codel" = "#af8dc3",
                               "G-CIMP-high" = "#af8dc3",
                               "G-CIMP-low" = "#af8dc3",
                               "Classic-like" = "#7fbf7b",
                               "Mesenchymal-like" = "#7fbf7b",
                               "LGm6-GBM" = "#7fbf7b",
                               "PA-like" = "#7fbf7b")) +
  labs(x="TCGA Pan-glioma DNAme classification", y="Bulk DNA disorder\n(high DNAme disorder genes in scRRBS)", fill="Sample\nclassification") +
  plot_theme + 
  stat_compare_means(method="kruskal.test") + 
  stat_n_text() +
  guides(fill=FALSE)
dev.off()

pdf(file = "results/Fig7/Fig7b-legend.pdf", height = 4, width = 6, useDingbats = FALSE, bg='transparent')
ggplot(tcga_clinical_filt_plot_encode, aes(x=`Supervised.DNA.Methylation.Cluster`, y=eITH, fill=Supervised.DNA.Methylation.Cluster)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Normal (ENCODE)" = "gray70",
                               "Codel" = "#af8dc3",
                               "G-CIMP-high" = "#af8dc3",
                               "G-CIMP-low" = "#af8dc3",
                               "Classic-like" = "#7fbf7b",
                               "Mesenchymal-like" = "#7fbf7b",
                               "LGm6-GBM" = "#7fbf7b",
                               "PA-like" = "#7fbf7b")) +
  labs(x="TCGA Pan-glioma DNAme classification", y="Bulk DNA disorder\n(high DNAme disorder genes in scRRBS)", fill="Sample\nclassification") +
  plot_theme +
  theme(legend.position="bottom") 
dev.off()


#### END ####