##################################
# Analyze HF GSC hypoxia biological replicate data
# Updated: 2020.08.13
# Author: Kevin J.
###################################

# Working directory for this analysis.
mybasedir = "/Users/johnsk/github/"
setwd(mybasedir)

###################################
# Necessary packages:
library(tidyverse)
library(openxlsx)
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


## Load data frame from qPCR experiments.
gsc_expression <- read.table("data/gsc-biological-replicates.txt", sep="\t", header=T, stringsAsFactors = F)

## Define factors.
oxygen_levels = c("21%", "2%", "1%")
gsc_expression$oxygen_concentration <- factor(x = gsc_expression$oxygen_concentration, levels = oxygen_levels)
gene_levels = c("SOX2", "POU5F1", "JUN", "HIF2A", "VEGFA")
gsc_expression$gene <- factor(x = gsc_expression$gene, levels = gene_levels)

## Create standard error of the mean function.
sem <- function(x) sd(x) / sqrt(length(x))

## Calculate the summary metrics.
gsc_expression_summary = gsc_expression %>% 
  group_by(cell_line, gene, oxygen_concentration) %>% 
  summarize(relative_ge_avg = mean(relative_gene_expression),
          relative_ge_sem = sem(relative_gene_expression))

## Error bar + bar plot.
ggplot(gsc_expression_summary, aes(x=gene, y=relative_ge_avg, fill= oxygen_concentration)) + 
  geom_bar(position= "dodge", stat="identity", colour="black") +
  geom_errorbar(aes(ymin=relative_ge_avg-relative_ge_sem, 
                    ymax=relative_ge_avg+relative_ge_sem), 
                width=.2, 
                position=position_dodge(.9)) +
  facet_wrap(~ cell_line, scales = "free_y") +
  labs(y = "Relative gene expression level", x = "Genes", fill = "Oxygen\nconc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                             "2%" = "#fdae61",  
                             "1%" = "#d7191c")) + 
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines"))



gsc_expression_summary = gsc_expression %>% 
  group_by(cell_line, gene, oxygen_concentration) %>% 
  summarize(relative_ge_avg = mean(relative_gene_expression),
            relative_ge_sem = sem(relative_gene_expression))

pdf("results/Fig4/SuppFig8a-qPCR-expression.pdf", width = 7.5, height = 4, useDingbats = FALSE)
ggplot(gsc_expression, aes(x=gene, y=relative_gene_expression, fill= oxygen_concentration)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cell_line, scales = "free_y") +
  labs(y = "Relative gene expression level", x = "Genes", fill = "Oxygen\nconc.") +
  scale_fill_manual(values=c("21%" = "#abd9e9",
                             "2%" = "#fdae61",  
                             "1%" = "#d7191c")) + 
  plot_theme +
  theme(panel.spacing.x = unit(1.5, "lines")) 
dev.off()


## Tukey's test:
########### HF2354 ############
gsc_expression_hf2354 = gsc_expression %>% 
  filter(cell_line == "HF2354")
gsc_expression_hf2354_sox2 = gsc_expression_hf2354 %>% 
  filter(gene == "SOX2")
gsc_expression_hf2354_pou5f1 = gsc_expression_hf2354 %>% 
  filter(gene == "POU5F1")
gsc_expression_hf2354_jun = gsc_expression_hf2354 %>% 
  filter(gene == "JUN")
gsc_expression_hf2354_hif2a = gsc_expression_hf2354 %>% 
  filter(gene == "HIF2A")
gsc_expression_hf2354_vegfa = gsc_expression_hf2354 %>% 
  filter(gene == "VEGFA")

SOX2.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf2354_sox2)
SOX2.av <- aov(SOX2.lm)
summary(SOX2.av)
SOX2.tukey.test <- TukeyHSD(SOX2.av)
SOX2.tukey.test

POU5F1.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf2354_pou5f1)
POU5F1.av <- aov(POU5F1.lm)
summary(POU5F1.av)
POU5F1.tukey.test <- TukeyHSD(POU5F1.av)
POU5F1.tukey.test

JUN.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf2354_jun)
JUN.av <- aov(JUN.lm)
summary(JUN.av)
JUN.tukey.test <- TukeyHSD(JUN.av)
JUN.tukey.test

HIF2A.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf2354_hif2a)
HIF2A.av <- aov(HIF2A.lm)
summary(HIF2A.av)
HIF2A.tukey.test <- TukeyHSD(HIF2A.av)
HIF2A.tukey.test

VEGFA.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf2354_vegfa)
VEGFA.av <- aov(VEGFA.lm)
summary(VEGFA.av)
VEGFA.tukey.test <- TukeyHSD(VEGFA.av)
VEGFA.tukey.test

########### HF3016 ############
gsc_expression_hf3016 = gsc_expression %>% 
  filter(cell_line == "HF3016")
gsc_expression_hf3016_sox2 = gsc_expression_hf3016 %>% 
  filter(gene == "SOX2")
gsc_expression_hf3016_pou5f1 = gsc_expression_hf3016 %>% 
  filter(gene == "POU5F1")
gsc_expression_hf3016_jun = gsc_expression_hf3016 %>% 
  filter(gene == "JUN")
gsc_expression_hf3016_hif2a = gsc_expression_hf3016 %>% 
  filter(gene == "HIF2A")
gsc_expression_hf3016_vegfa = gsc_expression_hf3016 %>% 
  filter(gene == "VEGFA")

SOX2.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf3016_sox2)
SOX2.av <- aov(SOX2.lm)
summary(SOX2.av)
SOX2.tukey.test <- TukeyHSD(SOX2.av)
SOX2.tukey.test

POU5F1.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf3016_pou5f1)
POU5F1.av <- aov(POU5F1.lm)
summary(POU5F1.av)
POU5F1.tukey.test <- TukeyHSD(POU5F1.av)
POU5F1.tukey.test

JUN.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf3016_jun)
JUN.av <- aov(JUN.lm)
summary(JUN.av)
JUN.tukey.test <- TukeyHSD(JUN.av)
JUN.tukey.test

HIF2A.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf3016_hif2a)
HIF2A.av <- aov(HIF2A.lm)
summary(HIF2A.av)
HIF2A.tukey.test <- TukeyHSD(HIF2A.av)
HIF2A.tukey.test

VEGFA.lm <- lm(relative_gene_expression ~ oxygen_concentration, data = gsc_expression_hf3016_vegfa)
VEGFA.av <- aov(VEGFA.lm)
summary(VEGFA.av)
VEGFA.tukey.test <- TukeyHSD(VEGFA.av)
VEGFA.tukey.test

### END ####