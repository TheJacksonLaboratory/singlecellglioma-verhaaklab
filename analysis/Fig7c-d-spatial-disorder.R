##############################################
# Assess DNAme disorder in spatially collected specimens
# Updated: 2021.05.15
# Author: Kevin J.
##################################################

## Objective: Use the FRONTIER data to determine relationship between disorder and spatial location.

########################################
# Necessary packages:
library(tidyverse)
library(minfi)
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

### Load relevant samples from FRONTIER data [entitled `all_data`].
load('/Users/johnsk/Documents/Single-Cell-DNAmethylation/450k/data/frontier/FRONTIER.QC.filtered.normalized.anno.final.Rdata')
epimut_high = readRDS("/Users/johnsk/Documents/Single-Cell-DNAmethylation/results/methylation/epimutation/high-epimutation-genebody-all-scRRBS.Rds")

### Relabel the tumor samples so that they can be grouped together.
all_data$Sample_Type[all_data$Sample_Type %in% c("Initial", "Recurrence", "Recurrence2", "Recurrence3")] = "Sample"

### Retrieve the already processed, normalized beta-values.
all_b = getBeta(all_data)
annot_b <- getAnnotation(all_data)

## Restrict to the high-epimutation genes only.
high_epimut_cpgs = which(annot_b$UCSC_RefGene_Name%in%epimut_high$Associated.Gene.Name)
all_b = all_b[high_epimut_cpgs, ]

## Plot density values:
densityPlot(all_b, sampGroups = all_data$Dataset)

## The goal is to reduce any intermediate DNA methylation value into the same identifier.
## Reducing the allele-specific methylation burden being between 0.25-0.75, we can estimate the DNA methylation ITH patterns.
all_b[] <- vapply(all_b, function(x) ifelse(x>0.25 & x<0.75, 1, 0), numeric(1))

## Generate an eITH metric that quantifies intermediate methylation 0.25 and 0.75 by summing the number of intermediate alleles
## divided by the total number of probes measured. The result is a metric where increasing values represent increasing ITH.
eITH <- as.data.frame(colSums(all_b)/dim(all_b)[1])
eITH$idat <- rownames(eITH)
colnames(eITH)[1] <- "eITH"

##############################
### Apply mDNAsi to FRONTIER
##############################
load("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/methylation/Files_OneClassModel_450K/pcbc-stemsig.p219.Rda")
w <- mm$w
w_filt = w[names(w)%in%rownames(all_b)]
X <- all_b[as.character(names(w_filt)),]

# Convert `X` into a matrix
X <- as.matrix(X)

# Score via linear model. The resulting variable `ss` will be row vector 9627 in length 
ss <- t(w_filt) %*% X

## Scale the scores into a ratio from 0 to 1 and store as data frame. 
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(t(ss))
colnames(ss) <- "mDNAsi"   
ss$Sentrix_Accession <- rownames(ss)

## There's a need to demonstrate that this method would provide higher values in tumor compared with normal tissues.
meta_data = pData(all_data) %>% 
  as.data.frame() %>% 
  left_join(eITH, by=c("Sentrix_Accession"="idat")) %>%  
  left_join(ss, by="Sentrix_Accession")

## Produce a boxplot representation of the eITH across DKFZ normal tissues and tumor samples.
ggplot(meta_data, aes(x=Dataset, y=eITH)) + geom_boxplot()
ggplot(meta_data, aes(x=Dataset, y=mDNAsi)) + geom_boxplot()
ggplot(meta_data, aes(x=eITH, y=mDNAsi)) + geom_point()

## Open coordinates for the spatial mapping
imaging_dat = read.table("/Users/johnsk/Documents/Single-Cell-DNAmethylation/data/clinical/FRONTIER.imaging-K2.081419.csv", sep = ",", header = T)


## IDH-wt tumors:
meta_data_wt = meta_data %>% 
  inner_join(imaging_dat, by="Sentrix_Accession") %>% 
  filter(IHD_IHC=="WT", !is.na(Dist_to_CE_surface)) 
ggplot(meta_data_wt, aes(x=Location, y=eITH)) + geom_boxplot() +
  facet_grid(~Patient)

## IDHwt:
## VUmc-02, VUmc-07, VUmc-08, Vumc-11, Vumc-13, and Vumc-17
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="VUmc-02"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="VUmc-02"], method="s")
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="VUmc-07"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="VUmc-07"], method="s")
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="VUmc-08"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="VUmc-08"], method="s")
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="Vumc-11"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="Vumc-11"], method="s")
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="Vumc-13"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="Vumc-13"], method="s")
cor.test(meta_data_wt$eITH[meta_data_wt$Patient=="Vumc-17"], meta_data_wt$Dist_to_CE_surface[meta_data_wt$Patient=="Vumc-17"], method="s")
cor.test(meta_data_wt$eITH, meta_data_wt$Dist_to_CE_surface, method="s")

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig7/Fig7c-IDHwt-spatial.pdf", width = 4, height = 3, useDingbats = FALSE, bg="transparent")
ggplot(meta_data_wt, aes(x=Dist_to_CE_surface, y=eITH, color=Patient)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  plot_theme +
  theme(legend.position="bottom") +
  labs(x="Dist. to contrast enhancement surface (mm)", y="DNAme disorder", color="IDHwt patients")
dev.off()
cor.test(meta_data_wt$eITH, meta_data_wt$Dist_to_CE_surface, method="s")

eITH_model <- glm(meta_data_wt$eITH~meta_data_wt$Dist_to_CE_surface+meta_data_wt$Patient, family = "gaussian")
summary(eITH_model)

ggplot(meta_data_wt, aes(x=Normalized_dist_to_CE_surface, y=eITH, color=Patient)) + geom_point() +
  geom_smooth(method = "lm") 
cor.test(meta_data_wt$eITH, meta_data_wt$Normalized_dist_to_CE_surface, method="s")


eITH_model <- glm(meta_data_wt$eITH~meta_data_wt$Normalized_dist_to_CE_surface+meta_data_wt$Patient, family = "gaussian")
summary(eITH_model)

## IDHmut tumors:
meta_data_mut = meta_data %>% 
  inner_join(imaging_dat, by="Sentrix_Accession") %>% 
  filter(IHD_IHC=="MT", !is.na(Dist_to_nCE_surface))
ggplot(meta_data_mut, aes(x=Location, y=eITH)) + geom_boxplot() +
  facet_grid(~Patient)

## VUmc-01 VUmc-04 VUmc-05 VUmc-06 VUmc-09 Vumc-10 Vumc-12 Vumc-15 
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="VUmc-01"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="VUmc-01"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="VUmc-04"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="VUmc-04"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="VUmc-05"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="VUmc-05"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="VUmc-06"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="VUmc-06"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="VUmc-09"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="VUmc-09"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="Vumc-10"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="Vumc-10"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="Vumc-12"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="Vumc-12"], method="s")
cor.test(meta_data_mut$eITH[meta_data_mut$Patient=="Vumc-15"], meta_data_mut$Dist_to_nCE_surface[meta_data_mut$Patient=="Vumc-15"], method="s")

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig7/Fig7d-IDHmut-spatial.pdf", width = 4, height = 3, useDingbats = FALSE, bg="transparent")
ggplot(meta_data_mut, aes(x=Dist_to_nCE_surface, y=eITH, color=Patient)) + 
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  plot_theme +
  theme(legend.position="bottom") +
  labs(x="Distance to non-enhancement surface (mm)", y="Bulk DNAme disorder", color="IDHmut patients")
dev.off()

pdf(file = "/Users/johnsk/Documents/Single-Cell-DNAmethylation/github/results/Fig7/Fig7d-IDHmut-spatial-legend.pdf", width = 6, height = 3, useDingbats = FALSE, bg="transparent")
ggplot(meta_data_mut, aes(x=Dist_to_nCE_surface, y=eITH, color=Patient)) + 
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  plot_theme +
  theme(legend.position="bottom") +
  labs(x="Distance to non-enhancement surface (mm)", y="Bulk DNAme disorder", color="IDHmut patients")
dev.off()
cor.test(meta_data_mut$eITH, meta_data_mut$Dist_to_nCE_surface, method="s")

ggplot(meta_data_mut, aes(x=Normalized_dist_to_nCE_surface, y=eITH, color=Patient)) + geom_point() +
  geom_smooth(method = "lm") 
cor.test(meta_data_mut$eITH, meta_data_mut$Normalized_dist_to_nCE_surface, method="s")

eITH_model <- glm(meta_data_mut$eITH~meta_data_mut$Dist_to_nCE_surface+meta_data_mut$Patient, family = "gaussian")
summary(eITH_model)

eITH_model <- glm(meta_data_mut$eITH~meta_data_mut$Normalized_dist_to_nCE_surface+meta_data_mut$Patient, family = "gaussian")
summary(eITH_model)

#### END ####








